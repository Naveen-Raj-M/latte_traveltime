!
! © 2024. Triad National Security, LLC. All rights reserved.
!
! Modified version with robust trial management
!

module inversion_step_size

    use libflit
    use parameters
    use vars
    use utility
    use gradient, only: PHASE_BASELINE, PHASE_TRIAL, PHASE_FINAL, compute_gradient_shots

    implicit none

#ifdef dim2
#define model_dimension dimension(:, :)
#endif

#ifdef dim3
#define model_dimension dimension(:, :, :)
#endif

    private

    ! Maximum number of step size search
    integer :: nsearch_max = 6
    logical :: step_suitable
    
    ! Trial tracking
    integer :: current_trial_number
    character(len=1024) :: dir_trial_base

    public :: compute_step_size

contains

    !
    !> Get directory for a specific trial
    !
    function dir_trial(trial_num) result(dir)
        integer, intent(in) :: trial_num
        character(len=:), allocatable :: dir
        character(len=1024) :: tmpdir
        
        tmpdir = tidy(dir_trial_base)//'/trial_'//num2str(trial_num)
        allocate(character(len=len_trim(tmpdir)) :: dir)
        dir = tidy(tmpdir)
    end function dir_trial

    !
    !> Initialize trial directory structure
    !
    subroutine init_trial_structure()
        
        dir_trial_base = tidy(dir_scratch)//'/iteration_'//num2str(iter)//'/step_size_trials'
        
        if (rankid == 0) then
            call make_directory(dir_trial_base)
            
            ! Create a log file for this iteration's trials
            open(99, file=tidy(dir_trial_base)//'/trial_log.txt', status='replace')
            write(99, '(a)') '# Trial log for iteration '//num2str(iter)
            write(99, '(a)') '# Format: trial_number step_size misfit status'
            close(99)
        end if
        
        call mpibarrier
        current_trial_number = 0
        
    end subroutine init_trial_structure

    !
    !> Log trial information
    !
    subroutine log_trial(trial_num, step_val, misfit_val, status)
        integer, intent(in) :: trial_num
        real, intent(in) :: step_val, misfit_val
        character(len=*), intent(in) :: status
        
        if (rankid == 0) then
            open(99, file=tidy(dir_trial_base)//'/trial_log.txt', &
                 status='old', position='append')
            write(99, '(i4,1x,es18.10,1x,es18.10,1x,a)') &
                trial_num, step_val, misfit_val, trim(status)
            close(99)
        end if
        
    end subroutine log_trial

    !
    !> Compute misfit for a trial step (writes to trial directory)
    !
    subroutine compute_trial_misfit(trial_num, step_scaling_factor, misfit)

        integer, intent(in) :: trial_num
        real, intent(in) :: step_scaling_factor
        real, intent(inout) :: misfit
        
        character(len=1024) :: trial_dir

        ! Create directory for this trial
        trial_dir = dir_trial(trial_num)
        if (rankid == 0) then
            call make_directory(trial_dir)
        end if
        call mpibarrier

        ! Update model with trial step
        call backup_current_model
        call update_model(step_scaling_factor)

        ! Forward modeling - redirect output to this trial directory
        yn_misfit_only = .true.

        ! Temporarily override dir_synthetic to point to the TRIAL directory
        dir_synthetic = trial_dir

        ! Phase-aware call - use PHASE_TRIAL
        call compute_gradient_shots(PHASE_TRIAL)
        call mpibarrier

        misfit = sum(step_misfit)

        ! ---- Restore globals that we temporarily overrode for the trial ----
        dir_synthetic = ''        ! so baseline/final in the next phases don’t inherit this path
        yn_misfit_only = .false.  ! restore default (if that’s your norm)

        ! Restore current model (before trial)
        call restore_current_model
        call mpibarrier
        
        ! Log this trial
        call log_trial(trial_num, step_scaling_factor, misfit, 'computed')

    end subroutine compute_trial_misfit

    !
    !> Copy accepted trial to final synthetic directory
    !
    subroutine finalize_accepted_trial(trial_num, final_misfit)
        
        integer, intent(in) :: trial_num
        real, intent(in) :: final_misfit
        
        character(len=1024) :: trial_dir, final_dir
        
        trial_dir = dir_trial(trial_num)
        final_dir = dir_iter_synthetic(iter)
        
        if (rankid == 0) then
            call warn(date_time_compact()//' Finalizing trial '//num2str(trial_num) &
                //' as accepted solution')
            call warn(date_time_compact()//' Copying from: '//tidy(trial_dir))
            call warn(date_time_compact()//' Copying to:   '//tidy(final_dir))
            
            ! Copy all synthetic data from trial to final
            call copy_directory(trial_dir, final_dir)
            
            ! Mark in log
            call log_trial(trial_num, -1.0, final_misfit, 'ACCEPTED')
        end if
        
        call mpibarrier
        
    end subroutine finalize_accepted_trial

    !
    !> Alternative: Do final forward model with accepted step
    !
    subroutine compute_final_forward_model(step_scaling_factor, final_misfit)
        
        real, intent(in) :: step_scaling_factor
        real, intent(inout) :: final_misfit
        
        if (rankid == 0) then
            call warn(date_time_compact()//' Computing final forward model with accepted step')
            call warn(date_time_compact()//' Step size: '//num2str(step_scaling_factor, '(es)'))
        end if
        
        ! Update model to final state
        call update_model(step_scaling_factor)
        
        ! Forward modeling - write to final iteration directory
        yn_misfit_only = .true.
        dir_synthetic = dir_iter_synthetic(iter)   ! gradient_tloc will force this anyway for "Final", but keep explicit

        call compute_gradient_shots(PHASE_FINAL)
        call mpibarrier

        
        final_misfit = sum(step_misfit)
        
        if (rankid == 0) then
            call warn(date_time_compact()//' Final misfit: '//num2str(final_misfit, '(es)'))
            
            ! Log final forward model
            open(99, file=tidy(dir_trial_base)//'/trial_log.txt', &
                 status='old', position='append')
            write(99, '(a,es18.10,a,es18.10)') '# FINAL_FORWARD: step=', &
                step_scaling_factor, ' misfit=', final_misfit
            close(99)
        end if
        
        call mpibarrier
        
    end subroutine compute_final_forward_model

    !
    !> Read synthetic data from a specific trial
    !
    function read_trial_synthetic(trial_num, shot_id, component) result(data)
        
        integer, intent(in) :: trial_num, shot_id
        character(len=*), intent(in) :: component
        real, allocatable, dimension(:, :) :: data
        
        character(len=1024) :: filename, trial_dir
        integer :: nr_local, ncol
        
        trial_dir = dir_trial(trial_num)
        
        if (yn_dd_no_st0) then
            ! For double-difference without st0
            filename = tidy(trial_dir)//'/t'//tidy(component)//'_all.bin'
            nr_local = nr_virtual
            ncol = ns
        else
            ! For standard case
            filename = tidy(trial_dir)//'/shot_'//num2str(shot_id) &
                //'_traveltime_'//tidy(component)//'.bin'
            nr_local = gmtr(ishot)%nr
            ncol = nrefl + 1
        end if
        
        if (.not. file_exists(filename)) then
            call warn(date_time_compact()//' ERROR: Trial synthetic file not found: ' &
                //tidy(filename))
            call mpibarrier
            call mpistop
        end if
        
        data = load(filename, nr_local, ncol)
        
    end function read_trial_synthetic

    !
    !> Check if search step size is within the desired range
    !
    subroutine check_step_range_single_model(step_model, &
            step_scale, srch, step_max)

        real, model_dimension, intent(in) :: srch
        real, intent(in) :: step_model, step_scale, step_max

        if (step_max == 0) then
            step_suitable = step_suitable .and. .true.
        else
            if (maxval(abs(step_model*srch*step_scale)) < step_max) then
                step_suitable = step_suitable .and. .true.
            else
                step_suitable = step_suitable .and. .false.
            end if
        end if

    end subroutine check_step_range_single_model

    !
    !> Compute the initial step size based on search direction values
    !
    subroutine initial_step_size(cij, srch, init_step, max_perturbation)

        real, model_dimension, intent(in) :: cij, srch
        real, intent(inout) :: init_step
        real, intent(in) :: max_perturbation

        real :: max_srch, max_cij

        max_cij = 0.0
        max_srch = 0.0

        max_srch = maxval(abs(srch))
        max_cij = maxval(abs(cij))

        if (max_srch == 0) then
            init_step = 0.0
            return
        end if

        init_step = max_perturbation/max_srch

    end subroutine initial_step_size

    !
    !> Get min max of model purturbation
    !
    subroutine min_max_step_size(step, srch, step_scaling_factor, mname)

        real, model_dimension, intent(in) :: srch
        real, intent(in) :: step, step_scaling_factor
        character(len=*) :: mname

        real, allocatable, model_dimension :: deltam
        real :: min_deltam, max_deltam

        deltam = srch*step*step_scaling_factor

        min_deltam = minval(deltam)
        max_deltam = maxval(deltam)

        if (rankid == 0) then
            call warn(date_time_compact()//' Perturbation '//tidy(mname)//' range: ' &
                //num2str(min_deltam, '(es12.4)')//' ' &
                //num2str(max_deltam, '(es12.4)'))
        end if

    end subroutine min_max_step_size

    !
    !> Output updated model
    !
    subroutine output_updated_model

        integer :: i

        do i = 1, nmodel
            call output_array(model_m(i)%array, &
                dir_iter_model(iter)//'/updated_'//tidy(model_name(i))//'.bin')
        end do

    end subroutine output_updated_model

    !
    !> Check if search step size is within the desired range
    !
    subroutine check_step_range(s)

        real, intent(in) :: s
        integer :: i

        step_suitable = .true.

        do i = 1, nmodel
            call check_step_range_single_model(model_step(i), &
                s, model_srch(i)%array, model_step_max(i))
        end do

    end subroutine check_step_range

    !
    !> Backup current model
    !
    subroutine backup_current_model

        integer :: i

        do i = 1, nmodel
            model_m_backup(i)%array = model_m(i)%array
        end do

    end subroutine backup_current_model

    !
    !> Restore current model
    !
    subroutine restore_current_model

        integer :: i

        do i = 1, nmodel
            model_m(i)%array = model_m_backup(i)%array
        end do

    end subroutine restore_current_model

    !
    !> Update model
    !
    subroutine update_model(step)

        real, intent(in) :: step
        integer :: i

        do i = 1, nmodel

            model_m(i)%array = model_m(i)%array + model_step(i)*model_srch(i)%array*step

            ! Box-clip model
            if (model_m(i)%name == 'vp' .or. model_m(i)%name == 'vs') then
                where (model_m(i)%array /= 0)
                    model_m(i)%array = clip(model_m(i)%array, model_min(i), model_max(i))
                end where
            else
                model_m(i)%array = clip(model_m(i)%array, model_min(i), model_max(i))
            end if

            if (model_m(i)%name == 'refl') then
                model_m(i)%array = model_srch(i)%array
            end if

        end do

        call clip_vpvsratio

    end subroutine update_model

    !
    !> Print perturbation information
    !
    subroutine print_perturb_info

        integer :: i

        do i = 1, nmodel
            call min_max_step_size(model_step(i), model_srch(i)%array, &
                step_scaling_factor, model_name(i))
        end do

    end subroutine print_perturb_info

    !
    !> Calcualte initial step size for all parametes
    !
    subroutine compute_initial_step_size

        integer :: i
        real :: valstep

        do i = 1, nmodel

            select case (remove_string_after(model_name(i), ['c', 'C']))
                case ('vp', 'vs', 'rho')
                    valstep = 100.0
                case ('epsilon', 'delta', 'gamma', 'eps', 'del', 'gam')
                    valstep = 0.1
                case ('c', 'C')
                    valstep = 1.0e9
                case ('refl')
                    valstep = 1.0
                case ('sx')
                    valstep = 0.1*(nx - 1)*dx
                case ('sy')
                    valstep = 0.1*(ny - 1)*dy
                case ('sz')
                    valstep = 0.1*(nz - 1)*dz
                case ('st0')
                    valstep = 0.1
            end select

            call readpar_xfloat(file_parameter, 'step_max_'//tidy(model_name(i)), &
                model_step_max(i), valstep, iter*1.0)
            model_step_max(i) = model_step_max(i)*step_max_scale_factor
            call initial_step_size(model_m(i)%array, model_srch(i)%array, &
                model_step(i), model_step_max(i))

        end do

    end subroutine compute_initial_step_size

    !
    !> Calculate optimal step size coefficient for some parameter
    !
    subroutine compute_step_coef_single_component(srcindex, component_name, &
        trial_num, sum1, sum2)

        integer, intent(in) :: srcindex, trial_num
        character(len=*) :: component_name
        real, intent(inout) :: sum1, sum2

        character(len=1024) :: file_recorded, file_synthetic_prev
        integer :: nr, i, j
        real, allocatable, dimension(:, :) :: seis_obs, seis_syn_trial, seis_syn_prev
        real, allocatable, dimension(:, :) :: weight
        real, allocatable, dimension(:, :) :: r1, r2, m1, m2
        real, allocatable, dimension(:, :) :: tobs_all, tsyn_trial_all, tsyn_prev_all

        if (yn_dd_no_st0) then

            ! For double-difference without st0
            file_recorded = dir_iter_record(iter)//'/t'//tidy(component_name)//'_all.bin'
            file_synthetic_prev = dir_iter_synthetic(iter)//'/t'//tidy(component_name)//'_all.bin'
            
            ! Read from trial directory
            tsyn_trial_all = read_trial_synthetic(trial_num, -1, component_name)
            tobs_all = load(file_recorded, nr_virtual, ns)
            tsyn_prev_all = load(file_synthetic_prev, nr_virtual, ns)

            nr = gmtr(srcindex)%nr
            r1 = zeros(nr, 1)
            r2 = zeros(nr, 1)
            
            !$omp parallel do private(i, j)
            do j = 1, nr
                if (maxval(tsyn_trial_all(j, :)) > 0 .and. &
                    maxval(tsyn_prev_all(j, :)) > 0 .and. &
                    maxval(tobs_all(j, :)) > 0) then
                    do i = 1, ns
                        if (tobs_all(j, srcindex) >= 0 .and. tobs_all(j, i) >= 0) then
                            r1(j, 1) = r1(j, 1) + &
                                (tsyn_prev_all(j, srcindex) - tsyn_prev_all(j, i)) - &
                                (tsyn_trial_all(j, srcindex) - tsyn_trial_all(j, i))
                            r2(j, 1) = r2(j, 1) + &
                                (tsyn_prev_all(j, srcindex) - tsyn_prev_all(j, i)) - &
                                (tobs_all(j, srcindex) - tobs_all(j, i))
                        end if
                    end do
                end if
            end do
            !$omp end parallel do

        else

            ! Standard case
            file_recorded = dir_iter_record(iter)//'/shot_'//num2str(gmtr(srcindex)%id) &
                //'_traveltime_'//tidy(component_name)//'.bin'
            file_synthetic_prev = dir_iter_synthetic(iter)//'/shot_' &
                //num2str(gmtr(srcindex)%id)//'_traveltime_'//tidy(component_name)//'.bin'

            ! Read data
            nr = gmtr(srcindex)%nr
            seis_obs = load(file_recorded, nr, nrefl + 1)
            seis_syn_prev = load(file_synthetic_prev, nr, nrefl + 1)
            
            ! Read from trial directory
            seis_syn_trial = read_trial_synthetic(trial_num, gmtr(srcindex)%id, component_name)

            call traveltime_residual(seis_syn_trial, seis_syn_prev, r1, m1)
            call traveltime_residual(seis_obs, seis_syn_prev, r2, m2)

        end if

        weight = zeros(nr, nrefl + 1)
        !$omp parallel do private(i, j)
        do i = 1, nr
            do j = 1, nrefl + 1

                weight(i, j) = gmtr(srcindex)%recr(i)%weight*misfit_weight(j)

                if (abs(r1(i, j)) > misfit_threshold .or. &
                    abs(r2(i, j)) > misfit_threshold) then
                    weight(i, j) = 0
                end if

                if (j > 1 .and. (gmtr(srcindex)%recr(i)%aoff < offset_min_refl &
                        .or. gmtr(srcindex)%recr(i)%aoff > offset_max_refl)) then
                    weight(i, j) = 0
                end if

            end do
        end do
        !$omp end parallel do

        sum1 = sum1 + sum(abs(r1)*abs(r2)*weight)
        sum2 = sum2 + sum(abs(r1)**2*weight)

    end subroutine compute_step_coef_single_component

    !
    !> Calculate optimal step size coefficient for all parameters
    !
    subroutine compute_step_coef(srcindex, trial_num, sum1, sum2)

        integer, intent(in) :: srcindex, trial_num
        real, intent(inout) :: sum1, sum2

        integer :: i

        do i = 1, ndata
            call compute_step_coef_single_component(srcindex, data_name(i), &
                trial_num, sum1, sum2)
        end do

    end subroutine compute_step_coef

    !
    !> Compute the optimal step size for a quasi-linear inversion
    !
    subroutine compute_step_size_linear

        real :: trial_misfit, sum1, sum2, data_misfit_current, best_misfit, best_step
        integer :: cnt, i
        logical :: use_final_forward_model
        integer :: best_trial_num

        if (rankid == 0) then
            call warn(date_time_compact()//' >>>>>>>>>> Computing optimal step size...')
        end if

        ! Initialize trial structure
        call init_trial_structure()

        ! Calculate the initial step size
        call compute_initial_step_size

        ! Default initial step scalar
        step_scaling_factor = 0.2

        ! Check if the initial step scalar is suitable
        step_suitable = .false.
        cnt = 1
        do while (.not. step_suitable .and. cnt < 20)
            call check_step_range(step_scaling_factor)
            step_scaling_factor = step_scaling_factor/2.0
            cnt = cnt + 1
        end do
        
        if (.not. step_suitable) then
            if (rankid == 0) then
                call warn(date_time_compact() &
                    //' Error: Cannot find a suitable initial step size. Exit. ')
            end if
            call mpibarrier
            call mpiend
        else
            if (rankid == 0) then
                call warn(date_time_compact()//' Initial step scaling factor = ' &
                    //num2str(step_scaling_factor, '(es)'))
            end if
        end if

        ! Step size coefficients
        sum1 = 0.0
        sum2 = 0.0

        ! Backup current model
        call backup_current_model

        ! Trial 1: Compute forward model with trial step
        current_trial_number = 1
        call compute_trial_misfit(current_trial_number, step_scaling_factor, trial_misfit)

        ! Compute optimal step coefficient using trial 1
        call mpibarrier
        do i = shot_in_rank(rankid, 1), shot_in_rank(rankid, 2)
            call compute_step_coef(i, current_trial_number, sum1, sum2)
        end do
        call mpibarrier

        ! Restore current iteration model
        call restore_current_model

        ! Calculate optimal step size
        call mpibarrier
        call allreduce(sum1)
        call allreduce(sum2)
        
        if (isnan(sum1) .or. isnan(sum2)) then
            call warn(date_time_compact()//' Step size coef sum 1 = '//num2str(sum1, '(es)'))
            call warn(date_time_compact()//' Step size coef sum 2 = '//num2str(sum2, '(es)'))
            call mpibarrier
            call mpiend
        else
            step_scaling_factor = sum1/(sum2 + float_tiny)
        end if
        
        if (rankid == 0) then
            call warn(date_time_compact()//' Step size coef sum 1 = '//num2str(sum1, '(es)'))
            call warn(date_time_compact()//' Step size coef sum 2 = '//num2str(sum2, '(es)'))
            call warn(date_time_compact()//' Optimal step size coef = ' &
                //num2str(step_scaling_factor, '(es)'))
        end if

        if (.not. (step_scaling_factor <= float_huge)) then
            if (rankid == 0) then
                call warn(date_time_compact() &
                    //' <compute_step_size_linear> Error: Step size is NaN. Exit. ')
                call mpi_abort(mpi_comm_world, mpi_err_other, mpi_ierr)
            end if
        end if

        ! Check step size with line search
        step_scaling_factor = step_scaling_factor*2.0
        data_misfit_current = data_misfit(iter)
        trial_misfit = data_misfit_current*2.0
        cnt = 2

        ! Read control flag
        call readpar_xlogical(file_parameter, 'yn_enforce_update', &
            yn_enforce_update, .false., iter*1.0)
        call readpar_logical(file_parameter, 'use_final_forward_model', &
            use_final_forward_model, .true.)

        best_trial_num = -1

        if (yn_enforce_update) then

            ! Enforce update regardless of misfit increase
            do while (cnt <= 2)

                step_scaling_factor = step_scaling_factor/2.0

                call check_step_range(step_scaling_factor)
                if (.not. step_suitable) then
                    cycle
                end if

                current_trial_number = cnt
                call compute_trial_misfit(current_trial_number, &
                    step_scaling_factor, trial_misfit)

                call print_step_info(current_trial_number, step_scaling_factor, trial_misfit)
                
                best_trial_num = current_trial_number

                cnt = cnt + 1

            end do

        else

            ! Conventional trial-error approach
            ! Track the BEST trial across all attempts
            
            best_misfit = data_misfit_current
            best_step = 0.0
            best_trial_num = -1
            
            do while (cnt < nsearch_max + 1)

                step_scaling_factor = step_scaling_factor/2.0

                call check_step_range(step_scaling_factor)
                if (.not. step_suitable) then
                    cycle
                end if

                current_trial_number = cnt
                call compute_trial_misfit(current_trial_number, &
                    step_scaling_factor, trial_misfit)

                call print_step_info(current_trial_number, step_scaling_factor, trial_misfit)

                ! Track BEST trial (not just first improvement)
                if (trial_misfit < best_misfit) then
                    best_trial_num = current_trial_number
                    best_misfit = trial_misfit
                    best_step = step_scaling_factor
                    if (rankid == 0) then
                        call warn(date_time_compact()//' >> New best trial: ' &
                            //num2str(best_trial_num)//' with misfit=' &
                            //num2str(best_misfit, '(es)'))
                    end if
                end if

                cnt = cnt + 1
                
                ! Optional: early exit if we found a good improvement
                if (trial_misfit < 0.99 * data_misfit_current) then
                    if (rankid == 0) then
                        call warn(date_time_compact() &
                            //' >> Sufficient improvement found, stopping trials')
                    end if
                    exit
                end if

            end do
            
            ! Use the BEST trial found
            trial_misfit = best_misfit
            step_scaling_factor = best_step

        end if

        ! Decision: accept or reject
        if (best_trial_num < 0 .or. best_misfit >= data_misfit_current) then
            ! REJECTED: No trial better than baseline
            step_scaling_factor = 0.0
            data_misfit(iter) = data_misfit_current
            
            if (rankid == 0) then
                call warn(date_time_compact()//' !! No trial accepted, step size = 0')
                call warn(date_time_compact() &
                    //' !! Keeping baseline synthetic in dir_synthetic(iter='//num2str(iter)//')')
                call log_trial(-1, 0.0, data_misfit_current, 'REJECTED_ALL')
            end if
            
        else
            ! ACCEPTED: Trial better than baseline
            
            if (use_final_forward_model) then
                ! Option 1B: Do final forward model with accepted step
                if (rankid == 0) then
                    call warn(date_time_compact() &
                        //' >> Using final forward model approach')
                end if
                
                call compute_final_forward_model(step_scaling_factor, data_misfit(iter))
                
            else
                ! Option 1A: Copy best trial to final directory
                if (rankid == 0) then
                    call warn(date_time_compact()//' >> Copying best trial to final')
                end if
                
                call finalize_accepted_trial(best_trial_num, trial_misfit)
                data_misfit(iter) = trial_misfit
            end if
            
        end if

        ! Print perturbation information
        call print_perturb_info

        call mpibarrier

        if (rankid == 0) then
            call warn(date_time_compact()//' >>>>>>>>>> Step size computation is done.')
        end if

        ! Update model and output
        if (step_scaling_factor > 0) then
            call update_model(step_scaling_factor)
        end if
        
        if (rankid == 0) then
            call output_updated_model
        end if

        ! Optional cleanup of trial directories
        if (clear_scratch) then
            if (rankid == 0) then
                call warn(date_time_compact()//' Cleaning trial directories in scratch...')
                call system('find ' // trim(dir_scratch) // ' -type d -name "trial_*" -exec rm -rf {} + 2>/dev/null')
            end if
            call mpibarrier
        end if


    end subroutine compute_step_size_linear

    !
    !> Calculate the optimal step size (quadratic search)
    !
    subroutine compute_step_size_quadratic

        integer :: cnt
        real :: step0, misfit0, step1, misfit1, step2, misfit2
        logical :: use_final_forward_model
        integer :: best_trial_num

        if (rankid == 0) then
            call warn(' >>>>>>>>>> Computing step size...')
        end if

        ! Initialize trial structure
        call init_trial_structure()

        ! Calculate the initial step size
        call compute_initial_step_size

        step_scaling_factor = 1.0

        ! Check if the initial step scalar is suitable
        step_suitable = .false.
        cnt = 1
        do while (.not. step_suitable .and. cnt < 20)
            call check_step_range(step_scaling_factor)
            step_scaling_factor = step_scaling_factor/2.0
            cnt = cnt + 1
        end do
        
        if (.not. step_suitable) then
            if (rankid == 0) then
                call warn(date_time_compact() &
                    //' Error: Cannot find a suitable initial step size. Exit. ')
            end if
            call mpibarrier
            call mpiend
        else
            if (rankid == 0) then
                call warn(date_time_compact()//' Initial step scaling factor = ' &
                    //num2str(step_scaling_factor, '(es)'))
            end if
        end if

        ! Trial 0: baseline (already computed)
        step0 = 0.0
        misfit0 = data_misfit(iter)
        call print_step_info(0, step0, misfit0)
        call mpibarrier

        ! Trial 1
        step1 = 0.1*step_scaling_factor
        current_trial_number = 1
        call compute_trial_misfit(current_trial_number, step1, misfit1)
        call print_step_info(1, step1, misfit1)
        call mpibarrier

        ! Trial 2
        step2 = step_scaling_factor
        current_trial_number = 2
        call compute_trial_misfit(current_trial_number, step2, misfit2)
        call print_step_info(2, step2, misfit2)
        call mpibarrier

        ! Quadratic fit
        step_scaling_factor = &
            0.5*((misfit0 - misfit2)*step1**2 + (-misfit0 + misfit1)*step2**2) &
            /(-(misfit2*step1) + misfit0*(step1 - step2) + misfit1*step2)
        
        ! Trial 3: optimal from quadratic fit
        current_trial_number = 3
        call compute_trial_misfit(current_trial_number, step_scaling_factor, data_misfit(iter))
        call print_step_info(3, step_scaling_factor, data_misfit(iter))
        call mpibarrier

        ! Read control flag
        call readpar_logical(file_parameter, 'use_final_forward_model', &
            use_final_forward_model, .true.)

        best_trial_num = 3

        ! Finalize
        if (use_final_forward_model) then
            if (rankid == 0) then
                call warn(date_time_compact()//' >> Using final forward model approach')
            end if
            call compute_final_forward_model(step_scaling_factor, data_misfit(iter))
        else
            if (rankid == 0) then
                call warn(date_time_compact()//' >> Copying best trial to final')
            end if
            call finalize_accepted_trial(best_trial_num, data_misfit(iter))
        end if

        call mpibarrier

        if (rankid == 0) then
            call warn(date_time_compact()//' >>>>>>>>>> Step size computation is done. ')
        end if

        ! Update model and output
        call update_model(step_scaling_factor)
        if (rankid == 0) then
            call output_updated_model
        end if

        ! Optional cleanup of trial directories
        if (clear_scratch) then
            if (rankid == 0) then
                call warn(date_time_compact()//' Cleaning trial directories in scratch...')
                call system('find ' // trim(dir_scratch) // ' -type d -name "trial_*" -exec rm -rf {} + 2>/dev/null')
            end if
            call mpibarrier
        end if


    end subroutine compute_step_size_quadratic

    !
    !> Calculate the optimal step size (line search with quadratic + bisection)
    !
    subroutine compute_step_size_linesearch

        integer :: cnt
        real :: al, ar, ac, fl, fr, fc, a0
        real :: acr, arl, alc, bcr, brl, blc, num, den, an, fn, a2, f2
        real :: misfit0, step1, misfit1, step2, misfit2
        logical :: use_final_forward_model
        integer :: best_trial_num
        real :: best_misfit, best_step

        if (rankid == 0) then
            call warn(' >>>>>>>>>> Computing step size... ')
        end if

        ! Initialize trial structure
        call init_trial_structure()

        ! Calculate the initial step size
        call compute_initial_step_size

        step_scaling_factor = 1.0

        ! Check if the initial step scalar is suitable
        step_suitable = .false.
        cnt = 1
        do while (.not. step_suitable .and. cnt < 20)
            call check_step_range(step_scaling_factor)
            step_scaling_factor = step_scaling_factor/2.0
            cnt = cnt + 1
        end do
        
        if (.not. step_suitable) then
            if (rankid == 0) then
                call warn(date_time_compact() &
                    //' Error: Cannot find a suitable initial step size. Exit. ')
            end if
            call mpibarrier
            call mpiend
        else
            if (rankid == 0) then
                call warn(date_time_compact()//' Initial step scaling factor = ' &
                    //num2str(step_scaling_factor, '(es)'))
            end if
        end if

        a0 = step_scaling_factor

        ! Trial 0: baseline
        al = 0.0
        fl = data_misfit(iter)
        misfit0 = fl
        call print_step_info(0, al, fl)
        call mpibarrier

        ! Trial 1
        ar = a0
        current_trial_number = 1
        call compute_trial_misfit(current_trial_number, ar, fr)
        call print_step_info(1, ar, fr)
        step1 = ar
        misfit1 = fr
        call mpibarrier

        ! If this step produces a smaller misfit, then skip trial
        if (fr < fl) then
            ac = ar
            fc = fr
            best_trial_num = 1
            best_misfit = fr
            best_step = ar
            call mpibarrier
            goto 123
        end if

        ! Trial 2
        ac = 0.5*(al + ar)
        current_trial_number = 2
        call compute_trial_misfit(current_trial_number, ac, fc)
        call print_step_info(2, ac, fc)
        step2 = ac
        misfit2 = fc
        call mpibarrier

        best_trial_num = 2
        best_misfit = fc
        best_step = ac

        ! If this step produces a smaller misfit, then skip trial
        if (fc < fl) then
            call mpibarrier
            goto 123
        end if

        cnt = 3

        call mpibarrier

        ! Quadratic fit + bisection line search
        do while (cnt < nsearch_max)

            if ((fc < fl) .and. (fc < fr)) then

                acr = ac - ar
                bcr = ac**2 - ar**2

                arl = ar - al
                brl = ar**2 - al**2

                alc = al - ac
                blc = al**2 - ac**2

                num = bcr*fl + brl*fc + blc*fr
                den = acr*fl + arl*fc + alc*fr

                if (den == 0.0) then
                    exit
                end if

                an = 0.5*num/den
                current_trial_number = cnt
                call compute_trial_misfit(current_trial_number, an, fn)
                call print_step_info(current_trial_number, an, fn)
                call mpibarrier

                if (fn < best_misfit) then
                    best_trial_num = current_trial_number
                    best_misfit = fn
                    best_step = an
                end if

                cnt = cnt + 1

                if (an > ac) then
                    if (fn >= fc) then
                        ar = an
                        fr = fn
                    else
                        al = ac
                        fl = fc
                        ac = an
                        fc = fn
                    end if
                else
                    if (fn >= fc) then
                        al = an
                        fl = fn
                    else
                        ar = ac
                        fr = fc
                        ac = an
                        fc = fn
                    end if
                end if

            else
                ! Bisection safe switchover

                a2 = ac + 1.0e-1*(ar - al)
                current_trial_number = cnt
                call compute_trial_misfit(current_trial_number, a2, f2)
                call print_step_info(current_trial_number, a2, f2)
                call mpibarrier

                if (f2 < best_misfit) then
                    best_trial_num = current_trial_number
                    best_misfit = f2
                    best_step = a2
                end if

                cnt = cnt + 1

                if (fc < f2) then
                    ar = a2
                    fr = f2
                else
                    al = ac
                    fl = fc
                end if

                ac = 0.5*(al + ar)
                current_trial_number = cnt
                call compute_trial_misfit(current_trial_number, ac, fc)
                call print_step_info(current_trial_number, ac, fc)
                call mpibarrier

                if (fc < best_misfit) then
                    best_trial_num = current_trial_number
                    best_misfit = fc
                    best_step = ac
                end if

                cnt = cnt + 1

            end if

            ! Stop search if range too small
            if (abs(ar - al) <= 0.05) then
                if (fc >= misfit1) then
                    ac = step1
                    fc = misfit1
                    best_trial_num = 1
                    best_misfit = misfit1
                    best_step = step1
                end if
                if (fc >= misfit2) then
                    ac = step2
                    fc = misfit2
                    best_trial_num = 2
                    best_misfit = misfit2
                    best_step = step2
                end if
                exit
            end if

        end do

        123 continue

        ! Read control flag
        call readpar_logical(file_parameter, 'use_final_forward_model', &
            use_final_forward_model, .true.)

        ! Select BEST trial overall (ac may not be the best in bisection cases)
        if (best_misfit >= misfit0 .or. best_trial_num < 0) then
            ! REJECTED: No trial better than baseline
            step_scaling_factor = 0.0
            data_misfit(iter) = misfit0
            
            if (rankid == 0) then
                call warn(date_time_compact()//' !! No trial accepted, step size = 0')
                call warn(date_time_compact() &
                    //' !! Keeping baseline synthetic in dir_synthetic(iter='//num2str(iter)//')')
                call log_trial(-1, 0.0, data_misfit(iter), 'REJECTED_ALL')
            end if
            
        else
            ! ACCEPTED: Use the BEST trial found during search
            step_scaling_factor = best_step   ! Use the actual best step
            
            if (use_final_forward_model) then
                if (rankid == 0) then
                    call warn(date_time_compact()//' >> Using final forward model approach')
                end if
                call compute_final_forward_model(step_scaling_factor, data_misfit(iter))
            else
                if (rankid == 0) then
                    call warn(date_time_compact()//' >> Copying best trial to final')
                end if
                call finalize_accepted_trial(best_trial_num, best_misfit)
                data_misfit(iter) = best_misfit
            end if
        end if

        call mpibarrier

        if (rankid == 0) then
            call warn(date_time_compact()//' >>>>>>>>>> Step size computation is done. ')
        end if

        ! Update model and output
        if (step_scaling_factor > 0) then
            call update_model(step_scaling_factor)
        end if
        
        if (rankid == 0) then
            call output_updated_model
        end if

        ! Optional cleanup of trial directories
        if (clear_scratch) then
            if (rankid == 0) then
                call warn(date_time_compact()//' Cleaning trial directories in scratch...')
                call system('find ' // trim(dir_scratch) // ' -type d -name "trial_*" -exec rm -rf {} + 2>/dev/null')
            end if
            call mpibarrier
        end if


    end subroutine compute_step_size_linesearch

    !
    !> Compute step size (main entry point)
    !
    subroutine compute_step_size

        select case (step_size_method)
            case ('linear')
                call compute_step_size_linear
            case ('quadratic')
                call compute_step_size_quadratic
            case ('line_search')
                call compute_step_size_linesearch
        end select

        call mpibarrier

    end subroutine compute_step_size

end module inversion_step_size
