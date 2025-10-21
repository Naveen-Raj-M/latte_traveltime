!
! © 2024. Triad National Security, LLC. All rights reserved.
!
! This program was produced under U.S. Government contract 89233218CNA000001
! for Los Alamos National Laboratory (LANL), which is operated by
! Triad National Security, LLC for the U.S. Department of Energy/National Nuclear
! Security Administration. All rights in the program are reserved by
! Triad National Security, LLC, and the U.S. Department of Energy/National
! Nuclear Security Administration. The Government is granted for itself and
! others acting on its behalf a nonexclusive, paid-up, irrevocable worldwide
! license in this material to reproduce, prepare. derivative works,
! distribute copies to the public, perform publicly and display publicly,
! and to permit others to do so.
!
! Author:
!    Kai Gao, kaigao@lanl.gov
!


module utility

    use libflit
    use parameters
    use vars

    implicit none

    private

    public :: divide_shots
    public :: dir_iter_model
    public :: dir_iter_synthetic
    public :: dir_iter_record
    public :: print_misfit
    public :: print_step_info
    public :: make_iter_dir
    public :: traveltime_residual
    public :: compute_shot_misfit
    public :: copy_directory

contains

    !
    !> Iteration directory of model
    !
    function dir_iter_model(iter) result(dir)

        !> The iteration number
        integer :: iter
        !> The directory to save gradients, search directions, and updated models
        character(len=:), allocatable :: dir
        character(len=1024) :: tmpdir

        tmpdir = tidy(dir_working)//'/iteration_'//num2str(iter)//'/model'

        allocate (character(len=len_trim(tmpdir)) :: dir)
        dir = tidy(tmpdir)

    end function dir_iter_model

    !
    !> Iteration directory of synthetic data
    !> MODIFIED: Never redirect to scratch during initial forward model
    !
    function dir_iter_synthetic(iter) result(dir)

        !> The iteration number
        integer :: iter
        !> The directory to save synthetic data
        character(len=:), allocatable :: dir
        character(len=1024) :: tmpdir

        ! ALWAYS use the proper iteration directory
        ! The trial management is now handled explicitly in step size routines
        tmpdir = tidy(dir_working)//'/iteration_'//num2str(iter)//'/synthetic'

        allocate (character(len=len_trim(tmpdir)) :: dir)
        dir = tidy(tmpdir)

    end function dir_iter_synthetic

    !
    !> Iteration directory of record data. Here, the directory is same for different iterations,
    !> as there is no operation like source encoding etc.
    !
    function dir_iter_record(iter) result(dir)

        !> The iteration number
        integer :: iter
        !> The directory to save processed observed data
        character(len=:), allocatable :: dir

        integer :: i

        i = iter
        dir = tidy(dir_record)

    end function dir_iter_record

    !
    !> Make necessary directories in each iteration
    !
    subroutine make_iter_dir
        character(len=1024) :: dir_iteration

        dir_iteration = tidy(dir_working)//'/iteration_'//num2str(iter)

        if (rankid == 0) then
            call make_directory(tidy(dir_iteration)//'/synthetic')
            call make_directory(tidy(dir_iteration)//'/model')
        end if

        call mpibarrier
    end subroutine make_iter_dir

    !
    !> Print misfit information
    !
    subroutine print_misfit
        integer :: i, tmpi, tmpr

        if (resume_from_iter == 1 .and. iter == 0) then
            misfit0 = 0.0
            call alloc_array(shot_misfit, [1, ns, 0, niter_max])
            call alloc_array(step_misfit, [1, ns])
            call alloc_array(data_misfit, [0, niter_max + 1])
            open (2, file=tidy(file_datamisfit), status='replace')
            close (2)
        end if

        if (resume_from_iter > 1 .and. iter == 0) then
            call alloc_array(shot_misfit, [1, ns, 0, niter_max])
            call alloc_array(step_misfit, [1, ns])
            call alloc_array(data_misfit, [0, niter_max + 1])
            open (2, file=tidy(file_datamisfit), status='old')
            do i = 0, resume_from_iter - 1
                read (2, *) tmpi, data_misfit(i), tmpr
            end do
            close (2)
            open (3, file=tidy(file_shotmisfit), status='old', access='stream', form='unformatted')
            do i = 1, resume_from_iter - 1
                read (3) shot_misfit(:, i)
            end do
            close (3)
            misfit0 = data_misfit(0)
        end if

        if (rankid == 0 .and. iter >= 1) then
            if (.not. (data_misfit(iter) <= float_huge) .or. isnan(data_misfit(iter))) then
                call warn(date_time_compact()//' Error: Misfit is NaN. ')
                call mpistop
            end if

            if (iter == 1 .and. misfit0 == 0) then
                misfit0 = data_misfit(0)
                open (2, file=tidy(file_datamisfit), position='append')
                write (2, '(i4,es18.10,es18.10)') 0, misfit0, 1.0
                close (2)
            end if

            call warn(date_time_compact()//' >>>>>>>>>> Data misfit: ' &
                //num2str(data_misfit(iter), '(es18.10)'))

            if (data_misfit(0) == 0) then
                call warn(date_time_compact()//' >>>>>>>>>> Normalized data misfit: ' &
                    //num2str(0.0, '(es18.10)'))
                stop
            else
                call warn(date_time_compact()//' >>>>>>>>>> Normalized data misfit: ' &
                    //num2str(data_misfit(iter)/data_misfit(0), '(es18.10)'))
            end if

            open (2, file=tidy(file_datamisfit), status='old', form='formatted', &
                access='direct', recl=41)
            write (2, '(i4,es18.10,es18.10,a)', rec=iter + 1) &
                iter, data_misfit(iter), data_misfit(iter)/data_misfit(0), char(10)
            close (2)

            shot_misfit(:, iter) = step_misfit
            call output_array(shot_misfit(:, 0:iter), tidy(file_shotmisfit))

        end if

        step_misfit = 0.0
        call mpibarrier

        if (iter >= 3) then
            if (data_misfit(iter) == data_misfit(iter - 1) .and. &
                data_misfit(iter - 1) == data_misfit(iter - 2)) then
                step_max_scale_factor = step_max_scale_factor/2.0
                if (yn_flat_stop) then
                    if (rankid == 0) then
                        call warn('')
                        call warn(date_time_compact()//' Inversion has three consecutive equal misfits. Exiting.')
                        call warn('')
                    end if
                    call mpibarrier
                    call mpiend
                end if
            end if
        end if
    end subroutine print_misfit

    !
    !> Print step number information
    !
    subroutine print_step_info(step_number, step, data_misfit)

        integer, intent(in) :: step_number
        real, intent(in) :: step, data_misfit

        if (rankid == 0) then
            if (.not. (data_misfit <= float_huge)) then
                call warn(date_time_compact()//' Error: Misfit is NaN. Inversion stops.')
                call mpistop
            end if
            call warn(' ')
            call warn(date_time_compact()//' ══════════════════════════════════════════════════')
            call warn(date_time_compact()//' Trial Number: '//num2str(step_number))
            call warn(date_time_compact()//' Step Size: '//num2str(step, '(es18.10)'))
            call warn(date_time_compact()//' Misfit: '//num2str(data_misfit, '(es18.10)'))
            call warn(date_time_compact()//' ══════════════════════════════════════════════════')
            call warn(' ')
        end if

    end subroutine print_step_info

    !
    !> Divide shots into different ranks
    !
    subroutine divide_shots

        real, allocatable, dimension(:, :) :: sxyz, rxyz
        type(meta_array2_real), allocatable, dimension(:) :: ttp_real, tts_real, ttp, tts
        integer :: i, j, l, nr
        real, allocatable, dimension(:, :, :) :: t_all
        integer :: nr_max
        
        ! Variables for mapping export
        integer :: mapping_unit, geom_unit, shot_mapping_unit
        character(len=256) :: mapping_filename, geom_filename, shot_mapping_filename
        integer :: matching_count, total_matches

        ! The number of sources might be >> the number of stations/receivers.
        ! In such a case, using reciprocity theorem, we exchange the source and receivers
        ! to reduce computational cost. The observed data must also be exchanged accordingly.
        !
        ! Also, for source location, because at a source location, a singularity exists,
        ! the sources and stations must be exchanged regardless of # of sources and receivers.
        if (yn_exchange_sr .or. which_program == 'tloc') then

            if (rankid == 0) then
                call warn('')
                call warn(date_time_compact()//' Exchanging source and receivers... ')
                
                ! Create comprehensive mapping files for Python analysis inside dir_working
                mapping_filename       = tidy(dir_working)//'/exchange_mapping.txt'
                geom_filename          = tidy(dir_working)//'/exchange_geometry.txt'
                shot_mapping_filename  = tidy(dir_working)//'/shot_mapping.txt'
                
                open(newunit=mapping_unit, file=mapping_filename, status='replace')
                open(newunit=geom_unit, file=geom_filename, status='replace')
                open(newunit=shot_mapping_unit, file=shot_mapping_filename, status='replace')
                
                ! Write headers
                write(mapping_unit, '(A)') '# Exchange Mapping File'
                write(mapping_unit, '(A)') '# Format: virtual_shot_id virtual_receiver_id original_shot_id original_receiver_id travel_time_p travel_time_s'
                write(mapping_unit, '(A)') '# This file maps every virtual shot-receiver pair to its original counterpart'
                
                write(geom_unit, '(A)') '# Exchange Geometry File'
                write(geom_unit, '(A)') '# Virtual Sources (Unique Receiver Positions)'
                write(geom_unit, '(A)') '# Format: virtual_shot_id x_coord y_coord z_coord'
                
                write(shot_mapping_unit, '(A)') '# Shot Index Mapping File'
                write(shot_mapping_unit, '(A)') '# Maps virtual receiver indices to original shot IDs'
                write(shot_mapping_unit, '(A)') '# Format: virtual_receiver_index original_shot_id'
            end if

            ! Find unique receivers
            nr = sum(gmtr(:)%nr)
            rxyz = zeros(nr, 3)
            sxyz = zeros(ns, 3)

            l = 0
            do i = 1, ns
                sxyz(i, :) = [gmtr(i)%srcr(1)%x, gmtr(i)%srcr(1)%y, gmtr(i)%srcr(1)%z]
                do j = 1, gmtr(i)%nr
                    rxyz(l + j, :) = [gmtr(i)%recr(j)%x, gmtr(i)%recr(j)%y, gmtr(i)%recr(j)%z]
                end do
                l = l + gmtr(i)%nr
            end do
            rxyz = unique(rxyz, cols=[1, 2, 3])
            nr = size(rxyz, 1)

            ! Write virtual source geometry (unique receiver positions)
            if (rankid == 0) then
                do i = 1, nr
                    write(geom_unit, '(I0,1X,F12.3,1X,F12.3,1X,F12.3)') i, rxyz(i, 1), rxyz(i, 2), rxyz(i, 3)
                end do
            end if

            allocate (ttp_real(1:ns))
            allocate (tts_real(1:ns))

            ! Every process reads some files, and share
            !$omp parallel do private(i)
            do i = 1, ns
                ttp_real(i)%array = zeros(gmtr(i)%nr, 1)
                tts_real(i)%array = zeros(gmtr(i)%nr, 1)
            end do
            !$omp end parallel do

            if (which_program /= 'eikonal') then

                nr_max = maxval(gmtr(:)%nr)
                t_all = zeros(nr_max, ns, ndata)

                call alloc_array(shot_in_rank, [0, nrank - 1, 1, 2])
                call cut(1, ns, nrank, shot_in_rank)
                !$omp parallel do private(i)
                do i = shot_in_rank(rankid, 1), shot_in_rank(rankid, 2)
                    select case (which_medium)
                        case ('acoustic-iso', 'acoustic-tti')
                            t_all(1:gmtr(i)%nr, i:i, 1) = load(tidy(dir_record)//'/shot_'//num2str(gmtr(i)%id)//'_traveltime_p.bin', gmtr(i)%nr, 1)
                        case ('elastic-iso', 'elastic-tti')
                            t_all(1:gmtr(i)%nr, i:i, 1) = load(tidy(dir_record)//'/shot_'//num2str(gmtr(i)%id)//'_traveltime_p.bin', gmtr(i)%nr, 1)
                            t_all(1:gmtr(i)%nr, i:i, 2) = load(tidy(dir_record)//'/shot_'//num2str(gmtr(i)%id)//'_traveltime_s.bin', gmtr(i)%nr, 1)
                    end select
                end do
                !$omp end parallel do

                call allreduce_array(t_all)

                !$omp parallel do private(i)
                do i = 1, ns
                    select case (which_medium)
                        case ('acoustic-iso', 'acoustic-tti')
                            ttp_real(i)%array = t_all(1:gmtr(i)%nr, i:i, 1)
                        case ('elastic-iso', 'elastic-tti')
                            ttp_real(i)%array = t_all(1:gmtr(i)%nr, i:i, 1)
                            tts_real(i)%array = t_all(1:gmtr(i)%nr, i:i, 2)
                    end select
                end do
                !$omp end parallel do

            end if
            call mpibarrier

            ! Exchange source and receivers
            gmtr_real = gmtr

            deallocate (gmtr)

            allocate (gmtr(1:nr))
            allocate (ttp(1:nr))
            allocate (tts(1:nr))

            ! Write shot mapping: virtual receiver index -> original shot ID
            if (rankid == 0) then
                do i = 1, ns
                    write(shot_mapping_unit, '(I0,1X,I0)') i-1, gmtr_real(i)%id  ! 0-based indexing for Python
                end do
            end if

            ! For each virtual source (i.e., real receiver)
            call alloc_array(shot_in_rank, [0, nrank - 1, 1, 2])
            call cut(1, nr, nrank, shot_in_rank)
            
            total_matches = 0
            
            do i = 1, nr

                gmtr(i)%ns = 1
                gmtr(i)%nr = ns
                gmtr(i)%id = i

                ! Virtual source location = unique real receiver location
                allocate (gmtr(i)%srcr(1:1))
                gmtr(i)%srcr(1)%x = rxyz(i, 1)
                gmtr(i)%srcr(1)%y = rxyz(i, 2)
                gmtr(i)%srcr(1)%z = rxyz(i, 3)
                gmtr(i)%srcr(1)%amp = 1.0

                allocate (gmtr(i)%recr(1:gmtr(i)%nr))

                select case (which_medium)
                    case ('acoustic-iso', 'acoustic-tti')
                        ttp(i)%array = zeros(gmtr(i)%nr, 1)
                    case ('elastic-iso', 'elastic-tti')
                        ttp(i)%array = zeros(gmtr(i)%nr, 1)
                        tts(i)%array = zeros(gmtr(i)%nr, 1)
                end select

                matching_count = 0

                ! For each virtual receiver (i.e., real source), nr_virtual = ns_real
                do j = 1, gmtr(i)%nr

                    ! Virtual receiver location = real source location
                    gmtr(i)%recr(j)%x = sxyz(j, 1)
                    gmtr(i)%recr(j)%y = sxyz(j, 2)
                    gmtr(i)%recr(j)%z = sxyz(j, 3)
                    gmtr(i)%recr(j)%t0 = gmtr_real(j)%srcr(1)%t0

                    ! Virtual receiver time - find reciprocal travel time
                    select case (which_medium)
                        case ('acoustic-iso', 'acoustic-tti')
                            do l = 1, gmtr_real(j)%nr
                                ! Check for position match using exact comparison
                                if (gmtr(i)%srcr(1)%x == gmtr_real(j)%recr(l)%x &
                                        .and. gmtr(i)%srcr(1)%y == gmtr_real(j)%recr(l)%y &
                                        .and. gmtr(i)%srcr(1)%z == gmtr_real(j)%recr(l)%z &
                                        .and. gmtr(i)%recr(j)%x == gmtr_real(j)%srcr(1)%x &
                                        .and. gmtr(i)%recr(j)%y == gmtr_real(j)%srcr(1)%y &
                                        .and. gmtr(i)%recr(j)%z == gmtr_real(j)%srcr(1)%z) then
                                    
                                    ttp(i)%array(j, 1) = ttp_real(j)%array(l, 1)
                                    gmtr(i)%recr(j)%weight = gmtr_real(j)%recr(l)%weight
                                    gmtr(i)%recr(j)%aoff = gmtr_real(j)%recr(l)%aoff
                                    
                                    ! Write mapping information
                                    if (rankid == 0) then
                                        write(mapping_unit, '(I0,1X,I0,1X,I0,1X,I0,1X,F12.6,1X,F12.6)') &
                                            i, j-1, gmtr_real(j)%id, l-1, ttp_real(j)%array(l, 1), 0.0
                                    end if
                                    
                                    matching_count = matching_count + 1
                                    cycle
                                end if
                            end do
                            
                        case ('elastic-iso', 'elastic-tti')
                            do l = 1, gmtr_real(j)%nr
                                if (gmtr(i)%srcr(1)%x == gmtr_real(j)%recr(l)%x &
                                        .and. gmtr(i)%srcr(1)%y == gmtr_real(j)%recr(l)%y &
                                        .and. gmtr(i)%srcr(1)%z == gmtr_real(j)%recr(l)%z &
                                        .and. gmtr(i)%recr(j)%x == gmtr_real(j)%srcr(1)%x &
                                        .and. gmtr(i)%recr(j)%y == gmtr_real(j)%srcr(1)%y &
                                        .and. gmtr(i)%recr(j)%z == gmtr_real(j)%srcr(1)%z) then
                                    
                                    ttp(i)%array(j, 1) = ttp_real(j)%array(l, 1)
                                    tts(i)%array(j, 1) = tts_real(j)%array(l, 1)
                                    gmtr(i)%recr(j)%weight = gmtr_real(j)%recr(l)%weight
                                    gmtr(i)%recr(j)%aoff = gmtr_real(j)%recr(l)%aoff
                                    
                                    ! Write mapping information
                                    if (rankid == 0) then
                                        write(mapping_unit, '(I0,1X,I0,1X,I0,1X,I0,1X,F12.6,1X,F12.6)') &
                                            i, j-1, gmtr_real(j)%id, l-1, ttp_real(j)%array(l, 1), tts_real(j)%array(l, 1)
                                    end if
                                    
                                    matching_count = matching_count + 1
                                    cycle
                                end if
                            end do
                    end select

                end do

                total_matches = total_matches + matching_count

                ! Save the exchanged data to processed observed data directory
                if (i >= shot_in_rank(rankid, 1) .and. i <= shot_in_rank(rankid, 2) .and. which_program /= 'eikonal') then
                    call make_directory(tidy(dir_working)//'/record_processed')
                    select case (which_medium)
                        case ('acoustic-iso', 'acoustic-tti')
                            call output_array(ttp(i)%array, tidy(dir_working)//'/record_processed/shot_'//num2str(i)//'_traveltime_p.bin')
                        case ('elastic-iso', 'elastic-tti')
                            call output_array(ttp(i)%array, tidy(dir_working)//'/record_processed/shot_'//num2str(i)//'_traveltime_p.bin')
                            call output_array(tts(i)%array, tidy(dir_working)//'/record_processed/shot_'//num2str(i)//'_traveltime_s.bin')
                    end select
                end if

            end do

            if (rankid == 0) then
                ! Write summary information
                write(mapping_unit, '(A)') '# Exchange Summary:'
                write(mapping_unit, '(A,I0)') '# Original shots: ', size(gmtr_real)
                write(mapping_unit, '(A,I0)') '# Virtual sources: ', nr
                write(mapping_unit, '(A,I0)') '# Virtual receivers per source: ', ns
                write(mapping_unit, '(A,I0)') '# Total reciprocal matches: ', total_matches
                write(mapping_unit, '(A,I0)') '# Total possible pairs: ', nr * ns
                write(mapping_unit, '(A,F8.3)') '# Match rate (%): ', real(total_matches) / real(nr * ns) * 100.0
                
                ! Close files
                close(mapping_unit)
                close(geom_unit)
                close(shot_mapping_unit)
                
                call warn(date_time_compact()//' Exchange mapping files created:')
                call warn('  - exchange_mapping.txt (detailed shot-receiver mappings)')
                call warn('  - exchange_geometry.txt (virtual source positions)')
                call warn('  - shot_mapping.txt (virtual receiver index to shot ID mapping)')
            end if

            call mpibarrier

            ! Change obs data directory to record_processed
            dir_record = tidy(dir_working)//'/record_processed'

            ! The # of virtual receivers is ns_real
            nr_virtual = ns

            ! The # of virtual sources is nr_real_unique
            ns = nr

            if (rankid == 0) then
                call warn(date_time_compact()//' Number of virtual sources = '//num2str(ns))
                call warn(date_time_compact()//' Number of virtual stations = '//num2str(nr_virtual))
            end if

        end if

        ! Divide effective sources into different MPI processess
        if (ns < nrank) then
            if (rankid == 0) then
                call warn(' <divide_shots> Error: # of shots '//num2str(ns)//' < # of MPI ranks '//num2str(nrank))
            end if
            call mpibarrier
            stop
        end if
        call alloc_array(shot_in_rank, [0, nrank - 1, 1, 2])
        call cut(1, ns, nrank, shot_in_rank)

    end subroutine divide_shots

    !
    !> Compute residual traveltime using absolute difference: $$\Delta T_i = T_{\text{syn}, i} - T_{\text{obs}, i}$$
    !> or double difference:
    !> $$\Delta T_i = \sum_{k = 1}^{N_r} (T_{\text{syn}, i} - T_{\text{syn}, k}) - (T_{\text{obs}, i} - T_{\text{obs}, k})$$
    !
    subroutine traveltime_residual(tsyn, tobs, tresidual, tmisfit)

        !> Synthetic traveltime, where each column represents a set of traveltime.
        !> For instance, fisrt-arrival, reflector 1, reflector 2, ...
        real, dimension(:, :), intent(in) :: tsyn
        !> Observed traveltime
        real, dimension(:, :), intent(inout) :: tobs
        !> Residual traveltime computed using either absolute difference or double difference
        real, allocatable, dimension(:, :), intent(inout) :: tresidual
        !> Misfit
        real, allocatable, dimension(:, :), intent(inout) :: tmisfit

        integer :: nr, nd, i, j, k
        real :: d

        nr = size(tsyn, 1)
        nd = size(tsyn, 2)
        tresidual = zeros_like(tsyn)
        tmisfit = zeros_like(tsyn)

        select case (misfit_type)

            case ('absolute-difference', 'ad')

                !$omp parallel do private(i, j, d)
                do i = 1, nd
                    ! Only compute for nonzero data
                    if (maxval(tsyn(:, i)) > 0 .and. maxval(tobs(:, i)) > 0) then
                        do j = 1, nr
                            ! If a value of field record is negative, then ignore
                            if (tobs(j, i) > 0) then
                                d = tsyn(j, i) - tobs(j, i)
                                tresidual(j, i) = d
                                tmisfit(j, i) = d**2
                            end if
                        end do
                    end if
                end do
                !$omp end parallel do

            case ('double-difference', 'dd')

                !$omp parallel do private(i, j, k, d)
                do i = 1, nd
                    ! Only compute for nonzero data
                    if (maxval(tsyn(:, i)) > 0 .and. maxval(tobs(:, i)) > 0) then
                        do j = 1, nr
                            ! If a value of field record is negative, then ignore
                            do k = 1, nr
                                if (tobs(j, i) >= 0 .and. tobs(k, i) >= 0) then
                                    d = (tsyn(j, i) - tsyn(k, i)) - (tobs(j, i) - tobs(k, i))
                                    tresidual(j, i) = tresidual(j, i) + d
                                    tmisfit(j, i) = tmisfit(j, i) + d**2
                                end if
                            end do
                        end do
                    end if
                end do
                !$omp end parallel do

        end select

    end subroutine traveltime_residual

    !
    !> Compute travletime misfit for a gather
    !
    subroutine compute_shot_misfit(ishot, dir_field, label, ttp_syn, ttp_obs, ttp_residual, m, tsyn_all, tobs_all)
        integer, intent(in) :: ishot
        character(len=*), intent(in) :: dir_field, label
        real, dimension(:, :), intent(in) :: ttp_syn
        real, allocatable, dimension(:, :), intent(inout) :: ttp_obs
        real, allocatable, dimension(:, :), intent(out) :: ttp_residual
        real, intent(inout) :: m
        real, dimension(:, :), intent(in), optional :: tsyn_all, tobs_all

        real, allocatable, dimension(:, :) :: weight, ttp_misfit
        integer :: i, j
        real :: d

        if (yn_dd_no_st0) then
            call assert(present(tsyn_all) .and. present(tobs_all), ' <compute_shot_misfit> Error: tloc-dd requires all data ')
            ttp_residual = zeros(gmtr(ishot)%nr, 1)
            ttp_misfit = zeros_like(ttp_residual)
            !$omp parallel do private(i, j, d)
            do j = 1, gmtr(ishot)%nr
                if (maxval(tsyn_all(j, :)) > 0 .and. maxval(tobs_all(j, :)) > 0) then
                    do i = 1, ns
                        if (tobs_all(j, ishot) >= 0 .and. tobs_all(j, i) >= 0) then
                            d = (tsyn_all(j, ishot) - tsyn_all(j, i)) - (tobs_all(j, ishot) - tobs_all(j, i))
                            ttp_residual(j, 1) = ttp_residual(j, 1) + d
                            ttp_misfit(j, 1) = ttp_misfit(j, 1) + d**2
                        end if
                    end do
                end if
            end do
            !$omp end parallel do
        else
            ttp_obs = load(tidy(dir_field)//'/shot_'//num2str(gmtr(ishot)%id)//'_traveltime_'//tidy(label)//'.bin', gmtr(ishot)%nr, nrefl + 1)
            call traveltime_residual(ttp_syn, ttp_obs, ttp_residual, ttp_misfit)
        end if

        weight = ones_like(ttp_residual)
        !$omp parallel do private(i, j)
        do i = 1, gmtr(ishot)%nr
            do j = 1, size(ttp_residual, 2)
                weight(i, j) = gmtr(ishot)%recr(i)%weight*misfit_weight(j)
                if (abs(ttp_residual(i, j)) > misfit_threshold) weight(i, j) = 0
                if (j > 1 .and. (gmtr(ishot)%recr(i)%aoff < offset_min_refl &
                    .or. gmtr(ishot)%recr(i)%aoff > offset_max_refl)) weight(i, j) = 0
            end do
        end do
        !$omp end parallel do

        m = m + sum(ttp_misfit*weight)
        ttp_residual = ttp_residual*weight

        do j = 2, size(ttp_residual, 2)
            where (weight(:, j) == 0)
                ttp_obs(:, j) = ttp_obs(:, j)*weight(:, j)
            end where
        end do
    end subroutine compute_shot_misfit

    !
    !> Copy entire directory contents from source to destination
    !
    subroutine copy_directory(source_dir, dest_dir)

        character(len=*), intent(in) :: source_dir, dest_dir
        character(len=2048) :: cmd
        integer :: stat

        ! Create destination directory if it doesn't exist
        call make_directory(dest_dir)

        ! Use cp -r to copy all files recursively
        ! The trailing /. ensures we copy the contents, not the directory itself
        cmd = 'cp -r '//tidy(source_dir)//'/. '//tidy(dest_dir)//'/'

        ! Execute the copy command
        call execute_command_line(trim(cmd), exitstat=stat)

        if (stat /= 0) then
            call warn(date_time_compact()//' Warning: copy_directory from ' &
                //tidy(source_dir)//' to '//tidy(dest_dir)//' may have failed (exit status: ' &
                //num2str(stat)//')')
        end if

    end subroutine copy_directory

end module utility