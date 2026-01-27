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
! license in this material to reproduce, prepare derivative works,
! distribute copies to the public, perform publicly and display publicly,
! and to permit others to do so.
!
! Author:
!    Kai Gao, kaigao@lanl.gov
!


program main

    use libflit
    use parameters
    use utility
    use vars
    use traveltime_iso
    use traveltime_iso_reflection
    use raypath

    implicit none

#ifdef dim2
    real, allocatable, dimension(:, :) :: ttp, tts, vp, vs, refl
    real, allocatable, dimension(:, :, :) :: ttp_all, tts_all
    ! Gradient field storage
    real, allocatable, dimension(:, :) :: pdx_p, pdz_p
    real, allocatable, dimension(:, :) :: pdx_s, pdz_s
#endif
#ifdef dim3
    real, allocatable, dimension(:, :, :) :: ttp, tts, vp, vs, refl
    real, allocatable, dimension(:, :, :, :) :: ttp_all, tts_all
    ! Gradient field storage
    real, allocatable, dimension(:, :, :) :: pdx_p, pdy_p, pdz_p
    real, allocatable, dimension(:, :, :) :: pdx_s, pdy_s, pdz_s
#endif
    real, allocatable, dimension(:, :) :: ttprecr, ttsrecr
    integer :: i, iz  ! TEMPORARY DEBUG: iz for depth loop

    ! Ray path storage variables
    integer, allocatable, dimension(:) :: raypath_npoints_p, raypath_npoints_s
    real, allocatable, dimension(:, :) :: raypath_coords_p, raypath_coords_s
    integer, allocatable, dimension(:) :: raypath_status_p, raypath_status_s
    logical :: need_gradients_for_raypath  ! Flag to track if gradients needed only for ray paths
    integer :: n_success_rays, n_total_rays

    which_program = 'eikonal'

    call mpistart

    if (command_argument_count() == 0) then
        if (rankid == 0) then
            call warn('')
            call warn(date_time_compact()//' Error: Parameter file not found. Exiting. ')
            call warn('')
        end if
        call mpibarrier
        call mpistop
    end if

    if (rankid == 0) then
        call warn('')
        call warn(tile('=', 80))
        call warn(center_substring('Eikonal Traveltime Computation Begins', 80))
        call warn('')
        call print_date_time
        call warn('')
    end if

    ! Read parameters
    call read_parameters

    ! Initialize ray path parameters
    call initialize_raypath_params

    ! Set dimensions
    call set_regular_space

    ! Load geometry
    call load_geometry

    ! Divide shots to MPI ranks
    call divide_shots

    ! Read models
    call prepare_model

    ! Get model
    do i = 1, nmodel
        select case (model_m(i)%name)
            case ('vp')
                vp = model_m(i)%array
            case ('vs')
                vs = model_m(i)%array
            case ('refl')
                refl = model_m(i)%array
        end select
    end do

    ! TEMPORARY DEBUG: Print 1D velocity profile at x=1, y=1 (or x=1 for 2D)
    if (rankid == 0) then
#ifdef dim3
        write(*, '(a)') '--- TEMPORARY DEBUG: VP and VS distribution at x=1, y=1 ---'
        if (allocated(vp)) then
            if (allocated(vs)) then
                write(*, '(a)') 'Depth(m)         VP(m/s)         VS(m/s)'
                do iz = 1, nz
                    write(*, '(f10.2, f16.2, f16.2)') oz + (iz-1)*dz, vp(iz, 1, 1), vs(iz, 1, 1)
                end do
            else
                write(*, '(a)') 'Depth(m)         VP(m/s)'
                do iz = 1, nz
                    write(*, '(f10.2, f16.2)') oz + (iz-1)*dz, vp(iz, 1, 1)
                end do
            end if
        end if
#endif
#ifdef dim2
        write(*, '(a)') '--- TEMPORARY DEBUG: VP and VS distribution at x=1 ---'
        if (allocated(vp)) then
            if (allocated(vs)) then
                write(*, '(a)') 'Depth(m)         VP(m/s)         VS(m/s)'
                do iz = 1, nz
                    write(*, '(f10.2, f16.2, f16.2)') oz + (iz-1)*dz, vp(iz, 1), vs(iz, 1)
                end do
            else
                write(*, '(a)') 'Depth(m)         VP(m/s)'
                do iz = 1, nz
                    write(*, '(f10.2, f16.2)') oz + (iz-1)*dz, vp(iz, 1)
                end do
            end if
        end if
#endif
        write(*, '(a)') '--- END TEMPORARY DEBUG ---'
    end if
    ! END TEMPORARY DEBUG

    ! Make directory
    call make_directory(dir_synthetic)

    call mpibarrier

    ! Show the medium parameter statistics
    if (rankid == 0) then
        do i = 1, nmodel
            call plot_histogram(model_m(i)%array, &
                label=date_time_compact()//' '//tidy(model_name(i))//' distribution ')
        end do
    end if
    call mpibarrier

    if (sum(snaps) > 0 .and. rankid == 0) then
        call make_directory(dir_snapshot)
    end if

    if (yn_save_traveltime_gradients .and. rankid == 0) then
        call make_directory(dir_gradient)
    end if

    ! Create ray path directory if needed
    if (raypath_save .and. rankid == 0) then
        call make_directory(raypath_dir)
        call warn(date_time_compact()//' Ray path output enabled. Directory: '//tidy(raypath_dir))
        if (raypath_adaptive_step) call warn(date_time_compact()//'   - Adaptive step size: enabled')
        if (raypath_hybrid_mode) call warn(date_time_compact()//'   - Hybrid mode: enabled')
        if (raypath_momentum_mode) call warn(date_time_compact()//'   - Momentum mode: enabled')
    end if

    call mpibarrier

    ! Traveltime computation source by source
    do ishot = shot_in_rank(rankid, 1), shot_in_rank(rankid, 2)

        ! Determine ranges of model selected for this source, if necessary
        call set_adaptive_model_range(gmtr(ishot))

        ! Solving eikional equation for traveltime
        select case (which_medium)

            case ('acoustic-iso')

                ! Determine if gradients are needed (for saving or ray path tracing)
                need_gradients_for_raypath = raypath_save .and. .not. yn_save_traveltime_gradients

#ifdef dim2
                if (allocated(refl)) then

                    call forward_iso_reflection( &
                        vp(shot_nzbeg:shot_nzend, shot_nxbeg:shot_nxend), &
                        [dx, dz], [shot_xbeg, shot_zbeg], gmtr(ishot), &
                        refl(shot_nzbeg:shot_nzend, shot_nxbeg:shot_nxend), ttp_all, ttprecr)

                else

                    ! Compute gradients if saving them OR if needed for ray paths
                    if (yn_save_traveltime_gradients .or. raypath_save) then
                        call forward_iso( &
                            vp(shot_nzbeg:shot_nzend, shot_nxbeg:shot_nxend), &
                            [dx, dz], [shot_xbeg, shot_zbeg], gmtr(ishot), &
                            ttp, ttprecr, pdx_p, pdz_p)
                    else
                        call forward_iso( &
                            vp(shot_nzbeg:shot_nzend, shot_nxbeg:shot_nxend), &
                            [dx, dz], [shot_xbeg, shot_zbeg], gmtr(ishot), ttp, ttprecr)
                    end if

                end if

#endif
#ifdef dim3
                if (allocated(refl)) then

                    call forward_iso_reflection( &
                        vp(shot_nzbeg:shot_nzend, shot_nybeg:shot_nyend, shot_nxbeg:shot_nxend), &
                        [dx, dy, dz], [shot_xbeg, shot_ybeg, shot_zbeg], gmtr(ishot), &
                        refl(shot_nzbeg:shot_nzend, shot_nybeg:shot_nyend, shot_nxbeg:shot_nxend), ttp_all, ttprecr)

                else

                    ! Compute gradients if saving them OR if needed for ray paths
                    if (yn_save_traveltime_gradients .or. raypath_save) then
                        call forward_iso( &
                            vp(shot_nzbeg:shot_nzend, shot_nybeg:shot_nyend, shot_nxbeg:shot_nxend), &
                            [dx, dy, dz], [shot_xbeg, shot_ybeg, shot_zbeg], gmtr(ishot), &
                            ttp, ttprecr, pdx_p, pdy_p, pdz_p)
                    else
                        call forward_iso( &
                            vp(shot_nzbeg:shot_nzend, shot_nybeg:shot_nyend, shot_nxbeg:shot_nxend), &
                            [dx, dy, dz], [shot_xbeg, shot_ybeg, shot_zbeg], gmtr(ishot), ttp, ttprecr)
                    end if

                end if
#endif
                call output_array(ttprecr, tidy(dir_synthetic)// &
                    '/shot_'//num2str(gmtr(ishot)%id)//'_traveltime_p.bin')
                if (sum(snaps) > 0) then
                    if (allocated(refl)) then
                        call output_array(ttp_all, tidy(dir_snapshot)// &
                            '/shot_'//num2str(gmtr(ishot)%id)//'_traveltime_p.bin')
                    else
                        call output_array(ttp, tidy(dir_snapshot)// &
                            '/shot_'//num2str(gmtr(ishot)%id)//'_traveltime_p.bin')
                    end if
                end if

                ! Trace and save ray paths for P-wave (acoustic)
                if (raypath_save .and. .not. allocated(refl)) then
#ifdef dim2
                    call trace_shot_ray_paths(pdx_p, pdz_p, [dx, dz], &
                        [shot_xbeg, shot_zbeg], gmtr(ishot), &
                        raypath_npoints_p, raypath_coords_p, raypath_status_p)
                    call write_ray_paths_binary(tidy(raypath_dir)// &
                        '/shot_'//num2str(gmtr(ishot)%id)//'_raypath_p.bin', &
                        raypath_npoints_p, raypath_coords_p)
#endif
#ifdef dim3
                    call trace_shot_ray_paths(pdx_p, pdy_p, pdz_p, &
                        [dx, dy, dz], [shot_xbeg, shot_ybeg, shot_zbeg], gmtr(ishot), &
                        raypath_npoints_p, raypath_coords_p, raypath_status_p)
                    call write_ray_paths_binary(tidy(raypath_dir)// &
                        '/shot_'//num2str(gmtr(ishot)%id)//'_raypath_p.bin', &
                        raypath_npoints_p, raypath_coords_p)
#endif
                    ! Report ray path statistics
                    n_total_rays = size(raypath_status_p)
                    n_success_rays = count(raypath_status_p == RAYSTATUS_SUCCESS)
                    if (verbose) then
                        call warn(date_time_compact()//' Shot '//num2str(gmtr(ishot)%id)// &
                            ' P-wave rays: '//num2str(n_success_rays)//'/'//num2str(n_total_rays)//' converged')
                    end if
                    ! Deallocate ray path arrays
                    if (allocated(raypath_npoints_p)) deallocate(raypath_npoints_p)
                    if (allocated(raypath_coords_p)) deallocate(raypath_coords_p)
                    if (allocated(raypath_status_p)) deallocate(raypath_status_p)
                end if

                if (yn_save_traveltime_gradients) then
                    call output_array(pdx_p, tidy(dir_gradient)// &
                        '/shot_'//num2str(gmtr(ishot)%id)//'_gradient_px.bin')
#ifdef dim3
                    call output_array(pdy_p, tidy(dir_gradient)// &
                        '/shot_'//num2str(gmtr(ishot)%id)//'_gradient_py.bin')
#endif
                    call output_array(pdz_p, tidy(dir_gradient)// &
                        '/shot_'//num2str(gmtr(ishot)%id)//'_gradient_pz.bin')
                    ! Write metadata for gradient field
#ifdef dim3
                    call write_gradient_metadata(gmtr(ishot)%id, dir_gradient, &
                        shot_xbeg, shot_ybeg, shot_zbeg, &
                        shot_nx, shot_ny, shot_nz, dx, dy, dz)
#else
                    call write_gradient_metadata(gmtr(ishot)%id, dir_gradient, &
                        shot_xbeg, 0.0, shot_zbeg, &
                        shot_nx, 1, shot_nz, dx, 1.0, dz)
#endif
                end if

                ! Deallocate gradients if they were only computed for ray paths
                if (need_gradients_for_raypath .and. .not. allocated(refl)) then
                    if (allocated(pdx_p)) deallocate(pdx_p)
                    if (allocated(pdz_p)) deallocate(pdz_p)
#ifdef dim3
                    if (allocated(pdy_p)) deallocate(pdy_p)
#endif
                end if

            case ('elastic-iso')

#ifdef dim2
                if (allocated(refl)) then

                    select case (incident_wave)
                        case ('p')
                            call forward_iso_reflection_elastic( &
                                vp(shot_nzbeg:shot_nzend, shot_nxbeg:shot_nxend), &
                                vs(shot_nzbeg:shot_nzend, shot_nxbeg:shot_nxend), &
                                [dx, dz], [shot_xbeg, shot_zbeg], gmtr(ishot), &
                                refl(shot_nzbeg:shot_nzend, shot_nxbeg:shot_nxend), &
                                ttp_all, tts_all, ttprecr, ttsrecr)
                        case ('s')
                            call forward_iso_reflection_elastic( &
                                vs(shot_nzbeg:shot_nzend, shot_nxbeg:shot_nxend), &
                                vp(shot_nzbeg:shot_nzend, shot_nxbeg:shot_nxend), &
                                [dx, dz], [shot_xbeg, shot_zbeg], gmtr(ishot), &
                                refl(shot_nzbeg:shot_nzend, shot_nxbeg:shot_nxend), &
                                tts_all, ttp_all, ttsrecr, ttprecr)
                    end select

                else

                    ! Compute gradients if saving them OR if needed for ray paths
                    if (yn_save_traveltime_gradients .or. raypath_save) then
                        call forward_iso( &
                            vp(shot_nzbeg:shot_nzend, shot_nxbeg:shot_nxend), &
                            [dx, dz], [shot_xbeg, shot_zbeg], gmtr(ishot), &
                            ttp, ttprecr, pdx_p, pdz_p)
                        call forward_iso( &
                            vs(shot_nzbeg:shot_nzend, shot_nxbeg:shot_nxend), &
                            [dx, dz], [shot_xbeg, shot_zbeg], gmtr(ishot), &
                            tts, ttsrecr, pdx_s, pdz_s)
                    else
                        call forward_iso( &
                            vp(shot_nzbeg:shot_nzend, shot_nxbeg:shot_nxend), &
                            [dx, dz], [shot_xbeg, shot_zbeg], gmtr(ishot), ttp, ttprecr)
                        call forward_iso( &
                            vs(shot_nzbeg:shot_nzend, shot_nxbeg:shot_nxend), &
                            [dx, dz], [shot_xbeg, shot_zbeg], gmtr(ishot), tts, ttsrecr)
                    end if

                end if

#endif
#ifdef dim3
                if (allocated(refl)) then

                    select case (incident_wave)
                        case ('p')
                            call forward_iso_reflection_elastic( &
                                vp(shot_nzbeg:shot_nzend, shot_nybeg:shot_nyend, shot_nxbeg:shot_nxend), &
                                vs(shot_nzbeg:shot_nzend, shot_nybeg:shot_nyend, shot_nxbeg:shot_nxend), &
                                [dx, dy, dz], [shot_xbeg, shot_ybeg, shot_zbeg], gmtr(ishot), &
                                refl(shot_nzbeg:shot_nzend, shot_nybeg:shot_nyend, shot_nxbeg:shot_nxend), &
                                ttp_all, tts_all, ttprecr, ttsrecr)
                        case ('s')
                            call forward_iso_reflection_elastic( &
                                vs(shot_nzbeg:shot_nzend, shot_nybeg:shot_nyend, shot_nxbeg:shot_nxend), &
                                vp(shot_nzbeg:shot_nzend, shot_nybeg:shot_nyend, shot_nxbeg:shot_nxend), &
                                [dx, dy, dz], [shot_xbeg, shot_ybeg, shot_zbeg], gmtr(ishot), &
                                refl(shot_nzbeg:shot_nzend, shot_nybeg:shot_nyend, shot_nxbeg:shot_nxend), &
                                tts_all, ttp_all, ttsrecr, ttprecr)
                    end select

                else

                    ! Compute gradients if saving them OR if needed for ray paths
                    if (yn_save_traveltime_gradients .or. raypath_save) then
                        call forward_iso( &
                            vp(shot_nzbeg:shot_nzend, shot_nybeg:shot_nyend, shot_nxbeg:shot_nxend), &
                            [dx, dy, dz], [shot_xbeg, shot_ybeg, shot_zbeg], gmtr(ishot), &
                            ttp, ttprecr, pdx_p, pdy_p, pdz_p)
                        call forward_iso( &
                            vs(shot_nzbeg:shot_nzend, shot_nybeg:shot_nyend, shot_nxbeg:shot_nxend), &
                            [dx, dy, dz], [shot_xbeg, shot_ybeg, shot_zbeg], gmtr(ishot), &
                            tts, ttsrecr, pdx_s, pdy_s, pdz_s)
                    else
                        call forward_iso( &
                            vp(shot_nzbeg:shot_nzend, shot_nybeg:shot_nyend, shot_nxbeg:shot_nxend), &
                            [dx, dy, dz], [shot_xbeg, shot_ybeg, shot_zbeg], gmtr(ishot), ttp, ttprecr)
                        call forward_iso( &
                            vs(shot_nzbeg:shot_nzend, shot_nybeg:shot_nyend, shot_nxbeg:shot_nxend), &
                            [dx, dy, dz], [shot_xbeg, shot_ybeg, shot_zbeg], gmtr(ishot), tts, ttsrecr)
                    end if

                end if
#endif
                call output_array(ttprecr, tidy(dir_synthetic)// &
                    '/shot_'//num2str(gmtr(ishot)%id)//'_traveltime_p.bin')
                call output_array(ttsrecr, tidy(dir_synthetic)// &
                    '/shot_'//num2str(gmtr(ishot)%id)//'_traveltime_s.bin')
                if (sum(snaps) > 0) then
                    if (allocated(refl)) then
                        call output_array(ttp_all, tidy(dir_snapshot)// &
                            '/shot_'//num2str(gmtr(ishot)%id)//'_traveltime_p.bin')
                        call output_array(tts_all, tidy(dir_snapshot)// &
                            '/shot_'//num2str(gmtr(ishot)%id)//'_traveltime_s.bin')
                    else
                        call output_array(ttp, tidy(dir_snapshot)// &
                            '/shot_'//num2str(gmtr(ishot)%id)//'_traveltime_p.bin')
                        call output_array(tts, tidy(dir_snapshot)// &
                            '/shot_'//num2str(gmtr(ishot)%id)//'_traveltime_s.bin')
                    end if
                end if

                ! Trace and save ray paths for P-wave and S-wave (elastic)
                if (raypath_save .and. .not. allocated(refl)) then
#ifdef dim2
                    ! P-wave ray paths
                    call trace_shot_ray_paths(pdx_p, pdz_p, [dx, dz], &
                        [shot_xbeg, shot_zbeg], gmtr(ishot), &
                        raypath_npoints_p, raypath_coords_p, raypath_status_p)
                    call write_ray_paths_binary(tidy(raypath_dir)// &
                        '/shot_'//num2str(gmtr(ishot)%id)//'_raypath_p.bin', &
                        raypath_npoints_p, raypath_coords_p)
                    ! S-wave ray paths
                    call trace_shot_ray_paths(pdx_s, pdz_s, [dx, dz], &
                        [shot_xbeg, shot_zbeg], gmtr(ishot), &
                        raypath_npoints_s, raypath_coords_s, raypath_status_s)
                    call write_ray_paths_binary(tidy(raypath_dir)// &
                        '/shot_'//num2str(gmtr(ishot)%id)//'_raypath_s.bin', &
                        raypath_npoints_s, raypath_coords_s)
#endif
#ifdef dim3
                    ! P-wave ray paths
                    call trace_shot_ray_paths(pdx_p, pdy_p, pdz_p, &
                        [dx, dy, dz], [shot_xbeg, shot_ybeg, shot_zbeg], gmtr(ishot), &
                        raypath_npoints_p, raypath_coords_p, raypath_status_p)
                    call write_ray_paths_binary(tidy(raypath_dir)// &
                        '/shot_'//num2str(gmtr(ishot)%id)//'_raypath_p.bin', &
                        raypath_npoints_p, raypath_coords_p)
                    ! S-wave ray paths
                    call trace_shot_ray_paths(pdx_s, pdy_s, pdz_s, &
                        [dx, dy, dz], [shot_xbeg, shot_ybeg, shot_zbeg], gmtr(ishot), &
                        raypath_npoints_s, raypath_coords_s, raypath_status_s)
                    call write_ray_paths_binary(tidy(raypath_dir)// &
                        '/shot_'//num2str(gmtr(ishot)%id)//'_raypath_s.bin', &
                        raypath_npoints_s, raypath_coords_s)
#endif
                    ! Report ray path statistics
                    if (verbose) then
                        n_total_rays = size(raypath_status_p)
                        n_success_rays = count(raypath_status_p == RAYSTATUS_SUCCESS)
                        call warn(date_time_compact()//' Shot '//num2str(gmtr(ishot)%id)// &
                            ' P-wave rays: '//num2str(n_success_rays)//'/'//num2str(n_total_rays)//' converged')
                        n_total_rays = size(raypath_status_s)
                        n_success_rays = count(raypath_status_s == RAYSTATUS_SUCCESS)
                        call warn(date_time_compact()//' Shot '//num2str(gmtr(ishot)%id)// &
                            ' S-wave rays: '//num2str(n_success_rays)//'/'//num2str(n_total_rays)//' converged')
                    end if
                    ! Deallocate ray path arrays
                    if (allocated(raypath_npoints_p)) deallocate(raypath_npoints_p)
                    if (allocated(raypath_coords_p)) deallocate(raypath_coords_p)
                    if (allocated(raypath_status_p)) deallocate(raypath_status_p)
                    if (allocated(raypath_npoints_s)) deallocate(raypath_npoints_s)
                    if (allocated(raypath_coords_s)) deallocate(raypath_coords_s)
                    if (allocated(raypath_status_s)) deallocate(raypath_status_s)
                end if

                if (yn_save_traveltime_gradients) then
                    call output_array(pdx_p, tidy(dir_gradient)// &
                        '/shot_'//num2str(gmtr(ishot)%id)//'_gradient_px.bin')
#ifdef dim3
                    call output_array(pdy_p, tidy(dir_gradient)// &
                        '/shot_'//num2str(gmtr(ishot)%id)//'_gradient_py.bin')
#endif
                    call output_array(pdz_p, tidy(dir_gradient)// &
                        '/shot_'//num2str(gmtr(ishot)%id)//'_gradient_pz.bin')
                    call output_array(pdx_s, tidy(dir_gradient)// &
                        '/shot_'//num2str(gmtr(ishot)%id)//'_gradient_sx.bin')
#ifdef dim3
                    call output_array(pdy_s, tidy(dir_gradient)// &
                        '/shot_'//num2str(gmtr(ishot)%id)//'_gradient_sy.bin')
#endif
                    call output_array(pdz_s, tidy(dir_gradient)// &
                        '/shot_'//num2str(gmtr(ishot)%id)//'_gradient_sz.bin')
                    ! Write metadata for gradient field
#ifdef dim3
                    call write_gradient_metadata(gmtr(ishot)%id, dir_gradient, &
                        shot_xbeg, shot_ybeg, shot_zbeg, &
                        shot_nx, shot_ny, shot_nz, dx, dy, dz)
#else
                    call write_gradient_metadata(gmtr(ishot)%id, dir_gradient, &
                        shot_xbeg, 0.0, shot_zbeg, &
                        shot_nx, 1, shot_nz, dx, 1.0, dz)
#endif
                end if

                ! Deallocate gradients if they were only computed for ray paths
                if (need_gradients_for_raypath .and. .not. allocated(refl)) then
                    if (allocated(pdx_p)) deallocate(pdx_p)
                    if (allocated(pdz_p)) deallocate(pdz_p)
                    if (allocated(pdx_s)) deallocate(pdx_s)
                    if (allocated(pdz_s)) deallocate(pdz_s)
#ifdef dim3
                    if (allocated(pdy_p)) deallocate(pdy_p)
                    if (allocated(pdy_s)) deallocate(pdy_s)
#endif
                end if

        end select

        call warn(date_time_compact()//' Shot '//num2str(gmtr(ishot)%id)//' traveltime computation completed. ')

    end do

    call mpibarrier

    if (rankid == 0) then
        call warn('')
        call print_date_time
        call warn('')
        call warn(center_substring('EIKONAL Completed', 80))
        call warn(tile('=', 80))
        call warn('')
    end if

    call mpiend

end program main
