!
! © 2024-2026. Triad National Security, LLC. All rights reserved.
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


module parameters

    use libflit
    use spherical

    implicit none

    character(len=1024) :: file_parameter

    character(len=32) :: which_program

    character(len=1024) :: dir_record
    character(len=1024) :: dir_working
    character(len=1024) :: dir_scratch
    character(len=1024) :: dir_synthetic
    character(len=1024) :: dir_snapshot
    character(len=1024) :: dir_adjoint

    character(len=1204) :: file_geometry

    integer :: ns = 1

    ! meditum type
    character(len=24) :: which_medium

    ! snapshot times
    real, allocatable, dimension(:) :: snaps


    ! shot index and division related
    integer, allocatable, dimension(:, :) :: shot_in_rank
    integer :: ishot
    integer, allocatable, dimension(:) :: shot_index, rec_index
    integer, allocatable, dimension(:) :: src_exclude, sid_exclude, rec_exclude

    real :: sxmin = -float_huge
    real :: sxmax = +float_huge
    real :: symin = -float_huge
    real :: symax = +float_huge
    real :: szmin = -float_huge
    real :: szmax = +float_huge
    real :: rxmin = -float_huge
    real :: rxmax = +float_huge
    real :: rymin = -float_huge
    real :: rymax = +float_huge
    real :: rzmin = -float_huge
    real :: rzmax = +float_huge

    ! only compute misfit
    logical :: yn_misfit_only
    real :: misfit
    real :: val_misfit

    ! search direction computation method in inversion
    character(len=32), allocatable, dimension(:) :: model_search_method
    character(len=32) :: search_method

    ! regularization for inversion
    character(len=32), allocatable, dimension(:) :: model_regularization_method
    logical :: yn_regularize_model = .false.

    ! inversion parameters
    integer :: iter
    integer :: niter_max
    character(len=32) :: misfit_type

    ! vp/vs ratio smoothing
    real :: vpvsratio_smoothx = 0
    real :: vpvsratio_smoothy = 0
    real :: vpvsratio_smoothz = 0

    ! method of computing step size: linear or line search
    character(len=32) :: step_size_method

    ! step
    real :: step_scaling_factor
    real :: step_max_scale_factor = 1.0

    ! data misfit arrays
    real, allocatable, dimension(:) :: step_misfit
    real, allocatable, dimension(:) :: data_misfit
    real, allocatable, dimension(:, :) :: shot_misfit
    real, allocatable, dimension(:) :: misfit_weight
    real :: jumpout_factor = 1.0

    ! maximum and minimum offset
    real :: offset_min = 0.0d0
    real :: offset_max = float_huge

    ! source id
    integer :: sid_min = 1
    integer :: sid_max = int4_huge

    ! maximum and minimum shot number
    integer :: shot_min = 1
    integer :: shot_max = 1
    integer :: shot_every = 1

    ! maximum and minimum receiver number
    integer :: rec_min = 1
    integer :: rec_max = int4_huge
    integer :: rec_every = 1

    ! necessity of adaptivity
    logical :: yn_adpx
    logical :: yn_adpy
    logical :: yn_adpz

    real :: adpextrax = 0.0
    real :: adpextray = 0.0
    real :: adpextraz = 0.0

    ! taper lengths at the sides of an adaptive range that lie inside the model
    real :: adptaperx = 0.0
    real :: adptapery = 0.0
    real :: adptaperz = 0.0

    integer, allocatable, dimension(:) :: sid_select, src_select

    character(len=24) :: aniso_param = 'thomsen'

    real :: min_vpvsratio, max_vpvsratio

    logical :: verbose
    integer :: sweep_niter_max = 10
    real :: sweep_stop_threshold = 1.0e-4
    character(len=12) :: adjoint_receiver = 'linear'
    real :: misfit_threshold


    character(len=1024) :: file_andfaux, file_andfcoh
    integer :: resume_from_iter

    character(len=1024) :: file_datamisfit, file_shotmisfit
    logical :: yn_continue_inv = .false.
    logical :: yn_enforce_update = .false.
    logical :: trigger_jumpout = .false.
    logical :: put_synthetic_in_scratch = .false.

    character(len=32), allocatable, dimension(:) :: process_shot_grad, process_grad

    character(len=24), allocatable, dimension(:) :: model_name, model_name_aux
    integer :: nmodel, nmodel_aux

    real, allocatable, dimension(:) :: model_step, model_step_max
    real, allocatable, dimension(:) :: model_min, model_max

    character(len=24), allocatable, dimension(:) :: data_name
    integer :: ndata

    logical :: require_model_interp = .false.
    logical :: uniform_processing = .true.

    logical :: yn_flat_stop = .false.

    character(len=32) :: data_prefix = 'traveltime'

    logical :: yn_precond = .true.
    logical :: yn_output_dws = .false.

    real :: precond_eps = 1.0e-4

    ! Dimension
    integer :: nx0, ny0, nz0
    real :: dx0, dy0, dz0
    real :: ox0, oy0, oz0

    integer :: nx, ny, nz
    real :: dx, dy, dz
    real :: ox, oy, oz

    real :: xmin = -float_large
    real :: ymin = -float_large
    real :: zmin = -float_large
    real :: xmax = +float_large
    real :: ymax = +float_large
    real :: zmax = +float_large

    integer :: nxbeg1, nxbeg2, nxend1, nxend2, nxbeg_e, nxend_e
    integer :: nybeg1, nybeg2, nyend1, nyend2, nybeg_e, nyend_e
    integer :: nzbeg1, nzbeg2, nzend1, nzend2, nzbeg_e, nzend_e

    real :: shot_xbeg, shot_ybeg, shot_zbeg
    real :: shot_xend, shot_yend, shot_zend
    integer :: shot_nxbeg, shot_nybeg, shot_nzbeg
    integer :: shot_nxend, shot_nyend, shot_nzend
    integer :: shot_nx, shot_ny, shot_nz

    integer :: regular_nx, regular_ny, regular_nz
    real :: regular_dx, regular_dy, regular_dz
    real :: regular_ox, regular_oy, regular_oz

    type source
        real :: x, y, z
        real :: t0 = 0.0
        real :: amp = 1.0
    end type source

    type receiver
        real :: x, y, z
        real :: aoff
        real :: weight = 1.0    ! master enable: non-zero if receiver is active for any wave
        real :: weight_p = 1.0  ! P-wave scalar weight (equals weight in 4-column format)
        real :: weight_s = 1.0  ! S-wave scalar weight (equals weight in 4-column format)
        real :: t0 = 0.0
    end type receiver

    type source_receiver_geometry
        integer :: id
        integer :: ns
        type(source), allocatable, dimension(:) :: srcr
        integer :: nr
        type(receiver), allocatable, dimension(:) :: recr
    end type source_receiver_geometry

    type(source_receiver_geometry), allocatable, dimension(:) :: gmtr

    ! For TLOC
    type(source_receiver_geometry), allocatable, dimension(:) :: gmtr_real
    logical :: yn_exchange_sr
    integer :: nr_virtual
    logical :: yn_dd_no_st0

    ! True when geometry files use the 5-column format (x y z weight_p weight_s)
    logical :: yn_separate_ps_weight = .false.

contains

    !
    ! Read global parameters for forward modeling or inversion
    !
    subroutine read_parameters

        character(len=1024) :: tmpparam
        integer :: i
        real :: valmin, valmax

        ! Get name of the parameter file
        call get_command_argument(1, file_parameter)

        ! Check existence of the parameter file
        if (.not. file_exists(file_parameter)) then
            call warn(date_time_compact()//' <read_parameters> Error: Parameter file = ' &
                //tidy(file_parameter)//' not found. Exiting. ')
            call mpibarrier
            call mpistop
        end if

        ! Copy parameter file to working directory and use it as the working parameter file
        select case (which_program)
            case ('eikonal')
                call readpar_string(file_parameter, 'dir_synthetic', dir_synthetic, './data_synthetic')
                call make_directory(dir_synthetic)
                if (rankid == 0) then
                    tmpparam = tidy(dir_synthetic)//'/parameters.'//tidy(which_program)//'.'//date_time_string()
                    call copy_file(file_parameter, tmpparam)
                    file_parameter = tmpparam
                    call warn(' Parameter file: '//tidy(file_parameter))
                end if
            case ('fatt', 'tloc')
                call readpar_string(file_parameter, 'dir_working', dir_working, './test')
                call make_directory(dir_working)
                if (rankid == 0) then
                    tmpparam = tidy(dir_working)//'/parameters.'//tidy(which_program)//'.'//date_time_string()
                    call copy_file(file_parameter, tmpparam)
                    file_parameter = tmpparam
                    call warn(' Parameter file: '//tidy(file_parameter))
                end if
        end select
        call mpibarrier

        call bcast(file_parameter)

        ! Stop on parameters of the Cartesian code or of unsupported features
        call check_spherical_parameter_file

        ! Radius of the spherical Earth, in the length unit of the model
        call readpar_double(file_parameter, 'earth_radius', earth_radius, 6371000.0d0)
        call assert(earth_radius > 0, ' <read_parameters> Error: earth_radius must be > 0')
        if (rankid == 0) then
            call warn(' Earth radius = '//num2str(real(earth_radius), '(es)'))
        end if

        ! which kind of medium
        call readpar_string(file_parameter, 'which_medium', which_medium, 'acoustic-iso')
        call readpar_float(file_parameter, 'min_vpvsratio', min_vpvsratio, 1.1)
        call readpar_float(file_parameter, 'max_vpvsratio', max_vpvsratio, 9.0)
        call readpar_float(file_parameter, 'vpvsratio_smooth_lon', vpvsratio_smoothx, 0.0)
        call readpar_float(file_parameter, 'vpvsratio_smooth_lat', vpvsratio_smoothy, 0.0)
        call readpar_float(file_parameter, 'vpvsratio_smooth_depth', vpvsratio_smoothz, 0.0)

        call readpar_float(file_parameter, 'misfit_threshold', misfit_threshold, float_huge)

        call readpar_int(file_parameter, 'ns', ns, 1)
        call readpar_string(file_parameter, 'file_geometry', file_geometry, '', required=.true.)
        call readpar_nint(file_parameter, 'src_index', shot_index, [1, 1, ns])
        shot_min = max(shot_index(1), 1)
        shot_every = max(shot_index(2), 1)
        shot_max = min(shot_index(3), ns)
        call readpar_nint(file_parameter, 'rec_index', rec_index, [1, 1, int4_huge])
        rec_min = max(rec_index(1), 1)
        rec_every = max(rec_index(2), 1)
        rec_max = min(rec_index(3), int4_huge)
        call readpar_nint(file_parameter, 'rec_exclude', rec_exclude, [0])
        call readpar_float(file_parameter, 'offset_min', offset_min, 0.0)
        call readpar_float(file_parameter, 'offset_max', offset_max, float_huge)
        call readpar_nint(file_parameter, 'sid_select', sid_select, [0])
        call readpar_nint(file_parameter, 'src_select', src_select, [0])
        call readpar_nint(file_parameter, 'sid_exclude', sid_exclude, [0])
        call readpar_nint(file_parameter, 'src_exclude', src_exclude, [0])
        call readpar_int(file_parameter, 'sid_min', sid_min, 0)
        call readpar_int(file_parameter, 'sid_max', sid_max, int4_huge)
        call readpar_float(file_parameter, 'slon_min', sxmin, -float_huge)
        call readpar_float(file_parameter, 'slon_max', sxmax, +float_huge)
        call readpar_float(file_parameter, 'slat_min', symin, -float_huge)
        call readpar_float(file_parameter, 'slat_max', symax, +float_huge)
        call readpar_float(file_parameter, 'sdepth_min', szmin, -float_huge)
        call readpar_float(file_parameter, 'sdepth_max', szmax, +float_huge)
        call readpar_float(file_parameter, 'rlon_min', rxmin, -float_huge)
        call readpar_float(file_parameter, 'rlon_max', rxmax, +float_huge)
        call readpar_float(file_parameter, 'rlat_min', rymin, -float_huge)
        call readpar_float(file_parameter, 'rlat_max', rymax, +float_huge)
        call readpar_float(file_parameter, 'rdepth_min', rzmin, -float_huge)
        call readpar_float(file_parameter, 'rdepth_max', rzmax, +float_huge)
        call readpar_logical(file_parameter, 'yn_exchange_sr', yn_exchange_sr, ifelse(which_program == 'tloc', .true., .false.))

        call readpar_string(file_parameter, 'dir_working', dir_working, './test')
        call readpar_string(file_parameter, 'dir_record', dir_record, './data')

        call readpar_nstring(file_parameter, 'model_regularization_method', model_regularization_method, [''])
        if (model_regularization_method(1) /= '') then
            yn_regularize_model = .true.
        end if
        call readpar_int(file_parameter, 'niter_max', niter_max, 100)

        call readpar_string(file_parameter, 'search_method', search_method, 'CG')

        call readpar_string(file_parameter, 'misfit_type', misfit_type, 'ad')
        call assert(any_in([misfit_type], ['absolute-difference', 'ad', 'double-difference', 'dd']), &
            ' <read_parameter> Error: misfit_type must be one of absolute-difference, ad, double-difference, dd')

        call readpar_logical(file_parameter, 'yn_adp_lon', yn_adpx, .false.)
        call readpar_logical(file_parameter, 'yn_adp_lat', yn_adpy, .false.)
        call readpar_logical(file_parameter, 'yn_adp_depth', yn_adpz, .false.)
        call readpar_float(file_parameter, 'adp_extra_lon', adpextrax, 0.0)
        call readpar_float(file_parameter, 'adp_extra_lat', adpextray, 0.0)
        call readpar_float(file_parameter, 'adp_extra_depth', adpextraz, 0.0)
        call readpar_float(file_parameter, 'adp_taper_lon', adptaperx, 0.0)
        call readpar_float(file_parameter, 'adp_taper_lat', adptapery, 0.0)
        call readpar_float(file_parameter, 'adp_taper_depth', adptaperz, 0.0)

        call readpar_string(file_parameter, 'step_size_method', step_size_method, 'linear')
        call readpar_logical(file_parameter, 'yn_precond', yn_precond, .true.)
        call readpar_float(file_parameter, 'precond_eps', precond_eps, 1.0e-4)
        call readpar_logical(file_parameter, 'yn_output_dws', yn_output_dws, .false.)

        call readpar_nfloat(file_parameter, 'snaps', snaps, [-1.0])
        call readpar_string(file_parameter, 'dir_snapshot', dir_snapshot, './snapshot')
        call readpar_string(file_parameter, 'dir_synthetic', dir_synthetic, './data_synthetic')
        call assert(dir_snapshot /= dir_synthetic, ' <read_parameters> Error: dir_snapshot and dir_synthetic cannot be the same. ')

        call readpar_logical(file_parameter, 'verbose', verbose, .false.)
        call readpar_float(file_parameter, 'sweep_stop_threshold', sweep_stop_threshold, 1.0e-4)
        call readpar_int(file_parameter, 'sweep_niter_max', sweep_niter_max, 10)
        ! How the traveltime residuals at the receivers enter the adjoint-state
        ! equation; the spherical branch supports only linear injection, for
        ! FATT and TLOC
        call readpar_string(file_parameter, 'adjoint_receiver', adjoint_receiver, 'linear')
        call assert(adjoint_receiver == 'linear', &
            ' <read_parameters> Error: the spherical branch supports only adjoint_receiver = linear')

        call readpar_string(file_parameter, 'file_data_misfit', file_datamisfit, tidy(dir_working)//'/data_misfit.txt')
        call readpar_string(file_parameter, 'file_shot_misfit', file_shotmisfit, tidy(dir_working)//'/shot_misfit.bin')

        call readpar_logical(file_parameter, 'yn_continue', yn_continue_inv, .false.)
        if (yn_continue_inv .and. file_exists(file_datamisfit)) then
            resume_from_iter = max(1, count_nonempty_lines(file_datamisfit))
        else
            call readpar_int(file_parameter, 'resume_from_iter', resume_from_iter, 1)
        end if

        call readpar_logical(file_parameter, 'uniform_processing', uniform_processing, .true.)
        call readpar_logical(file_parameter, 'yn_flat_stop', yn_flat_stop, .false.)

        if (which_program == 'eikonal') then

            call readpar_nstring(file_parameter, 'model_name', model_name, [''], required=.true.)
            nmodel = size(model_name)

        else

            call readpar_nstring(file_parameter, 'model_update', model_name, ['vp'])
            call assert(model_name(1) /= '', ' <read_parameter> Error: model_update cannot be empty or start with null')
            nmodel = size(model_name)
            call readpar_nstring(file_parameter, 'model_aux', model_name_aux, [''])
            if (model_name_aux(1) == '') then
                nmodel_aux = 0
            else
                nmodel_aux = size(model_name_aux)
            end if
            model_min = zeros(nmodel)
            model_max = zeros(nmodel)
            model_step_max = zeros(nmodel)
            model_step = zeros(nmodel)
            allocate(model_search_method(1:nmodel))
            do i = 1, nmodel
                select case (remove_string_after(model_name(i), ['c', 'C']))
                    case ('vp', 'vs')
                        valmin = 0.0
                        valmax = 1.0e5
                    case ('epsilon', 'delta', 'gamma', 'eps', 'del', 'gam')
                        valmin = 0.0
                        valmax = 0.5
                    case ('c', 'C')
                        valmin = 0.0
                        valmax = 1.0e9
                end select
                call readpar_float(file_parameter, 'min_'//tidy(model_name(i)), model_min(i), valmin)
                call readpar_float(file_parameter, 'max_'//tidy(model_name(i)), model_max(i), valmax)
                ! The program allows using different search method for different parameters
                call readpar_string(file_parameter, 'search_method_'//tidy(model_name(i)), model_search_method(i), search_method)
            end do

            if (which_program == 'tloc' .and. .not.any_in(['st0'], model_name) .and. (misfit_type == 'dd' .or. misfit_type == 'double-difference')) then
                yn_dd_no_st0 = .true.
            else
                yn_dd_no_st0 = .false.
            end if

        end if

        select case (which_medium)
            case ('acoustic-iso', 'acoustic-tti')
                call readpar_nstring(file_parameter, 'data_name', data_name, ['p'])
            case('elastic-iso', 'elastic-tti')
                call readpar_nstring(file_parameter, 'data_name', data_name, ['p', 's'])
        end select
        ndata = size(data_name)

        call readpar_nfloat(file_parameter, 'misfit_weight', misfit_weight, ones(1))

    end subroutine read_parameters

    subroutine read_additional_parameters_tloc

        integer :: i
        real :: valmin, valmax

        do i = 1, nmodel
            if (any(model_name(i) == [character(len=6) :: 'slon', 'slat', 'sdepth', 'st0'])) then
                select case (model_name(i))
                    case ('slon')
                        valmin = ox
                        valmax = ox + (nx - 1)*dx
                    case ('slat')
                        valmin = oy
                        valmax = oy + (ny - 1)*dy
                    case ('sdepth')
                        valmin = oz
                        valmax = oz + (nz - 1)*dz
                    case ('st0')
                        valmin = 0.0
                        valmax = 2.0e9
                end select
                call readpar_float(file_parameter, 'min_'//tidy(model_name(i)), model_min(i), valmin)
                call readpar_float(file_parameter, 'max_'//tidy(model_name(i)), model_max(i), valmax)
                ! The program allows using different search method for different parameters
                call readpar_string(file_parameter, 'search_method_'//tidy(model_name(i)), model_search_method(i), search_method)
            end if
        end do

    end subroutine read_additional_parameters_tloc

    !
    !> Load geometry
    !
    subroutine load_geometry

        integer :: i, j, l
        character(len=1024) :: shot_file_geometry, dir_geometry
        logical, allocatable, dimension(:) :: qs
        integer, allocatable, dimension(:) :: usid
        type(source_receiver_geometry), allocatable, dimension(:) :: g
        character(len=512) :: geom_line
        real :: wp, ws
        integer :: ios5
        double precision, dimension(1:3) :: xc, xr

        if (rankid == 0) then
            call warn(date_time_compact()//' Loading geometry... ')
        end if

        if (allocated(gmtr)) then
            deallocate (gmtr)
        end if
        allocate (gmtr(1:ns))
        if (.not. file_exists(file_geometry)) then
            call warn(date_time_compact()//' @'//get_hostname() //' Error: Geometry file not found. ')
            stop
        else
            dir_geometry = get_file_directory(file_geometry)
        end if

        open (3, file=tidy(file_geometry), action='read', status='old')
        do i = 1, ns

            read(3, *) shot_file_geometry

            open (4, file=tidy(dir_geometry)//tidy(shot_file_geometry), action='read', status='old')

            ! Read geometry ID
            read(4, *) gmtr(i)%id

            ! Read source parameters
            read(4, *) gmtr(i)%ns
            allocate(gmtr(i)%srcr(1:gmtr(i)%ns))
            do j = 1, gmtr(i)%ns
                read (4, *) gmtr(i)%srcr(j)%x, gmtr(i)%srcr(j)%y, gmtr(i)%srcr(j)%z, gmtr(i)%srcr(j)%t0
            end do

            ! Read receiver parameters
            read(4, *) gmtr(i)%nr
            allocate(gmtr(i)%recr(1:gmtr(i)%nr))
            do j = 1, gmtr(i)%nr
                read(4, '(A)') geom_line
                ! Try reading 5 values (x y z weight_p weight_s); ios5 /= 0 means only 4 values
                read(geom_line, *, iostat=ios5) &
                    gmtr(i)%recr(j)%x, gmtr(i)%recr(j)%y, gmtr(i)%recr(j)%z, wp, ws
                if (ios5 == 0) then
                    ! 5-column format: separate P- and S-wave weights
                    gmtr(i)%recr(j)%weight_p = wp
                    gmtr(i)%recr(j)%weight_s = ws
                    ! master enable is non-zero when either wave has data
                    gmtr(i)%recr(j)%weight = max(wp, ws)
                    yn_separate_ps_weight = .true.
                else
                    ! 4-column format: single unified weight (backward-compatible)
                    read(geom_line, *) &
                        gmtr(i)%recr(j)%x, gmtr(i)%recr(j)%y, gmtr(i)%recr(j)%z, &
                        gmtr(i)%recr(j)%weight
                    gmtr(i)%recr(j)%weight_p = gmtr(i)%recr(j)%weight
                    gmtr(i)%recr(j)%weight_s = gmtr(i)%recr(j)%weight
                end if
            end do

            close(4)

            if ((mod(i, max(nint(ns/10.0), 1)) == 0 .or. i == ns .or. i == 1) .and. rankid == 0) then
                call warn(date_time_compact()//' Loading geometry '//num2str(i)//' of '//num2str(ns))
            end if

        end do
        close(3)

        ! For acoustic media, S-wave weight is not applicable
        if (which_medium == 'acoustic-iso' .or. which_medium == 'acoustic-tti') then
            do i = 1, ns
                gmtr(i)%recr(:)%weight_s = 0.0
            end do
        end if

        ! Check validity of geometry
        do i = 1, ns

            ! Check source
            l = 0
            do j =  1, gmtr(i)%ns
                if ( &
                        gmtr(i)%srcr(j)%x >= xmin .and. gmtr(i)%srcr(j)%x <= xmax .and. &
                        gmtr(i)%srcr(j)%y >= ymin .and. gmtr(i)%srcr(j)%y <= ymax .and. &
                        gmtr(i)%srcr(j)%z >= zmin .and. gmtr(i)%srcr(j)%z <= zmax .and. &
                        gmtr(i)%srcr(j)%x >= sxmin .and. gmtr(i)%srcr(j)%x <= sxmax .and. &
                        gmtr(i)%srcr(j)%y >= symin .and. gmtr(i)%srcr(j)%y <= symax .and. &
                        gmtr(i)%srcr(j)%z >= szmin .and. gmtr(i)%srcr(j)%z <= szmax .and. &
                        gmtr(i)%id >= sid_min .and. gmtr(i)%id <= sid_max) then
                    l = l + 1
                else
                    gmtr(i)%srcr(j)%amp = 0
                end if
            end do

            if (l == 0) then
                gmtr(i)%ns = 0
                if (rankid == 0) then
                    call warn(date_time_compact()//' Warning: Shot '//num2str(i)//' (ID: '//num2str(gmtr(i)%id)//') removed - no valid sources')
                end if
            end if

            ! Centroid of the source points of this shot, averaged in
            ! Earth-centred Cartesian coordinates
            xc = 0.0d0
            do j = 1, gmtr(i)%ns
                xc = xc + lonlatdepth_to_ecef(dble(gmtr(i)%srcr(j)%x), dble(gmtr(i)%srcr(j)%y), dble(gmtr(i)%srcr(j)%z))
            end do
            xc = xc/max(gmtr(i)%ns, 1)

            ! Check receiver
            l = 0
            do j = 1, gmtr(i)%nr
                ! Absolute offset: the straight-line distance from the receiver
                ! to the centroid of the source points
                xr = lonlatdepth_to_ecef(dble(gmtr(i)%recr(j)%x), dble(gmtr(i)%recr(j)%y), dble(gmtr(i)%recr(j)%z))
                gmtr(i)%recr(j)%aoff = real(norm2(xr - xc))

                ! Check if receiver in the range
                if ( &
                        gmtr(i)%recr(j)%x >= xmin .and. gmtr(i)%recr(j)%x <= xmax .and. &
                        gmtr(i)%recr(j)%y >= ymin .and. gmtr(i)%recr(j)%y <= ymax .and. &
                        gmtr(i)%recr(j)%z >= zmin .and. gmtr(i)%recr(j)%z <= zmax .and. &
                        gmtr(i)%recr(j)%x >= rxmin .and. gmtr(i)%recr(j)%x <= rxmax .and. &
                        gmtr(i)%recr(j)%y >= rymin .and. gmtr(i)%recr(j)%y <= rymax .and. &
                        gmtr(i)%recr(j)%z >= rzmin .and. gmtr(i)%recr(j)%z <= rzmax .and. &
                        gmtr(i)%recr(j)%aoff >= offset_min .and. gmtr(i)%recr(j)%aoff <= offset_max .and. &
                        j >= rec_min .and. j <= rec_max .and. mod(j - rec_min, rec_every) == 0 .and. &
                        .not. any(rec_exclude == j) .and. &
                        gmtr(i)%recr(j)%weight /= 0) then
                    l = l + 1
                else
                    gmtr(i)%recr(j)%weight   = 0.0
                    gmtr(i)%recr(j)%weight_p = 0.0
                    gmtr(i)%recr(j)%weight_s = 0.0
                end if

            end do
            close (4)

            ! check if any receiver in the model
            if (l == 0) then
                gmtr(i)%nr = 0
                if (rankid == 0) then
                    call warn(date_time_compact()//' Warning: Shot '//num2str(i)//' (ID: '//num2str(gmtr(i)%id)//') removed - no valid receivers')
                end if
            end if

        end do

        ! Remove common-shot gathers that are not qualified
        qs = falses(ns)
        ! Select sources based on offset, gain and representational id restrictions
        do i = shot_min, shot_max, shot_every
            if (gmtr(i)%ns /= 0 .and. gmtr(i)%nr /= 0) then
                qs(i) = .true.
            end if
        end do
        if (allocated(src_select)) then
            if (any(src_select /= 0)) then
                ! Select sources based on their representional id
                qs = .false.
                do i = 1, ns
                    if (gmtr(i)%ns /= 0 .and. gmtr(i)%nr /= 0 .and. any(src_select == i)) then
                        qs(i) = .true.
                    end if
                end do
            end if
        end if
        if (allocated(sid_select)) then
            if (any(sid_select /= 0)) then
                ! Select source based on their field id
                qs = .false.
                do i = 1, ns
                    if (gmtr(i)%ns /= 0 .and. gmtr(i)%nr /= 0 .and. any(sid_select == gmtr(i)%id)) then
                        qs(i) = .true.
                    end if
                end do
            end if
        end if
        if (allocated(src_exclude)) then
            if (any(src_exclude /= 0)) then
                ! Exclude sources based on their representional id
                do i = 1, ns
                    if (gmtr(i)%ns /= 0 .and. gmtr(i)%nr /= 0 .and. qs(i) .and. any(src_exclude == i)) then
                        qs(i) = .false.
                    end if
                end do
            end if
        end if
        if (allocated(sid_exclude)) then
            if (any(sid_exclude /= 0)) then
                ! Exclude source based on their field id
                do i = 1, ns
                    if (gmtr(i)%ns /= 0 .and. gmtr(i)%nr /= 0 .and. qs(i) .and. any(sid_exclude == gmtr(i)%id)) then
                        qs(i) = .false.
                    end if
                end do
            end if
        end if

        ! If no qualified shots, then stop
        if (.not. any(qs)) then
            call warn(date_time_compact()//' @'//get_hostname()//' Error: No source or receiver in the model. ')
            stop
        end if

        ! Select all qualified shots
        allocate(g(1:ns))
        l = 1
        do i = 1, ns
            if (qs(i)) then
                g(l) = gmtr(i)
                l = l + 1
            end if
        end do

        gmtr = g(1:l - 1)
        deallocate(g)

        ns = size(gmtr)

        ! If any duplicate source id, then stop
        if (ns >= 2) then
            call alloc_array(sid_select, [1, ns])
            do i = 1, ns
                sid_select(i) = gmtr(i)%id
            end do
            usid = unique(sid_select)
            if (size(usid) < ns) then
                call warn(date_time_compact()//' <load_geometry> Erorr: Duplicate source ID found. ')
                stop
            end if
        end if

        ! Print info
        if (rankid == 0) then
            call warn(date_time_compact()//' Geometry is loaded with '//num2str(ns)//' sources. ')
        end if
        call mpibarrier

    end subroutine load_geometry

    !
    !> Get and set model dimensions
    !
    subroutine set_regular_space

        ! Original model dimensions
        call readpar_int(file_parameter, 'nlon', nx0, 0, required=.true.)
        call readpar_int(file_parameter, 'nlat', ny0, 0, required=.true.)
        call readpar_int(file_parameter, 'ndepth', nz0, 0, required=.true.)
        call readpar_float(file_parameter, 'dlon', dx0, 1.0, required=.true.)
        call readpar_float(file_parameter, 'dlat', dy0, 1.0, required=.true.)
        call readpar_float(file_parameter, 'ddepth', dz0, 1.0, required=.true.)
        call readpar_float(file_parameter, 'olon', ox0, 0.0, required=.true.)
        call readpar_float(file_parameter, 'olat', oy0, 0.0, required=.true.)
        call readpar_float(file_parameter, 'odepth', oz0, 0.0)

        ! Target model dimensions
        call readpar_int(file_parameter, 'nnlon', nx, nx0)
        call readpar_int(file_parameter, 'nnlat', ny, ny0)
        call readpar_int(file_parameter, 'nndepth', nz, nz0)
        call readpar_float(file_parameter, 'ddlon', dx, dx0)
        call readpar_float(file_parameter, 'ddlat', dy, dy0)
        call readpar_float(file_parameter, 'dddepth', dz, dz0)
        call readpar_float(file_parameter, 'oolon', ox, ox0)
        call readpar_float(file_parameter, 'oolat', oy, oy0)
        call readpar_float(file_parameter, 'oodepth', oz, oz0)

        ! Target model ranges
        call readpar_float(file_parameter, 'lonmin', xmin, ox)
        call readpar_float(file_parameter, 'lonmax', xmax, ox + (nx - 1)*dx)
        call assert(xmax >= xmin, ' <set_regular_space> Error: lonmax must be >= lonmin')
        call readpar_float(file_parameter, 'latmin', ymin, oy)
        call readpar_float(file_parameter, 'latmax', ymax, oy + (ny - 1)*dy)
        call assert(ymax >= ymin, ' <set_regular_space> Error: latmax must be >= latmin')
        call readpar_float(file_parameter, 'depthmin', zmin, oz)
        call readpar_float(file_parameter, 'depthmax', zmax, oz + (nz - 1)*dz)
        call assert(zmax >= zmin, ' <set_regular_space> Error: depthmax must be >= depthmin')

        ! Target model origins
        ! Note: n# are computed with nint rather than ceiling, because ceiling is
        ! sensitive to floating-point round-off. For instance, with ox = 0,
        ! dx = 11.2220842 and nx = 16, the default xmax = ox + (nx - 1)*dx gives
        ! (xmax - xmin)/dx = 15.000001 instead of 15, and ceiling would add a
        ! spurious grid point. The price is that in rare cases, when xmax - xmin
        ! is not a multiple of dx, the model will be shortened by up to half a
        ! grid interval, e.g., xmax - xmin = 120.4 with dx = 10 gives a model
        ! ending at xmin + 120 rather than xmin + 130. A rigorous solution, which
        ! always covers xmax and still ignores round-off, is
        !     nx = ceiling((xmax - xmin - 4*spacing(max(abs(xmin), abs(xmax))))/dx) + 1
        ! where spacing(x) is the gap between adjacent floats near x. The round-off
        ! in xmax - xmin is at most about 2 such gaps, so 4 gives a safety margin.
        ox = xmin
        oy = ymin
        oz = zmin
        nx = nint((xmax - xmin)/dx) + 1
        ny = nint((ymax - ymin)/dy) + 1
        nz = nint((zmax - zmin)/dz) + 1

        ! Check if requires model interpolation
        if (nx /= nx0 .or. ny /= ny0 .or. nz /= nz0 .or. &
                dx /= dx0 .or. dy /= dy0 .or. dz /= dz0 .or. &
                ox /= ox0 .or. oy /= oy0 .or. oz /= oz0) then
            require_model_interp = .true.
        end if

        ! Check non-negativeness
        call assert(nx >= 1 .and. ny >= 1 .and. nz >= 1, ' <set_regular_space> Error: nlon, nlat and ndepth must be >= 1')
        call assert(dx > 0 .and. dy > 0 .and. dz > 0, ' <set_regular_space> Error: dlon, dlat and ddepth must be > 0')

        ! The exact geometry of the spherical grid, set once here
        call init_spherical_grid(nx, ny, nz, dble(ox), dble(oy), dble(oz), dble(dx), dble(dy), dble(dz))

        ! Print original and target model dimensions
        if (rankid == 0) then
            call warn(' ')
            call warn(' Original model dimensions (lon, lat, depth): ')
            call warn('      nlon, nlat, ndepth = '//num2str(nx0)//', '//num2str(ny0)//', '//num2str(nz0))
            call warn('      dlon, dlat, ddepth = '//num2str(dx0, '(es)')//', '//num2str(dy0, '(es)')//', '//num2str(dz0, '(es)'))
            call warn('      olon, olat, odepth = '//num2str(ox0, '(es)')//', '//num2str(oy0, '(es)')//', '//num2str(oz0, '(es)'))
            call warn('      elon, elat, edepth = ' &
                //num2str(ox0 + (nx0 - 1)*dx0, '(es)')//', ' &
                //num2str(oy0 + (ny0 - 1)*dy0, '(es)')//', ' &
                //num2str(oz0 + (nz0 - 1)*dz0, '(es)'))
            call warn(' ')
            call warn(' Target model dimensions (lon, lat, depth): ')
            call warn('      nlon, nlat, ndepth = '//num2str(nx)//', '//num2str(ny)//', '//num2str(nz))
            call warn('      dlon, dlat, ddepth = '//num2str(dx, '(es)')//', '//num2str(dy, '(es)')//', '//num2str(dz, '(es)'))
            call warn('      olon, olat, odepth = '//num2str(ox, '(es)')//', '//num2str(oy, '(es)')//', '//num2str(oz, '(es)'))
            call warn('      elon, elat, edepth = ' &
                //num2str(ox + (nx - 1)*dx, '(es)')//', ' &
                //num2str(oy + (ny - 1)*dy, '(es)')//', ' &
                //num2str(oz + (nz - 1)*dz, '(es)'))
            call warn(' ')
        end if

    end subroutine set_regular_space

    !
    !> Set adaptive range
    !
    subroutine set_adaptive_model_range(geom)

        type(source_receiver_geometry), intent(inout) :: geom

        double precision, dimension(1:2) :: mpd
        real :: extra_lon, extra_lat

        ! The margins adp_extra_lon and adp_extra_lat are lengths; in degrees,
        ! they use the smallest length of a degree anywhere in the grid, so
        ! that the margin is never shorter than asked for
        mpd = min_length_per_degree()
        extra_lon = real(adpextrax/mpd(1))
        extra_lat = real(adpextray/mpd(2))

        ! range selection
        if (yn_adpx) then
            shot_xbeg = max(xmin, min( &
                minval(geom%recr(:)%x, mask=geom%recr(:)%weight /= 0), &
                minval(geom%srcr(:)%x)) - extra_lon)
            shot_xend = min(xmax, max( &
                maxval(geom%recr(:)%x, mask=geom%recr(:)%weight /= 0), &
                maxval(geom%srcr(:)%x)) + extra_lon)
            shot_xbeg = floor((shot_xbeg - ox)/dx)*dx + ox
            shot_xend = ceiling((shot_xend - ox)/dx)*dx + ox
        else
            shot_xbeg = ox
            shot_xend = ox + (nx - 1)*dx
        end if
        where (geom%srcr(:)%x < shot_xbeg) geom%srcr(:)%amp = 0.0
        where (geom%srcr(:)%x > shot_xend) geom%srcr(:)%amp = 0.0
        where (geom%recr(:)%x < shot_xbeg) geom%recr(:)%weight   = 0.0
        where (geom%recr(:)%x < shot_xbeg) geom%recr(:)%weight_p = 0.0
        where (geom%recr(:)%x < shot_xbeg) geom%recr(:)%weight_s = 0.0
        where (geom%recr(:)%x > shot_xend) geom%recr(:)%weight   = 0.0
        where (geom%recr(:)%x > shot_xend) geom%recr(:)%weight_p = 0.0
        where (geom%recr(:)%x > shot_xend) geom%recr(:)%weight_s = 0.0

        if (yn_adpy) then
            shot_ybeg = max(ymin, min( &
                minval(geom%recr(:)%y, mask=geom%recr(:)%weight /= 0), &
                minval(geom%srcr(:)%y)) - extra_lat)
            shot_yend = min(ymax, max( &
                maxval(geom%recr(:)%y, mask=geom%recr(:)%weight /= 0), &
                maxval(geom%srcr(:)%y)) + extra_lat)
            shot_ybeg = floor((shot_ybeg - oy)/dy)*dy + oy
            shot_yend = ceiling((shot_yend - oy)/dy)*dy + oy
        else
            shot_ybeg = oy
            shot_yend = oy + (ny - 1)*dy
        end if
        where (geom%srcr(:)%y < shot_ybeg) geom%srcr(:)%amp = 0.0
        where (geom%srcr(:)%y > shot_yend) geom%srcr(:)%amp = 0.0
        where (geom%recr(:)%y < shot_ybeg) geom%recr(:)%weight   = 0.0
        where (geom%recr(:)%y < shot_ybeg) geom%recr(:)%weight_p = 0.0
        where (geom%recr(:)%y < shot_ybeg) geom%recr(:)%weight_s = 0.0
        where (geom%recr(:)%y > shot_yend) geom%recr(:)%weight   = 0.0
        where (geom%recr(:)%y > shot_yend) geom%recr(:)%weight_p = 0.0
        where (geom%recr(:)%y > shot_yend) geom%recr(:)%weight_s = 0.0

        if (yn_adpz) then
            shot_zbeg = max(zmin, min( &
                minval(geom%recr(:)%z, mask=geom%recr(:)%weight /= 0), &
                minval(geom%srcr(:)%z)) - adpextraz)
            shot_zend = min(zmax, max( &
                maxval(geom%recr(:)%z, mask=geom%recr(:)%weight /= 0), &
                maxval(geom%srcr(:)%z)) + adpextraz)
            shot_zbeg = floor((shot_zbeg - oz)/dz)*dz + oz
            shot_zend = ceiling((shot_zend - oz)/dz)*dz + oz
        else
            shot_zbeg = oz
            shot_zend = oz + (nz - 1)*dz
        end if
        where (geom%srcr(:)%z < shot_zbeg) geom%srcr(:)%amp = 0.0
        where (geom%srcr(:)%z > shot_zend) geom%srcr(:)%amp = 0.0
        where (geom%recr(:)%z < shot_zbeg) geom%recr(:)%weight   = 0.0
        where (geom%recr(:)%z < shot_zbeg) geom%recr(:)%weight_p = 0.0
        where (geom%recr(:)%z < shot_zbeg) geom%recr(:)%weight_s = 0.0
        where (geom%recr(:)%z > shot_zend) geom%recr(:)%weight   = 0.0
        where (geom%recr(:)%z > shot_zend) geom%recr(:)%weight_p = 0.0
        where (geom%recr(:)%z > shot_zend) geom%recr(:)%weight_s = 0.0

        ! Range for each shot
        shot_nxbeg = clip(nint((shot_xbeg - ox)/dx) + 1, 1, nx)
        shot_nxend = clip(nint((shot_xend - ox)/dx) + 1, 1, nx)
        shot_nybeg = clip(nint((shot_ybeg - oy)/dy) + 1, 1, ny)
        shot_nyend = clip(nint((shot_yend - oy)/dy) + 1, 1, ny)
        shot_nzbeg = clip(nint((shot_zbeg - oz)/dz) + 1, 1, nz)
        shot_nzend = clip(nint((shot_zend - oz)/dz) + 1, 1, nz)

        shot_nx = shot_nxend - shot_nxbeg + 1
        shot_ny = shot_nyend - shot_nybeg + 1
        shot_nz = shot_nzend - shot_nzbeg + 1

        if (verbose) then
            call warn(date_time_compact()//' Shot '//num2str(gmtr(ishot)%id)//' dimensions:')
            call warn(date_time_compact()//'      lonmin, lonmax = ' &
                //num2str(shot_xbeg, '(es)')//', '//num2str(shot_xend, '(es)'))
            if (shot_ny > 1) then
                call warn(date_time_compact()//'      latmin, latmax = ' &
                    //num2str(shot_ybeg, '(es)')//', '//num2str(shot_yend, '(es)'))
            end if
            call warn(date_time_compact()//'      depthmin, depthmax = ' &
                //num2str(shot_zbeg, '(es)')//', '//num2str(shot_zend, '(es)'))
            if (shot_ny > 1) then
                call warn(date_time_compact()//'      nlon, nlat, ndepth = ' &
                    //num2str(shot_nx)//', '//num2str(shot_ny)//', '//num2str(shot_nz))
            else
                call warn(date_time_compact()//'      nlon, ndepth = ' &
                    //num2str(shot_nx)//', '//num2str(shot_nz))
            end if
            call warn(date_time_compact()//'      ilonmin, ilonmax = ' &
                //num2str(shot_nxbeg)//', '//num2str(shot_nxend))
            if (shot_ny > 1) then
                call warn(date_time_compact()//'      ilatmin, ilatmax = ' &
                    //num2str(shot_nybeg)//', '//num2str(shot_nyend))
            end if
            call warn(date_time_compact()//'      idepthmin, idepthmax = ' &
                //num2str(shot_nzbeg)//', '//num2str(shot_nzend))
            call warn(date_time_compact()//'      ns, nr = ' &
                //num2str(count(geom%srcr(:)%amp /= 0))//', '//num2str(count(geom%recr(:)%weight /= 0)))
        end if

    end subroutine set_adaptive_model_range

    !
    !> Stop if the parameter file uses a parameter name, model name or option
    !> value of the Cartesian LATTE that the spherical branch renamed, or a
    !> parameter of a feature that the spherical branch does not support. The
    !> file is read with the same rules as the FLIT readpar routines: one
    !> "name = value" per line, names are case-insensitive, a line with an
    !> empty value is not a parameter, and reading stops at a line "exit".
    !> Every offending line is reported before the program stops.
    !
    subroutine check_spherical_parameter_file

        character(len=1024) :: line, key, val
        integer :: funit, eq, ioerr, nbad

        nbad = 0
        open (newunit=funit, file=tidy(file_parameter), status='old', action='read')
        do
            read (funit, '(a)', iostat=ioerr) line
            if (ioerr /= 0) exit
            if (tidy(line) == 'exit') exit
            eq = index(line, '=')
            if (eq < 2) cycle
            key = to_lower(tidy(line(1:eq - 1)))
            val = to_lower(tidy(line(eq + 1:)))
            if (len_trim(val) == 0) cycle
            call check_spherical_parameter(key, val, nbad)
        end do
        close (funit)

        if (nbad > 0) then
            if (rankid == 0) then
                call warn(date_time_compact()//' <check_spherical_parameter_file> Error: ' &
                    //num2str(nbad)//' parameter(s) above are not valid in the spherical branch. Exiting. ')
            end if
            call mpibarrier
            call mpistop
        end if

    end subroutine check_spherical_parameter_file

    !
    !> Check one parameter of the parameter file; see check_spherical_parameter_file
    !
    subroutine check_spherical_parameter(key, val, nbad)

        character(len=*), intent(in) :: key, val
        integer, intent(inout) :: nbad

        character(len=64), dimension(10), parameter :: unsupported_names = [character(len=64) :: &
            'source_regularization_method', 'yn_save_traveltime_gradients', 'dir_gradient', &
            'nrefl', 'offset_min_refl', 'offset_max_refl', 'offset_max_refl_depth', &
            'step_max_refl', 'yn_refl_image', 'incident_wave']
        character(len=32), dimension(6), parameter :: unsupported_prefixes = [character(len=32) :: &
            'clustering_fit_', 'reg_ml_', 'raypath_', 'reg_andf_', 'refl_', 'reflector_']
        character(len=1024) :: newkey
        character(len=64), allocatable, dimension(:) :: items
        integer :: i

        ! Parameters of unsupported features; wavenumber filtering along
        ! latitude is not supported, because the spacing along a latitude
        ! line is not even
        if (any(key == unsupported_names) .or. index(key, 'dipfilt') > 0 .or. par_has_token(key, 'refl') &
                .or. par_ends_with(key, '_wavenum_lat') .or. par_ends_with(key, '_wavenum_lat_amps')) then
            call report_unsupported('parameter '//tidy(key))
            return
        end if
        do i = 1, size(unsupported_prefixes)
            if (par_starts_with(key, unsupported_prefixes(i))) then
                call report_unsupported('parameter '//tidy(key))
                return
            end if
        end do

        ! Values that hold model names
        if (any(key == [character(len=16) :: 'model_update', 'model_aux', 'model_name'])) then
            items = par_split_list(val)
            do i = 1, size(items)
                select case (items(i))
                    case ('sx')
                        call report_renamed('model name sx in '//tidy(key), 'slon')
                    case ('sy')
                        call report_renamed('model name sy in '//tidy(key), 'slat')
                    case ('sz')
                        call report_renamed('model name sz in '//tidy(key), 'sdepth')
                    case ('refl')
                        call report_unsupported('model name refl in '//tidy(key))
                end select
            end do
        end if

        ! Values that hold gradient processing options
        if (par_starts_with(key, 'process_')) then
            items = par_split_list(val)
            do i = 1, size(items)
                select case (items(i))
                    case ('rmsbalx')
                        call report_renamed('option rmsbalx in '//tidy(key), 'rmsbal_lon')
                    case ('rmsbaly')
                        call report_renamed('option rmsbaly in '//tidy(key), 'rmsbal_lat')
                    case ('rmsbalxy')
                        call report_renamed('option rmsbalxy in '//tidy(key), 'rmsbal_lonlat')
                    case ('rmsbalz')
                        call report_renamed('option rmsbalz in '//tidy(key), 'rmsbal_depth')
                    case ('dipfilt')
                        call report_unsupported('option dipfilt in '//tidy(key))
                end select
            end do
        end if

        ! Values that hold model regularization methods
        if (key == 'model_regularization_method') then
            items = par_split_list(val)
            do i = 1, size(items)
                if (items(i) == 'structure') then
                    call report_unsupported('option structure in '//tidy(key))
                end if
            end do
        end if

        ! Renamed parameter names
        newkey = spherical_parameter_name(key)
        if (newkey /= '') then
            call report_renamed('parameter '//tidy(key), tidy(newkey))
        end if

    contains

        subroutine report_unsupported(what)
            character(len=*), intent(in) :: what
            nbad = nbad + 1
            if (rankid == 0) then
                call warn(' Error: '//what//' belongs to a feature that the spherical branch does not support.')
            end if
        end subroutine report_unsupported

        subroutine report_renamed(what, replacement)
            character(len=*), intent(in) :: what, replacement
            nbad = nbad + 1
            if (rankid == 0) then
                call warn(' Error: '//what//' is a Cartesian name; the spherical branch uses '//replacement//'.')
            end if
        end subroutine report_renamed

    end subroutine check_spherical_parameter

    !
    !> The spherical name of a parameter name of the Cartesian LATTE, or an
    !> empty string if the name was not renamed. x, y and z become lon, lat
    !> and depth; the TLOC source parameters sx, sy and sz become slon, slat
    !> and sdepth wherever they appear as a part of a name, e.g., file_sx.
    !
    function spherical_parameter_name(key) result(newkey)

        character(len=*), intent(in) :: key
        character(len=1024) :: newkey

        character(len=1), dimension(3), parameter :: axes = ['x', 'y', 'z']
        character(len=5), dimension(3), parameter :: names = [character(len=5) :: 'lon', 'lat', 'depth']
        character(len=32), dimension(13) :: words
        character(len=64), allocatable, dimension(:) :: tokens
        character(len=1) :: a
        character(len=5) :: n
        character(len=32) :: w
        integer :: i, l
        logical :: changed

        newkey = ''

        do i = 1, 3

            a = axes(i)
            n = names(i)

            ! Grid, ranges, source and receiver selection, adaptive range,
            ! regularization and MPI decomposition
            if (key == 'n'//a) newkey = 'n'//tidy(n)
            if (key == 'd'//a) newkey = 'd'//tidy(n)
            if (key == 'o'//a) newkey = 'o'//tidy(n)
            if (key == 'nn'//a) newkey = 'nn'//tidy(n)
            if (key == 'dd'//a) newkey = 'dd'//tidy(n)
            if (key == 'oo'//a) newkey = 'oo'//tidy(n)
            if (key == a//'min') newkey = tidy(n)//'min'
            if (key == a//'max') newkey = tidy(n)//'max'
            if (key == 's'//a//'min') newkey = 's'//tidy(n)//'_min'
            if (key == 's'//a//'max') newkey = 's'//tidy(n)//'_max'
            if (key == 'r'//a//'min') newkey = 'r'//tidy(n)//'_min'
            if (key == 'r'//a//'max') newkey = 'r'//tidy(n)//'_max'
            if (key == 'rank'//a) newkey = 'rank_'//tidy(n)
            if (key == 'reg_smooth'//a) newkey = 'reg_smooth_'//tidy(n)
            if (par_starts_with(key, 'reg_smooth'//a//'_')) then
                newkey = 'reg_smooth_'//tidy(n)//'_'//key(len('reg_smooth'//a//'_') + 1:)
            end if

            ! Names that end with a word followed by the axis
            words = [character(len=32) :: 'smooth', 'taper', 'movingbal', 'medianfilt', 'rmsbal', &
                'wavenum', 'andf_rank', 'extra', 'adp', 'conemute', '', '', '']
            do l = 1, 10
                w = words(l)
                if (w == 'conemute' .and. a /= 'z') then
                    cycle
                end if
                if (par_ends_with(key, '_'//tidy(w)//a)) then
                    newkey = key(1:len_trim(key) - len_trim(w) - 1)//tidy(w)//'_'//tidy(n)
                end if
                if (w == 'wavenum' .and. par_ends_with(key, '_wavenum'//a//'_amps')) then
                    newkey = key(1:len_trim(key) - len('wavenum'//a//'_amps'))//'wavenum_'//tidy(n)//'_amps'
                end if
            end do

        end do

        ! The TLOC source parameters as parts of a name
        if (newkey == '') then
            tokens = par_split_tokens(key)
            changed = .false.
            do i = 1, size(tokens)
                select case (tokens(i))
                    case ('sx')
                        tokens(i) = 'slon'
                        changed = .true.
                    case ('sy')
                        tokens(i) = 'slat'
                        changed = .true.
                    case ('sz')
                        tokens(i) = 'sdepth'
                        changed = .true.
                end select
            end do
            if (changed) then
                newkey = tidy(tokens(1))
                do i = 2, size(tokens)
                    newkey = tidy(newkey)//'_'//tidy(tokens(i))
                end do
            end if
        end if

    end function spherical_parameter_name

    !
    !> Whether a string starts with a prefix
    !
    pure logical function par_starts_with(s, prefix)
        character(len=*), intent(in) :: s, prefix
        par_starts_with = .false.
        if (len_trim(s) >= len_trim(prefix)) then
            par_starts_with = s(1:len_trim(prefix)) == prefix(1:len_trim(prefix))
        end if
    end function par_starts_with

    !
    !> Whether a string ends with a suffix
    !
    pure logical function par_ends_with(s, suffix)
        character(len=*), intent(in) :: s, suffix
        integer :: ls, lf
        ls = len_trim(s)
        lf = len_trim(suffix)
        par_ends_with = .false.
        if (ls >= lf .and. lf > 0) then
            par_ends_with = s(ls - lf + 1:ls) == suffix(1:lf)
        end if
    end function par_ends_with

    !
    !> The parts of a name separated by underscores
    !
    function par_split_tokens(s) result(tokens)
        character(len=*), intent(in) :: s
        character(len=64), allocatable, dimension(:) :: tokens
        integer :: i, b, n
        n = 1
        do i = 1, len_trim(s)
            if (s(i:i) == '_') n = n + 1
        end do
        allocate (tokens(1:n))
        n = 0
        b = 1
        do i = 1, len_trim(s) + 1
            if (i > len_trim(s)) then
                n = n + 1
                tokens(n) = s(b:len_trim(s))
            else if (s(i:i) == '_') then
                n = n + 1
                tokens(n) = s(b:i - 1)
                b = i + 1
            end if
        end do
    end function par_split_tokens

    !
    !> Whether one of the underscore-separated parts of a name is a token
    !
    logical function par_has_token(s, token)
        character(len=*), intent(in) :: s, token
        character(len=64), allocatable, dimension(:) :: tokens
        tokens = par_split_tokens(s)
        par_has_token = any(tokens == token)
    end function par_has_token

    !
    !> The items of a comma-separated value, without blanks and without a
    !> leading iteration schedule "n:" or "n~m:" of the readpar_x routines
    !
    function par_split_list(s) result(items)
        character(len=*), intent(in) :: s
        character(len=64), allocatable, dimension(:) :: items
        integer :: i, b, n, c
        character(len=1024) :: t
        n = 1
        do i = 1, len_trim(s)
            if (s(i:i) == ',') n = n + 1
        end do
        allocate (items(1:n))
        n = 0
        b = 1
        do i = 1, len_trim(s) + 1
            if (i > len_trim(s) .or. s(min(i, len(s)):min(i, len(s))) == ',') then
                n = n + 1
                t = adjustl(s(b:i - 1))
                c = index(t, ':')
                if (c > 0) then
                    t = adjustl(t(c + 1:))
                end if
                items(n) = trim(t)
                b = i + 1
            end if
        end do
    end function par_split_list

end module parameters
