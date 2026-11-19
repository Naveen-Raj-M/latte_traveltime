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


module gradient

    use libflit
    use parameters
    use vars
    use utility
    use spherical_filter

    implicit none


    interface
        module subroutine compute_gradient_shots_fatt
        end subroutine
        module subroutine compute_gradient_shots_tloc
        end subroutine
        module subroutine compute_gradient_shots_tloc_dd
        end subroutine
    end interface

    private

    public :: zero_gradient
    public :: compute_gradient_shots
    public :: process_gradient
    public :: output_gradient
    public :: output_dws

contains

    !
    !> Initialize gradient arrays
    !
    subroutine zero_gradient

        integer :: i

        do i = 1, nmodel
            select case (model_name(i))
                case default
                    model_grad(i)%array = zeros(nz, ny, nx)
                case ('slon', 'slat', 'sdepth', 'st0')
                    model_grad(i)%array = zeros(nr_virtual, 1, 1)
            end select
        end do

        ! Initialize DWS arrays
        if (yn_output_dws .and. yn_precond) then
            do i = 1, nmodel
                select case (model_name(i))
                    case default
                        model_dws(i)%array = zeros(nz, ny, nx)
                    case ('slon', 'slat', 'sdepth', 'st0')
                        ! DWS not applicable for source parameters
                        model_dws(i)%array = zeros(nr_virtual, 1, 1)
                end select
            end do
        end if

    end subroutine zero_gradient

    !
    !> Compute gradients shot by shot
    !
    subroutine compute_gradient_shots

        select case (which_program)
            case ('fatt')
                call compute_gradient_shots_fatt
            case ('tloc')
                if (yn_dd_no_st0) then
                    call compute_gradient_shots_tloc_dd
                else
                    call compute_gradient_shots_tloc
                end if
        end select

    end subroutine compute_gradient_shots

    !
    !> Process gradient
    !
    subroutine process_gradient

        integer :: i

        do i = 1, nmodel

            if (.not. any(model_m(i)%name == [character(len=6) :: 'slon', 'slat', 'sdepth', 'st0'])) then

                if (uniform_processing) then
                    call process_model_single_parameter(model_grad(i)%array, 'grad', param_name=model_m(i)%name)
                else
                    call process_model_single_parameter(model_grad(i)%array, 'grad_'//tidy(model_name(i)), param_name=model_m(i)%name)
                end if

            end if

        end do

    end subroutine process_gradient

    !
    !> Output gradient
    !
    subroutine output_gradient

        integer :: i

        if (rankid == 0) then

            do i = 1, nmodel
                call output_array(model_grad(i)%array, &
                    dir_iter_model(iter)//'/grad_'//tidy(model_name(i))//'.bin')
            end do

            call warn(date_time_compact()//' >>>>>>>>>> Gradient saved. ')

        end if

    end subroutine output_gradient

    !
    !> Output DWS (Derivative Weighted Sum)
    !
    subroutine output_dws

        integer :: i

        if (yn_output_dws .and. yn_precond) then

            if (rankid == 0) then

                do i = 1, nmodel
                    ! Only output DWS for velocity model parameters (not source parameters)
                    if (.not. any(model_name(i) == [character(len=6) :: 'slon', 'slat', 'sdepth', 'st0'])) then
                        call output_array(model_dws(i)%array, &
                            dir_iter_model(iter)//'/dws_'//tidy(model_name(i))//'.bin')
                    end if
                end do

                call warn(date_time_compact()//' >>>>>>>>>> DWS saved. ')

            end if

        end if

    end subroutine output_dws


    !
    !> Process gradient of one shot
    !>
    !> Spherical branch: lengths are in the model length unit (meters).
    !> Smoothing, tapers, the cone mute and the adaptive-range taper use exact
    !> distances on the spherical grid (module spherical_filter). The window
    !> filters (movingbal, medianfilt) work on a fixed number of cells, which
    !> is the length over the reference cell size at the centre of the grid,
    !> top level (reference_spacing); laplacefilt works in grid-index space.
    !> Default lengths are the Cartesian defaults in cells times the
    !> reference cell size.
    !
    subroutine process_model_single_shot(ishot, w, ws, name)

        integer, intent(in) :: ishot
        real, dimension(:, :, :), intent(inout) :: w, ws
        character(len=*), intent(in) :: name

        integer :: i, j, k, ix, iy, jg, kg
        character(len=1024) :: dir_mask
        real :: shot_w_movingbalx, shot_w_movingbaly, shot_w_movingbalz
        real, allocatable, dimension(:) :: shot_w_taperx, shot_w_tapery, shot_w_taperz
        real :: shot_w_medianfiltx, shot_w_medianfilty, shot_w_medianfiltz
        real :: shot_w_smoothx, shot_w_smoothy, shot_w_smoothz
        real, allocatable, dimension(:) :: wavenums, wamps
        real, allocatable, dimension(:, :, :) :: w_mask
        character(len=32), allocatable, dimension(:) :: process_shot_grad_w
        real :: conemuter, conemutez, conemutepower, conemutetaper
        character(len=24) :: conemuteorigin
        character(len=1024) :: file_mask
        real, dimension(1:3) :: href
        real, dimension(1:6) :: tlen
        double precision :: depth, srcz, radius, dist, tp
        double precision, dimension(1:3) :: xc, uc, un

        ! Reference cell size: east-west, north-south, vertical
        href = real(reference_spacing())

        call readpar_nstring(file_parameter, 'process_shot_'//tidy(name)//'', process_shot_grad_w, [''])

        ws = return_normal(ws)

        do i = 1, size(process_shot_grad_w)

            select case (process_shot_grad_w(i))

                case ('smooth')
                    ! Gaussian smooth, exact distances; the lengths are the
                    ! standard deviations
                    call readpar_float(file_parameter, 'shot_'//tidy(name)//'_smooth_lon', shot_w_smoothx, 3*href(1))
                    call readpar_float(file_parameter, 'shot_'//tidy(name)//'_smooth_lat', shot_w_smoothy, 3*href(2))
                    call readpar_float(file_parameter, 'shot_'//tidy(name)//'_smooth_depth', shot_w_smoothz, 3*href(3))
                    ws = sph_gauss_filt(ws, shot_nybeg, shot_nzbeg, [shot_w_smoothx, shot_w_smoothy, shot_w_smoothz])

                case ('maxbal')
                    if (maxval(abs(ws)) /= 0) then
                        ws = ws/maxval(abs(ws))
                    end if

                case ('rmsbal')
                    ! Normalize with shot image energy
                    if (mean(ws, 2) /= 0) then
                        ws = ws/mean(ws, 2)
                    end if

                case ('movingbal')
                    ! Moving balance; window in cells of the reference size
                    call readpar_float(file_parameter, 'shot_'//tidy(name)//'_movingbal_lon', shot_w_movingbalx, 6*href(1))
                    call readpar_float(file_parameter, 'shot_'//tidy(name)//'_movingbal_lat', shot_w_movingbaly, 6*href(2))
                    call readpar_float(file_parameter, 'shot_'//tidy(name)//'_movingbal_depth', shot_w_movingbalz, 6*href(3))
                    ws = balance_filt(ws, nint([0.5*shot_w_movingbalz/href(3), 0.5*shot_w_movingbaly/href(2), &
                        0.5*shot_w_movingbalx/href(1)]), 0.01)

                case ('medianfilt')
                    ! Median filtering; window in cells of the reference size
                    call readpar_float(file_parameter, 'shot_'//tidy(name)//'_medianfilt_lon', shot_w_medianfiltx, href(1))
                    call readpar_float(file_parameter, 'shot_'//tidy(name)//'_medianfilt_lat', shot_w_medianfilty, href(2))
                    call readpar_float(file_parameter, 'shot_'//tidy(name)//'_medianfilt_depth', shot_w_medianfiltz, href(3))
                    ws = median_filt(ws, nint([shot_w_medianfiltz/href(3), shot_w_medianfilty/href(2), &
                        shot_w_medianfiltx/href(1)]))

                case ('removenan')
                    ! Remove NaN
                    ws = return_normal(ws)

                case ('laplacefilt')
                    ! Laplace filter in grid-index space
                    ws = laplace_filt(ws)

                case ('wavenumfilt')
                    ! Wavenumber-domain filtering along longitude, with the
                    ! exact spacing of each line
                    call readpar_nfloat(file_parameter, 'shot_'//tidy(name)//'_wavenum_lon', wavenums, [-1.0])
                    call readpar_nfloat(file_parameter, 'shot_'//tidy(name)//'_wavenum_lon_amps', wamps, [-1.0])
                    if (wavenums(1) >= 0) then
                        ws = sph_fourier_filt_lon(ws, shot_nybeg, shot_nzbeg, wavenums, wamps)
                    end if
                    ! Wavenumber-domain filtering along depth; along latitude
                    ! the spacing is not even, and the spherical branch does
                    ! not support it (check_spherical_parameter_file)
                    call readpar_nfloat(file_parameter, 'shot_'//tidy(name)//'_wavenum_depth', wavenums, [-1.0])
                    call readpar_nfloat(file_parameter, 'shot_'//tidy(name)//'_wavenum_depth_amps', wamps, [-1.0])
                    if (wavenums(1) >= 0) then
                        ws = sph_fourier_filt_depth(ws, wavenums, wamps)
                    end if

                case ('mask')
                    ! Masking
                    call readpar_string(file_parameter, 'dir_shot_'//tidy(name)//'_mask', dir_mask, '')
                    if (dir_mask == '') then
                        call readpar_string(file_parameter, 'shot_'//tidy(name)//'_mask', file_mask, '')
                    else
                        file_mask = tidy(dir_mask)//'/shot_'//num2str(gmtr(ishot)%id)//'_mask.bin'
                    end if
                    if (file_mask /= '') then
                        call prepare_model_single_parameter(w_mask, 'mask', file_mask, update=.false.)
                    else
                        w_mask = ones(nz, ny, nx)
                    end if
                    call alloc_array(w_mask, [1, shot_nz, 1, shot_ny, 1, shot_nx], &
                        source=w_mask(shot_nzbeg:shot_nzend, shot_nybeg:shot_nyend, shot_nxbeg:shot_nxend))
                    ws = ws*w_mask

                case ('taper')
                    ! Tapering, exact distances from the sides of the shot range
                    call readpar_nfloat(file_parameter, 'shot_'//tidy(name)//'_taper_lon', shot_w_taperx, [0.0, 0.0])
                    call readpar_nfloat(file_parameter, 'shot_'//tidy(name)//'_taper_lat', shot_w_tapery, [0.0, 0.0])
                    call readpar_nfloat(file_parameter, 'shot_'//tidy(name)//'_taper_depth', shot_w_taperz, [0.0, 0.0])
                    if (size(shot_w_taperx) == 1) then
                        call alloc_array(shot_w_taperx, [1, 2], source=[shot_w_taperx(1), shot_w_taperx(1)])
                    end if
                    if (size(shot_w_tapery) == 1) then
                        call alloc_array(shot_w_tapery, [1, 2], source=[shot_w_tapery(1), shot_w_tapery(1)])
                    end if
                    if (size(shot_w_taperz) == 1) then
                        call alloc_array(shot_w_taperz, [1, 2], source=[shot_w_taperz(1), shot_w_taperz(1)])
                    end if
                    ws = sph_taper(ws, shot_nybeg, shot_nzbeg, [shot_w_taperx(1:2), shot_w_tapery(1:2), shot_w_taperz(1:2)])

                case ('conemute')
                    ! Mute outside a cone below the source (or receivers): full
                    ! weight within the radius, a linear ramp to zero over the
                    ! taper length, and zero above the depth of the origin. The
                    ! horizontal distance of a node is the great-circle angle
                    ! between its column and the origin, times the node's radius.
                    call readpar_string(file_parameter, 'shot_'//tidy(name)//'_conemuteorigin', conemuteorigin, 'source')
                    call readpar_float(file_parameter, 'shot_'//tidy(name)//'_conemuter', conemuter, -1.0)
                    call readpar_float(file_parameter, 'shot_'//tidy(name)//'_conemute_depth', conemutez, -1.0)
                    call readpar_float(file_parameter, 'shot_'//tidy(name)//'_conemutepower', conemutepower, 2.0)
                    call readpar_float(file_parameter, 'shot_'//tidy(name)//'_conemutetaper', conemutetaper, 10*mean(href))
                    if (conemuter >= 0 .and. conemutez >= 0) then

                        conemutepower = max(1.0, conemutepower)

                        ! Origin: the direction of the Earth-centred centroid of
                        ! the points, and the depth of reference
                        xc = 0.0d0
                        select case (conemuteorigin)
                            case ('source')
                                do j = 1, gmtr(ishot)%ns
                                    xc = xc + lonlatdepth_to_ecef(dble(gmtr(ishot)%srcr(j)%x), &
                                        dble(gmtr(ishot)%srcr(j)%y), dble(gmtr(ishot)%srcr(j)%z))
                                end do
                                srcz = mean(gmtr(ishot)%srcr(:)%z)
                            case ('receiver')
                                do j = 1, gmtr(ishot)%nr
                                    xc = xc + lonlatdepth_to_ecef(dble(gmtr(ishot)%recr(j)%x), &
                                        dble(gmtr(ishot)%recr(j)%y), dble(gmtr(ishot)%recr(j)%z))
                                end do
                                srcz = minval(gmtr(ishot)%recr(:)%z)
                        end select
                        uc = xc/norm2(xc)

                        do k = 1, shot_nz

                            kg = shot_nzbeg + k - 1
                            depth = shot_zbeg + (k - 1)*dble(dz)

                            if (depth < srcz) then
                                ! When the depth is shallower than source then set zero
                                ws(k, :, :) = 0.0
                            else
                                ! The radius of the full weight at this depth
                                radius = conemuter*(min(depth - srcz, dble(conemutez))/conemutez)**(1.0d0/conemutepower)
                                !$omp parallel do private(ix, iy, jg, un, dist, tp)
                                do ix = 1, shot_nx
                                    do iy = 1, shot_ny
                                        jg = shot_nybeg + iy - 1
                                        un = [sph_coslat(jg)*sph_coslon(shot_nxbeg + ix - 1), &
                                            sph_coslat(jg)*sph_sinlon(shot_nxbeg + ix - 1), sph_sinlat(jg)]
                                        dist = sph_r(kg)*atan2(norm2(cross_product(un, uc)), dot_product(un, uc))
                                        if (dist <= radius) then
                                            tp = 1.0d0
                                        else if (dist >= radius + conemutetaper .or. conemutetaper <= 0) then
                                            tp = 0.0d0
                                        else
                                            tp = 1.0d0 - (dist - radius)/conemutetaper
                                        end if
                                        ws(k, iy, ix) = real(ws(k, iy, ix)*tp)
                                    end do
                                end do
                                !$omp end parallel do
                            end if

                        end do

                    end if

            end select

            if (process_shot_grad_w(i) /= '') then
                call warn(' Model value range = '//num2str(minval(ws), '(es)')//', '//num2str(maxval(ws), '(es)'))
                call warn(date_time_compact()//' shot '//num2str(gmtr(ishot)%id)//' '//tidy(name) &
                    //' processing ('//tidy(process_shot_grad_w(i))//') finished. ')
            end if

        end do

        ! Taper the sides of the shot range that lie inside the model, so that
        ! the merged gradient has no sharp edges along the adaptive range;
        ! an axis without adaptive range is never tapered. Exact distances,
        ! lengths adp_taper_lon, adp_taper_lat, adp_taper_depth.
        tlen = [ &
            merge(adptaperx, 0.0, yn_adpx .and. shot_nxbeg > 1), &
            merge(adptaperx, 0.0, yn_adpx .and. shot_nxend < nx), &
            merge(adptapery, 0.0, yn_adpy .and. shot_nybeg > 1), &
            merge(adptapery, 0.0, yn_adpy .and. shot_nyend < ny), &
            merge(adptaperz, 0.0, yn_adpz .and. shot_nzbeg > 1), &
            merge(adptaperz, 0.0, yn_adpz .and. shot_nzend < nz)]
        if (any(tlen > 0)) then
            ws = sph_taper(ws, shot_nybeg, shot_nzbeg, tlen)
        end if

        ! Merge image
        w(shot_nzbeg:shot_nzend, shot_nybeg:shot_nyend, shot_nxbeg:shot_nxend) = w(shot_nzbeg:shot_nzend, shot_nybeg:shot_nyend, shot_nxbeg:shot_nxend) + ws

        call warn(date_time_compact()//' Shot '//num2str(gmtr(ishot)%id)//' '//tidy(name)//' is merged. ')

    contains

        pure function cross_product(a, b) result(c)
            double precision, dimension(1:3), intent(in) :: a, b
            double precision, dimension(1:3) :: c
            c = [a(2)*b(3) - a(3)*b(2), a(3)*b(1) - a(1)*b(3), a(1)*b(2) - a(2)*b(1)]
        end function cross_product

    end subroutine process_model_single_shot

    !
    !> Processing gradients for a single model parameter
    !
    !> Smoothing gradients can remove high-wavenumber noises
    !> generated during gradient calculations. This functiionality is
    !> not as fancy as anisotropic diffusion denoising, but
    !> is much more less computational expensive.
    !
    subroutine process_model_single_parameter(w, name, param_name)

        real, dimension(:, :, :), intent(inout) :: w
        character(len=*), intent(in) :: name
        character(len=*), intent(in), optional :: param_name

        integer :: i, j, k
        real, allocatable, dimension(:) :: w_taperx, w_tapery, w_taperz
        real :: w_movingbalx, w_movingbaly, w_movingbalz
        real :: w_medianfiltx, w_medianfilty, w_medianfiltz
        real :: w_smoothx, w_smoothy, w_smoothz
        type(andf_param) :: param
        real, allocatable, dimension(:, :, :) :: andfaux, andfcoh
        character(len=32), allocatable, dimension(:) :: process_grad_w
        real :: w_signed_power, w_scalar
        real, allocatable, dimension(:, :, :) :: wt, w_mask
        integer :: wrx, wry, wrz
        real :: w_rmsbalx, w_rmsbaly, w_rmsbalz
        integer, allocatable, dimension(:) :: update_iter
        character(len=1024) :: file_mask
        real, dimension(1:3) :: href

        ! Reference cell size: east-west, north-south, vertical
        href = real(reference_spacing())

        call readpar_nstring(file_parameter, 'process_'//tidy(name), process_grad_w, [''])

        if (present(param_name)) then
            call readpar_nint(file_parameter, tidy(param_name)//'_update_iter', update_iter, [1, niter_max])
            if (size(update_iter) == 1) then
                update_iter = [update_iter(1), niter_max]
            end if
            if (iter < update_iter(1) .or. iter > update_iter(2)) then
                w = 0
                return
            end if
        end if

        ! Process gradient
        do i = 1, size(process_grad_w)

            select case (process_grad_w(i))

                case ('scale')
                    call readpar_xfloat(file_parameter, tidy(name)//'_scale', w_scalar, 1.0, iter*1.0)
                    w = w*w_scalar

                case ('signed_power')
                    call readpar_xfloat(file_parameter, tidy(name)//'_signed_power', w_signed_power, 0.5, iter*1.0)
                    w = sign(1.0, w)*(abs(w))**w_signed_power

                case ('movingbal')
                    call readpar_xfloat(file_parameter, tidy(name)//'_movingbal_lon', w_movingbalx, 3*href(1), iter*1.0)
                    call readpar_xfloat(file_parameter, tidy(name)//'_movingbal_lat', w_movingbaly, 3*href(2), iter*1.0)
                    call readpar_xfloat(file_parameter, tidy(name)//'_movingbal_depth', w_movingbalz, 3*href(3), iter*1.0)
                    w = balance_filt(w, nint([0.5*w_movingbalz/href(3), 0.5*w_movingbaly/href(2), 0.5*w_movingbalx/href(1)]), 0.01)

                case ('taper')
                    call readpar_nfloat(file_parameter, tidy(name)//'_taper_lon', w_taperx, [0.0, 0.0])
                    call readpar_nfloat(file_parameter, tidy(name)//'_taper_lat', w_tapery, [0.0, 0.0])
                    call readpar_nfloat(file_parameter, tidy(name)//'_taper_depth', w_taperz, [0.0, 0.0])
                    if (size(w_taperx) == 1) then
                        call alloc_array(w_taperx, [1, 2], source=[w_taperx(1), w_taperx(1)])
                    end if
                    if (size(w_tapery) == 1) then
                        call alloc_array(w_tapery, [1, 2], source=[w_tapery(1), w_tapery(1)])
                    end if
                    if (size(w_taperz) == 1) then
                        call alloc_array(w_taperz, [1, 2], source=[w_taperz(1), w_taperz(1)])
                    end if
                    w = sph_taper(w, 1, 1, [w_taperx(1:2), w_tapery(1:2), w_taperz(1:2)])

                case ('smooth')
                    call readpar_xfloat(file_parameter, tidy(name)//'_smooth_lon', w_smoothx, 3*href(1), iter*1.0)
                    call readpar_xfloat(file_parameter, tidy(name)//'_smooth_lat', w_smoothy, 3*href(2), iter*1.0)
                    call readpar_xfloat(file_parameter, tidy(name)//'_smooth_depth', w_smoothz, 3*href(3), iter*1.0)
                    w = sph_gauss_filt(w, 1, 1, [w_smoothx, w_smoothy, w_smoothz])

                case ('andf')
                    call readpar_xfloat(file_parameter, tidy(name)//'_andf_smooth_lon', param%smooth3, 2*href(1), iter*1.0)
                    call readpar_xfloat(file_parameter, tidy(name)//'_andf_smooth_lat', param%smooth2, 2*href(2), iter*1.0)
                    call readpar_xfloat(file_parameter, tidy(name)//'_andf_smooth_depth', param%smooth1, 8*href(3), iter*1.0)
                    param%smooth3 = param%smooth3/href(1)
                    param%smooth2 = param%smooth2/href(2)
                    param%smooth1 = param%smooth1/href(3)
                    call readpar_xfloat(file_parameter, tidy(name)//'_andf_powerm', param%powerm, 1.0, iter*1.0)
                    call readpar_xint(file_parameter, tidy(name)//'_andf_t', param%niter, 5, iter*1.0)
                    call readpar_xfloat(file_parameter, tidy(name)//'_andf_sigma', param%sigma, 6*maxval(href), iter*1.0)
                    param%sigma = param%sigma/maxval(href)
                    call readpar_xfloat(file_parameter, tidy(name)//'_andf_alpha', param%lambda1, 1.0e-3, iter*1.0)
                    call readpar_xfloat(file_parameter, tidy(name)//'_andf_beta', param%lambda2, 1.0, iter*1.0)
                    call readpar_xfloat(file_parameter, tidy(name)//'_andf_gamma', param%lambda3, 1.0, iter*1.0)
                    call readpar_xstring(file_parameter, tidy(name)//'_andf_aux', file_andfaux, '', iter*1.0)
                    call readpar_xstring(file_parameter, tidy(name)//'_andf_coh', file_andfcoh, '', iter*1.0)
                    call readpar_int(file_parameter, tidy(name)//'_andf_rank_lon', rank3, 1)
                    call readpar_int(file_parameter, tidy(name)//'_andf_rank_lat', rank2, 1)
                    call readpar_int(file_parameter, tidy(name)//'_andf_rank_depth', rank1, 1)
                    if (file_andfaux == '' .and. file_andfcoh == '') then
                        w = andf_filt_mpi(w, param)
                    else if (file_andfaux /= '' .and. file_andfcoh == '') then
                        call prepare_model_single_parameter(andfaux, 'andfaux', file_andfaux, update=.false.)
                        w = andf_filt_mpi(w, param, aux=andfaux)
                    else if (file_andfaux == '' .and. file_andfcoh /= '') then
                        call prepare_model_single_parameter(andfcoh, 'andfcoh', file_andfcoh, update=.false.)
                        w = andf_filt_mpi(w, param, acoh=andfcoh)
                    else if (file_andfaux /= '' .and. file_andfcoh /= '') then
                        call prepare_model_single_parameter(andfaux, 'andfaux', file_andfaux, update=.false.)
                        call prepare_model_single_parameter(andfcoh, 'andfcoh', file_andfcoh, update=.false.)
                        w = andf_filt_mpi(w, param, aux=andfaux, acoh=andfcoh)
                    end if

                case ('medianfilt')
                    call readpar_xfloat(file_parameter, tidy(name)//'_medianfilt_lon', w_medianfiltx, href(1), iter*1.0)
                    call readpar_xfloat(file_parameter, tidy(name)//'_medianfilt_lat', w_medianfilty, href(2), iter*1.0)
                    call readpar_xfloat(file_parameter, tidy(name)//'_medianfilt_depth', w_medianfiltz, href(3), iter*1.0)
                    w = median_filt(w, nint([w_medianfiltz/href(3), w_medianfilty/href(2), w_medianfiltx/href(1)]))

                case ('maxbal')
                    if (maxval(abs(w)) /= 0) then
                        w = w/maxval(abs(w))
                    end if

                case ('rmsbal')
                    if (mean(w, 2) /= 0) then
                        w = w/mean(w, 2)
                    end if

                case ('rmsbal_lon')
                    call readpar_xfloat(file_parameter, tidy(name)//'_rmsbal_lon', w_rmsbalx, 1.0*href(1), iter*1.0)
                    wrx = nint(w_rmsbalx/href(1))
                    if (mod(wrx, 2) == 0) then
                        wrx = wrx - 1
                    end if
                    wt = w
                    call pad_array(wt, [0, 0, 0, 0, wrx, wrx])
                    !$omp parallel do private(j)
                    do j = 1, size(w, 3)
                        w(:, :, j) = w(:, :, j)/norm2(wt(:, :, j - (wrx - 1)/2:j + (wrx - 1)/2))
                    end do
                    !$omp end parallel do
                    w = return_normal(w)

                case ('rmsbal_lat')
                    call readpar_xfloat(file_parameter, tidy(name)//'_rmsbal_lat', w_rmsbaly, 1.0*href(2), iter*1.0)
                    wry = nint(w_rmsbaly/href(2))
                    if (mod(wry, 2) == 0) then
                        wry = wry - 1
                    end if
                    wt = w
                    call pad_array(wt, [0, 0, wry, wry, 0, 0])
                    !$omp parallel do private(j)
                    do j = 1, size(w, 2)
                        w(:, j, :) = w(:, j, :)/norm2(wt(:, j - (wry - 1)/2:j + (wry - 1)/2, :))
                    end do
                    !$omp end parallel do
                    w = return_normal(w)

                case ('rmsbal_lonlat')
                    call readpar_xfloat(file_parameter, tidy(name)//'_rmsbal_lon', w_rmsbalx, 1.0*href(1), iter*1.0)
                    call readpar_xfloat(file_parameter, tidy(name)//'_rmsbal_lat', w_rmsbaly, 1.0*href(2), iter*1.0)
                    wrx = nint(w_rmsbalx/href(1))
                    if (mod(wrx, 2) == 0) then
                        wrx = wrx - 1
                    end if
                    wry = nint(w_rmsbaly/href(2))
                    if (mod(wry, 2) == 0) then
                        wry = wry - 1
                    end if
                    wt = w
                    call pad_array(wt, [0, 0, wry, wry, wrx, wrx])
                    !$omp parallel do private(j, k)
                    do k = 1, size(w, 3)
                        do j = 1, size(w, 2)
                            w(:, j, k) = w(:, j, k)/norm2(wt(:, j - (wry - 1)/2:j + (wry - 1)/2, k - (wrx - 1)/2:k + (wrx - 1)/2))
                        end do
                    end do
                    !$omp end parallel do
                    w = return_normal(w)

                case ('rmsbal_depth')
                    call readpar_xfloat(file_parameter, tidy(name)//'_rmsbal_depth', w_rmsbalz, 1.0*href(3), iter*1.0)
                    wrz = nint(w_rmsbalz/href(3))
                    if (mod(wrz, 2) == 0) then
                        wrz = wrz - 1
                    end if
                    wt = w
                    call pad_array(wt, [wrz, wrz, 0, 0, 0, 0])
                    !$omp parallel do private(j)
                    do j = 1, size(w, 1)
                        w(j, :, :) = w(j, :, :)/norm2(wt(j - (wrz - 1)/2:j + (wrz - 1)/2, :, :))
                    end do
                    !$omp end parallel do
                    w = return_normal(w)

                case ('mask')
                    call readpar_xstring(file_parameter, tidy(name)//'_mask', file_mask, '', iter*1.0)
                    if (file_mask /= '') then
                        call prepare_model_single_parameter(w_mask, 'mask', file_mask, update=.false.)
                    else
                        w_mask = ones_like(w)
                    end if
                    w = w*w_mask

            end select

            if (rankid == 0) then
                if (process_grad_w(i) /= '') then
                    call warn(' Model value range = '//num2str(minval(w), '(es)')//', '//num2str(maxval(w), '(es)'))
                    call warn(date_time_compact()//' '//tidy(name)//' processing ('//tidy(process_grad_w(i))//') finished. ')
                end if
            end if

        end do

    end subroutine process_model_single_parameter

end module gradient
