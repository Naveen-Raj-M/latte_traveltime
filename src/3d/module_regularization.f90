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


module regularization

    use libflit
    use parameters
    use vars
    use utility
    use spherical_filter

    implicit none

contains

    !
    !> Solve the Tikhonov regularization in a hard way:
    !>        Decompose the inversion using ADMM, and then solve
    !>             min_u  lambda ||grad u||^2 + ||u - m||^2,
    !>        both norms integrated over the volume. The first-order
    !>        optimality condition is
    !>             -lambda div(grad u) + u = m,
    !>        which is solved using the Gauss-Seidel method.
    !>
    !> Spherical branch: the Laplacian is the exact finite-volume one of the
    !> spherical cells. For a node with cell volume V, its neighbours n, the
    !> area A_n of the face shared with each and the arc length h_n to it,
    !> one Gauss-Seidel update is
    !>
    !>     u = (V m + lambda sum_n c_n u_n) / (V + lambda sum_n c_n),  c_n = A_n/h_n,
    !>
    !> with no neighbour across the boundary of the grid. lambda is a squared
    !> length (model length unit squared). In a Cartesian grid of unit cells
    !> this is the update of the Cartesian LATTE.
    !
    function tikhonov_denoise(m, tikhonov_lambda) result(u)

        real, dimension(:, :, :) :: m
        real :: tikhonov_lambda
        real, allocatable, dimension(:, :, :) :: u

        real, allocatable, dimension(:, :, :) :: pu
        double precision, allocatable, dimension(:, :) :: vol, clon, csouth, cnorth, cup, cdown
        integer :: i, j, k, iter
        integer :: n1, n2, n3, niter
        double precision :: num, den, lam

        u = m
        niter = 200
        lam = dble(tikhonov_lambda)

        ! Model order: i along depth, j along latitude, k along longitude
        n1 = size(u, 1)
        n2 = size(u, 2)
        n3 = size(u, 3)

        ! Cell volumes and face coefficients A/h, which do not depend on longitude
        allocate (vol(n1, n2), clon(n1, n2), csouth(n1, n2), cnorth(n1, n2), cup(n1, n2), cdown(n1, n2))
        do j = 1, n2
            do i = 1, n1
                vol(i, j) = cell_volume(j, i)
                clon(i, j) = face_weight_lon(j, i)*vol(i, j)/arc_lon(j, i)
                csouth(i, j) = 0.0d0
                cnorth(i, j) = 0.0d0
                cup(i, j) = 0.0d0
                cdown(i, j) = 0.0d0
                if (j > 1) then
                    csouth(i, j) = face_weight_lat(j, j, i)*vol(i, j)/arc_lat(j - 1, j, i)
                end if
                if (j < n2) then
                    cnorth(i, j) = face_weight_lat(j + 1, j, i)*vol(i, j)/arc_lat(j, j + 1, i)
                end if
                if (i > 1) then
                    cup(i, j) = face_weight_radial(i, j, i)*vol(i, j)/sph_ddepth
                end if
                if (i < n1) then
                    cdown(i, j) = face_weight_radial(i + 1, j, i)*vol(i, j)/sph_ddepth
                end if
            end do
        end do

        ! Gauss-Seidel iterations
        do iter = 1, niter

            pu = u

            do k = 1, n3
                do j = 1, n2
                    do i = 1, n1
                        num = vol(i, j)*m(i, j, k)
                        den = vol(i, j)
                        if (k > 1) then
                            num = num + lam*clon(i, j)*u(i, j, k - 1)
                            den = den + lam*clon(i, j)
                        end if
                        if (k < n3) then
                            num = num + lam*clon(i, j)*u(i, j, k + 1)
                            den = den + lam*clon(i, j)
                        end if
                        if (j > 1) then
                            num = num + lam*csouth(i, j)*u(i, j - 1, k)
                            den = den + lam*csouth(i, j)
                        end if
                        if (j < n2) then
                            num = num + lam*cnorth(i, j)*u(i, j + 1, k)
                            den = den + lam*cnorth(i, j)
                        end if
                        if (i > 1) then
                            num = num + lam*cup(i, j)*u(i - 1, j, k)
                            den = den + lam*cup(i, j)
                        end if
                        if (i < n1) then
                            num = num + lam*cdown(i, j)*u(i + 1, j, k)
                            den = den + lam*cdown(i, j)
                        end if
                        u(i, j, k) = real(num/den)
                    end do
                end do
            end do

            ! progress
            if (rankid == 0 .and. (mod(iter, max(nint(niter/10.0), 1)) == 0 .or. iter == 1)) then
                call warn(date_time_compact()//'>> Tikhonov iteration '//tidy(num2str(iter, '(i)')) &
                    //' of '//tidy(num2str(niter, '(i)')// &
                    ' relative norm2 diff = '//tidy(num2str(norm2(u - pu), '(es)'))))
            end if

        end do

    end function tikhonov_denoise

    !
    !> Apply regularization to a single model parameter
    !
    subroutine regularize_single_parameter(m, mr, name, default_value, smooth_inverse)

        real, dimension(:, :, :), intent(in) :: m
        real, dimension(:, :, :), intent(inout) :: mr
        character(len=*), intent(in) :: name
        real, intent(in) :: default_value
        logical, intent(in), optional :: smooth_inverse

        integer :: i
        real, allocatable, dimension(:, :, :) :: mt
        real :: reg_smoothx, reg_smoothy, reg_smoothz
        real :: tv_mu
        real :: tv_lambda1
        real :: tv_lambda2
        real :: tv_norm
        integer :: tv_niter
        real :: tikhonov_lambda
        real :: sim_smoothx, sim_smoothy, sim_smoothz
        real, allocatable, dimension(:, :, :) :: vp, vs, r
        real :: rmin, rmax
        real, dimension(1:3) :: href

        ! Reference cell size: east-west, north-south, vertical
        href = real(reference_spacing())

        mt = m

        call readpar_int(file_parameter, 'rank_lon', rank3, floor(nrank**0.33333333d0))
        call readpar_int(file_parameter, 'rank_lat', rank2, floor(nrank**0.33333333d0))
        call readpar_int(file_parameter, 'rank_depth', rank1, floor(nrank**0.33333333d0))

        do i = 1, size(model_regularization_method)

            select case (model_regularization_method(i))

                case ('Tikhonov', 'tikhonov')
                    ! A squared length; the default is 10 squared cells of the
                    ! mean size of the reference cell
                    call readpar_xfloat(file_parameter, 'reg_tikhonov_lambda', tikhonov_lambda, &
                        10.0*(sum(href)/3.0)**2, iter*1.0)
                    if (rankid == 0) then
                        call warn(date_time_compact()//' Tikhonov lambda = '//num2str(tikhonov_lambda, '(es)'))
                    end if
                    mt = tikhonov_denoise(mt, tikhonov_lambda)

                    if (rankid == 0) then
                        call warn(date_time_compact()//' >>>>>>>>>> Tikhonov regularization finished. ')
                    end if

                case ('smooth')
                    call readpar_xfloat(file_parameter, 'reg_smooth_lon', reg_smoothx, -1.0, iter*1.0)
                    if (reg_smoothx < 0) then
                        call readpar_xfloat(file_parameter, 'reg_smooth_lon_'//tidy(name), reg_smoothx, 1.0*href(1), iter*1.0)
                    end if
                    call readpar_xfloat(file_parameter, 'reg_smooth_lat', reg_smoothy, -1.0, iter*1.0)
                    if (reg_smoothy < 0) then
                        call readpar_xfloat(file_parameter, 'reg_smooth_lat_'//tidy(name), reg_smoothy, 1.0*href(2), iter*1.0)
                    end if
                    call readpar_xfloat(file_parameter, 'reg_smooth_depth', reg_smoothz, -1.0, iter*1.0)
                    if (reg_smoothz < 0) then
                        call readpar_xfloat(file_parameter, 'reg_smooth_depth_'//tidy(name), reg_smoothz, 1.0*href(3), iter*1.0)
                    end if
                    if (present(smooth_inverse) .and. smooth_inverse) then
                        mt = 1.0/sph_gauss_filt(1.0/mt, 1, 1, [reg_smoothx, reg_smoothy, reg_smoothz])
                    else
                        mt = sph_gauss_filt(mt, 1, 1, [reg_smoothx, reg_smoothy, reg_smoothz])
                    end if
                    if (rankid == 0) then
                        call warn(date_time_compact()//' >>>>>>>>>> Smoothing regularization finished. ')
                    end if

                case ('TGpV', 'tgpv')
                    call readpar_xfloat(file_parameter, 'reg_tv_mu_'//tidy(name), tv_mu, default_value, iter*1.0)
                    call readpar_xfloat(file_parameter, 'reg_tv_lambda1', tv_lambda1, 1.0, iter*1.0)
                    call readpar_xfloat(file_parameter, 'reg_tv_lambda2', tv_lambda2, 1.0, iter*1.0)
                    call readpar_xfloat(file_parameter, 'reg_tv_norm', tv_norm, 0.5, iter*1.0)
                    call readpar_xint(file_parameter, 'reg_tv_niter', tv_niter, 50, iter*1.0)
                    mt = tgpv_filt_mpi(mt, tv_mu, tv_lambda1, tv_lambda2, tv_niter, tv_norm)

                    if (rankid == 0) then
                        call warn(date_time_compact()//' >>>>>>>>>> TGpV regularization finished. ')
                    end if

                case ('vpvs_similarity')
                    if (name == 'vs') then

                        call readpar_xfloat(file_parameter, 'reg_similarity_smooth_lon', sim_smoothx, 5*href(1), iter*1.0)
                        call readpar_xfloat(file_parameter, 'reg_similarity_smooth_lat', sim_smoothy, 5*href(2), iter*1.0)
                        call readpar_xfloat(file_parameter, 'reg_similarity_smooth_depth', sim_smoothz, 5*href(3), iter*1.0)
                        call readpar_xfloat(file_parameter, 'reg_similarity_vpvs_ratio_min', rmin, min_vpvsratio, iter*1.0)
                        call readpar_xfloat(file_parameter, 'reg_similarity_vpvs_ratio_max', rmax, max_vpvsratio, iter*1.0)

                        if (any('vp' == model_name)) then
                            vp = get_meta_array_core(model_m, 'vp')
                        else if (any('vp' == model_name_aux)) then
                            vp = get_meta_array_core(model_aux, 'vp')
                        end if
                        vs = mt

                        mt = ones_like(vs)
                        where (vs == 0)
                            mt = 0
                        end where

                        r = zeros_like(vp)
                        where (mt /= 0)
                            r = vp/vs
                        end where
                        r = clip(r, rmin, rmax)
                        ! Median over the 3 x 3 x 3 cells around each node, then
                        ! the exact Gaussian
                        r = median_filt(r)
                        r = sph_gauss_filt(r, 1, 1, [sim_smoothx, sim_smoothy, sim_smoothz])

                        where (mt /= 0)
                            vs = vp/r
                        end where
                        mt = vs

                    end if

            end select

        end do

        call mpibarrier
        mr = mt

        ! Output regularized model
        if (rankid == 0) then
            call output_array(mr, dir_iter_model(iter)//'/reg_'//tidy(name)//'.bin')
        end if

    end subroutine regularize_single_parameter

    !
    !> Apply regularization to model parameters
    !
    subroutine model_regularization

        integer :: i

        do i = 1, nmodel

            if (.not. any(model_m(i)%name == [character(len=6) :: 'slon', 'slat', 'sdepth', 'st0'])) then

                if (model_m(i)%name == 'vp' .or. model_m(i)%name == 'vs') then
                    call regularize_single_parameter(model_m(i)%array, &
                        model_reg(i)%array, model_name(i), 100.0/(maxval(abs(model_m(i)%array)) + float_tiny), &
                        smooth_inverse=.true.)
                else
                    call regularize_single_parameter(model_m(i)%array, &
                        model_reg(i)%array, model_name(i), 100.0/(maxval(abs(model_m(i)%array)) + float_tiny))
                end if

            end if

        end do

    end subroutine model_regularization

end module regularization
