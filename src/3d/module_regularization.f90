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

    implicit none

contains

    !
    !> Solve the Tikhonov regularization in a hard way:
    !>        Decompose the inversion using ADMM, and then solve
    !>             ||\nabla u||_2^2 = ||u-m^{k+1}||_2^2.
    !>        The first-order optimality condition of the above
    !>        equation is
    !>             \nabla^T \nabla u = u - m^{k+1},
    !>        which can be solved using the Gauss-Seidel method.
    !
    function tikhonov_denoise(m, tikhonov_lambda) result(u)

        real, dimension(:, :, :) :: m
        real :: tikhonov_lambda
        real, allocatable, dimension(:, :, :) :: u

        real, allocatable, dimension(:, :, :) :: pu
        integer :: i, j, k, iter
        real :: sumu
        integer :: n1, n2, n3, niter
        integer :: i1, i2, j1, j2, k1, k2
        integer :: ii1, ii2, jj1, jj2, kk1, kk2
        real :: mu

        u = m
        mu = 1.0
        niter = 200
        n1 = size(u, 1)
        n2 = size(u, 2)
        n3 = size(u, 3)

        ! TGpV iteration
        do iter = 1, niter

            pu = u

            ! minimization u
            do k = 1, n3
                do j = 1, n2
                    do i = 1, n1

                        if (i == 1) then
                            i1 = 0
                            ii1 = i
                        else
                            i1 = 1
                            ii1 = i-1
                        end if

                        if (i == n1) then
                            i2 = 0
                            ii2 = i
                        else
                            i2 = 1
                            ii2 = i+1
                        end if

                        if (j == 1) then
                            j1 = 0
                            jj1 = j
                        else
                            j1 = 1
                            jj1 = j-1
                        end if
                        if (j == n2) then
                            j2 = 0
                            jj2 = j
                        else
                            j2 = 1
                            jj2 = j+1
                        end if

                        if (k == 1) then
                            k1 = 0
                            kk1 = k
                        else
                            k1 = 1
                            kk1 = k-1
                        end if
                        if (k == n3) then
                            k2 = 0
                            kk2 = k
                        else
                            k2 = 1
                            kk2 = k+1
                        end if

                        sumu = tikhonov_lambda*( &
                            + i2*u(ii2, j, k) + i1*u(ii1, j, k) &
                            + j2*u(i, jj2, k) + j1*u(i, jj1, k) &
                            + k2*u(i, j, kk2) + k1*u(i, j, kk1)) &
                            + mu*m(i, j, k)

                        u(i, j, k) = sumu/(mu + (6.0 - (1-i1) - (1-i2) - (1-j1) - (1-j2) - (1-k1) - (1-k2))*tikhonov_lambda)

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

        mt = m

        call readpar_int(file_parameter, 'rank_lon', rank3, floor(nrank**0.33333333d0))
        call readpar_int(file_parameter, 'rank_lat', rank2, floor(nrank**0.33333333d0))
        call readpar_int(file_parameter, 'rank_depth', rank1, floor(nrank**0.33333333d0))

        do i = 1, size(model_regularization_method)

            select case (model_regularization_method(i))

                case ('Tikhonov', 'tikhonov')
                    call readpar_xfloat(file_parameter, 'reg_tikhonov_lambda', tikhonov_lambda, 10.0, iter*1.0)
                    mt = tikhonov_denoise(mt, tikhonov_lambda)

                    if (rankid == 0) then
                        call warn(date_time_compact()//' >>>>>>>>>> Tikhonov regularization finished. ')
                    end if

                case ('smooth')
                    call readpar_xfloat(file_parameter, 'reg_smooth_lon', reg_smoothx, -1.0, iter*1.0)
                    if (reg_smoothx < 0) then
                        call readpar_xfloat(file_parameter, 'reg_smooth_lon_'//tidy(name), reg_smoothx, 1.0*dx, iter*1.0)
                    end if
                    call readpar_xfloat(file_parameter, 'reg_smooth_lat', reg_smoothy, -1.0, iter*1.0)
                    if (reg_smoothy < 0) then
                        call readpar_xfloat(file_parameter, 'reg_smooth_lat_'//tidy(name), reg_smoothy, 1.0*dy, iter*1.0)
                    end if
                    call readpar_xfloat(file_parameter, 'reg_smooth_depth', reg_smoothz, -1.0, iter*1.0)
                    if (reg_smoothz < 0) then
                        call readpar_xfloat(file_parameter, 'reg_smooth_depth_'//tidy(name), reg_smoothz, 1.0*dz, iter*1.0)
                    end if
                    if (present(smooth_inverse) .and. smooth_inverse) then
                        mt = 1.0/gauss_filt(1.0/mt, [reg_smoothz/dz, reg_smoothy/dy, reg_smoothx/dx])
                    else
                        mt = gauss_filt(mt, [reg_smoothz/dz, reg_smoothy/dy, reg_smoothx/dx])
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

                        call readpar_xfloat(file_parameter, 'reg_similarity_smooth_lon', sim_smoothx, 5*dx, iter*1.0)
                        call readpar_xfloat(file_parameter, 'reg_similarity_smooth_lat', sim_smoothy, 5*dy, iter*1.0)
                        call readpar_xfloat(file_parameter, 'reg_similarity_smooth_depth', sim_smoothz, 5*dz, iter*1.0)
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
                        r = median_filt(r)
                        r = gauss_filt(r, [sim_smoothz, sim_smoothy, sim_smoothx]/[dz, dy, dx])

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
