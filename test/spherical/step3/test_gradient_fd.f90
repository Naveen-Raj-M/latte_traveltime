!
! Finite-difference check of the adjoint-state gradient in spherical
! coordinates (step 3).
!
! Functional: J(v) = sum over receivers of T_r(v), the sum of the receiver
! traveltimes of one source. Its adjoint-state gradient is the FATT gradient
! of x_fatt3 with every residual equal to 1:
!
!     g = -lambda * V / v^3,   lambda from adjoint_iso with residuals 1,
!
! where V is the cell volume. For a smooth velocity perturbation dv, the
! first-order change of J is sum over nodes of g*dv. The test compares it
! with the central finite difference
!
!     (J(v + eps dv) - J(v - eps dv)) / (2 eps)
!
! for Gaussian perturbations at different latitudes and depths. With the
! cell volume, the ratio of the two is the same for every perturbation, up
! to discretization error. Without it (V replaced by a constant), the ratio
! changes with latitude and depth by the variation of the cell volume.
!
! Perturbations 5 and 6 overlap the source (at 300 km depth) and the
! surface receivers. The adjoint field is singular at both, and the
! adjoint-state gradient is less accurate there, as in the Cartesian code of
! the main branch.
!
! Grid: longitude 0 to 20 degrees, latitude 10 to 60 degrees, depth 0 to
! 600 km. Velocity 6000 + 0.01 d m/s at depth d. Source at (10, 35, 300 km);
! receivers on the surface every 0.5 degree.
!
program test_gradient_fd

    use libflit
    use parameters
    use spherical
    use traveltime_iso

    implicit none

    integer, parameter :: nlon = 81, nlat = 201, ndep = 61
    double precision, parameter :: olon = 0.0d0, olat = 10.0d0, odep = 0.0d0
    double precision, parameter :: dlon = 0.25d0, dlat = 0.25d0, ddep = 10000.0d0
    real, allocatable, dimension(:, :, :) :: v, vpert, t, lam, g, gnov
    real, allocatable, dimension(:, :) :: trec, ones_r
    type(source_receiver_geometry) :: geom
    integer :: i, j, k, l, nr, ib
    double precision :: eps, jp, jm, fd, pred, prednov, lon, lat, dep, d2
    double precision, dimension(1:7, 1:4) :: blobs
    double precision, dimension(1:3) :: xb, xn

    ! Tight sweep convergence for the finite differences
    sweep_niter_max = 50
    sweep_stop_threshold = 1.0e-9

    call init_spherical_grid(nlon, nlat, ndep, olon, olat, odep, dlon, dlat, ddep)

    ! Model, in the model order (depth, latitude, longitude)
    allocate (v(ndep, nlat, nlon))
    do k = 1, ndep
        v(k, :, :) = real(6000.0d0 + 0.01d0*(odep + (k - 1)*ddep))
    end do

    ! Geometry: one source, surface receivers every 0.5 degree
    geom%id = 1
    geom%ns = 1
    allocate (geom%srcr(1))
    geom%srcr(1)%x = 10.0
    geom%srcr(1)%y = 35.0
    geom%srcr(1)%z = 300000.0
    nr = 0
    geom%nr = 39*99
    allocate (geom%recr(geom%nr))
    do i = 1, 39
        do j = 1, 99
            nr = nr + 1
            geom%recr(nr)%x = real(0.5d0*i)
            geom%recr(nr)%y = real(10.0d0 + 0.5d0*j)
            geom%recr(nr)%z = 0.0
        end do
    end do
    ones_r = ones(geom%nr, 1)

    ! Adjoint-state gradient with every residual equal to 1
    call forward_iso(v, real([dlon, dlat, ddep]), real([olon, olat, odep]), geom, t, trec)
    call adjoint_iso(v, real([dlon, dlat, ddep]), real([olon, olat, odep]), geom, t, ones_r, lam)
    allocate (g(ndep, nlat, nlon), gnov(ndep, nlat, nlon))
    do j = 1, nlat
        do k = 1, ndep
            g(k, j, :) = real(-lam(k, j, :)*cell_volume(j, k)/v(k, j, :)**3)
            gnov(k, j, :) = real(-lam(k, j, :)*sph_vmean/v(k, j, :)**3)
        end do
    end do

    ! Gaussian perturbations: lon, lat, depth (m), radius (m)
    blobs(1, :) = [10.0d0, 15.0d0, 80000.0d0, 150000.0d0]
    blobs(2, :) = [10.0d0, 25.0d0, 80000.0d0, 150000.0d0]
    blobs(3, :) = [10.0d0, 45.0d0, 80000.0d0, 150000.0d0]
    blobs(4, :) = [10.0d0, 55.0d0, 80000.0d0, 150000.0d0]
    blobs(5, :) = [10.0d0, 35.0d0, 150000.0d0, 150000.0d0]
    blobs(6, :) = [12.0d0, 30.0d0, 30000.0d0, 150000.0d0]
    blobs(7, :) = [10.0d0, 29.0d0, 220000.0d0, 150000.0d0]

    eps = 1.0d0
    print '(a)', ' lon    lat    depth(km)   finite difference   adjoint with V    ratio    adjoint without V   ratio'
    do ib = 1, 7
        xb = lonlatdepth_to_ecef(blobs(ib, 1), blobs(ib, 2), blobs(ib, 3))
        allocate (vpert(ndep, nlat, nlon))
        do i = 1, nlon
            do j = 1, nlat
                do k = 1, ndep
                    xn = node_ecef(i, j, k)
                    d2 = sum((xn - xb)**2)
                    vpert(k, j, i) = real(100.0d0*exp(-d2/blobs(ib, 4)**2))
                end do
            end do
        end do
        call forward_iso(v + real(eps)*vpert, real([dlon, dlat, ddep]), real([olon, olat, odep]), geom, t, trec)
        jp = sum(dble(trec))
        call forward_iso(v - real(eps)*vpert, real([dlon, dlat, ddep]), real([olon, olat, odep]), geom, t, trec)
        jm = sum(dble(trec))
        fd = (jp - jm)/(2.0d0*eps)
        pred = sum(dble(g)*dble(vpert))
        prednov = sum(dble(gnov)*dble(vpert))
        print '(f5.1, f7.1, f10.0, 4x, es14.6, 4x, es14.6, f9.4, 4x, es14.6, f9.4)', &
            blobs(ib, 1), blobs(ib, 2), blobs(ib, 3)/1000.0d0, fd, pred, pred/fd, prednov, prednov/fd
        deallocate (vpert)
    end do

end program test_gradient_fd
