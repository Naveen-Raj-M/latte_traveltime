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

module traveltime_iso

    ! To ensure accuracy, computations in the internal subroutines/functions of this module
    ! are double-precision, but the input/output variables are single-precision.

    use libflit
    use parameters
    use spherical
    use utility
    use omp_lib

    implicit none

    real, parameter :: huge_value = sqrt(float_huge)

    double precision, allocatable, dimension(:, :, :) :: t0, pdxt0, pdyt0, pdzt0, tt, lambda, adjsrc
    ! Offsets of the local (per-shot) grid within the global grid: local node
    ! (i, j, k) is global node (i + ioff, j + joff, k + koff); the spherical
    ! geometry arrays of module spherical are indexed by global nodes
    integer :: ioff = 0, joff = 0, koff = 0

    private
    public :: forward_iso
    public :: adjoint_iso

contains


    !
    ! Fast sweeping local solver
    !
    function local_solver_2d(t0c, t0x, t0z, pdxt0c, pdzt0c, taux, tauz, dx, dz, signx, signz, v) result(txz)

        double precision :: t0c, t0x, t0z, pdxt0c, pdzt0c, taux, tauz
        double precision :: dx, dz, v
        integer :: signx, signz
        double precision :: txz

        double precision :: a1, a2, b1, b2
        double precision :: a, b, c
        double precision :: taucx, taucz
        double precision :: tau1, tau2
        logical :: causality1, causality2

        a1 = pdxt0c + t0c/dx*signx
        a2 = pdzt0c + t0c/dz*signz
        b1 = -taux*t0c/dx*signx
        b2 = -tauz*t0c/dz*signz

        a = a1**2 + a2**2
        b = 2*(a1*b1 + a2*b2)
        c = b1**2 + b2**2 - 1.0/v**2

        if (b**2 - 4*a*c > 0) then

            tau1 = (-b + sqrt(b**2 - 4*a*c))/(2*a)
            tau2 = (-b - sqrt(b**2 - 4*a*c))/(2*a)

            causality1 = tau1*t0c >= taux*t0x .and. tau1*t0c >= tauz*t0z
            causality2 = tau2*t0c >= taux*t0x .and. tau2*t0c >= tauz*t0z

            if (causality1 .and. causality2) then
                txz = min(tau1, tau2)
            else if (causality1 .and. .not. causality2) then
                txz = tau1
            else if (.not. causality1 .and. causality2) then
                txz = tau2
            else
                taucx = max((t0c*taux + dx/v)/(t0c + pdxt0c*dx*signx), taux*t0x/t0c)
                taucz = max((t0c*tauz + dz/v)/(t0c + pdzt0c*dz*signz), tauz*t0z/t0c)
                txz = min(taucx, taucz)
            end if

        else

            taucx = max((t0c*taux + dx/v)/(t0c + pdxt0c*dx*signx), taux*t0x/t0c)
            taucz = max((t0c*tauz + dz/v)/(t0c + pdzt0c*dz*signz), tauz*t0z/t0c)
            txz = min(taucx, taucz)

        end if

    end function local_solver_2d

    function local_solver_3d(t0c, t0x, t0y, t0z, pdxt0c, pdyt0c, pdzt0c, taux, tauy, tauz, dx, dy, dz, signx, signy, signz, v) result(txyz)

        double precision :: t0c, t0x, t0y, t0z, pdxt0c, pdyt0c, pdzt0c, taux, tauy, tauz
        double precision :: dx, dy, dz, v
        integer :: signx, signy, signz
        double precision :: txyz

        double precision :: a1, a2, a3, b1, b2, b3
        double precision :: a, b, c
        double precision :: tau1, tau2
        logical :: causality1, causality2

        a1 = pdxt0c + t0c/dx*signx
        a2 = pdyt0c + t0c/dy*signy
        a3 = pdzt0c + t0c/dz*signz
        b1 = -taux*t0c/dx*signx
        b2 = -tauy*t0c/dy*signy
        b3 = -tauz*t0c/dz*signz

        a = a1**2 + a2**2 + a3**2
        b = 2*(a1*b1 + a2*b2 + a3*b3)
        c = b1**2 + b2**2 + b3**2 - 1.0/v**2

        if (b**2 - 4*a*c > 0) then

            tau1 = (-b + sqrt(b**2 - 4*a*c))/(2*a)
            tau2 = (-b - sqrt(b**2 - 4*a*c))/(2*a)

            causality1 = tau1*t0c >= max(taux*t0x, tauy*t0y, tauz*t0z)
            causality2 = tau2*t0c >= max(taux*t0x, tauy*t0y, tauz*t0z)

            if (causality1 .and. causality2) then
                txyz = min(tau1, tau2)
            else if (causality1 .and. .not. causality2) then
                txyz = tau1
            else if (.not. causality1 .and. causality2) then
                txyz = tau2
            else
                txyz = min( &
                    local_solver_2d(t0c, t0x, t0y, pdxt0c, pdyt0c, taux, tauy, dx, dy, signx, signy, v), &
                    local_solver_2d(t0c, t0x, t0z, pdxt0c, pdzt0c, taux, tauz, dx, dz, signx, signz, v), &
                    local_solver_2d(t0c, t0y, t0z, pdyt0c, pdzt0c, tauy, tauz, dy, dz, signy, signz, v))
            end if

        else

            txyz = min( &
                local_solver_2d(t0c, t0x, t0y, pdxt0c, pdyt0c, taux, tauy, dx, dy, signx, signy, v), &
                local_solver_2d(t0c, t0x, t0z, pdxt0c, pdzt0c, taux, tauz, dx, dz, signx, signz, v), &
                local_solver_2d(t0c, t0y, t0z, pdyt0c, pdzt0c, tauy, tauz, dy, dz, signy, signz, v))

        end if

    end function local_solver_3d

    !
    !> Parallel fast sweeping
    !>
    !> The global scheme is based on
    !>      Detrixhe et al., 2013, JCP,
    !>      A parallel fast sweeping method for the eikonal equation
    !>      doi: 10.1016/j.jcp.2012.11.042
    !>
    !> The local scheme is based on
    !>      Fomel et al., 2009, JCP,
    !>      Fast sweeping method for the factored eikonal equation
    !>      doi: 10.1016/j.jcp.2009.05.029
    !>
    !> I made modifications to these schemes for efficiency
    !>
    !> Spherical coordinates: the axes are longitude (i, eastward), geodetic
    !> latitude (j, northward) and depth (k, downward). The eikonal equation
    !> in the local orthonormal frame of east, north and down is
    !>
    !>     (dT/de)^2 + (dT/dn)^2 + (dT/dd)^2 = 1/v^2,
    !>     dT/de = 1/(r cos(theta)) dT/dlambda, dT/dn = 1/r dT/dtheta, dT/dd = dT/dd,
    !>
    !> and the one-sided differences of the local solver use the arc length to
    !> the upwind neighbour: h_lon = r_k cos(theta_j) dlambda toward i - 1 or
    !> i + 1, r_k (theta_j - theta_j-1) toward j - 1, r_k (theta_j+1 - theta_j)
    !> toward j + 1, and ddepth toward k - 1 or k + 1. px0, py0 and pz0 are the
    !> east, north and down components of grad(T0).
    !
    subroutine fast_sweep_forward(nx, ny, nz, nxd, nyd, nzd, vp, t0, px0, py0, pz0, tau)

        integer, intent(in) :: nx, ny, nz, nxd, nyd, nzd
        double precision, dimension(:, :, :), intent(in) :: vp, t0, px0, py0, pz0
        double precision, dimension(:, :, :), intent(inout) :: tau

        integer :: nxbeg, nxend, nybeg, nyend, nzbeg, nzend
        integer :: i, j, k
        integer :: level, jbeg, jend, kbeg, kend
        double precision :: txyz1, txyz2, txyz3, txyz4, txyz5, txyz6, txyz7, txyz8
        logical :: valid1, valid2, valid3, valid4, valid5, valid6, valid7, valid8
        double precision :: hx, hym, hyp, hz

        if (nxd > 0) then
            nxbeg = 1
            nxend = nx
        else
            nxbeg = nx
            nxend = 1
        end if

        if (nyd > 0) then
            nybeg = 1
            nyend = ny
        else
            nybeg = ny
            nyend = 1
        end if

        if (nzd > 0) then
            nzbeg = 1
            nzend = nz
        else
            nzbeg = nz
            nzend = 1
        end if

        ! In 3D, the maximum level is nz + (ny - 1) + (nx - 1)
        do level = 1, nx + ny + nz - 2

            kbeg = ifelse(level <= nx + ny - 1, nzbeg, nzbeg + (level - (nx + ny - 1))*nzd)
            kend = ifelse(level <= nz, nzbeg + (level - 1)*nzd, nzend)

            jbeg = ifelse(level <= nx + nz - 1, nybeg, nybeg + (level - (nx + nz - 1))*nyd)
            jend = ifelse(level <= ny, nybeg + (level - 1)*nyd, nyend)

            !$omp parallel do private(i, j, k, txyz1, txyz2, txyz3, txyz4, txyz5, txyz6, txyz7, txyz8, &
                !$omp valid1, valid2, valid3, valid4, valid5, valid6, valid7, valid8, hx, hym, hyp, hz) schedule(dynamic)
            do k = kbeg, kend, nzd
                do j = jbeg, jend, nyd

                    ! The condition is
                    ! abs(i - nxbeg) + abs(j - nybeg) + abs(k - nzbeg) + 1 = level
                    i = ((level - 1) - abs(k - nzbeg) - abs(j - nybeg))/nxd + nxbeg

                    ! Different from 2D, in 3D, the above condition can satisfy with out-of-range i values.
                    ! This can potentially cause load unbalance, but is the best I can do at this moment.
                    if (i < 1 .or. i > nx) then
                        cycle
                    end if

                    txyz1 = huge_value
                    txyz2 = huge_value
                    txyz3 = huge_value
                    txyz4 = huge_value
                    txyz5 = huge_value
                    txyz6 = huge_value
                    txyz7 = huge_value
                    txyz8 = huge_value

                    valid1 = .true.
                    valid2 = .true.
                    valid3 = .true.
                    valid4 = .true.
                    valid5 = .true.
                    valid6 = .true.
                    valid7 = .true.
                    valid8 = .true.

                    if (i == 1) then
                        valid1 = .false.
                        valid4 = .false.
                        valid5 = .false.
                        valid8 = .false.
                    else if (i == nx) then
                        valid2 = .false.
                        valid3 = .false.
                        valid6 = .false.
                        valid7 = .false.
                    end if
                    if (j == 1) then
                        valid1 = .false.
                        valid2 = .false.
                        valid5 = .false.
                        valid6 = .false.
                    else if (j == ny) then
                        valid3 = .false.
                        valid4 = .false.
                        valid7 = .false.
                        valid8 = .false.
                    end if
                    if (k == 1) then
                        valid1 = .false.
                        valid2 = .false.
                        valid3 = .false.
                        valid4 = .false.
                    else if (k == nz) then
                        valid5 = .false.
                        valid6 = .false.
                        valid7 = .false.
                        valid8 = .false.
                    end if

                    ! Arc lengths to the neighbours of this node
                    hx = arc_lon(j + joff, k + koff)
                    hym = 0.0d0
                    hyp = 0.0d0
                    if (j > 1) then
                        hym = arc_lat(j + joff - 1, j + joff, k + koff)
                    end if
                    if (j < ny) then
                        hyp = arc_lat(j + joff, j + joff + 1, k + koff)
                    end if
                    hz = sph_ddepth

                    if (valid1) then
                        txyz1 = local_solver_3d(t0(i, j, k), t0(i - 1, j, k), t0(i, j - 1, k), t0(i, j, k - 1), &
                            px0(i, j, k), py0(i, j, k), pz0(i, j, k), tau(i - 1, j, k), tau(i, j - 1, k), tau(i, j, k - 1), &
                            hx, hym, hz, +1, +1, +1, vp(i, j, k))
                    end if

                    if (valid2) then
                        txyz2 = local_solver_3d(t0(i, j, k), t0(i + 1, j, k), t0(i, j - 1, k), t0(i, j, k - 1), &
                            px0(i, j, k), py0(i, j, k), pz0(i, j, k), tau(i + 1, j, k), tau(i, j - 1, k), tau(i, j, k - 1), &
                            hx, hym, hz, -1, +1, +1, vp(i, j, k))
                    end if

                    if (valid3) then
                        txyz3 = local_solver_3d(t0(i, j, k), t0(i + 1, j, k), t0(i, j + 1, k), t0(i, j, k - 1), &
                            px0(i, j, k), py0(i, j, k), pz0(i, j, k), tau(i + 1, j, k), tau(i, j + 1, k), tau(i, j, k - 1), &
                            hx, hyp, hz, -1, -1, +1, vp(i, j, k))
                    end if

                    if (valid4) then
                        txyz4 = local_solver_3d(t0(i, j, k), t0(i - 1, j, k), t0(i, j + 1, k), t0(i, j, k - 1), &
                            px0(i, j, k), py0(i, j, k), pz0(i, j, k), tau(i - 1, j, k), tau(i, j + 1, k), tau(i, j, k - 1), &
                            hx, hyp, hz, +1, -1, +1, vp(i, j, k))
                    end if

                    if (valid5) then
                        txyz5 = local_solver_3d(t0(i, j, k), t0(i - 1, j, k), t0(i, j - 1, k), t0(i, j, k + 1), &
                            px0(i, j, k), py0(i, j, k), pz0(i, j, k), tau(i - 1, j, k), tau(i, j - 1, k), tau(i, j, k + 1), &
                            hx, hym, hz, +1, +1, -1, vp(i, j, k))
                    end if

                    if (valid6) then
                        txyz6 = local_solver_3d(t0(i, j, k), t0(i + 1, j, k), t0(i, j - 1, k), t0(i, j, k + 1), &
                            px0(i, j, k), py0(i, j, k), pz0(i, j, k), tau(i + 1, j, k), tau(i, j - 1, k), tau(i, j, k + 1), &
                            hx, hym, hz, -1, +1, -1, vp(i, j, k))
                    end if

                    if (valid7) then
                        txyz7 = local_solver_3d(t0(i, j, k), t0(i + 1, j, k), t0(i, j + 1, k), t0(i, j, k + 1), &
                            px0(i, j, k), py0(i, j, k), pz0(i, j, k), tau(i + 1, j, k), tau(i, j + 1, k), tau(i, j, k + 1), &
                            hx, hyp, hz, -1, -1, -1, vp(i, j, k))
                    end if

                    if (valid8) then
                        txyz8 = local_solver_3d(t0(i, j, k), t0(i - 1, j, k), t0(i, j + 1, k), t0(i, j, k + 1), &
                            px0(i, j, k), py0(i, j, k), pz0(i, j, k), tau(i - 1, j, k), tau(i, j + 1, k), tau(i, j, k + 1), &
                            hx, hyp, hz, +1, -1, -1, vp(i, j, k))
                    end if

                    tau(i, j, k) = min(tau(i, j, k), txyz1, txyz2, txyz3, txyz4, txyz5, txyz6, txyz7, txyz8)

                end do
            end do
            !$omp end parallel do

        end do

    end subroutine fast_sweep_forward


    !
    !> Fast sweeping to solve the adjoint state equation
    !>
    !>     div(lambda grad T) = sum_r dt_r delta(x - x_r)
    !>
    !> in the conservative upwind form of the Cartesian LATTE (Leung and
    !> Qian, 2006; equation A-9 of Taillandier et al., 2009), written as a
    !> finite-volume balance of the cell of each node: the six face fluxes
    !> times the face areas over the cell volume. For each axis, with the
    !> traveltime gradients a+ = (T_+ - T_c)/h_+ and a- = (T_c - T_-)/h_- on
    !> the two faces, the face weights w+ and w- (face area over cell
    !> volume), pos(a) = (a + |a|)/2 and neg(a) = (a - |a|)/2,
    !>
    !>     lhs = sum over axes of ( w+ neg(a+) - w- pos(a-) )
    !>     rhs = sum over axes of ( w- neg(a-) lambda_- - w+ pos(a+) lambda_+ )
    !>     lambda_c = (rhs - s_c)/lhs
    !>
    !> where s_c is the point source density at the node. In Cartesian
    !> coordinates every w is 1/h, which is the original scheme. Here h is the
    !> arc length between the nodes and w the exact face area of the spherical
    !> cell over its exact volume (module spherical). A face on the boundary
    !> of the grid has a = 0, i.e., no flux, as in the Cartesian code.
    !
    subroutine fast_sweep_adjoint(n1, n2, n3, n1d, n2d, n3d, lambda)

        integer, intent(in) :: n1, n2, n3, n1d, n2d, n3d
        double precision, dimension(:, :, :), intent(inout) :: lambda

        integer :: i, j, k, jg, kg
        double precision :: app, amp, apm, amm
        double precision :: bpp, bmp, bpm, bmm
        double precision :: cpp, cmp, cpm, cmm
        double precision :: ap, am, bp, bm, cp, cm
        double precision :: wxp, wxm, wyp, wym, wzp, wzm
        double precision :: hx, hym, hyp, hz
        double precision :: lhs, rhs, t
        integer :: n1beg, n1end, n2beg, n2end, n3beg, n3end
        integer :: i1, i2, j1, j2, k1, k2
        integer :: level, kbeg, kend, jbeg, jend

        if (n1d > 0) then
            n1beg = 1
            n1end = n1
        else
            n1beg = n1
            n1end = 1
        end if

        if (n2d > 0) then
            n2beg = 1
            n2end = n2
        else
            n2beg = n2
            n2end = 1
        end if

        if (n3d > 0) then
            n3beg = 1
            n3end = n3
        else
            n3beg = n3
            n3end = 1
        end if

        hz = sph_ddepth

        ! In 3D, the maximum level is n3 + (n2 - 1) + (n1 - 1)
        do level = 1, n1 + n2 + n3 - 2

            kbeg = ifelse(level <= n1 + n2 - 1, n3beg, n3beg + (level - (n1 + n2 - 1))*n3d)
            kend = ifelse(level <= n3, n3beg + (level - 1)*n3d, n3end)

            jbeg = ifelse(level <= n1 + n3 - 1, n2beg, n2beg + (level - (n1 + n3 - 1))*n2d)
            jend = ifelse(level <= n2, n2beg + (level - 1)*n2d, n2end)

            !$omp parallel do private(i, j, k, jg, kg, i1, i2, j1, j2, k1, k2, &
                !$omp app, amp, apm, amm, &
                !$omp bpp, bmp, bpm, bmm, &
                !$omp cpp, cmp, cpm, cmm, &
                !$omp ap, am, bp, bm, cp, cm, &
                !$omp wxp, wxm, wyp, wym, wzp, wzm, hx, hym, hyp, &
                !$omp lhs, rhs, t) collapse(2)
            do k = kbeg, kend, n3d
                do j = jbeg, jend, n2d

                    ! The condition is
                    ! abs(i - n1beg) + abs(j - n2beg) + abs(k - n3beg) + 1 = level
                    i = ((level - 1) - abs(k - n3beg) - abs(j - n2beg))/n1d + n1beg

                    ! Different from 2D, in 3D, the above condition can satisfy with out-of-range i values.
                    ! This can potentially cause load unbalance, but is the best I can do at this moment.
                    if (i < 1 .or. i > n1) then
                        cycle
                    end if

                    if (i == 1) then
                        i1 = i
                    else
                        i1 = i - 1
                    end if
                    if (i == n1) then
                        i2 = i
                    else
                        i2 = i + 1
                    end if

                    if (j == 1) then
                        j1 = j
                    else
                        j1 = j - 1
                    end if
                    if (j == n2) then
                        j2 = j
                    else
                        j2 = j + 1
                    end if

                    if (k == 1) then
                        k1 = k
                    else
                        k1 = k - 1
                    end if
                    if (k == n3) then
                        k2 = k
                    else
                        k2 = k + 1
                    end if

                    ! Global indices of the node for the spherical geometry
                    jg = j + joff
                    kg = k + koff

                    ! Arc lengths between the node and its neighbours; on a
                    ! boundary face the difference below is zero, and the
                    ! arc length there only needs to be nonzero
                    hx = arc_lon(jg, kg)
                    hym = 1.0d0
                    hyp = 1.0d0
                    if (j > 1) then
                        hym = arc_lat(jg - 1, jg, kg)
                    end if
                    if (j < n2) then
                        hyp = arc_lat(jg, jg + 1, kg)
                    end if

                    ! Face areas over the cell volume
                    wxm = face_weight_lon(jg, kg)
                    wxp = wxm
                    wym = face_weight_lat(jg, jg, kg)
                    wyp = face_weight_lat(jg + 1, jg, kg)
                    wzm = face_weight_radial(kg, jg, kg)
                    wzp = face_weight_radial(kg + 1, jg, kg)

                    ! Traveltime gradients on the faces: east, north, down
                    ap = (tt(i2, j, k) - tt(i, j, k))/hx
                    am = (tt(i, j, k) - tt(i1, j, k))/hx

                    bp = (tt(i, j2, k) - tt(i, j, k))/hyp
                    bm = (tt(i, j, k) - tt(i, j1, k))/hym

                    cp = (tt(i, j, k2) - tt(i, j, k))/hz
                    cm = (tt(i, j, k) - tt(i, j, k1))/hz

                    app = (ap + abs(ap))/2.0
                    apm = (ap - abs(ap))/2.0

                    amp = (am + abs(am))/2.0
                    amm = (am - abs(am))/2.0

                    bpp = (bp + abs(bp))/2.0
                    bpm = (bp - abs(bp))/2.0

                    bmp = (bm + abs(bm))/2.0
                    bmm = (bm - abs(bm))/2.0

                    cpp = (cp + abs(cp))/2.0
                    cpm = (cp - abs(cp))/2.0

                    cmp = (cm + abs(cm))/2.0
                    cmm = (cm - abs(cm))/2.0

                    lhs = (wxp*apm - wxm*amp) + (wyp*bpm - wym*bmp) + (wzp*cpm - wzm*cmp)
                    rhs = (wxm*amm*lambda(i1, j, k) - wxp*app*lambda(i2, j, k)) &
                        + (wym*bmm*lambda(i, j1, k) - wyp*bpp*lambda(i, j2, k)) &
                        + (wzm*cmm*lambda(i, j, k1) - wzp*cpp*lambda(i, j, k2))

                    if (lhs == 0) then
                        t = 0
                    else
                        t = (rhs - adjsrc(i, j, k))/lhs
                    end if

                    lambda(i, j, k) = min(lambda(i, j, k), t)

                end do
            end do
        end do

    end subroutine fast_sweep_adjoint

    !
    !> Background time field of a source of one or many points: at every
    !> grid point the analytic time from the source point that arrives
    !> first, t0 = min_l (|x - s_l|/v_l + t0_l), and its derivatives, into
    !> the module arrays t0, pdxt0, pdyt0, pdzt0.
    !>
    !> Spherical coordinates: x and s_l are Earth-centred Cartesian positions,
    !> so |x - s_l| is the exact straight-line distance and t0 is the exact
    !> traveltime of a homogeneous sphere of velocity v_l. pdxt0, pdyt0 and
    !> pdzt0 are the components of grad(t0) = (x - s_l)/(v_l |x - s_l|) along
    !> the local east, north and down unit vectors of the grid node. The
    !> binning below works on the Earth-centred positions, and the distance
    !> from a node to a bin's bounding box is still a lower bound of its
    !> distance to any point of the bin, so the pruning stays exact.
    !
    !> The minimum used to re-interpolate every source point's velocity at
    !> every grid point and fill four ns-long arrays per grid point:
    !> nx*ny*nz*ns interpolations, with the source points taken one after
    !> another. A reflector ensemble (forward_iso_reflection) has tens of
    !> thousands of them. Now each source point's velocity is interpolated
    !> once, and an ensemble larger than background_direct is binned in
    !> space: a bin whose lower bound -- distance to its box over its
    !> fastest velocity, plus its earliest start time -- is above the best
    !> time found so far cannot hold the minimum, and is skipped. The
    !> pruning is exact: the same source point wins, by minloc's rule for
    !> ties (the lowest index), and its time and derivatives are the same
    !> expressions, so the field is identical to evaluating every point.
    !> The loop runs over grid rows in parallel; along a row, each grid
    !> point starts from its neighbour's winner, which is usually its own.
    !
    subroutine background_time(geom, source_inside, nx, ny, nz, ox, oy, oz, dx, dy, dz, vp)

        type(source_receiver_geometry), intent(in) :: geom
        logical, dimension(:), intent(in) :: source_inside
        integer, intent(in) :: nx, ny, nz
        double precision, intent(in) :: ox, oy, oz, dx, dy, dz
        double precision, dimension(:, :, :), intent(in) :: vp

        ! Up to this many source points, all are evaluated at every grid point
        integer, parameter :: background_direct = 16
        ! sx, sy, sz: Earth-centred Cartesian positions of the source points
        double precision, allocatable, dimension(:) :: sx, sy, sz, sv, st
        double precision, allocatable, dimension(:) :: bxl, bxh, byl, byh, bzl, bzh, bvmax, btmin
        integer, allocatable, dimension(:) :: member, bin_of, bin_count, bin_start, order, bfirst, blast
        integer :: ns, nb, nbin, nbx, nby, nbz, npos, i, j, k, l, m, b, h, win, prev, first
        double precision :: ex, ey, ez, cell, best, tval, lb, lbmin, tol
        double precision, dimension(1:3) :: xs

        t0 = zeros(nx, ny, nz) + huge_value
        pdxt0 = zeros(nx, ny, nz)
        pdyt0 = zeros(nx, ny, nz)
        pdzt0 = zeros(nx, ny, nz)

        ! The source points that take part, in index order, with what the
        ! minimum needs of each
        ns = count(source_inside)
        if (ns == 0) then
            return
        end if
        member = pack([(l, l = 1, geom%ns)], source_inside)
        allocate (sx(1:ns), sy(1:ns), sz(1:ns), sv(1:ns), st(1:ns))
        !$omp parallel do private(m, l, xs)
        do m = 1, ns
            l = member(m)
            xs = lonlatdepth_to_ecef(dble(geom%srcr(l)%x), dble(geom%srcr(l)%y), dble(geom%srcr(l)%z))
            sx(m) = xs(1)
            sy(m) = xs(2)
            sz(m) = xs(3)
            ! Linear interpolation to get velocity at the source position
            sv(m) = get_point_value_inside(dble([geom%srcr(l)%x, geom%srcr(l)%y, geom%srcr(l)%z]), &
                [nx, ny, nz], [ox, oy, oz], [dx, dy, dz], vp)
            st(m) = geom%srcr(l)%t0
            !                ! In comparison, the nearest grid point apporach
            !                isx = nint((geom%srcr(l)%x - ox)/dx) + 1
            !                isy = nint((geom%srcr(l)%y - oy)/dy) + 1
            !                isz = nint((geom%srcr(l)%z - oz)/dz) + 1
            !                vsource = vp(isx, isy, isz)
            !                dsx = (i - isx)*dx
            !                dsy = (j - isy)*dy
            !                dsz = (k - isz)*dz
        end do
        !$omp end parallel do

        if (ns <= background_direct) then

            !$omp parallel do collapse(2) private(i, j, k, m, best, win, tval)
            do k = 1, nz
                do j = 1, ny
                    do i = 1, nx
                        best = huge_value
                        win = 0
                        do m = 1, ns
                            tval = source_time(m, i, j, k)
                            if (tval < best) then
                                best = tval
                                win = m
                            end if
                        end do
                        call set_background(win, i, j, k)
                    end do
                end do
            end do
            !$omp end parallel do
            return

        end if

        ! Bins: a uniform grid over the box of the source points, about
        ! 2*sqrt(ns) cells -- per grid point, one bound per non-empty bin
        ! plus the points of the few bins the bound does not rule out
        ex = maxval(sx) - minval(sx)
        ey = maxval(sy) - minval(sy)
        ez = maxval(sz) - minval(sz)
        npos = count([ex, ey, ez] > 0)
        nbx = 1
        nby = 1
        nbz = 1
        if (npos > 0) then
            cell = (product(pack([ex, ey, ez], [ex, ey, ez] > 0))/(2.0d0*sqrt(dble(ns))))**(1.0d0/npos)
            if (ex > 0) nbx = max(1, min(ns, ceiling(ex/cell)))
            if (ey > 0) nby = max(1, min(ns, ceiling(ey/cell)))
            if (ez > 0) nbz = max(1, min(ns, ceiling(ez/cell)))
        end if
        nbin = nbx*nby*nbz
        allocate (bin_of(1:ns))
        do m = 1, ns
            bin_of(m) = bin_index(sx(m), minval(sx), ex, nbx) &
                + nbx*(bin_index(sy(m), minval(sy), ey, nby) - 1) &
                + nbx*nby*(bin_index(sz(m), minval(sz), ez, nbz) - 1)
        end do

        ! A counting sort by bin, stable, so each bin lists its points in
        ! index order
        allocate (bin_count(1:nbin), bin_start(1:nbin + 1))
        bin_count = 0
        do m = 1, ns
            bin_count(bin_of(m)) = bin_count(bin_of(m)) + 1
        end do
        bin_start(1) = 1
        do b = 1, nbin
            bin_start(b + 1) = bin_start(b) + bin_count(b)
        end do
        allocate (order(1:ns))
        bin_count = 0
        do m = 1, ns
            b = bin_of(m)
            order(bin_start(b) + bin_count(b)) = m
            bin_count(b) = bin_count(b) + 1
        end do

        ! The non-empty bins: their members, tight box, fastest velocity
        ! and earliest start time
        nb = count(bin_count > 0)
        allocate (bfirst(1:nb), blast(1:nb), bxl(1:nb), bxh(1:nb), byl(1:nb), byh(1:nb), &
            bzl(1:nb), bzh(1:nb), bvmax(1:nb), btmin(1:nb))
        h = 0
        do b = 1, nbin
            if (bin_count(b) == 0) then
                cycle
            end if
            h = h + 1
            bfirst(h) = bin_start(b)
            blast(h) = bin_start(b + 1) - 1
            associate (mm => order(bin_start(b):bin_start(b + 1) - 1))
                bxl(h) = minval(sx(mm))
                bxh(h) = maxval(sx(mm))
                byl(h) = minval(sy(mm))
                byh(h) = maxval(sy(mm))
                bzl(h) = minval(sz(mm))
                bzh(h) = maxval(sz(mm))
                bvmax(h) = maxval(sv(mm))
                btmin(h) = minval(st(mm))
            end associate
        end do

        !$omp parallel do collapse(2) schedule(dynamic) &
            !$omp private(i, j, k, m, b, h, best, win, prev, first, tval, lb, lbmin, tol)
        do k = 1, nz
            do j = 1, ny
                prev = 0
                do i = 1, nx

                    best = huge_value
                    win = 0
                    if (prev > 0) then
                        ! the neighbour's winner, which is usually this
                        ! point's as well: a tight bound from the start
                        best = source_time(prev, i, j, k)
                        win = prev
                    else
                        ! otherwise the bin with the lowest bound goes first
                        lbmin = huge_value
                        first = 1
                        do b = 1, nb
                            lb = bin_bound(b, i, j, k)
                            if (lb < lbmin) then
                                lbmin = lb
                                first = b
                            end if
                        end do
                        do h = bfirst(first), blast(first)
                            m = order(h)
                            tval = source_time(m, i, j, k)
                            if (tval < best .or. (tval == best .and. m < win)) then
                                best = tval
                                win = m
                            end if
                        end do
                    end if

                    ! Every bin that could hold a time at or below the best
                    ! so far; the margin keeps rounding in the bound from
                    ! ever ruling out the true minimum
                    tol = 1.0d-12*abs(best)
                    do b = 1, nb
                        if (bin_bound(b, i, j, k) > best + tol) then
                            cycle
                        end if
                        do h = bfirst(b), blast(b)
                            m = order(h)
                            tval = source_time(m, i, j, k)
                            if (tval < best .or. (tval == best .and. m < win)) then
                                best = tval
                                win = m
                                tol = 1.0d-12*abs(best)
                            end if
                        end do
                    end do

                    call set_background(win, i, j, k)
                    prev = win

                end do
            end do
        end do
        !$omp end parallel do

    contains

        ! The loop indices come in as arguments: inside a parallel loop a
        ! contained procedure would see the original i, j, k through host
        ! association, not the thread's private copies

        !> The time from source point m to the local grid point (i, j, k)
        double precision function source_time(m, i, j, k)
            integer, intent(in) :: m, i, j, k
            double precision, dimension(1:3) :: xn
            xn = node_ecef(i + ioff, j + joff, k + koff)
            source_time = sqrt((xn(1) - sx(m))**2 + (xn(2) - sy(m))**2 + (xn(3) - sz(m))**2)/sv(m)
            source_time = source_time + st(m)
        end function source_time

        !> A lower bound of the time from any source point of bin b: the
        !> distance from the node to the bin's box over the fastest velocity
        !> of the bin, plus its earliest start time
        double precision function bin_bound(b, i, j, k)
            integer, intent(in) :: b, i, j, k
            double precision, dimension(1:3) :: xn
            double precision :: cx, cy, cz
            xn = node_ecef(i + ioff, j + joff, k + koff)
            cx = max(bxl(b) - xn(1), 0.0d0, xn(1) - bxh(b))
            cy = max(byl(b) - xn(2), 0.0d0, xn(2) - byh(b))
            cz = max(bzl(b) - xn(3), 0.0d0, xn(3) - bzh(b))
            bin_bound = sqrt(cx**2 + cy**2 + cz**2)/bvmax(b) + btmin(b)
        end function bin_bound

        !> The winner's time at the local grid point (i, j, k), and the east,
        !> north and down components of its gradient
        subroutine set_background(m, i, j, k)
            integer, intent(in) :: m, i, j, k
            double precision, dimension(1:3) :: xn, ds, ue, un, ud
            double precision :: dist, t1
            integer :: ig, jg
            if (m == 0) then
                return
            end if
            ig = i + ioff
            jg = j + joff
            xn = node_ecef(ig, jg, k + koff)
            ds = xn - [sx(m), sy(m), sz(m)]
            dist = norm2(ds)
            t1 = dist/sv(m)
            if (t1 == 0) then
                pdxt0(i, j, k) = 0
                pdyt0(i, j, k) = 0
                pdzt0(i, j, k) = 0
            else
                ! Local unit vectors of the node: east, north and down
                ue = [-sph_sinlon(ig), sph_coslon(ig), 0.0d0]
                un = [-sph_sinlat(jg)*sph_coslon(ig), -sph_sinlat(jg)*sph_sinlon(ig), sph_coslat(jg)]
                ud = -[sph_coslat(jg)*sph_coslon(ig), sph_coslat(jg)*sph_sinlon(ig), sph_sinlat(jg)]
                pdxt0(i, j, k) = dot_product(ds, ue)/(sv(m)*dist)
                pdyt0(i, j, k) = dot_product(ds, un)/(sv(m)*dist)
                pdzt0(i, j, k) = dot_product(ds, ud)/(sv(m)*dist)
            end if
            t0(i, j, k) = t1 + st(m)
        end subroutine set_background

    end subroutine background_time

    !
    !> Which of n bins a coordinate falls in, over [low, low + extent]
    !
    pure integer function bin_index(x, low, extent, n)
        double precision, intent(in) :: x, low, extent
        integer, intent(in) :: n
        if (extent > 0 .and. n > 1) then
            bin_index = min(n, int((x - low)/extent*n) + 1)
        else
            bin_index = 1
        end if
    end function bin_index

    !
    !> A point's grid cell as one integer: the floor and ceiling indices of
    !> each coordinate, exactly what within_the_same_grid compares
    !
    function cell_key(p, n, o, d) result(key)
        double precision, dimension(:), intent(in) :: p, o, d
        integer, dimension(:), intent(in) :: n
        integer(kind=8) :: key
        integer(kind=8) :: c1, c2, c3
        double precision :: p1, p2, p3
        p1 = p(1) - o(1)
        p2 = p(2) - o(2)
        p3 = p(3) - o(3)
        c1 = 2*int(floor(p1/d(1)), 8) + (ceiling(p1/d(1)) - floor(p1/d(1))) + 4
        c2 = 2*int(floor(p2/d(2)), 8) + (ceiling(p2/d(2)) - floor(p2/d(2))) + 4
        c3 = 2*int(floor(p3/d(3)), 8) + (ceiling(p3/d(3)) - floor(p3/d(3))) + 4
        key = c1 + (2*int(n(1), 8) + 8)*(c2 + (2*int(n(2), 8) + 8)*c3)
    end function cell_key

    !
    !> The two grid points ii around a point at the 0-based fractional grid
    !> coordinate p, 0 <= p <= n - 1, and their linear interpolation weights
    !> w; for n = 1, both are the single grid point
    !
    pure subroutine linear_weights(p, n, ii, w)
        double precision, intent(in) :: p
        integer, intent(in) :: n
        integer, dimension(1:2), intent(out) :: ii
        double precision, dimension(1:2), intent(out) :: w
        if (n == 1) then
            ii = 1
            w = [1.0d0, 0.0d0]
        else
            ii(1) = min(max(floor(p), 0), n - 2) + 1
            ii(2) = ii(1) + 1
            w(2) = min(max(p - (ii(1) - 1), 0.0d0), 1.0d0)
            w(1) = 1 - w(2)
        end if
    end subroutine linear_weights

    !
    !> Set the offsets ioff, joff, koff of a per-shot grid of nx x ny x nz
    !> nodes whose first node lies at o = (lon, lat, depth) within the global
    !> grid of module spherical
    !
    subroutine set_grid_offsets(o, nx, ny, nz)

        real, dimension(1:3), intent(in) :: o
        integer, intent(in) :: nx, ny, nz

        ioff = nint((dble(o(1)) - sph_olon)/sph_dlon)
        joff = nint((dble(o(2)) - sph_olat)/sph_dlat)
        koff = nint((dble(o(3)) - sph_odepth)/sph_ddepth)

        call assert(ioff >= 0 .and. ioff + nx <= sph_nlon .and. &
            joff >= 0 .and. joff + ny <= sph_nlat .and. &
            koff >= 0 .and. koff + nz <= sph_ndepth, &
            ' <set_grid_offsets> Error: the shot grid lies outside the global grid')

    end subroutine set_grid_offsets

    !
    !> Fast sweeping factorized eikonal solver
    !
    subroutine forward_iso(v, d, o, geom, t, trec)

        real, dimension(:, :, :), intent(in) :: v
        real, dimension(1:3), intent(in) :: d, o
        type(source_receiver_geometry), intent(in) :: geom
        real, allocatable, dimension(:, :, :), intent(out) :: t
        real, allocatable, dimension(:, :), intent(out) :: trec

        double precision, allocatable, dimension(:, :, :) :: tt, tt_prev, vp
        double precision :: dx, dy, dz, ox, oy, oz
        integer :: nx, ny, nz, niter, i, j, k, l, isx, isy, isz, irx, iry, irz, itx, ity, itz
        double precision :: ttdiff, vsource
        double precision, dimension(1:3) :: xr, xs
        integer :: sw
        integer(kind=8) :: key
        integer(kind=8), allocatable, dimension(:) :: cellkey
        logical, allocatable, dimension(:) :: source_inside, receiver_inside
        double precision :: time1, time2

        dx = d(1)
        dy = d(2)
        dz = d(3)
        ox = o(1)
        oy = o(2)
        oz = o(3)

        vp = permute(v, 321)
        nx = size(vp, 1)
        ny = size(vp, 2)
        nz = size(vp, 3)

        ! Position of this (per-shot) grid within the global grid, for the
        ! spherical geometry arrays
        call set_grid_offsets(o, nx, ny, nz)

        ! Check if sources are at zero-velocity points
        source_inside = falses(geom%ns)
        !$omp parallel do private(l)
        do l = 1, geom%ns
            if (geom%srcr(l)%amp == 1) then
                source_inside(l) = point_in_domain(dble([geom%srcr(l)%x, geom%srcr(l)%y, geom%srcr(l)%z]), &
                    [nx, ny, nz], [ox, oy, oz], [dx, dy, dz], vp)
            end if
        end do
        !$omp end parallel do

        ! Check if receivers are at zero-velocity points
        receiver_inside = falses(geom%nr)
        !$omp parallel do private(i)
        do i = 1, geom%nr
            if (geom%recr(i)%weight /= 0) then
                receiver_inside(i) = point_in_domain(dble([geom%recr(i)%x, geom%recr(i)%y, geom%recr(i)%z]), &
                    [nx, ny, nz], [ox, oy, oz], [dx, dy, dz], vp)
            end if
        end do
        !$omp end parallel do

        ! Background time field: the analytic time from the source point
        ! that arrives first, and its derivatives
        call background_time(geom, source_inside, nx, ny, nz, ox, oy, oz, dx, dy, dz, vp)

        ! Initialize the multiplicative time field
        tt_prev = zeros(nx, ny, nz)
        tt = zeros(nx, ny, nz) + huge_value
        ! Increase the width of the multiplicative field around the source
        ! in the initialization seems to significantly improve the resulting accuracy,
        ! especially when the source point is not on an integer grid point
        sw = 3
        !$omp parallel do private(i, isx, isy, isz, itx, ity, itz)
        do i = 1, geom%ns
            if (source_inside(i)) then
                isx = max(floor((geom%srcr(i)%x - ox)/dx) + 1 - sw, 1)
                isy = max(floor((geom%srcr(i)%y - oy)/dy) + 1 - sw, 1)
                isz = max(floor((geom%srcr(i)%z - oz)/dz) + 1 - sw, 1)
                itx = min(ceiling((geom%srcr(i)%x - ox)/dx) + 1 + sw, nx)
                ity = min(ceiling((geom%srcr(i)%y - oy)/dy) + 1 + sw, ny)
                itz = min(ceiling((geom%srcr(i)%z - oz)/dz) + 1 + sw, nz)
                tt(isx:itx, isy:ity, isz:itz) = 1.0
            end if
        end do
        !$omp end parallel do

        !        ! As a comparison, the nearest grid point approach is
        !        !$omp parallel do private(i, isx, isy, isz, itx, ity, itz)
        !        do i = 1, geom%ns
        !            if (source_inside(i)) then
        !                isx = nint((geom%srcr(i)%x - ox)/dx) + 1
        !                isy = nint((geom%srcr(i)%y - oy)/dy) + 1
        !                isz = nint((geom%srcr(i)%z - oz)/dz) + 1
        !                tt(isx, isy, isz) = 1.0
        !            end if
        !        end do
        !        !$omp end parallel do

        ! Fast sweeping via Gauss-Seidel iterations
        niter = 0
        ttdiff = huge_value
        where (vp == 0)
            vp = float_tiny
        end where

        ! For timing the modeling, commented out for clarity
        time1 = omp_get_wtime()

        do while (ttdiff >= sweep_stop_threshold .and. niter < sweep_niter_max)

            tt_prev = tt
            call fast_sweep_forward(nx, ny, nz, 1, 1, 1, vp, t0, pdxt0, pdyt0, pdzt0, tt)
            call fast_sweep_forward(nx, ny, nz, 1, 1, -1, vp, t0, pdxt0, pdyt0, pdzt0, tt)
            call fast_sweep_forward(nx, ny, nz, 1, -1, 1, vp, t0, pdxt0, pdyt0, pdzt0, tt)
            call fast_sweep_forward(nx, ny, nz, 1, -1, -1, vp, t0, pdxt0, pdyt0, pdzt0, tt)
            call fast_sweep_forward(nx, ny, nz, -1, 1, 1, vp, t0, pdxt0, pdyt0, pdzt0, tt)
            call fast_sweep_forward(nx, ny, nz, -1, 1, -1, vp, t0, pdxt0, pdyt0, pdzt0, tt)
            call fast_sweep_forward(nx, ny, nz, -1, -1, 1, vp, t0, pdxt0, pdyt0, pdzt0, tt)
            call fast_sweep_forward(nx, ny, nz, -1, -1, -1, vp, t0, pdxt0, pdyt0, pdzt0, tt)

            ttdiff = mean(abs(tt - tt_prev))
            niter = niter + 1

        end do

        ! For timing the modeling, commented out for clarity
        time2 = omp_get_wtime()
        call warn(' Wall-clock time = '//num2str(time2 - time1, '(es)'))

        tt = t0*tt
        where (vp == float_tiny)
            tt = 0
            vp = 0
        end where

        ! Get the traveltime values at the receivers, if necessary
        trec = zeros(geom%nr, 1)
        ! A receiver in the same grid cell as a source point takes the
        ! analytic time from it -- from the last such source point, as the
        ! pairwise test within_the_same_grid always did. That test ran for
        ! every receiver against every source point; each source point's
        ! cell is now one integer (cell_key), computed once.
        allocate (cellkey(1:geom%ns))
        cellkey = -1
        !$omp parallel do private(l)
        do l = 1, geom%ns
            if (source_inside(l)) then
                cellkey(l) = cell_key(dble([geom%srcr(l)%x, geom%srcr(l)%y, geom%srcr(l)%z]), &
                    [nx, ny, nz], [ox, oy, oz], [dx, dy, dz])
            end if
        end do
        !$omp end parallel do
        !$omp parallel do private(i, irx, iry, irz, l, key, vsource, xr, xs)
        do i = 1, geom%nr
            if (receiver_inside(i)) then

                ! Linear interpolation; same as for the source points,
                ! here the implementation automatically handles the case of
                ! a receiver falling on integer grid points.
                trec(i, 1) = get_point_value_inside(dble([geom%recr(i)%x, geom%recr(i)%y, geom%recr(i)%z]), &
                    [nx, ny, nz], [ox, oy, oz], [dx, dy, dz], tt)
                key = cell_key(dble([geom%recr(i)%x, geom%recr(i)%y, geom%recr(i)%z]), &
                    [nx, ny, nz], [ox, oy, oz], [dx, dy, dz])
                do l = geom%ns, 1, -1
                    if (cellkey(l) == key) then
                        vsource = get_point_value_inside(dble([geom%srcr(l)%x, geom%srcr(l)%y, geom%srcr(l)%z]), &
                            [nx, ny, nz], [ox, oy, oz], [dx, dy, dz], vp)
                        ! Straight-line distance between the Earth-centred positions
                        xr = lonlatdepth_to_ecef(dble(geom%recr(i)%x), dble(geom%recr(i)%y), dble(geom%recr(i)%z))
                        xs = lonlatdepth_to_ecef(dble(geom%srcr(l)%x), dble(geom%srcr(l)%y), dble(geom%srcr(l)%z))
                        trec(i, 1) = norm2(xr - xs)/vsource + geom%srcr(l)%t0
                        exit
                    end if
                end do

                ! Add virtual receiver time for TLOC
                trec(i, 1) = trec(i, 1) + geom%recr(i)%t0

                !                ! As a comparison, the nearest grid point approach is
                !                irx = nint((geom%recr(i)%x - ox)/dx) + 1
                !                iry = nint((geom%recr(i)%y - oy)/dy) + 1
                !                irz = nint((geom%recr(i)%z - oz)/dz) + 1
                !                trec(i, 1) = t(irx, iry, irz) + geom%recr(i)%t0

            end if
        end do
        !$omp end parallel do

        t = permute(tt, 321)

        call warn(date_time_compact()//' Fast sweeping eikonal niter = '//num2str(niter)// &
            ', relative diff = '//num2str(ttdiff, '(es)'))

    end subroutine forward_iso

    !
    !> Compute the adjoint-state field of a shot, for FATT and TLOC
    !>
    !> The traveltime residuals at the receivers, dt_r in tresidual, enter
    !> the adjoint-state equation as true point sources,
    !>
    !>     -div(lambda grad T) = sum_r dt_r delta(x - x_r),
    !>
    !> the linear injection of the Cartesian LATTE: the source of a receiver
    !> goes to the eight corners of the receiver's grid cell with the
    !> trilinear weights by which forward_iso interpolates the receiver
    !> traveltime, i.e., the transpose of that interpolation, divided by the
    !> exact spherical cell volume of each corner node. No grid point is
    !> fixed, so the sources of nearby receivers add up. The spherical branch
    !> supports only this injection.
    !
    subroutine adjoint_iso(v, d, o, geom, t, tresidual, tadj)

        real, dimension(:, :, :), intent(in) :: v, t
        real, dimension(1:3), intent(in) :: d, o
        type(source_receiver_geometry), intent(in) :: geom
        real, dimension(:, :), intent(in) :: tresidual
        real, dimension(:, :, :), allocatable, intent(out) :: tadj

        double precision, allocatable, dimension(:, :, :) :: vp, lambda_prev
        integer :: nx, ny, nz, iter, i, a, b, c
        integer, dimension(1:2) :: ii, jj, kk
        double precision :: lambda_diff, dx, dy, dz, ox, oy, oz
        double precision, dimension(1:2) :: wx, wy, wz

        dx = d(1)
        dy = d(2)
        dz = d(3)
        ox = o(1)
        oy = o(2)
        oz = o(3)

        tt = permute(t, 321)
        vp = permute(v, 321)
        nx = size(tt, 1)
        ny = size(tt, 2)
        nz = size(tt, 3)

        ! Position of this (per-shot) grid within the global grid, for the
        ! spherical geometry arrays
        call set_grid_offsets(o, nx, ny, nz)

        lambda_prev = zeros(nx, ny, nz)
        lambda_diff = huge_value

        ! The point sources of the receivers inside the medium, the same
        ! receivers as in forward_iso; serial because receivers may share
        ! grid points
        adjsrc = zeros(nx, ny, nz)
        do i = 1, geom%nr
            if (geom%recr(i)%weight /= 0) then
                if (point_in_domain(dble([geom%recr(i)%x, geom%recr(i)%y, geom%recr(i)%z]), &
                        [nx, ny, nz], [ox, oy, oz], [dx, dy, dz], vp)) then
                    call linear_weights((geom%recr(i)%x - ox)/dx, nx, ii, wx)
                    call linear_weights((geom%recr(i)%y - oy)/dy, ny, jj, wy)
                    call linear_weights((geom%recr(i)%z - oz)/dz, nz, kk, wz)
                    do c = 1, 2
                        do b = 1, 2
                            do a = 1, 2
                                adjsrc(ii(a), jj(b), kk(c)) = adjsrc(ii(a), jj(b), kk(c)) &
                                    + tresidual(i, 1)*wx(a)*wy(b)*wz(c) &
                                    /cell_volume(jj(b) + joff, kk(c) + koff)
                            end do
                        end do
                    end do
                end if
            end if
        end do
        lambda = zeros(nx, ny, nz) + huge_value

        ! Fast sweep to solve the adjoint-state equation
        iter = 0

        do while (lambda_diff >= sweep_stop_threshold .and. iter < sweep_niter_max)

            ! Save previous step
            lambda_prev = lambda

            ! Fast sweeps
            call fast_sweep_adjoint(nx, ny, nz, 1, 1, 1, lambda)
            call fast_sweep_adjoint(nx, ny, nz, 1, 1, -1, lambda)
            call fast_sweep_adjoint(nx, ny, nz, 1, -1, 1, lambda)
            call fast_sweep_adjoint(nx, ny, nz, 1, -1, -1, lambda)
            call fast_sweep_adjoint(nx, ny, nz, -1, 1, 1, lambda)
            call fast_sweep_adjoint(nx, ny, nz, -1, 1, -1, lambda)
            call fast_sweep_adjoint(nx, ny, nz, -1, -1, 1, lambda)
            call fast_sweep_adjoint(nx, ny, nz, -1, -1, -1, lambda)

            ! Check threshold
            lambda_diff = mean(abs(lambda - lambda_prev))

            iter = iter + 1

        end do

        tadj = permute(lambda, 321)

        call warn(date_time_compact()//' Fast sweeping adjoint niter = '//num2str(iter)// &
            ', relative diff = '//num2str(lambda_diff, '(es)'))

    end subroutine adjoint_iso

end module traveltime_iso
