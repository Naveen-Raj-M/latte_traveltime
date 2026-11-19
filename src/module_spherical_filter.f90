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

!
!> Filters of the spherical branch that use exact distances on the
!> spherical grid of module spherical.
!>
!> Arrays are blocks of the global grid in the model order (depth,
!> latitude, longitude): w(k, j, i) is global node (i + ioff, jbeg + j - 1,
!> kbeg + k - 1). Nothing depends on the longitude offset, because the
!> geometry of the grid does not depend on longitude.
!>
!> Distances along the three kinds of grid lines are
!>
!>     longitude line at row j, level k:  s_i = (i - 1) r_k cos(theta_j) dlambda
!>     latitude line at level k:          s_j = r_k (theta_j - theta_jbeg)
!>     depth line:                        s_k = (k - 1) ddepth
!>
!> i.e., arc lengths along the line itself, where theta is the geocentric
!> latitude of a row.
!
module spherical_filter

    use libflit
    use spherical

    implicit none

    private
    public :: sph_gauss_filt
    public :: sph_taper
    public :: sph_fourier_filt_lon
    public :: sph_fourier_filt_depth

contains

    !
    !> Positions along the grid lines of a block; see the module header
    !
    pure function line_lon(n, j, k) result(s)
        integer, intent(in) :: n, j, k
        double precision, dimension(1:n) :: s
        integer :: i
        do i = 1, n
            s(i) = (i - 1)*arc_lon(j, k)
        end do
    end function line_lon

    pure function line_lat(n, jbeg, k) result(s)
        integer, intent(in) :: n, jbeg, k
        double precision, dimension(1:n) :: s
        integer :: j
        do j = 1, n
            s(j) = sph_r(k)*(sph_lat(jbeg + j - 1) - sph_lat(jbeg))
        end do
    end function line_lat

    pure function line_depth(n) result(s)
        integer, intent(in) :: n
        double precision, dimension(1:n) :: s
        integer :: k
        do k = 1, n
            s(k) = (k - 1)*sph_ddepth
        end do
    end function line_depth

    !
    !> Gaussian filter of standard deviation sigma (a length) along a line of
    !> samples f at the positions s, increasing. The filter is the direct
    !> convolution with the Gaussian of the exact distances, truncated at
    !> 6 sigma like the sampled Gaussian filter of FLIT, with the boundary
    !> treated as FLIT does: the line is extended by repeated mirror
    !> reflection about its end samples. The weights are normalized by their
    !> sum.
    !
    subroutine gauss_line(f, s, sigma)

        double precision, dimension(:), intent(inout) :: f
        double precision, dimension(:), intent(in) :: s
        double precision, intent(in) :: sigma

        double precision, allocatable, dimension(:) :: g
        double precision :: period, p, p0, wsum, fsum, wgt, cut
        integer :: n, m, e, dirn, idx

        n = size(f)
        if (n <= 1 .or. sigma <= 0) then
            return
        end if

        period = 2.0d0*(s(n) - s(1))
        cut = 6.0d0*sigma
        g = f

        do m = 1, n
            p0 = s(m) - s(1)
            wsum = 0.0d0
            fsum = 0.0d0
            ! The sample itself, then outward in both directions through the
            ! mirrored extension until the distance exceeds the cut
            do dirn = -1, 1, 2
                e = m - 1
                if (dirn == 1) then
                    e = e + 1
                end if
                do
                    if (dirn == -1) then
                        call unfold(e, idx, p)
                        e = e - 1
                    else
                        call unfold(e, idx, p)
                        e = e + 1
                    end if
                    if (abs(p - p0) > cut) then
                        exit
                    end if
                    wgt = exp(-0.5d0*((p - p0)/sigma)**2)
                    wsum = wsum + wgt
                    fsum = fsum + wgt*g(idx)
                end do
            end do
            f(m) = fsum/wsum
        end do

    contains

        !> Sample index and unfolded position of extended index e (0-based,
        !> any integer) of the mirrored periodic extension of the line
        subroutine unfold(e, idx, p)
            integer, intent(in) :: e
            integer, intent(out) :: idx
            double precision, intent(out) :: p
            integer :: q, r
            q = floor(dble(e)/dble(2*(n - 1)))
            r = e - q*2*(n - 1)
            if (r <= n - 1) then
                idx = r + 1
                p = q*period + (s(idx) - s(1))
            else
                idx = 2*(n - 1) - r + 1
                p = q*period + period - (s(idx) - s(1))
            end if
        end subroutine unfold

    end subroutine gauss_line

    !
    !> Separable Gaussian filter of a block w (depth, latitude, longitude),
    !> whose first row and level are the global row jbeg and level kbeg,
    !> with standard deviations sigma = [lon, lat, depth] as lengths along the
    !> grid lines; a nonpositive sigma leaves that axis unfiltered. The passes
    !> run along longitude, then latitude, then depth.
    !
    function sph_gauss_filt(w, jbeg, kbeg, sigma) result(wr)

        real, dimension(:, :, :), intent(in) :: w
        integer, intent(in) :: jbeg, kbeg
        real, dimension(1:3), intent(in) :: sigma
        real, allocatable, dimension(:, :, :) :: wr

        double precision, allocatable, dimension(:) :: f, s
        integer :: nk, nj, ni, i, j, k

        nk = size(w, 1)
        nj = size(w, 2)
        ni = size(w, 3)
        wr = w

        if (sigma(1) > 0 .and. ni > 1) then
            !$omp parallel do collapse(2) private(j, k, f, s)
            do j = 1, nj
                do k = 1, nk
                    s = line_lon(ni, jbeg + j - 1, kbeg + k - 1)
                    f = dble(wr(k, j, :))
                    call gauss_line(f, s, dble(sigma(1)))
                    wr(k, j, :) = real(f)
                end do
            end do
            !$omp end parallel do
        end if

        if (sigma(2) > 0 .and. nj > 1) then
            !$omp parallel do collapse(2) private(i, k, f, s)
            do i = 1, ni
                do k = 1, nk
                    s = line_lat(nj, jbeg, kbeg + k - 1)
                    f = dble(wr(k, :, i))
                    call gauss_line(f, s, dble(sigma(2)))
                    wr(k, :, i) = real(f)
                end do
            end do
            !$omp end parallel do
        end if

        if (sigma(3) > 0 .and. nk > 1) then
            s = line_depth(nk)
            !$omp parallel do collapse(2) private(i, j, f)
            do i = 1, ni
                do j = 1, nj
                    f = dble(wr(:, j, i))
                    call gauss_line(f, s, dble(sigma(3)))
                    wr(:, j, i) = real(f)
                end do
            end do
            !$omp end parallel do
        end if

    end function sph_gauss_filt

    !
    !> Blackman ramp of FLIT's taper at the fraction x of the taper length:
    !> 0.42 - 0.5 cos(pi x) + 0.08 cos(2 pi x), 0 at x = 0 and 1 at x >= 1
    !
    elemental function blackman_ramp(x) result(b)
        double precision, intent(in) :: x
        double precision :: b
        double precision :: y
        y = min(max(x, 0.0d0), 1.0d0)
        b = 0.42d0 - 0.5d0*cos(sph_pi*y) + 0.08d0*cos(2.0d0*sph_pi*y)
    end function blackman_ramp

    !
    !> Taper the sides of a block w (depth, latitude, longitude), whose first
    !> row and level are the global row jbeg and level kbeg. len = [west,
    !> east, south, north, top, bottom] are taper lengths; a nonpositive
    !> length leaves that side untapered. A node at distance d from a side,
    !> measured along its grid line from the edge node, is multiplied by
    !> blackman_ramp(d/len), which is 0 at the edge and 1 from the length on.
    !
    function sph_taper(w, jbeg, kbeg, len) result(wr)

        real, dimension(:, :, :), intent(in) :: w
        integer, intent(in) :: jbeg, kbeg
        real, dimension(1:6), intent(in) :: len
        real, allocatable, dimension(:, :, :) :: wr

        integer :: nk, nj, ni, i, j, k, jg, kg
        double precision :: f, d, h

        nk = size(w, 1)
        nj = size(w, 2)
        ni = size(w, 3)
        wr = w

        !$omp parallel do collapse(2) private(i, j, k, jg, kg, f, d, h)
        do i = 1, ni
            do j = 1, nj
                do k = 1, nk
                    jg = jbeg + j - 1
                    kg = kbeg + k - 1
                    f = 1.0d0
                    h = arc_lon(jg, kg)
                    if (len(1) > 0) then
                        d = (i - 1)*h
                        f = f*blackman_ramp(d/len(1))
                    end if
                    if (len(2) > 0) then
                        d = (ni - i)*h
                        f = f*blackman_ramp(d/len(2))
                    end if
                    if (len(3) > 0) then
                        d = sph_r(kg)*(sph_lat(jg) - sph_lat(jbeg))
                        f = f*blackman_ramp(d/len(3))
                    end if
                    if (len(4) > 0) then
                        d = sph_r(kg)*(sph_lat(jbeg + nj - 1) - sph_lat(jg))
                        f = f*blackman_ramp(d/len(4))
                    end if
                    if (len(5) > 0) then
                        d = (k - 1)*sph_ddepth
                        f = f*blackman_ramp(d/len(5))
                    end if
                    if (len(6) > 0) then
                        d = (nk - k)*sph_ddepth
                        f = f*blackman_ramp(d/len(6))
                    end if
                    wr(k, j, i) = real(w(k, j, i)*f)
                end do
            end do
        end do
        !$omp end parallel do

    end function sph_taper

    !
    !> Wavenumber-domain filter of FLIT along the longitude lines of a block,
    !> each with its own exact spacing r_k cos(theta_j) dlambda
    !
    function sph_fourier_filt_lon(w, jbeg, kbeg, wavenums, amps) result(wr)

        real, dimension(:, :, :), intent(in) :: w
        integer, intent(in) :: jbeg, kbeg
        real, dimension(:), intent(in) :: wavenums, amps
        real, allocatable, dimension(:, :, :) :: wr
        integer :: j, k

        wr = w
        !$omp parallel do collapse(2) private(j, k)
        do j = 1, size(w, 2)
            do k = 1, size(w, 1)
                wr(k, j, :) = fourier_filt(w(k, j, :), real(arc_lon(jbeg + j - 1, kbeg + k - 1)), wavenums, amps)
            end do
        end do
        !$omp end parallel do

    end function sph_fourier_filt_lon

    !
    !> Wavenumber-domain filter of FLIT along the depth lines of a block,
    !> with the spacing ddepth
    !
    function sph_fourier_filt_depth(w, wavenums, amps) result(wr)

        real, dimension(:, :, :), intent(in) :: w
        real, dimension(:), intent(in) :: wavenums, amps
        real, allocatable, dimension(:, :, :) :: wr
        integer :: i, j

        wr = w
        !$omp parallel do collapse(2) private(i, j)
        do i = 1, size(w, 3)
            do j = 1, size(w, 2)
                wr(:, j, i) = fourier_filt(w(:, j, i), real(sph_ddepth), wavenums, amps)
            end do
        end do
        !$omp end parallel do

    end function sph_fourier_filt_depth

end module spherical_filter
