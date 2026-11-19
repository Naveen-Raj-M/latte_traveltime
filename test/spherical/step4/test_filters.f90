!
! Unit test of the exact-distance filters of step 4 (module spherical_filter
! and the Tikhonov solver of module regularization).
!
! Grid: longitude 0 to 20 degrees, geodetic latitude 10 to 70 degrees, depth
! 0 to 100 km; 0.2 degree x 0.2 degree x 5 km cells. A degree of longitude is
! about half as long at 60 degrees as at 20 degrees latitude, so filters
! that work in cells give different physical widths there, and exact ones
! do not.
!
! Checks:
!   1. Gaussian smoothing of a spike: the standard deviation, in meters, of
!      the smoothed spike along a parallel at 20 and at 60 degrees and along
!      a meridian equals the requested 30 km, and in depth the requested
!      10 km (a smaller one, so that the 100 km depth range does not fold
!      the tails back at its edges).
!   2. Gaussian smoothing of a constant: stays constant, also at the edges.
!   3. Taper: the weights equal the Blackman ramp at the exact distances.
!   4. Tikhonov: a spike spreads over the same physical east-west distance
!      at 20 and 60 degrees latitude, and a constant stays constant.
!   5. Wavenumber filtering along longitude: a low pass (amplitude 1 up to
!      1/100 km, 0 from 1/60 km) removes a 40 km wavelength and keeps a
!      300 km wavelength, both in meters along the line, at 20 and at 60
!      degrees latitude. Measured in the middle 60% of each line, away from
!      the edge effects of the FFT.
!
program test_filters

    use libflit
    use parameters
    use spherical
    use spherical_filter
    use regularization

    implicit none

    integer, parameter :: mlon = 101, mlat = 301, mdep = 21
    real, allocatable, dimension(:, :, :) :: w, ws
    double precision :: s0, s1, s2, d, expect, maxerr
    integer :: i, j, k, jj, kk, ic
    integer, dimension(1:2) :: rows

    call init_spherical_grid(mlon, mlat, mdep, 0.0d0, 10.0d0, 0.0d0, 0.2d0, 0.2d0, 5000.0d0)
    ic = 51
    kk = 11
    rows = [51, 251]

    print '(a)', ' 1. Gaussian smoothing of a spike, requested standard deviations 30000 m (lon, lat), 10000 m (depth)'
    do jj = 1, 2
        allocate (w(mdep, mlat, mlon))
        w = 0.0
        w(kk, rows(jj), ic) = 1.0
        ws = sph_gauss_filt(w, 1, 1, [30000.0, 30000.0, 10000.0])
        ! East-west: along the parallel through the spike
        s0 = 0; s1 = 0; s2 = 0
        do i = 1, mlon
            d = (i - ic)*arc_lon(rows(jj), kk)
            s0 = s0 + ws(kk, rows(jj), i)
            s1 = s1 + ws(kk, rows(jj), i)*d
            s2 = s2 + ws(kk, rows(jj), i)*d**2
        end do
        print '(a, f6.2, a, f10.1, a, f8.2, a)', '    latitude ', 10.0 + 0.2*(rows(jj) - 1), &
            ': east-west standard deviation = ', sqrt(s2/s0 - (s1/s0)**2), ' m (', &
            sqrt(s2/s0 - (s1/s0)**2)/arc_lon(rows(jj), kk), ' cells)'
        if (jj == 1) then
            ! North-south and depth at the first spike
            s0 = 0; s1 = 0; s2 = 0
            do j = 1, mlat
                d = sph_r(kk)*(sph_lat(j) - sph_lat(rows(jj)))
                s0 = s0 + ws(kk, j, ic)
                s1 = s1 + ws(kk, j, ic)*d
                s2 = s2 + ws(kk, j, ic)*d**2
            end do
            print '(a, f10.1, a)', '    north-south standard deviation = ', sqrt(s2/s0 - (s1/s0)**2), ' m'
            s0 = 0; s1 = 0; s2 = 0
            do k = 1, mdep
                d = (k - kk)*sph_ddepth
                s0 = s0 + ws(k, rows(jj), ic)
                s1 = s1 + ws(k, rows(jj), ic)*d
                s2 = s2 + ws(k, rows(jj), ic)*d**2
            end do
            print '(a, f10.1, a)', '    vertical standard deviation    = ', sqrt(s2/s0 - (s1/s0)**2), ' m'
        end if
        deallocate (w)
    end do

    print '(a)', ' 2. Gaussian smoothing of a constant field'
    allocate (w(mdep, mlat, mlon))
    w = 3.0
    ws = sph_gauss_filt(w, 1, 1, [80000.0, 80000.0, 30000.0])
    print '(a, es10.3)', '    max |smoothed - 3| = ', maxval(abs(ws - 3.0))

    print '(a)', ' 3. Taper of a constant field: west 200 km, north 300 km, bottom 40 km'
    w = 1.0
    ws = sph_taper(w, 1, 1, [200000.0, 0.0, 0.0, 300000.0, 0.0, 40000.0])
    maxerr = 0
    do i = 1, mlon
        do j = 1, mlat
            do k = 1, mdep
                expect = ramp((i - 1)*arc_lon(j, k)/200000.0d0) &
                    *ramp(sph_r(k)*(sph_lat(mlat) - sph_lat(j))/300000.0d0) &
                    *ramp((mdep - k)*sph_ddepth/40000.0d0)
                maxerr = max(maxerr, abs(ws(k, j, i) - expect))
            end do
        end do
    end do
    print '(a, es10.3)', '    max |weight - Blackman ramp at exact distance| = ', maxerr
    i = 1 + nint(100000.0d0/arc_lon(251, kk))
    print '(a, f9.1, a, f8.5, a, f8.5)', '    node ', (i - 1)*arc_lon(251, kk), &
        ' m from the west edge at 60 deg: weight ', ws(kk, 251, i), &
        ', Blackman ramp at that distance ', ramp((i - 1)*arc_lon(251, kk)/200000.0d0)

    print '(a)', ' 4. Tikhonov, lambda = (20 km)^2'
    do jj = 1, 2
        w = 0.0
        w(kk, rows(jj), ic) = 1.0
        ws = tikhonov_denoise(w, 4.0e8)
        s0 = 0; s1 = 0; s2 = 0
        do i = 1, mlon
            d = (i - ic)*arc_lon(rows(jj), kk)
            s0 = s0 + ws(kk, rows(jj), i)
            s1 = s1 + ws(kk, rows(jj), i)*d
            s2 = s2 + ws(kk, rows(jj), i)*d**2
        end do
        print '(a, f6.2, a, f10.1, a)', '    latitude ', 10.0 + 0.2*(rows(jj) - 1), &
            ': east-west spread of the spike (standard deviation) = ', sqrt(s2/s0 - (s1/s0)**2), ' m'
    end do
    w = 3.0
    ws = tikhonov_denoise(w, 4.0e8)
    print '(a, es10.3)', '    constant field: max |result - 3| = ', maxval(abs(ws - 3.0))

    print '(a)', ' 5. Wavenumber low pass along longitude: kept fraction of the RMS amplitude'
    do jj = 1, 2
        do kk = 1, 2
            do i = 1, mlon
                d = (i - 1)*arc_lon(rows(jj), 11)
                w(11, rows(jj), i) = real(sin(2.0d0*sph_pi*d/merge(40000.0d0, 300000.0d0, kk == 1)))
            end do
            ws = sph_fourier_filt_lon(w(11:11, rows(jj):rows(jj), :), rows(jj), 11, &
                [0.0, 1.0/100000.0, 1.0/60000.0], [1.0, 1.0, 0.0])
            print '(a, f6.2, a, f8.0, a, f8.4)', '    latitude ', 10.0 + 0.2*(rows(jj) - 1), &
                ', wavelength ', merge(40000.0, 300000.0, kk == 1), ' m: ', &
                norm2(ws(1, 1, 21:81))/norm2(w(11, rows(jj), 21:81))
        end do
    end do

contains

    pure double precision function ramp(x)
        double precision, intent(in) :: x
        double precision :: y
        y = min(max(x, 0.0d0), 1.0d0)
        ramp = 0.42d0 - 0.5d0*cos(sph_pi*y) + 0.08d0*cos(2.0d0*sph_pi*y)
    end function ramp

end program test_filters
