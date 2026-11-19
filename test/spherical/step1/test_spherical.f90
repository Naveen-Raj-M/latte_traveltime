program test_spherical
    use libflit
    use spherical
    implicit none
    double precision :: v, vexact, latc, x(3), m(2), h(3), r1, r2, s1, s2
    integer :: j, k
    ! 1. Geodetic -> geocentric at 45 deg, and round trip
    latc = geodetic_to_geocentric(45.0d0)
    print '(a,f12.8)', ' geocentric latitude of geodetic 45 deg (deg)   = ', latc/deg2rad
    print '(a,es12.4)', ' round trip error (deg)                          = ', geocentric_to_geodetic(latc) - 45.0d0
    print '(a,es12.4)', ' round trip error at 89.9 deg                    = ', geocentric_to_geodetic(geodetic_to_geocentric(89.9d0)) - 89.9d0
    ! 2. Earth-centred position of (0, 0, 0) and (90 E, 0, 1000 m)
    x = lonlatdepth_to_ecef(0.0d0, 0.0d0, 0.0d0)
    print '(a,3f14.3)', ' ECEF of (0,0,0)                                  = ', x
    x = lonlatdepth_to_ecef(90.0d0, 0.0d0, 1000.0d0)
    print '(a,3f14.3)', ' ECEF of (90E,0,1000 m)                           = ', x
    ! 3. Grid: 101 x 81 x 41, lon 100..110, lat 20..28, depth -2 km..78 km
    call init_spherical_grid(101, 81, 41, 100.0d0, 20.0d0, -2000.0d0, 0.1d0, 0.1d0, 2000.0d0)
    v = 0
    do k = 1, sph_ndepth
        do j = 1, sph_nlat
            v = v + cell_volume(j, k)*sph_nlon
        end do
    end do
    r1 = earth_radius - (-2000.0d0 - 1000.0d0)
    r2 = earth_radius - (-2000.0d0 + 40.5d0*2000.0d0)
    s1 = sin(geodetic_to_geocentric(20.0d0 - 0.05d0))
    s2 = sin(geodetic_to_geocentric(28.0d0 + 0.05d0))
    vexact = (r1**3 - r2**3)/3.0d0*(s2 - s1)*(101*0.1d0*deg2rad)
    print '(a,es20.12)', ' sum of cell volumes                              = ', v
    print '(a,es20.12)', ' exact volume of the region they tile             = ', vexact
    print '(a,es12.4)',  ' relative difference                              = ', (v - vexact)/vexact
    ! 4. Arc lengths against chord-based check: two neighbouring nodes along a parallel
    x = node_ecef(2, 41, 1) - node_ecef(1, 41, 1)
    print '(a,2f16.6)', ' arc along parallel vs chord (m)                  = ', arc_lon(41, 1), norm2(x)
    x = node_ecef(1, 42, 1) - node_ecef(1, 41, 1)
    print '(a,2f16.6)', ' arc along meridian vs chord (m)                  = ', arc_lat(41, 42, 1), norm2(x)
    h = reference_spacing()
    print '(a,3f14.4)', ' reference spacing at centre, top (m)             = ', h
    m = min_length_per_degree()
    print '(a,2f14.4)', ' smallest length of 1 deg lon, lat in grid (m)     = ', m
end program
