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
!> Spherical-Earth geometry of the spherical branch of LATTE.
!>
!> The model grid is regular in longitude (degrees), WGS84 geodetic latitude
!> (degrees) and depth below mean sea level (length unit of the model, by
!> default meters). The Earth is a sphere of radius earth_radius. A point
!> (lon, lat, depth) lies on the sphere at
!>
!>     longitude  lambda = lon
!>     latitude   theta  = geocentric latitude of the WGS84 geodetic latitude lat
!>     radius     r      = earth_radius - depth
!>
!> where the geodetic-to-geocentric conversion is the exact relation on the
!> WGS84 ellipsoid surface, tan(theta) = (1 - e^2) tan(lat). The same
!> conversion is used at every depth, so a grid column of constant (lon, lat)
!> is a radial line of the sphere. Every geometric quantity below is the exact
!> one of this sphere; the only approximation is the spherical Earth itself.
!>
!> The arrays of the grid geometry are allocated and filled once, by
!> init_spherical_grid during initialization, in double precision, and do
!> not change during a run.
!
module spherical

    use libflit

    implicit none

    ! WGS84 ellipsoid, used only to convert geodetic latitude
    double precision, parameter :: wgs84_a = 6378137.0d0
    double precision, parameter :: wgs84_f = 1.0d0/298.257223563d0
    double precision, parameter :: wgs84_e2 = wgs84_f*(2.0d0 - wgs84_f)

    double precision, parameter :: sph_pi = 3.141592653589793238462643383279502884d0
    double precision, parameter :: deg2rad = sph_pi/180.0d0

    !> Radius of the spherical Earth, in the length unit of the model
    double precision :: earth_radius = 6371000.0d0

    !> The global grid, as set by init_spherical_grid
    integer :: sph_nlon = 0, sph_nlat = 0, sph_ndepth = 0
    double precision :: sph_olon, sph_dlon, sph_olat, sph_dlat, sph_odepth, sph_ddepth
    !> Longitude interval in radians
    double precision :: sph_dlon_rad

    !> Cosine and sine of the longitude of each grid column, 1:nlon
    double precision, allocatable, dimension(:) :: sph_coslon, sph_sinlon
    !> Geocentric latitude (radians) of each grid row, with its cosine and sine, 1:nlat
    double precision, allocatable, dimension(:) :: sph_lat, sph_coslat, sph_sinlat
    !> Geocentric latitude (radians) of the cell faces between rows, 1:nlat + 1;
    !> face j lies at the geodetic latitude olat + (j - 1.5)*dlat, so row j owns
    !> the cell between faces j and j + 1
    double precision, allocatable, dimension(:) :: sph_lat_face, sph_coslat_face, sph_sinlat_face
    !> Radius of each grid level, 1:ndepth
    double precision, allocatable, dimension(:) :: sph_r
    !> Radius of the cell faces between levels, 1:ndepth + 1; face k lies at the
    !> depth odepth + (k - 1.5)*ddepth, so level k owns the cell between faces
    !> k (shallower) and k + 1 (deeper)
    double precision, allocatable, dimension(:) :: sph_r_face

    private :: sph_check

contains

    !
    !> Geocentric latitude (radians) of a WGS84 geodetic latitude (degrees)
    !
    elemental function geodetic_to_geocentric(latg) result(latc)

        double precision, intent(in) :: latg
        double precision :: latc

        latc = atan2((1.0d0 - wgs84_e2)*sin(latg*deg2rad), cos(latg*deg2rad))

    end function geodetic_to_geocentric

    !
    !> WGS84 geodetic latitude (degrees) of a geocentric latitude (radians)
    !
    elemental function geocentric_to_geodetic(latc) result(latg)

        double precision, intent(in) :: latc
        double precision :: latg

        latg = atan2(sin(latc), (1.0d0 - wgs84_e2)*cos(latc))/deg2rad

    end function geocentric_to_geodetic

    !
    !> Derivative of the geocentric latitude (radians) with respect to the
    !> geodetic latitude (radians), d(theta)/d(lat), at a geodetic latitude
    !> (degrees)
    !
    elemental function dgeocentric_dgeodetic(latg) result(dd)

        double precision, intent(in) :: latg
        double precision :: dd
        double precision :: c, s

        c = cos(latg*deg2rad)
        s = sin(latg*deg2rad)
        dd = (1.0d0 - wgs84_e2)/(c**2 + (1.0d0 - wgs84_e2)**2*s**2)

    end function dgeocentric_dgeodetic

    !
    !> Earth-centred Cartesian position of a point given by its longitude
    !> (degrees), geodetic latitude (degrees) and depth below mean sea level
    !
    pure function lonlatdepth_to_ecef(lon, latg, depth) result(x)

        double precision, intent(in) :: lon, latg, depth
        double precision, dimension(1:3) :: x
        double precision :: r, latc

        r = earth_radius - depth
        latc = geodetic_to_geocentric(latg)
        x(1) = r*cos(latc)*cos(lon*deg2rad)
        x(2) = r*cos(latc)*sin(lon*deg2rad)
        x(3) = r*sin(latc)

    end function lonlatdepth_to_ecef

    !
    !> Earth-centred Cartesian position of the global grid node (i, j, k),
    !> i along longitude, j along latitude, k along depth
    !
    pure function node_ecef(i, j, k) result(x)

        integer, intent(in) :: i, j, k
        double precision, dimension(1:3) :: x

        x(1) = sph_r(k)*sph_coslat(j)*sph_coslon(i)
        x(2) = sph_r(k)*sph_coslat(j)*sph_sinlon(i)
        x(3) = sph_r(k)*sph_sinlat(j)

    end function node_ecef

    !
    !> Volume of the cell of the global grid node (i, j, k), which does not
    !> depend on i:
    !>
    !>     V = (r_k-1/2^3 - r_k+1/2^3)/3 * (sin(theta_j+1/2) - sin(theta_j-1/2)) * dlambda
    !
    pure function cell_volume(j, k) result(v)

        integer, intent(in) :: j, k
        double precision :: v

        v = (sph_r_face(k)**3 - sph_r_face(k + 1)**3)/3.0d0 &
            *(sph_sinlat_face(j + 1) - sph_sinlat_face(j))*sph_dlon_rad

    end function cell_volume

    !
    !> Arc length between neighbouring nodes along a parallel, at row j and
    !> level k: r_k cos(theta_j) dlambda
    !
    pure function arc_lon(j, k) result(h)

        integer, intent(in) :: j, k
        double precision :: h

        h = sph_r(k)*sph_coslat(j)*sph_dlon_rad

    end function arc_lon

    !
    !> Arc length along a meridian between rows j1 and j2 at level k:
    !> r_k |theta_j2 - theta_j1|
    !
    pure function arc_lat(j1, j2, k) result(h)

        integer, intent(in) :: j1, j2, k
        double precision :: h

        h = sph_r(k)*abs(sph_lat(j2) - sph_lat(j1))

    end function arc_lat

    !
    !> Set up the geometry of the global grid. Called once at initialization;
    !> the arguments are the longitude, geodetic latitude and depth axes.
    !
    subroutine init_spherical_grid(nlon, nlat, ndepth, olon, olat, odepth, dlon, dlat, ddepth)

        integer, intent(in) :: nlon, nlat, ndepth
        double precision, intent(in) :: olon, olat, odepth, dlon, dlat, ddepth

        integer :: i, j, k
        double precision :: lon

        call sph_check(earth_radius > 0, 'earth_radius must be > 0')
        call sph_check(nlon >= 2 .and. nlat >= 2 .and. ndepth >= 2, &
            'nlon, nlat and ndepth must each be >= 2')
        call sph_check(dlon > 0 .and. dlat > 0 .and. ddepth > 0, &
            'dlon, dlat and ddepth must be > 0')
        call sph_check(nlon*dlon <= 360.0d0, &
            'the longitude cells cover more than 360 degrees: nlon*dlon must be <= 360')
        call sph_check(olat - 0.5d0*dlat > -90.0d0 .and. olat + (nlat - 0.5d0)*dlat < 90.0d0, &
            'the latitude cells, which extend half a grid interval beyond the first and last rows, ' &
            //'must lie strictly between -90 and 90 degrees')
        call sph_check(odepth + (ndepth - 0.5d0)*ddepth < earth_radius, &
            'the deepest cell, which extends half a grid interval below the last level, ' &
            //'must lie above the centre of the Earth')

        sph_nlon = nlon
        sph_nlat = nlat
        sph_ndepth = ndepth
        sph_olon = olon
        sph_olat = olat
        sph_odepth = odepth
        sph_dlon = dlon
        sph_dlat = dlat
        sph_ddepth = ddepth
        sph_dlon_rad = dlon*deg2rad

        allocate (sph_coslon(1:nlon), sph_sinlon(1:nlon))
        do i = 1, nlon
            lon = olon + (i - 1)*dlon
            sph_coslon(i) = cos(lon*deg2rad)
            sph_sinlon(i) = sin(lon*deg2rad)
        end do

        allocate (sph_lat(1:nlat), sph_coslat(1:nlat), sph_sinlat(1:nlat))
        do j = 1, nlat
            sph_lat(j) = geodetic_to_geocentric(olat + (j - 1)*dlat)
        end do
        sph_coslat = cos(sph_lat)
        sph_sinlat = sin(sph_lat)

        allocate (sph_lat_face(1:nlat + 1), sph_coslat_face(1:nlat + 1), sph_sinlat_face(1:nlat + 1))
        do j = 1, nlat + 1
            sph_lat_face(j) = geodetic_to_geocentric(olat + (j - 1.5d0)*dlat)
        end do
        sph_coslat_face = cos(sph_lat_face)
        sph_sinlat_face = sin(sph_lat_face)

        allocate (sph_r(1:ndepth), sph_r_face(1:ndepth + 1))
        do k = 1, ndepth
            sph_r(k) = earth_radius - (odepth + (k - 1)*ddepth)
        end do
        do k = 1, ndepth + 1
            sph_r_face(k) = earth_radius - (odepth + (k - 1.5d0)*ddepth)
        end do

    end subroutine init_spherical_grid

    !
    !> The physical cell size at the reference point of the grid: the centre of
    !> the grid in longitude and latitude, at the top level. Returns the
    !> east-west, north-south and vertical sizes, in the length unit of the
    !> model, of one grid interval there. Used to convert the window lengths of
    !> the filters that work on a fixed number of cells.
    !
    function reference_spacing() result(h)

        double precision, dimension(1:3) :: h
        double precision :: latg, r

        latg = sph_olat + 0.5d0*(sph_nlat - 1)*sph_dlat
        r = sph_r(1)
        h(1) = r*cos(geodetic_to_geocentric(latg))*sph_dlon_rad
        h(2) = r*dgeocentric_dgeodetic(latg)*sph_dlat*deg2rad
        h(3) = sph_ddepth

    end function reference_spacing

    !
    !> The smallest length, in the length unit of the model, of one degree of
    !> longitude and of one degree of geodetic latitude anywhere in the grid,
    !> i.e., at the deepest level and, for longitude, at the row farthest from
    !> the equator. Converting a length to degrees with these never gives fewer
    !> degrees than the length needs anywhere in the grid.
    !
    function min_length_per_degree() result(m)

        double precision, dimension(1:2) :: m
        double precision :: latg_lo, latg_hi, rmin

        rmin = sph_r(sph_ndepth)
        latg_lo = sph_olat
        latg_hi = sph_olat + (sph_nlat - 1)*sph_dlat
        ! Longitude: r cos(theta) is smallest at the row with the largest |theta|
        m(1) = rmin*min(cos(geodetic_to_geocentric(latg_lo)), cos(geodetic_to_geocentric(latg_hi)))*deg2rad
        ! Latitude: d(theta)/d(lat) is smallest at the row nearest the equator
        if (latg_lo <= 0 .and. latg_hi >= 0) then
            m(2) = rmin*dgeocentric_dgeodetic(0.0d0)*deg2rad
        else
            m(2) = rmin*min(dgeocentric_dgeodetic(latg_lo), dgeocentric_dgeodetic(latg_hi))*deg2rad
        end if

    end function min_length_per_degree

    !
    !> Stop with an error message unless a condition holds
    !
    subroutine sph_check(condition, message)

        logical, intent(in) :: condition
        character(len=*), intent(in) :: message

        if (.not. condition) then
            call warn(' <spherical> Error: '//message)
            stop
        end if

    end subroutine sph_check

end module spherical
