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
! Ray Path Tracing Module for LATTE Traveltime Tomography
!
! This module implements ray path tracing using steepest descent on the
! traveltime gradient field. Rays are traced backward from receivers to sources
! following the negative gradient direction.
!
! Features:
!   - Trilinear and nearest-neighbor gradient interpolation
!   - Adaptive step sizing (optional)
!   - Hybrid mode: blend gradient and direct directions (optional)
!   - Momentum mode: smooth oscillations (optional)
!   - OpenMP parallelization for multiple rays
!   - MPI-compatible (each rank traces rays for its assigned shots)
!
! Author:
!    Claude Code Assistant (based on Python postprocess implementation)
!

module raypath

    use libflit
    use parameters

    implicit none

    private

    ! Public subroutines
    public :: trace_shot_ray_paths
    public :: write_ray_paths_binary
    public :: initialize_raypath_params
    public :: raypath_needs_gradients

    ! Module-level raypath parameters (set via initialize_raypath_params)
    ! These use raypath_ prefix to avoid clashes with existing variables
    logical, public :: raypath_save = .false.
    character(len=1024), public :: raypath_dir = './ray_paths'
    real, public :: raypath_step_size = -1.0     ! Negative means auto-compute
    real, public :: raypath_tolerance = -1.0     ! Negative means auto-compute
    integer, public :: raypath_max_steps = 50000
    character(len=24), public :: raypath_interpolation = 'trilinear'
    logical, public :: raypath_adaptive_step = .false.
    logical, public :: raypath_hybrid_mode = .false.
    logical, public :: raypath_momentum_mode = .false.
    real, public :: raypath_min_gradient = 1.0e-10
    real, public :: raypath_hybrid_alpha = 0.7
    real, public :: raypath_momentum_beta = 0.9

    ! Constants for ray tracing
    real, parameter :: RAYPATH_MIN_STEP_FACTOR = 0.1
    real, parameter :: RAYPATH_MAX_STEP_FACTOR = 2.0

    ! Ray path termination status codes (use RAYSTATUS_ prefix to avoid name collision)
    integer, parameter, public :: RAYSTATUS_SUCCESS = 0
    integer, parameter, public :: RAYSTATUS_MAX_STEPS = 1
    integer, parameter, public :: RAYSTATUS_OUT_OF_BOUNDS = 2
    integer, parameter, public :: RAYSTATUS_GRADIENT_TOO_SMALL = 3
    integer, parameter, public :: RAYSTATUS_OSCILLATING = 4

contains

    !
    !> Initialize raypath parameters from parameter file
    !> Must be called after file_parameter is set
    !
    subroutine initialize_raypath_params

        call readpar_logical(file_parameter, 'raypath_save', raypath_save, .false.)
        call readpar_string(file_parameter, 'raypath_dir', raypath_dir, './ray_paths')
        call readpar_float(file_parameter, 'raypath_step_size', raypath_step_size, -1.0)
        call readpar_float(file_parameter, 'raypath_tolerance', raypath_tolerance, -1.0)
        call readpar_int(file_parameter, 'raypath_max_steps', raypath_max_steps, 50000)
        call readpar_string(file_parameter, 'raypath_interpolation', raypath_interpolation, 'trilinear')
        call readpar_logical(file_parameter, 'raypath_adaptive_step', raypath_adaptive_step, .false.)
        call readpar_logical(file_parameter, 'raypath_hybrid_mode', raypath_hybrid_mode, .false.)
        call readpar_logical(file_parameter, 'raypath_momentum_mode', raypath_momentum_mode, .false.)
        call readpar_float(file_parameter, 'raypath_min_gradient', raypath_min_gradient, 1.0e-10)
        call readpar_float(file_parameter, 'raypath_hybrid_alpha', raypath_hybrid_alpha, 0.7)
        call readpar_float(file_parameter, 'raypath_momentum_beta', raypath_momentum_beta, 0.9)

        ! Validate interpolation method
        if (raypath_interpolation /= 'trilinear' .and. raypath_interpolation /= 'nearest') then
            call warn(' <initialize_raypath_params> Warning: Invalid raypath_interpolation = '// &
                tidy(raypath_interpolation)//'. Using trilinear.')
            raypath_interpolation = 'trilinear'
        end if

    end subroutine initialize_raypath_params

    !
    !> Check if gradients are needed for ray path computation
    !> Returns true if raypath_save is true (we need gradients to trace rays)
    !
    function raypath_needs_gradients() result(needs)
        logical :: needs
        needs = raypath_save
    end function raypath_needs_gradients

#ifdef dim2

    !
    !> Trace ray paths for all receivers in a shot (2D version)
    !> Returns arrays of ray path coordinates
    !
    !> @param[in] pdx       Gradient field in x direction (nz, nx)
    !> @param[in] pdz       Gradient field in z direction (nz, nx)
    !> @param[in] d         Grid spacing [dx, dz]
    !> @param[in] o         Grid origin [ox, oz]
    !> @param[in] geom      Source-receiver geometry
    !> @param[out] ray_npoints  Number of points in each ray (nrays)
    !> @param[out] ray_coords   All ray coordinates concatenated (total_points, 2)
    !> @param[out] ray_status   Termination status for each ray (nrays)
    !
    subroutine trace_shot_ray_paths(pdx, pdz, d, o, geom, ray_npoints, ray_coords, ray_status)

        real, dimension(:, :), intent(in) :: pdx, pdz
        real, dimension(1:2), intent(in) :: d, o
        type(source_receiver_geometry), intent(in) :: geom
        integer, allocatable, dimension(:), intent(out) :: ray_npoints
        real, allocatable, dimension(:, :), intent(out) :: ray_coords
        integer, allocatable, dimension(:), intent(out) :: ray_status

        ! Local variables
        integer :: nr, ir, npts_total, max_pts_per_ray
        integer :: local_nz, local_nx
        real :: step_size, tolerance
        real :: src_x, src_z
        real, allocatable, dimension(:, :, :) :: temp_paths  ! (max_pts, 2, nr)
        integer, allocatable, dimension(:) :: temp_npoints
        integer, allocatable, dimension(:) :: temp_status
        integer :: ipt, total_pts

        ! Get grid dimensions
        local_nz = size(pdx, 1)
        local_nx = size(pdx, 2)

        ! Compute step size and tolerance if not specified
        step_size = raypath_step_size
        if (step_size < 0.0) then
            step_size = 0.5 * min(d(1), d(2))
        end if

        tolerance = raypath_tolerance
        if (tolerance < 0.0) then
            tolerance = 1.0 * max(d(1), d(2))
        end if

        ! Get source position (use first source if multiple)
        src_x = geom%srcr(1)%x
        src_z = geom%srcr(1)%z

        ! Number of receivers
        nr = geom%nr

        ! Allocate temporary storage for ray paths
        max_pts_per_ray = raypath_max_steps + 1
        allocate(temp_paths(max_pts_per_ray, 2, nr))
        allocate(temp_npoints(nr))
        allocate(temp_status(nr))

        temp_paths = 0.0
        temp_npoints = 0
        temp_status = RAYSTATUS_MAX_STEPS

        ! Trace rays for each receiver (OpenMP parallel)
        !$omp parallel do private(ir) schedule(dynamic)
        do ir = 1, nr
            ! Skip receivers with zero weight
            if (geom%recr(ir)%weight /= 0.0) then
                call trace_single_ray_2d(pdx, pdz, d, o, local_nx, local_nz, &
                    geom%recr(ir)%x, geom%recr(ir)%z, src_x, src_z, &
                    step_size, tolerance, &
                    temp_paths(:, :, ir), temp_npoints(ir), temp_status(ir))
            end if
        end do
        !$omp end parallel do

        ! Count total points and allocate output arrays
        npts_total = sum(temp_npoints)

        allocate(ray_npoints(nr))
        allocate(ray_status(nr))
        allocate(ray_coords(npts_total, 2))

        ray_npoints = temp_npoints
        ray_status = temp_status

        ! Copy ray coordinates to output array
        total_pts = 0
        do ir = 1, nr
            do ipt = 1, temp_npoints(ir)
                total_pts = total_pts + 1
                ray_coords(total_pts, 1) = temp_paths(ipt, 1, ir)
                ray_coords(total_pts, 2) = temp_paths(ipt, 2, ir)
            end do
        end do

        ! Clean up
        deallocate(temp_paths, temp_npoints, temp_status)

    end subroutine trace_shot_ray_paths

    !
    !> Trace a single ray from receiver to source (2D version)
    !
    subroutine trace_single_ray_2d(pdx, pdz, d, o, local_nx, local_nz, &
            recv_x, recv_z, src_x, src_z, step_size, tolerance, &
            ray_path, npoints, status)

        real, dimension(:, :), intent(in) :: pdx, pdz
        real, dimension(1:2), intent(in) :: d, o
        integer, intent(in) :: local_nx, local_nz
        real, intent(in) :: recv_x, recv_z, src_x, src_z
        real, intent(in) :: step_size, tolerance
        real, dimension(:, :), intent(out) :: ray_path
        integer, intent(out) :: npoints
        integer, intent(out) :: status

        ! Local variables
        real :: curr_x, curr_z
        real :: gx, gz, grad_mag
        real :: dir_x, dir_z
        real :: dist_to_src, prev_dist
        real :: actual_step
        real :: vel_x, vel_z  ! For momentum mode
        integer :: istep, oscillation_count
        integer :: max_pts

        max_pts = size(ray_path, 1)

        ! Initialize
        curr_x = recv_x
        curr_z = recv_z
        npoints = 0
        status = RAYSTATUS_MAX_STEPS
        oscillation_count = 0
        prev_dist = huge(1.0)
        vel_x = 0.0
        vel_z = 0.0

        ! Main ray tracing loop
        do istep = 1, raypath_max_steps

            ! Store current position
            npoints = npoints + 1
            if (npoints > max_pts) then
                npoints = max_pts
                exit
            end if
            ray_path(npoints, 1) = curr_x
            ray_path(npoints, 2) = curr_z

            ! Check distance to source
            dist_to_src = sqrt((curr_x - src_x)**2 + (curr_z - src_z)**2)

            if (dist_to_src < tolerance) then
                status = RAYSTATUS_SUCCESS
                exit
            end if

            ! Check for oscillation (not making progress)
            if (dist_to_src >= prev_dist) then
                oscillation_count = oscillation_count + 1
                if (oscillation_count > 100) then
                    status = RAYSTATUS_OSCILLATING
                    exit
                end if
            else
                oscillation_count = 0
            end if
            prev_dist = dist_to_src

            ! Interpolate gradient at current position
            call interpolate_gradient_2d(pdx, pdz, d, o, local_nx, local_nz, &
                curr_x, curr_z, gx, gz, status)

            if (status == RAYSTATUS_OUT_OF_BOUNDS) then
                exit
            end if

            ! Compute gradient magnitude
            grad_mag = sqrt(gx**2 + gz**2)

            ! Check for very small gradient
            if (grad_mag < raypath_min_gradient) then
                ! Use direct fallback to source
                gx = src_x - curr_x
                gz = src_z - curr_z
                grad_mag = sqrt(gx**2 + gz**2)
                if (grad_mag < 1.0e-15) then
                    status = RAYSTATUS_GRADIENT_TOO_SMALL
                    exit
                end if
            end if

            ! Compute ray direction (negative gradient direction)
            dir_x = -gx / grad_mag
            dir_z = -gz / grad_mag

            ! Apply hybrid mode if enabled (blend with direct direction)
            if (raypath_hybrid_mode) then
                call apply_hybrid_direction_2d(dir_x, dir_z, curr_x, curr_z, &
                    src_x, src_z, dist_to_src, grad_mag)
            end if

            ! Compute actual step size
            actual_step = step_size
            if (raypath_adaptive_step) then
                actual_step = compute_adaptive_step(step_size, grad_mag, d(1), d(2))
            end if

            ! Apply momentum mode if enabled
            if (raypath_momentum_mode) then
                vel_x = raypath_momentum_beta * vel_x + (1.0 - raypath_momentum_beta) * dir_x
                vel_z = raypath_momentum_beta * vel_z + (1.0 - raypath_momentum_beta) * dir_z
                grad_mag = sqrt(vel_x**2 + vel_z**2)
                if (grad_mag > 1.0e-15) then
                    dir_x = vel_x / grad_mag
                    dir_z = vel_z / grad_mag
                end if
            end if

            ! Update position
            curr_x = curr_x + actual_step * dir_x
            curr_z = curr_z + actual_step * dir_z

        end do

    end subroutine trace_single_ray_2d

    !
    !> Interpolate gradient at a given position (2D version)
    !
    subroutine interpolate_gradient_2d(pdx, pdz, d, o, local_nx, local_nz, &
            px, pz, gx, gz, status)

        real, dimension(:, :), intent(in) :: pdx, pdz
        real, dimension(1:2), intent(in) :: d, o
        integer, intent(in) :: local_nx, local_nz
        real, intent(in) :: px, pz
        real, intent(out) :: gx, gz
        integer, intent(inout) :: status

        real :: x_coord, z_coord
        integer :: i0, i1, k0, k1
        real :: wx, wz
        real :: w00, w01, w10, w11

        ! Convert to grid coordinates
        x_coord = (px - o(1)) / d(1)
        z_coord = (pz - o(2)) / d(2)

        ! Check bounds
        if (x_coord < -0.5 .or. x_coord > local_nx - 0.5 .or. &
            z_coord < -0.5 .or. z_coord > local_nz - 0.5) then
            status = RAYSTATUS_OUT_OF_BOUNDS
            gx = 0.0
            gz = 0.0
            return
        end if

        if (raypath_interpolation == 'nearest') then
            ! Nearest neighbor interpolation
            i0 = max(1, min(local_nx, nint(x_coord) + 1))
            k0 = max(1, min(local_nz, nint(z_coord) + 1))
            gx = pdx(k0, i0)
            gz = pdz(k0, i0)
        else
            ! Trilinear (bilinear in 2D) interpolation
            i0 = max(1, min(local_nx - 1, int(floor(x_coord)) + 1))
            k0 = max(1, min(local_nz - 1, int(floor(z_coord)) + 1))
            i1 = min(local_nx, i0 + 1)
            k1 = min(local_nz, k0 + 1)

            wx = x_coord - (i0 - 1)
            wz = z_coord - (k0 - 1)

            wx = max(0.0, min(1.0, wx))
            wz = max(0.0, min(1.0, wz))

            ! Bilinear weights
            w00 = (1.0 - wx) * (1.0 - wz)
            w01 = (1.0 - wx) * wz
            w10 = wx * (1.0 - wz)
            w11 = wx * wz

            gx = w00 * pdx(k0, i0) + w01 * pdx(k1, i0) + &
                 w10 * pdx(k0, i1) + w11 * pdx(k1, i1)
            gz = w00 * pdz(k0, i0) + w01 * pdz(k1, i0) + &
                 w10 * pdz(k0, i1) + w11 * pdz(k1, i1)
        end if

    end subroutine interpolate_gradient_2d

    !
    !> Apply hybrid direction blending (2D version)
    !
    subroutine apply_hybrid_direction_2d(dir_x, dir_z, curr_x, curr_z, &
            src_x, src_z, dist_to_src, grad_mag)

        real, intent(inout) :: dir_x, dir_z
        real, intent(in) :: curr_x, curr_z, src_x, src_z
        real, intent(in) :: dist_to_src, grad_mag

        real :: direct_x, direct_z, direct_mag
        real :: alpha, combined_x, combined_z, combined_mag

        ! Compute direct direction to source
        direct_x = src_x - curr_x
        direct_z = src_z - curr_z
        direct_mag = sqrt(direct_x**2 + direct_z**2)

        if (direct_mag < 1.0e-15) return

        direct_x = direct_x / direct_mag
        direct_z = direct_z / direct_mag

        ! Adaptive alpha based on distance and gradient strength
        alpha = raypath_hybrid_alpha
        if (dist_to_src > 100.0 * max(dx, dz)) then
            alpha = alpha * 0.5
        end if
        if (grad_mag < raypath_min_gradient * 10.0) then
            alpha = alpha * 0.3
        end if

        ! Blend directions
        combined_x = alpha * dir_x + (1.0 - alpha) * direct_x
        combined_z = alpha * dir_z + (1.0 - alpha) * direct_z
        combined_mag = sqrt(combined_x**2 + combined_z**2)

        if (combined_mag > 1.0e-15) then
            dir_x = combined_x / combined_mag
            dir_z = combined_z / combined_mag
        end if

    end subroutine apply_hybrid_direction_2d

    !
    !> Write ray paths to binary file (2D version)
    !> Format: [nrays][npoints_array][concatenated x,z coordinates]
    !
    subroutine write_ray_paths_binary(filename, ray_npoints, ray_coords)

        character(len=*), intent(in) :: filename
        integer, dimension(:), intent(in) :: ray_npoints
        real, dimension(:, :), intent(in) :: ray_coords

        integer :: funit, nrays, total_pts, ir, ipt, coord_idx

        nrays = size(ray_npoints)
        total_pts = size(ray_coords, 1)

        ! Open file for binary writing
        open(newunit=funit, file=filename, status='replace', &
             access='stream', form='unformatted')

        ! Write header: number of rays
        write(funit) nrays

        ! Write number of points per ray
        write(funit) ray_npoints

        ! Write coordinates (x, z for each point)
        coord_idx = 0
        do ir = 1, nrays
            do ipt = 1, ray_npoints(ir)
                coord_idx = coord_idx + 1
                write(funit) ray_coords(coord_idx, 1)  ! x
                write(funit) ray_coords(coord_idx, 2)  ! z
            end do
        end do

        close(funit)

    end subroutine write_ray_paths_binary

#endif

#ifdef dim3

    !
    !> Trace ray paths for all receivers in a shot (3D version)
    !
    subroutine trace_shot_ray_paths(pdx, pdy, pdz, d, o, geom, ray_npoints, ray_coords, ray_status)

        real, dimension(:, :, :), intent(in) :: pdx, pdy, pdz
        real, dimension(1:3), intent(in) :: d, o
        type(source_receiver_geometry), intent(in) :: geom
        integer, allocatable, dimension(:), intent(out) :: ray_npoints
        real, allocatable, dimension(:, :), intent(out) :: ray_coords
        integer, allocatable, dimension(:), intent(out) :: ray_status

        ! Local variables
        integer :: nr, ir, npts_total, max_pts_per_ray
        integer :: local_nz, local_ny, local_nx
        real :: step_size, tolerance
        real :: src_x, src_y, src_z
        real, allocatable, dimension(:, :, :) :: temp_paths  ! (max_pts, 3, nr)
        integer, allocatable, dimension(:) :: temp_npoints
        integer, allocatable, dimension(:) :: temp_status
        integer :: ipt, total_pts

        ! Get grid dimensions (pdx is nz, ny, nx)
        local_nz = size(pdx, 1)
        local_ny = size(pdx, 2)
        local_nx = size(pdx, 3)

        ! Compute step size and tolerance if not specified
        step_size = raypath_step_size
        if (step_size < 0.0) then
            step_size = 0.5 * min(d(1), d(2), d(3))
        end if

        tolerance = raypath_tolerance
        if (tolerance < 0.0) then
            tolerance = 1.0 * max(d(1), d(2), d(3))
        end if

        ! Get source position
        src_x = geom%srcr(1)%x
        src_y = geom%srcr(1)%y
        src_z = geom%srcr(1)%z

        ! Number of receivers
        nr = geom%nr

        ! Allocate temporary storage
        max_pts_per_ray = raypath_max_steps + 1
        allocate(temp_paths(max_pts_per_ray, 3, nr))
        allocate(temp_npoints(nr))
        allocate(temp_status(nr))

        temp_paths = 0.0
        temp_npoints = 0
        temp_status = RAYSTATUS_MAX_STEPS

        ! Trace rays for each receiver (OpenMP parallel)
        !$omp parallel do private(ir) schedule(dynamic)
        do ir = 1, nr
            if (geom%recr(ir)%weight /= 0.0) then
                call trace_single_ray_3d(pdx, pdy, pdz, d, o, &
                    local_nx, local_ny, local_nz, &
                    geom%recr(ir)%x, geom%recr(ir)%y, geom%recr(ir)%z, &
                    src_x, src_y, src_z, step_size, tolerance, &
                    temp_paths(:, :, ir), temp_npoints(ir), temp_status(ir))
            end if
        end do
        !$omp end parallel do

        ! Count total points
        npts_total = sum(temp_npoints)

        allocate(ray_npoints(nr))
        allocate(ray_status(nr))
        allocate(ray_coords(npts_total, 3))

        ray_npoints = temp_npoints
        ray_status = temp_status

        ! Copy ray coordinates
        total_pts = 0
        do ir = 1, nr
            do ipt = 1, temp_npoints(ir)
                total_pts = total_pts + 1
                ray_coords(total_pts, 1) = temp_paths(ipt, 1, ir)
                ray_coords(total_pts, 2) = temp_paths(ipt, 2, ir)
                ray_coords(total_pts, 3) = temp_paths(ipt, 3, ir)
            end do
        end do

        deallocate(temp_paths, temp_npoints, temp_status)

    end subroutine trace_shot_ray_paths

    !
    !> Trace a single ray from receiver to source (3D version)
    !
    subroutine trace_single_ray_3d(pdx, pdy, pdz, d, o, &
            local_nx, local_ny, local_nz, &
            recv_x, recv_y, recv_z, src_x, src_y, src_z, &
            step_size, tolerance, ray_path, npoints, status)

        real, dimension(:, :, :), intent(in) :: pdx, pdy, pdz
        real, dimension(1:3), intent(in) :: d, o
        integer, intent(in) :: local_nx, local_ny, local_nz
        real, intent(in) :: recv_x, recv_y, recv_z, src_x, src_y, src_z
        real, intent(in) :: step_size, tolerance
        real, dimension(:, :), intent(out) :: ray_path
        integer, intent(out) :: npoints
        integer, intent(out) :: status

        real :: curr_x, curr_y, curr_z
        real :: gx, gy, gz, grad_mag
        real :: dir_x, dir_y, dir_z
        real :: dist_to_src, prev_dist
        real :: actual_step
        real :: vel_x, vel_y, vel_z
        integer :: istep, oscillation_count
        integer :: max_pts

        max_pts = size(ray_path, 1)

        curr_x = recv_x
        curr_y = recv_y
        curr_z = recv_z
        npoints = 0
        status = RAYSTATUS_MAX_STEPS
        oscillation_count = 0
        prev_dist = huge(1.0)
        vel_x = 0.0
        vel_y = 0.0
        vel_z = 0.0

        do istep = 1, raypath_max_steps

            npoints = npoints + 1
            if (npoints > max_pts) then
                npoints = max_pts
                exit
            end if
            ray_path(npoints, 1) = curr_x
            ray_path(npoints, 2) = curr_y
            ray_path(npoints, 3) = curr_z

            dist_to_src = sqrt((curr_x - src_x)**2 + (curr_y - src_y)**2 + (curr_z - src_z)**2)

            if (dist_to_src < tolerance) then
                status = RAYSTATUS_SUCCESS
                exit
            end if

            if (dist_to_src >= prev_dist) then
                oscillation_count = oscillation_count + 1
                if (oscillation_count > 100) then
                    status = RAYSTATUS_OSCILLATING
                    exit
                end if
            else
                oscillation_count = 0
            end if
            prev_dist = dist_to_src

            call interpolate_gradient_3d(pdx, pdy, pdz, d, o, &
                local_nx, local_ny, local_nz, &
                curr_x, curr_y, curr_z, gx, gy, gz, status)

            if (status == RAYSTATUS_OUT_OF_BOUNDS) exit

            grad_mag = sqrt(gx**2 + gy**2 + gz**2)

            if (grad_mag < raypath_min_gradient) then
                gx = src_x - curr_x
                gy = src_y - curr_y
                gz = src_z - curr_z
                grad_mag = sqrt(gx**2 + gy**2 + gz**2)
                if (grad_mag < 1.0e-15) then
                    status = RAYSTATUS_GRADIENT_TOO_SMALL
                    exit
                end if
            end if

            dir_x = -gx / grad_mag
            dir_y = -gy / grad_mag
            dir_z = -gz / grad_mag

            if (raypath_hybrid_mode) then
                call apply_hybrid_direction_3d(dir_x, dir_y, dir_z, &
                    curr_x, curr_y, curr_z, src_x, src_y, src_z, &
                    dist_to_src, grad_mag)
            end if

            actual_step = step_size
            if (raypath_adaptive_step) then
                actual_step = compute_adaptive_step(step_size, grad_mag, d(1), d(2))
            end if

            if (raypath_momentum_mode) then
                vel_x = raypath_momentum_beta * vel_x + (1.0 - raypath_momentum_beta) * dir_x
                vel_y = raypath_momentum_beta * vel_y + (1.0 - raypath_momentum_beta) * dir_y
                vel_z = raypath_momentum_beta * vel_z + (1.0 - raypath_momentum_beta) * dir_z
                grad_mag = sqrt(vel_x**2 + vel_y**2 + vel_z**2)
                if (grad_mag > 1.0e-15) then
                    dir_x = vel_x / grad_mag
                    dir_y = vel_y / grad_mag
                    dir_z = vel_z / grad_mag
                end if
            end if

            curr_x = curr_x + actual_step * dir_x
            curr_y = curr_y + actual_step * dir_y
            curr_z = curr_z + actual_step * dir_z

        end do

    end subroutine trace_single_ray_3d

    !
    !> Interpolate gradient at a given position (3D version)
    !
    subroutine interpolate_gradient_3d(pdx, pdy, pdz, d, o, &
            local_nx, local_ny, local_nz, &
            px, py, pz, gx, gy, gz, status)

        real, dimension(:, :, :), intent(in) :: pdx, pdy, pdz
        real, dimension(1:3), intent(in) :: d, o
        integer, intent(in) :: local_nx, local_ny, local_nz
        real, intent(in) :: px, py, pz
        real, intent(out) :: gx, gy, gz
        integer, intent(inout) :: status

        real :: x_coord, y_coord, z_coord
        integer :: i0, i1, j0, j1, k0, k1
        real :: wx, wy, wz
        real :: w000, w001, w010, w011, w100, w101, w110, w111

        x_coord = (px - o(1)) / d(1)
        y_coord = (py - o(2)) / d(2)
        z_coord = (pz - o(3)) / d(3)

        if (x_coord < -0.5 .or. x_coord > local_nx - 0.5 .or. &
            y_coord < -0.5 .or. y_coord > local_ny - 0.5 .or. &
            z_coord < -0.5 .or. z_coord > local_nz - 0.5) then
            status = RAYSTATUS_OUT_OF_BOUNDS
            gx = 0.0
            gy = 0.0
            gz = 0.0
            return
        end if

        if (raypath_interpolation == 'nearest') then
            i0 = max(1, min(local_nx, nint(x_coord) + 1))
            j0 = max(1, min(local_ny, nint(y_coord) + 1))
            k0 = max(1, min(local_nz, nint(z_coord) + 1))
            gx = pdx(k0, j0, i0)
            gy = pdy(k0, j0, i0)
            gz = pdz(k0, j0, i0)
        else
            ! Trilinear interpolation
            i0 = max(1, min(local_nx - 1, int(floor(x_coord)) + 1))
            j0 = max(1, min(local_ny - 1, int(floor(y_coord)) + 1))
            k0 = max(1, min(local_nz - 1, int(floor(z_coord)) + 1))
            i1 = min(local_nx, i0 + 1)
            j1 = min(local_ny, j0 + 1)
            k1 = min(local_nz, k0 + 1)

            wx = x_coord - (i0 - 1)
            wy = y_coord - (j0 - 1)
            wz = z_coord - (k0 - 1)

            wx = max(0.0, min(1.0, wx))
            wy = max(0.0, min(1.0, wy))
            wz = max(0.0, min(1.0, wz))

            w000 = (1.0 - wx) * (1.0 - wy) * (1.0 - wz)
            w001 = (1.0 - wx) * (1.0 - wy) * wz
            w010 = (1.0 - wx) * wy * (1.0 - wz)
            w011 = (1.0 - wx) * wy * wz
            w100 = wx * (1.0 - wy) * (1.0 - wz)
            w101 = wx * (1.0 - wy) * wz
            w110 = wx * wy * (1.0 - wz)
            w111 = wx * wy * wz

            gx = w000 * pdx(k0, j0, i0) + w001 * pdx(k1, j0, i0) + &
                 w010 * pdx(k0, j1, i0) + w011 * pdx(k1, j1, i0) + &
                 w100 * pdx(k0, j0, i1) + w101 * pdx(k1, j0, i1) + &
                 w110 * pdx(k0, j1, i1) + w111 * pdx(k1, j1, i1)
            gy = w000 * pdy(k0, j0, i0) + w001 * pdy(k1, j0, i0) + &
                 w010 * pdy(k0, j1, i0) + w011 * pdy(k1, j1, i0) + &
                 w100 * pdy(k0, j0, i1) + w101 * pdy(k1, j0, i1) + &
                 w110 * pdy(k0, j1, i1) + w111 * pdy(k1, j1, i1)
            gz = w000 * pdz(k0, j0, i0) + w001 * pdz(k1, j0, i0) + &
                 w010 * pdz(k0, j1, i0) + w011 * pdz(k1, j1, i0) + &
                 w100 * pdz(k0, j0, i1) + w101 * pdz(k1, j0, i1) + &
                 w110 * pdz(k0, j1, i1) + w111 * pdz(k1, j1, i1)
        end if

    end subroutine interpolate_gradient_3d

    !
    !> Apply hybrid direction blending (3D version)
    !
    subroutine apply_hybrid_direction_3d(dir_x, dir_y, dir_z, &
            curr_x, curr_y, curr_z, src_x, src_y, src_z, &
            dist_to_src, grad_mag)

        real, intent(inout) :: dir_x, dir_y, dir_z
        real, intent(in) :: curr_x, curr_y, curr_z, src_x, src_y, src_z
        real, intent(in) :: dist_to_src, grad_mag

        real :: direct_x, direct_y, direct_z, direct_mag
        real :: alpha, combined_x, combined_y, combined_z, combined_mag

        direct_x = src_x - curr_x
        direct_y = src_y - curr_y
        direct_z = src_z - curr_z
        direct_mag = sqrt(direct_x**2 + direct_y**2 + direct_z**2)

        if (direct_mag < 1.0e-15) return

        direct_x = direct_x / direct_mag
        direct_y = direct_y / direct_mag
        direct_z = direct_z / direct_mag

        alpha = raypath_hybrid_alpha
        if (dist_to_src > 100.0 * max(dx, dy, dz)) then
            alpha = alpha * 0.5
        end if
        if (grad_mag < raypath_min_gradient * 10.0) then
            alpha = alpha * 0.3
        end if

        combined_x = alpha * dir_x + (1.0 - alpha) * direct_x
        combined_y = alpha * dir_y + (1.0 - alpha) * direct_y
        combined_z = alpha * dir_z + (1.0 - alpha) * direct_z
        combined_mag = sqrt(combined_x**2 + combined_y**2 + combined_z**2)

        if (combined_mag > 1.0e-15) then
            dir_x = combined_x / combined_mag
            dir_y = combined_y / combined_mag
            dir_z = combined_z / combined_mag
        end if

    end subroutine apply_hybrid_direction_3d

    !
    !> Write ray paths to binary file (3D version)
    !
    subroutine write_ray_paths_binary(filename, ray_npoints, ray_coords)

        character(len=*), intent(in) :: filename
        integer, dimension(:), intent(in) :: ray_npoints
        real, dimension(:, :), intent(in) :: ray_coords

        integer :: funit, nrays, total_pts, ir, ipt, coord_idx

        nrays = size(ray_npoints)
        total_pts = size(ray_coords, 1)

        open(newunit=funit, file=filename, status='replace', &
             access='stream', form='unformatted')

        write(funit) nrays
        write(funit) ray_npoints

        coord_idx = 0
        do ir = 1, nrays
            do ipt = 1, ray_npoints(ir)
                coord_idx = coord_idx + 1
                write(funit) ray_coords(coord_idx, 1)  ! x
                write(funit) ray_coords(coord_idx, 2)  ! y
                write(funit) ray_coords(coord_idx, 3)  ! z
            end do
        end do

        close(funit)

    end subroutine write_ray_paths_binary

#endif

    !
    !> Compute adaptive step size based on gradient magnitude
    !> Step is inversely proportional to gradient magnitude
    !
    function compute_adaptive_step(base_step, grad_mag, d1, d2) result(step)

        real, intent(in) :: base_step, grad_mag, d1, d2
        real :: step

        real :: min_step, max_step

        min_step = RAYPATH_MIN_STEP_FACTOR * min(d1, d2)
        max_step = RAYPATH_MAX_STEP_FACTOR * max(d1, d2)

        if (grad_mag > 1.0e-6) then
            step = base_step / grad_mag
        else
            step = max_step
        end if

        ! Clamp to valid range
        step = max(min_step, min(max_step, step))

    end function compute_adaptive_step

end module raypath
