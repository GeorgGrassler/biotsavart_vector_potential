module biotsavart
use, intrinsic :: iso_fortran_env, only: dp => real64
implicit none

real(dp), parameter :: clight = 2.99792458d10


!> Filamentary coils as a polygonal chain. Current zero marks end of one coil.
type coils_t
    real(dp), dimension(:), allocatable :: x, y, z, current
end type coils_t


contains


subroutine init_coils(x, y, z, current, coils)
    real(dp), intent(in) :: x(:), y(:), z(:), current(:)
    type(coils_t), intent(out) :: coils

    integer :: n_points

    n_points = size(x)
    call allocate_coils_data(coils, n_points)

    coils%x = x
    coils%y = y
    coils%z = z
    coils%current = current
end subroutine init_coils


subroutine deinit_coils(coils)
    type(coils_t), intent(inout) :: coils

    call deallocate_coils_data(coils)
end subroutine deinit_coils


subroutine load_coils_from_file(filename, coils)
    character(*), intent(in) :: filename
    type(coils_t), intent(out) :: coils

    integer :: unit
    integer :: i, n_points

    open(newunit=unit, file=filename, status="old", action="read")
    read(unit, *) n_points
    call allocate_coils_data(coils, n_points)
    do i = 1, n_points
        read(unit, *) coils%x(i), coils%y(i), coils%z(i), coils%current(i)
    end do
    close(unit)
end subroutine load_coils_from_file


subroutine save_coils_to_file(filename, coils)
    character(*), intent(in) :: filename
    type(coils_t), intent(in) :: coils

    integer :: unit
    integer :: i, n_points

    n_points = size(coils%x)
    open(newunit=unit, file=filename, status="replace", action="write")
    write(unit, *) n_points
    do i = 1, n_points
        write(unit, *) coils%x(i), coils%y(i), coils%z(i), coils%current(i)
    end do
    close(unit)
end subroutine save_coils_to_file


subroutine allocate_coils_data(coils, n_points)
    type(coils_t), intent(out) :: coils
    integer, intent(in) :: n_points

    allocate(coils%x(n_points), coils%y(n_points), coils%z(n_points))
    allocate(coils%current(n_points))
end subroutine allocate_coils_data


subroutine deallocate_coils_data(coils)
    type(coils_t), intent(inout) :: coils

    deallocate(coils%x, coils%y, coils%z, coils%current)
end subroutine deallocate_coils_data


subroutine equalize_coils_segments_lengths(coils)
    type(coils_t), intent(inout) :: coils

    real(dp), dimension((size(coils%x) - 1)) :: lengths
    real(dp) :: min_length
    integer, dimension((size(coils%x) - 1)) :: cuts_per_knot

    lengths = compute_coils_segments_lengths(coils)
    min_length = minval(lengths)
    cuts_per_knot = max(floor(lengths / min_length - 1), 0)
    call cut_coils_segments(coils, cuts_per_knot)

end subroutine equalize_coils_segments_lengths


function compute_coils_segments_lengths(coils) result(lenghts)
    type(coils_t), intent(in) :: coils

    real(dp), dimension(size(coils%x) - 1) :: lenghts
    integer :: segment

    do segment = 1, size(coils%x) - 1

        lenghts(segment) = sqrt((coils%x(segment+1) - coils%x(segment))**2 + &
                                (coils%y(segment+1) - coils%y(segment))**2 + &
                                (coils%z(segment+1) - coils%z(segment))**2)
    end do
end function compute_coils_segments_lengths


subroutine cut_coils_segments(coils, cuts_per_knot)
    type(coils_t), intent(inout) :: coils
    integer, dimension(:), intent(in) :: cuts_per_knot
    
    integer :: total_knots
    real(dp), dimension(:), allocatable :: x, y, z, current
    integer :: knot, subknot, global_knot_index
    real(dp), dimension(3) :: subknot_xyz

    total_knots = sum(cuts_per_knot) + size(coils%x)
    allocate(x(total_knots), y(total_knots), &
             z(total_knots), current(total_knots))
    global_knot_index = 0
    do knot = 1, size(cuts_per_knot)
        do subknot = 0, cuts_per_knot(knot)
            subknot_xyz = calc_subknot_xyz(coils, knot, subknot, cuts_per_knot(knot))
            global_knot_index = global_knot_index + 1
            x(global_knot_index) = subknot_xyz(1)
            y(global_knot_index) = subknot_xyz(2)
            z(global_knot_index) = subknot_xyz(3)
            current(global_knot_index) = coils%current(knot)
        end do
    end do
    
    x(total_knots) = coils%x(size(coils%x))
    y(total_knots) = coils%y(size(coils%y))
    z(total_knots) = coils%z(size(coils%z))
    current(total_knots) = 0.0d0
    call deinit_coils(coils)
    call init_coils(x, y, z, current, coils)
    deallocate(x, y, z, current)
end subroutine cut_coils_segments

function calc_subknot_xyz(coils, knot, subknot, cuts_per_knot) result(xyz)
    type(coils_t), intent(in) :: coils
    integer, intent(in) :: knot, subknot
    integer, intent(in) :: cuts_per_knot

    real(dp) :: ratio
    real(dp) :: xyz(3)

    ratio = 1.0d0 / (cuts_per_knot + 1) * subknot
    xyz(1) = coils%x(knot) + (coils%x(knot+1) - coils%x(knot)) * ratio
    xyz(2) = coils%y(knot) + (coils%y(knot+1) - coils%y(knot)) * ratio
    xyz(3) = coils%z(knot) + (coils%z(knot+1) - coils%z(knot)) * ratio
end function calc_subknot_xyz

!> Formula of Hanson and Hirshman (2002)
function compute_vector_potential(coils, x) result(A)
    type(coils_t), intent(in) :: coils
    real(dp), intent(in) :: x(3)

    real(dp) :: A(3), dx_i(3), dx_f(3), dl(3), R_i, R_f, L, eps, log_term
    integer :: i

    A = 0.0d0
    do i = 1, size(coils%x) - 1
        dl(1) = coils%x(i+1) - coils%x(i)
        dl(2) = coils%y(i+1) - coils%y(i)
        dl(3) = coils%z(i+1) - coils%z(i)
        dx_i(1) = x(1) - coils%x(i)
        dx_i(2) = x(2) - coils%y(i)
        dx_i(3) = x(3) - coils%z(i)
        dx_f(1) = x(1) - coils%x(i+1)
        dx_f(2) = x(2) - coils%y(i+1)
        dx_f(3) = x(3) - coils%z(i+1)
        R_i = sqrt(dx_i(1)**2 + dx_i(2)**2 + dx_i(3)**2)
        R_f = sqrt(dx_f(1)**2 + dx_f(2)**2 + dx_f(3)**2)
        L = sqrt(dl(1)**2 + dl(2)**2 + dl(3)**2)
        eps = L / (R_i + R_f)
        log_term = log((1.0d0 + eps) / (1.0d0 - eps))
        A = A + (dl / (clight*L)) * log_term
    end do
end function compute_vector_potential


end module biotsavart
