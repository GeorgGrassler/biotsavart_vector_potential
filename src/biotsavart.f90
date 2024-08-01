module biotsavart
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none

    real(dp), parameter :: clight = 2.99792458d10

    type CoilsData
        real(dp), dimension(:), allocatable :: x, y, z, current
    end type CoilsData

    contains


    subroutine init_coils_data(x, y, z, current, coils)
        real(dp), intent(in) :: x(:), y(:), z(:), current(:)
        type(CoilsData), intent(out) :: coils

        integer :: n_points

        n_points = size(x)
        call allocate_coils_data(coils, n_points)

        coils%x = x
        coils%y = y
        coils%z = z
        coils%current = current
    end subroutine init_coils_data


    subroutine deinit_coils_data(coils)
        type(CoilsData), intent(inout) :: coils

        call deallocate_coils_data(coils)
    end subroutine deinit_coils_data


    subroutine load_coils_file(filename, coils)
        character(*), intent(in) :: filename
        type(CoilsData), intent(out) :: coils

        integer :: unit
        integer :: i, n_points

        open(newunit=unit, file=filename, status="old", action="read")
        read(unit, *) n_points
        call allocate_coils_data(coils, n_points)
        do i = 1, n_points
            read(unit, *) coils%x(i), coils%y(i), coils%z(i), coils%current(i)
        end do
        close(unit)

    end subroutine load_coils_file


    subroutine save_coils_file(filename, coils)
        character(*), intent(in) :: filename
        type(CoilsData), intent(in) :: coils

        integer :: unit
        integer :: i, n_points

        n_points = size(coils%x)
        open(newunit=unit, file=filename, status="replace", action="write")
        write(unit, *) n_points
        do i = 1, n_points
            write(unit, *) coils%x(i), coils%y(i), coils%z(i), coils%current(i)
        end do
        close(unit)
    end subroutine save_coils_file


    subroutine allocate_coils_data(coils, n_points)
        type(CoilsData), intent(out) :: coils
        integer, intent(in) :: n_points

        allocate(coils%x(n_points), coils%y(n_points), coils%z(n_points))
        allocate(coils%current(n_points))
    end subroutine allocate_coils_data


    subroutine deallocate_coils_data(coils)
        type(CoilsData), intent(inout) :: coils

        deallocate(coils%x, coils%y, coils%z, coils%current)
    end subroutine deallocate_coils_data


    function compute_vector_potential(coils, x) result(A)
        ! Formula of Hanson and Hirshman (2002)
        type(CoilsData), intent(in) :: coils
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