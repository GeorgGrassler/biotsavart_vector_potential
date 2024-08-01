program test_biotsavart
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use test_util, only: print_test, print_ok, print_fail

    implicit none

    character(*), parameter :: test_coils_file = "coils.test"

    real(dp), parameter :: large_distance = 1.0d3

    call test_load_coils_file
    call test_compute_coils_segment_lengths
    call test_equalize_coils_lenghts
    call test_compute_vector_potential

    contains

    subroutine test_load_coils_file
        use biotsavart, only: CoilsData, load_coils_from_file, deinit_coils

        type(CoilsData) :: coils

        call print_test("load_coils_file")

        call create_test_coils_file
        call load_coils_from_file(test_coils_file, coils)
        call remove_test_coils_file

        if (size(coils%x) /= 4) then
            print *, "Coil length mismatch"
            print *, "len(coils%x) = ", size(coils%x)
            call print_fail
            error stop
        end if

        call deinit_coils(coils)

        call print_ok
    end subroutine test_load_coils_file


    subroutine test_compute_coils_segment_lengths
        use biotsavart, only: CoilsData, init_coils, & 
                              compute_coils_segment_lengths, deinit_coils

        real(dp), parameter :: tol = 1.0e-9
        integer, parameter :: N_KNOTS = 5

        type(CoilsData) :: coils
        real(dp), dimension(:), allocatable :: lengths
        real(dp), dimension(N_KNOTS) :: x, y, z, current, expected_lengths
        integer :: i

        call print_test("compute_coils_segment_lengths")

        x = [-1.0d0, 0.0d0, 1.0d0, 0.0d0, -1.0d0]
        y = [1.0d0, 0.0d0, -1.0d0, 0.0d0, 1.0d0]
        z = [0.0d0, 1.0d0, 0.0d0, -1.0d0, 0.0d0]
        current = [0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0]
        call init_coils(x, y, z, current, coils)

        lengths = compute_coils_segment_lengths(coils)
        expected_lengths = [sqrt(3.0d0), sqrt(3.0d0), sqrt(3.0d0), sqrt(3.0d0), sqrt(3.0d0)]

        do i = 1, 3
            if (abs(lengths(i) - expected_lengths(i)) > tol) then
                print *, "lengths(i) = ", lengths(i)
                print *, "expected_lengths(i) = ", expected_lengths(i)
                call print_fail
                error stop
            end if
        end do

        call deinit_coils(coils)

        call print_ok
    end subroutine test_compute_coils_segment_lengths

    subroutine test_equalize_coils_lenghts
        use biotsavart, only: CoilsData, init_coils, equalize_coils_lenghts, &
                              deinit_coils

        type(CoilsData) :: coils

        call print_test("equalize_coils_lenghts")

        call init_straight_wire_coils(coils)
        call equalize_coils_lenghts(coils)

        call deinit_coils(coils)

        call print_ok
    end subroutine test_equalize_coils_lenghts

    subroutine test_compute_vector_potential
        use biotsavart, only: CoilsData, compute_vector_potential, &
            deinit_coils, clight

        real(dp), parameter :: tol = 1.0e-9
        integer, parameter :: N_TEST = 3

        type(CoilsData) :: coils
        real(dp) :: x_test(3, N_TEST)
        real(dp), dimension(3) :: x, A, A_analytic, x_too_close
        integer :: i

        call print_test("compute_vector_potential")

        x_test(:, 1) = [0.4, 0.3, 0.8]
        x_test(:, 2) = [0.0, 0.2, -0.3]
        x_test(:, 3) = [1.0d2, -1.0d2, 1.0d2]

        call init_straight_wire_coils(coils)

        do i = 1, N_TEST
            x = x_test(:, i)
            A_analytic = vector_potential_straight_wire(x, large_distance, 1.0d0)
            A = compute_vector_potential(coils, x)
            if (any(abs(A - A_analytic)*clight > tol)) then
                print *, "A = ", A(3)*clight
                print *, "A_analytic = ", A_analytic(3)*clight
                print *, "Ratio = ", A(3) / A_analytic(3)
                call print_fail
                error stop
            end if
        end do

        x_too_close = [1.0e-5, 0.0, 0.0]
        A_analytic = vector_potential_straight_wire(x_too_close, large_distance, 1.0d0)
        A = compute_vector_potential(coils, x_too_close)
        if (all(abs(A - A_analytic)*clight < tol)) then
            call print_fail
            error stop
        end if

        call deinit_coils(coils)

        call print_ok
    end subroutine test_compute_vector_potential

    function vector_potential_straight_wire(x, L, current) result(A)
        use biotsavart, only: clight

        real(dp), dimension(3), intent(in) :: x
        real(dp), intent(in) :: L, current

        real(dp), dimension(3) :: A
        real(dp) :: R, z, L_half, A_z

        R = Rcyl(x)
        z = x(3)
        L_half = L/2.0d0
        A_z = current/clight*log((L_half - z + sqrt((L_half - z)**2 + R**2)) / &
                                 (-L_half - z + sqrt((-L_half - z)**2 + R**2)))
        A = [0.0d0, 0.0d0, A_z]
    end function vector_potential_straight_wire

    function Rcyl(x)
        real(dp), dimension(3), intent(in) :: x
        real(dp) :: Rcyl

        Rcyl = sqrt(x(1)**2 + x(2)**2)
    end function Rcyl

    subroutine create_test_coils_file
        use biotsavart, only: CoilsData, save_coils_to_file

        type(CoilsData) :: coils

        call init_straight_wire_coils(coils)
        call save_coils_to_file(test_coils_file, coils)
    end subroutine create_test_coils_file

    subroutine init_straight_wire_coils(coils)
        use biotsavart, only: CoilsData, init_coils

        type(CoilsData), intent(out) :: coils

        real(dp), dimension(4) :: x, y, z, current
        real(dp)               :: L

        x = [0.0d0, 0.0d0, 0.0d0, 0.0d0]
        y = [0.0d0, 0.0d0, 0.0d0, 0.0d0]
        L = large_distance
        z = [-L/2.0d0, -L/4.0d0, 0.0d0, L/2.0d0]
        current = [1.0d0, 1.0d0, 1.0d0, 0.0d0]

        call init_coils(x, y, z, current, coils)
    end subroutine init_straight_wire_coils

    subroutine remove_test_coils_file
        call system("rm -f " // test_coils_file)
    end subroutine remove_test_coils_file

end program test_biotsavart