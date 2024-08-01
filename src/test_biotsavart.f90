program test_biotsavart
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use test_util, only: print_test, print_ok, print_fail

    implicit none

    character(*), parameter :: test_coils_file = "coils.test"
    character(*), parameter :: real_coils_file = "/proj/plasma/data/W7X/COILS/coils.w7x"

    real(dp), parameter :: large_distance = 1.0d3

    call test_load_coils_file
    call test_compute_vector_potential

    contains

    subroutine test_load_coils_file
        use biotsavart, only: CoilsData, load_coils_file, deinit_coils_data

        type(CoilsData) :: coils

        call print_test("load_coils_file")

        call create_test_coils_file
        call load_coils_file(test_coils_file, coils)
        call remove_test_coils_file

        if (size(coils%x) /= 2) then
            print *, "Coil length mismatch"
            print *, "len(coils%x) = ", size(coils%x)
            call print_fail
            error stop
        end if

        call deinit_coils_data(coils)

        call print_ok
    end subroutine test_load_coils_file


    subroutine test_compute_vector_potential
        use biotsavart, only: CoilsData, compute_vector_potential, &
            deinit_coils_data, clight

        real(dp), parameter :: tol = 1.0e-9
        integer, parameter :: N_TEST = 3

        type(CoilsData) :: coils
        real(dp) :: x_test(3, N_TEST)
        real(dp), dimension(3) :: x, A, A_expected
        real(dp) :: sqrt_term, L, R, AZ_expected
        integer :: i

        call print_test("compute_vector_potential")

        x_test(:, 1) = [0.4, 0.3, 0.8]
        x_test(:, 2) = [0.0, 0.1, 0.3]
        x_test(:, 3) = [1e-5, 0.0, 0.0]

        call init_test_coils_data(coils)

        do i = 1, N_TEST
            x = x_test(:, i)
            L = large_distance
            R = Rcyl(x)
            sqrt_term = sqrt(L**2 + R**2)
            AZ_expected = 1.0d0/clight*log((L + sqrt_term)/(-L + sqrt_term))
            A_expected = [0.0d0, 0.0d0, AZ_expected]
            A = compute_vector_potential(coils, x)
            print *, "A = ", A(3)*clight
            print *, "A_expected = ", A_expected(3)*clight
            print *, "Ratio = ", A(3) / A_expected(3)
            if (any(abs(A - A_expected)*clight > tol)) then
                print *, "A = ", A(3)*clight
                print *, "A_expected = ", A_expected(3)*clight
                print *, "Ratio = ", A(3) / A_expected(3)
                call print_fail
                error stop
            end if
        end do

        call deinit_coils_data(coils)

        call print_ok
    end subroutine test_compute_vector_potential


    function Rcyl(x)
        real(dp), dimension(3), intent(in) :: x
        real(dp) :: Rcyl

        Rcyl = sqrt(x(1)**2 + x(2)**2)
    end function Rcyl


    subroutine init_test_coils_data(coils)
        use biotsavart, only: CoilsData, init_coils_data

        type(CoilsData), intent(out) :: coils

        real(dp), dimension(2) :: x, y, z, current

        x = [0.0d0, 0.0d0]
        y = [0.0d0, 0.0d0]
        z = [-large_distance/2.0d0, large_distance/2.0d0]
        current = [1.0d0, 0.0d0]

        call init_coils_data(x, y, z, current, coils)
    end subroutine init_test_coils_data


    subroutine create_test_coils_file
        use biotsavart, only: CoilsData, init_coils_data, save_coils_file

        type(CoilsData) :: coils

        call init_test_coils_data(coils)
        call save_coils_file(test_coils_file, coils)
    end subroutine create_test_coils_file


    subroutine remove_test_coils_file
        call system("rm -f " // test_coils_file)
    end subroutine remove_test_coils_file

end program test_biotsavart