program test_biotsavart
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use test_util, only: print_test, print_ok, print_fail

    implicit none

    character(*), parameter :: test_coils_file = "coils.test"

    real(dp), parameter :: large_distance = 1.0d3

    call test_load_coils_file
    call test_compute_coils_segment_lengths
    call test_calc_subknot_xyz
    call test_cut_coils_segments
    call test_equalize_coils_lenghts
    call test_compute_vector_potential

    contains

    subroutine test_load_coils_file
        use biotsavart, only: coils_t, load_coils_from_file, deinit_coils

        type(coils_t) :: coils

        call print_test("load_coils_file")

        call create_straight_wire_coils_file
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
        use biotsavart, only: coils_t, init_coils, &
                              compute_coils_segments_lengths, deinit_coils

        real(dp), parameter :: tol = 1.0e-9

        type(coils_t) :: coils
        real(dp), dimension(5) :: lengths
        real(dp), dimension(5) :: expected_lengths
        integer :: i

        call print_test("compute_coils_segment_lengths")

        call init_diamond_wire_coils(coils)

        lengths = compute_coils_segments_lengths(coils)
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


    subroutine test_calc_subknot_xyz
        use biotsavart, only: coils_t, init_coils, calc_subknot_xyz, deinit_coils

        real(dp), parameter :: tol = 1.0e-9

        type(coils_t) :: coils
        integer :: knot
        real(dp), dimension(3) :: xyz, expected_xyz

        call print_test("calc_subknot_xyz")

        call init_diamond_wire_coils(coils)

        do knot = 1, size(coils%x) - 1
            xyz = calc_subknot_xyz(coils, knot, subknot=2, cuts_per_knot=3)
            expected_xyz = [(coils%x(knot) + coils%x(knot+1))/2.0d0, &
                            (coils%y(knot) + coils%y(knot+1))/2.0d0, &
                            (coils%z(knot) + coils%z(knot+1))/2.0d0]
            if (any(abs(xyz - expected_xyz) > tol)) then
                print *, "xyz = ", xyz
                print *, "expected_xyz = ", expected_xyz
                call print_fail
                error stop
            end if
        end do

        call deinit_coils(coils)
        call print_ok
    end subroutine test_calc_subknot_xyz


    subroutine test_cut_coils_segments
        use biotsavart, only: coils_t, init_coils, cut_coils_segments, &
                              compute_coils_segments_lengths, deinit_coils

        real(dp), parameter :: tol = 1.0e-9

        type(coils_t) :: coils
        integer, dimension(4) :: cuts_per_knot
        real(dp) :: expected_x = 0.5d0, expected_y = 1.0d0, expected_current = 1.0d0

        call print_test("cut_coils_segments")

        call init_diamond_wire_coils(coils)

        cuts_per_knot = [0, 1, 0, 1]

        call cut_coils_segments(coils, cuts_per_knot)
        if (size(coils%x) /= 7) then
            print *, "Coil length mismatch"
            print *, "len(coils%x) = ", size(coils%x)
            call print_fail
            error stop
        end if
        if ((coils%x(3) - expected_x) > tol) then
            print *, "Coil x mismatch"
            print *, "coils%x(3) = ", coils%x(3)
            print *, "expected = ", expected_x
            call print_fail
            error stop
        end if
        if ((coils%y(7)-expected_y) > tol) then
            print *, "Coil y mismatch"
            print *, "coils%y(7) = ", coils%y(3)
            print *, "expected = ", expected_y
            call print_fail
            error stop
        end if
        if ((coils%current(6)-expected_current) > tol) then
            print *, "Coil current mismatch"
            print *, "coils%current(6) = ", coils%current(6)
            print *, "expected = ", expected_current
            call print_fail
            error stop
        end if

        call deinit_coils(coils)
        call print_ok
    end subroutine test_cut_coils_segments


    subroutine test_equalize_coils_lenghts
        use biotsavart, only: coils_t, init_coils, equalize_coils_segments_lengths, &
                              deinit_coils, compute_coils_segments_lengths

        type(coils_t) :: coils, old_coils
        real(dp), dimension(4) :: lengths
        real(dp) :: min_length

        call print_test("equalize_coils_lenghts")

        call init_diamond_wire_coils(coils)
        call equalize_coils_segments_lengths(coils)
        call init_diamond_wire_coils(old_coils)
        if (.not.(are_coils_equal(coils, old_coils))) then
            print *, "Equalized coil is not unchanged"
            call print_fail
            error stop
        end if
        call deinit_coils(coils)
        call deinit_coils(old_coils)

        call init_straight_wire_coils(coils)
        min_length = minval(compute_coils_segments_lengths(coils))
        call equalize_coils_segments_lengths(coils)
        lengths = compute_coils_segments_lengths(coils)
        if (any(abs(lengths - min_length) > min_length)) then
            print *, "Coil segments lengths mismatch"
            print *, "lengths = ", lengths
            call print_fail
            error stop
        end if
        print *, compute_coils_segments_lengths(coils)

        call deinit_coils(coils)
        call print_ok
    end subroutine test_equalize_coils_lenghts

    function are_coils_equal(coil1, coil2)
        use biotsavart, only: coils_t

        real(dp), parameter :: tol = 1.0e-9

        type(coils_t), intent(in) :: coil1, coil2
        logical :: are_coils_equal
        integer :: i

        are_coils_equal = .true.
        if (size(coil1%x) /= size(coil2%x)) then
            are_coils_equal = .false.
            return
        end if
        do i = 1, size(coil1%x)
            if (abs(coil1%x(i) - coil2%x(i)) > tol) then
                are_coils_equal = .false.
                return
            end if
            if (abs(coil1%y(i) - coil2%y(i)) > tol) then
                are_coils_equal = .false.
                return
            end if
            if (abs(coil1%z(i) - coil2%z(i)) > tol) then
                are_coils_equal = .false.
                return
            end if
            if (abs(coil1%current(i) - coil2%current(i)) > tol) then
                are_coils_equal = .false.
                return
            end if
        end do
    end function are_coils_equal

    subroutine test_compute_vector_potential
        use biotsavart, only: coils_t, compute_vector_potential, &
            deinit_coils, clight

        real(dp), parameter :: tol = 1.0e-9
        integer, parameter :: N_TEST = 3

        type(coils_t) :: coils
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


    subroutine create_straight_wire_coils_file
        use biotsavart, only: coils_t, save_coils_to_file

        type(coils_t) :: coils

        call init_straight_wire_coils(coils)
        call save_coils_to_file(test_coils_file, coils)
    end subroutine create_straight_wire_coils_file


    subroutine init_straight_wire_coils(coils)
        use biotsavart, only: coils_t, init_coils

        type(coils_t), intent(out) :: coils

        real(dp), dimension(4) :: x, y, z, current
        real(dp)               :: L

        x = [0.0d0, 0.0d0, 0.0d0, 0.0d0]
        y = [0.0d0, 0.0d0, 0.0d0, 0.0d0]
        L = large_distance
        z = [-L/2.0d0, -L/2.5d0, 0.0d0, L/2.0d0]
        current = [1.0d0, 1.0d0, 1.0d0, 0.0d0]

        call init_coils(x, y, z, current, coils)
    end subroutine init_straight_wire_coils


    subroutine init_diamond_wire_coils(coils)
        use biotsavart, only: coils_t, init_coils

        type(coils_t), intent(out) :: coils

        real(dp), dimension(5) :: x, y, z, current

        x = [-1.0d0, 0.0d0, 1.0d0, 0.0d0, -1.0d0]
        y = [1.0d0, 0.0d0, -1.0d0, 0.0d0, 1.0d0]
        z = [0.0d0, 1.0d0, 0.0d0, -1.0d0, 0.0d0]
        current = [1.0d0, 1.0d0, 1.0d0, 1.0d0, 0.0d0]

        call init_coils(x, y, z, current, coils)
    end subroutine init_diamond_wire_coils


    subroutine remove_test_coils_file
        call system("rm -f " // test_coils_file)
    end subroutine remove_test_coils_file

    
end program test_biotsavart
