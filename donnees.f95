MODULE donnees

    use math

    implicit none

contains

    ! second membre du problème
    function f(x)
        ! paramètres
        real(rp), intent(in) :: x

        ! return
        real(rp) :: f

        f = 1.0_rp
    end function



    subroutine buildA(x, A)
        ! paramètres
        real(rp), dimension(:), intent(in) :: x
        real(rp), dimension(:, :), intent(out) :: A

        ! variables locales
        integer :: i, n
        real(rp) :: t1, t2

        n = size(x) - 2

        do i = 1, n - 1
            t1 = 1.0_rp / (x(i + 1) - x(i))
            t2 = 1.0_rp / (x(i + 2) - x(i + 1))
            A(i, i) = t1 + t2
            A(i, i + 1) = -t1
            A(i + 1, i) = -t2
        end do
        A(n, n) = 1.0_rp / (x(n + 1) - x(n)) + 1.0_rp / (x(n + 2) - x(n + 1))
    end subroutine



    subroutine build_b(x, fc, ul, ur, b)
        ! paramètres
        real(rp), dimension(:), intent(in) :: x
        real(rp), dimension(:), intent(out) :: b
        real(rp), external :: fc
        real(rp), intent(in) :: ul, ur

        ! variables locales
        integer :: i, n
        real(rp) :: t

        n = size(x) - 2

        do i = 1, n
            b(i) = hi(x, i + 1) * fc(x(i + 1))
        end do

        t = 1.0_rp / (x(i + 1) - x(i))
        b(1) = b(1) + ul / t
        t = 1.0_rp / (x(i + 2) - x(i + 1))
    end subroutine

END MODULE donnees
