MODULE donnees

    use math
    use maillage

    implicit none

contains

    ! -------------------------------------------------------------------------------------------------------
    ! second membre du problème
    ! -------------------------------------------------------------------------------------------------------
    function f(x)
        ! paramètres
        real(rp), intent(in) :: x

        ! return
        real(rp) :: f

        f = 1.0_rp
    end function



    ! -------------------------------------------------------------------------------------------------------
    ! Construction matrice A pour résolution A x = b dans le cas de l'équation de Poisson
    ! -------------------------------------------------------------------------------------------------------
    subroutine build_A(m, A)
        ! paramètres
        type(Mesh), intent(in) :: m
        real(rp), dimension(:, :), intent(out) :: A

        ! variables locales
        integer :: i

        do i = 1, m%l - 1
            A(i, i) = 1.0_rp / m%h2(i) - 1.0_rp / m%h2(i + 1)
            A(i, i + 1) = -1.0_rp / m%h2(i + 1)
            A(i + 1, i) = A(i, i + 1)
        end do
        A(m%l, m%l) = 1.0_rp / m%h2(m%l) + 1.0_rp / m%h2(m%l + 1)
    end subroutine build_A



    ! -------------------------------------------------------------------------------------------------------
    ! Construction vecteur b pour résolution A x = b dans le cas de l'équation de Poisson
    ! -------------------------------------------------------------------------------------------------------
    subroutine build_b(m, fc, ul, ur, b)
        ! paramètres
        type(Mesh), intent(in) :: m
        real(rp), dimension(:), intent(out) :: b
        real(rp), external :: fc
        real(rp), intent(in) :: ul, ur

        ! variables locales
        integer :: i

        do i = 1, m%l
            b(i) = m%h(i) * f(m%x(i + 1))
        end do

        b(1) = b(1) + ul / m%h2(1)
        b(m%l) = b(m%l) + ur / m%h2(m%l + 1)
    end subroutine build_b

END MODULE donnees
