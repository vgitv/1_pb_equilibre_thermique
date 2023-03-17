MODULE donnees

    use math
    use maillage

    implicit none

    real(rp), parameter :: delta_n = 1.0_rp
    real(rp), parameter :: delta_p = 1.0_rp

contains

    ! =======================================================================================================
    ! PROBLÈME DE POISSON
    ! =======================================================================================================

    ! -------------------------------------------------------------------------------------------------------
    ! second membre du problème de Poisson
    ! -------------------------------------------------------------------------------------------------------
    function f(x)
        ! paramètres
        real(rp), intent(in) :: x

        ! return
        real(rp) :: f

        f = 1.0_rp
    end function



    ! -------------------------------------------------------------------------------------------------------
    ! Solution de l'équation de Poisson avec f = 1 et ul = ur = 0
    ! -------------------------------------------------------------------------------------------------------
    function u(x)
        ! paramètres
        real(rp), intent(in) :: x

        ! return
        real(rp) :: u

        u = x * (1.0_rp - x) / 2.0_rp
    end function



    ! -------------------------------------------------------------------------------------------------------
    ! Construction matrice A pour résolution A x = b dans le cas de l'équation de Poisson (idem éq thermique)
    ! -------------------------------------------------------------------------------------------------------
    subroutine build_A(m, A)
        ! paramètres
        type(Mesh), intent(in) :: m
        real(rp), dimension(:, :), intent(out) :: A

        ! variables locales
        integer :: i

        A = 0.0_rp
        do i = 1, m%l - 1
            A(i, i) = 1.0_rp / m%h2(i) + 1.0_rp / m%h2(i + 1)
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



    ! =======================================================================================================
    ! NEWTON 1D
    ! =======================================================================================================

    ! -------------------------------------------------------------------------------------------------------
    ! fonction test pour algo newton 1D : une fonction g et sa dérivée gp ( = gprime)
    ! -------------------------------------------------------------------------------------------------------
    function g(x)
        ! paramètres
        real(rp), intent(in) :: x

        ! return
        real(rp) :: g

        g = (x - 1.0_rp) * (x + 1.0_rp) * (x - 3.0_rp)
    end function

    function gp(x)
        ! paramètres
        real(rp), intent(in) :: x

        ! return
        real(rp) :: gp

        gp = (x + 1.0_rp) * (x - 3.0_rp) + (x - 1.0_rp) * (x - 3.0_rp) + (x - 1.0_rp) * (x + 1.0_rp)
    end function



    ! =======================================================================================================
    ! PROBLÈME À L'ÉQUILIBRE THERMIQUE
    ! =======================================================================================================

    ! -------------------------------------------------------------------------------------------------------
    ! dopage
    ! -------------------------------------------------------------------------------------------------------
    function c_dopage(x)
        ! paramètres
        real(rp), intent(in) :: x

        ! return
        real(rp) :: c_dopage

        c_dopage = 0.0_rp
    end function



    ! -------------------------------------------------------------------------------------------------------
    ! fonction b du système A.psi + fc_b(psi) + bd = 0   (même matrice A que pour l'équation de Poisson)
    ! -------------------------------------------------------------------------------------------------------
    subroutine fc_b(m, psi, s)
        ! paramètres
        type(Mesh), intent(in) :: m
        real(rp), dimension(:), intent(in) :: psi
        real(rp), dimension(m%l), intent(out) :: s

        ! variables locales
        integer :: i

        do i = 1, m%l
            s(i) = m%h(i) * (delta_n * exp(psi(i)) - delta_p * exp(-psi(i)) - c_dopage(m%x(i + 1)))
        end do
    end subroutine



    ! -------------------------------------------------------------------------------------------------------
    ! terme constant bd du système A.psi + fc_b(psi) + bd = 0 (conditions de bord)
    ! -------------------------------------------------------------------------------------------------------
    subroutine cste_bd(m, psi_l, psi_r, s)
        ! paramètres
        type(Mesh), intent(in) :: m
        real(rp), intent(in) :: psi_l, psi_r
        real(rp), dimension(m%l), intent(out) :: s

        s = 0.0_rp
        s(1) = -psi_l / m%h2(1)
        s(m%l) = -psi_r / m%h2(m%l + 1)
    end subroutine



    ! -------------------------------------------------------------------------------------------------------
    ! f(psi) = A*psi + fc_b(psi) + bd
    ! -------------------------------------------------------------------------------------------------------
    subroutine fmain(m, A, psi, bd, s)
        ! paramètres
        type(Mesh), intent(in) :: m
        real(rp), dimension(m%l, m%l), intent(in) :: A
        real(rp), dimension(m%l), intent(in) :: psi, bd
        real(rp), dimension(m%l), intent(out) :: s

        ! variables locales
        real(rp), dimension(m%l) :: temp

        call fc_b(m, psi, temp)
        s = matmul(A, psi) + temp + bd
    end subroutine



    ! -------------------------------------------------------------------------------------------------------
    ! Jacobienne du système non linéaire pour le pb de dérive diffusion stationnaire à l'équilibre thermique
    ! -------------------------------------------------------------------------------------------------------
    subroutine Jfmain(m, A, psi, jacobienne)
        ! paramètres
        type(Mesh), intent(in) :: m
        real(rp), dimension(m%l, m%l), intent(in) :: A
        real(rp), dimension(m%l), intent(in) :: psi
        real(rp), dimension(m%l, m%l), intent(out) :: jacobienne

        ! variables locales
        integer :: i

        jacobienne = A
        do i = 1, m%l
            jacobienne(i, i) = jacobienne(i, i) + m%h(i) * (delta_n * exp(psi(i)) + delta_p * exp(-psi(i)))
        end do
    end subroutine

END MODULE donnees
