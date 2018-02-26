MODULE variables

    use maillage

    implicit none

    real(rp), dimension(:, :), allocatable :: A
    real(rp), parameter :: delta_n = 1.0_rp
    real(rp), parameter :: delta_p = 1.0_rp
    real(rp), save :: psi_l = 1.0_rp
    real(rp), save :: psi_r = 2.0_rp
    type(Mesh), save :: maill

contains

    ! -------------------------------------------------------------------------------------------------------
    ! dopage
    ! -------------------------------------------------------------------------------------------------------
    function c_dopage(x)
        ! paramètres
        real(rp), intent(in) :: x

        ! return
        real(rp) :: c_dopage

        c_dopage = 3.0_rp
    end function

END MODULE variables
