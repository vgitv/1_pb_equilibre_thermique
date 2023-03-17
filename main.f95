! ===========================================================================================================
! Approximation du pb de dérive diffusion stationnaire à l'équilibre thermique. Contenu :
!
! Rappel équation de Poisson.
! Rappel méthode de Newton 1D.
! Résolution approchée équilibre thermique.
! ===========================================================================================================

PROGRAM main

    use math
    use maillage
    use donnees

    implicit none

    integer :: n, i
    real(rp) :: ul, ur
    real(rp) :: xmin, xmax
    real(rp), dimension(:, :), allocatable :: A
    real(rp), dimension(:), allocatable :: b, x, approx, sol
    type(Mesh) :: maill
    real(rp) :: zero

    ! lecture des données
    open(unit = 0, file = "entrees/constantes")
    read (0, *) n
    read (0, *) ul
    read (0, *) ur
    read (0, *) xmin
    read (0, *) xmax
    close(0)

    allocate(A(n, n), b(n), x(n + 2), approx(n), sol(n))

    print *, "ÉQUATION DE POISSON"
    write (*, '("n = ",1I4,"   ul = ",1F4.1,"   ur = ",1F4.1,/,"I = [",1F4.1,","1F4.1," ]")') &
        n, ul, ur, xmin, xmax

    x = linspace(xmin, xmax, n + 2)
    call newMesh(x, maill)



    ! -------------------------------------------------------------------------------------------------------
    ! Équation de Poisson
    ! -------------------------------------------------------------------------------------------------------
    ! construction système linéaire
    call build_A(maill, A)
    call build_b(maill, f, ul, ur, b)

    ! calcul approximation et solution exacte
    call linSolve("plu", maill%l, A, b, approx)
    sol = evaluate(u, maill%x)

    ! écriture des données
    call saveSol(maill%x, (/ul, approx, ur/), "sorties/approx.dat")
    call saveSol(maill%x, sol, "sorties/sol.dat")



    ! -------------------------------------------------------------------------------------------------------
    ! Méthode de Newton 1D
    ! -------------------------------------------------------------------------------------------------------
    print *
    print *, "MÉTODE DE NEWTON"
    call newton1D(-7.0_rp, g, gp, 0.0001_rp, 1000, zero)
    write (*, '("zero 1 :",1F5.2)') zero
    call newton1D(0.5_rp, g, gp, 0.0001_rp, 1000, zero)
    write (*, '("zero 2 :",1F5.2)') zero
    call newton1D(5.0_rp, g, gp, 0.0001_rp, 1000, zero)
    write (*, '("zero 3 :",1F5.2)') zero



    ! -------------------------------------------------------------------------------------------------------
    ! Calcul d'une approximation de l'équilibre thermique
    ! -------------------------------------------------------------------------------------------------------
    print *
    print *, "ÉQUATION THERMIQUE"



    !********************************************************************************************************
    ! désallocations finales
    deallocate(A, b, x, approx, sol)
    call rmMesh(maill)

END PROGRAM main
