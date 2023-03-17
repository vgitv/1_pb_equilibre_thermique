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

    ! lecture des données
    open(unit = 0, file = "entrees/constantes")
    read (0, *) n
    read (0, *) ul
    read (0, *) ur
    read (0, *) xmin
    read (0, *) xmax
    close(0)

    allocate(A(n, n), b(n), x(n + 2), approx(n), sol(n))

    write (*, '("n = ",1I4,"   ul = ",1F4.1,"   ur = ",1F4.1,/,"I = [",1F4.1,","1F4.1," ]")') &
        n, ul, ur, xmin, xmax

    x = linspace(xmin, xmax, n + 2)
    call newMesh(x, maill)

    ! construction système linéaire
    call build_A(maill, A)
    call build_b(maill, f, ul, ur, b)

    ! calcul approximation et solution exacte
    call linSolve("plu", maill%l, A, b, approx)
    sol = evaluate(u, maill%x)

    ! écriture des données
    call saveSol(maill%x, (/ul, approx, ur/), "sorties/approx.dat")
    call saveSol(maill%x, sol, "sorties/sol.dat")


    !********************************************************************************************************
    ! désallocations finales
    deallocate(A, b, x, approx, sol)
    call rmMesh(maill)

END PROGRAM main
