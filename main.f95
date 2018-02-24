PROGRAM main

    use math
    use maillage
    use donnees

    implicit none

    integer :: n, i
    real(rp) :: ul, ur
    real(rp) :: xmin, xmax
    real(rp), dimension(:, :), allocatable :: A
    real(rp), dimension(:), allocatable :: b, x, u
    type(Mesh) :: maill

    open(unit = 0, file = "entrees/constantes")
    read (0, *) n
    read (0, *) ul
    read (0, *) ur
    read (0, *) xmin
    read (0, *) xmax
    close(0)

    allocate(A(n, n), b(n), x(n + 2), u(n))

    write (*, '("n = ",1I4,"   ul = ",1F4.1,"   ur = ",1F4.1,/,"I = [",1F4.1,","1F4.1," ]")') n, ul, ur, &
        xmin, xmax

    x = linspace(xmin, xmax, n + 2)
    call newMesh(x, maill)

    print *, "### l"
    print *, maill%l
    print *, "### x"
    print *, maill%x
    print *, "### h2"
    print *, maill%h2
    print *, "### x2"
    print *, maill%x2
    print *, "### h"
    print *, maill%h

    call buildA(x, A)


    deallocate(A, b, x, u)

END PROGRAM main
