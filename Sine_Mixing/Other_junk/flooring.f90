program flooring
    implicit none
    integer(8)          :: j, N
    real(8)             :: a, w1
    real(8), parameter  :: pi = 4.D0*atan(1.D0)

    N = 4
    a = 300.0 

    do j = 1, N
        w1 = a*cos((j-1)*pi*2/N) - floor(a*cos((j-1)*pi*2/N))
        print *, 'w1 =', w1, "j-1 = ", j-1
    end do
end program flooring