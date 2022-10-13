program rev_test
    implicit None

    integer(8)              :: i, j, l, nn, k
    integer(8), parameter   :: N = 512
    real(8)                 :: a, amo, w1, w2
    real(8), parameter      :: pi = 4.D0*atan(1.D0)
    real(8), dimension(N,N) :: c1, c2, c_new, c_miin

    a = 500
    amo = 0.7
    c1 = 0.0

    do i = 1, N
        do j = 1, N
            c1(i,j) = -cos(j * pi*2/N)
        end do
    end do
    
    !do i = 1, N
    !    print *, c1(i,:)
    !end do

    c2 = 0.0
    do i = 1, N ! y
        do j = 1, N ! x
            l = a*sin(i*2*pi/N)
            !print *, "l =", l
            c2(i,j) = c1(i, iand((j+l-1), (N-1))+1)
            !print *, (j+l-1), "AND", N-1, "  =  ", iand((j+l-1), (N-1))+1
        end do
    end do
    !print *, 499, "AND", 3, "  =  ", iand((j+l+1), (N-1))
    !do i = 1, N
    !    print *, c2(i,:)
    !end do

    c_miin = 0.0
    do i = 1, N
        do j = 1, N
            w1 = a*cos((j-1)*pi*2/N) - floor(a*cos((j-1)*pi*2/N))
            !print *, ' eqv w1 = ', w1
            w2 = 1 - w1
            !print *, "w2 =", w2
            if (j == 1) then
                k = 4
                c_miin(i,j) = w1*c2(i,k)
            else
                k = j-1
                c_miin(i,j) = w1*c2(i,k)
            end if
            c2(i,j) =  c_miin(i,j) + w2*c2(i,j)
        end do
    end do

    !do i = 1, N
    !    print *, c2(i,:)
    !end do

    c_new = 0.0
    do i = 1, N
        do j = 1, N
            nn = n+1
            c_new(i,j) = c2(nn-i,nn-j)
        end do
    end do

    !do i = 1, N
    !    print *, c_new(i, :)
    !end do
    c_miin = 0.0
    do i = 1, N
        do j = 1, N
            w1 = a*cos((j-1)*pi*2/N) - floor(a*cos((j-1)*pi*2/N)) 
            !print *, ' eqv w1 new= ', w1
            w2 = w1 - 1
            !print *, "w2 =", w2
            if (j == 1) then
                k = 4
                c_miin(i,j) = w1*c2(i,k)
            else
                k = j-1
                c_miin(i,j) = w1*c2(i,k)
            end if
            c2(i,j) =  c_miin(i,j) + w2*c2(i,j)
        end do
    end do

    !do i = 1, N
    !    print *, c_new(i, :)
    !end do

    
        ! Formatting for CSV
    101 format(1x, *(g0, ", "))
    ! Open connection (i.e. create file where to write)
    open(unit = 10, access = "sequential", action = "write", &
        status = "replace", file = "data1.csv", form = "formatted")
    ! Loop across rows
    do i= 1, n
        write(10, 101) c_new(i,:)
    end do
    ! Close connection
    close(10)

end program