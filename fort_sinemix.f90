program sine_flow

    implicit none
    integer(8)                  :: i, j, k, l
    integer(8), parameter       :: n = 2**8     ! n = 2**8
    !integer(4), parameter       :: dim = 2, nx = 100, ny = 100, niter = 1000, amp = 2
    real(8), parameter          :: pi = 4.0*atan(1.0), amp = 250.0 !amp = 300.0
    real(8), dimension(n)       :: x, y
    real(8), dimension(n, n)    :: tracer, tracer_new
    real(8)                     :: Lx, Ly, hx, hy

    Lx = 1.0
    Ly = 1.0
    hx = Lx/(float(n)-1)
    hy = Ly/(float(n)-1)

    do i = 1, n
        do j = 1, n
            tracer(i,j) = sin(i * pi * 2.0/n)
        end do
    end do
    tracer_new = 0.0

    open(unit=9,file='x_y.csv')
    do i = 1, n
        x(i) = 0. + (i-1)*hx
        y(i) = 0. + (i-1)*hy
    end do
    ! Printing in a separate file x_y.csv
    do i = 1, n
        write ( 9, * ) x(i), ',',y(i)
    end do
    ! Close connection
    close(9)

    !do n = 1, 512 ! cycles
        do i = 1, n ! y
            do j = 1, n ! x
                k = amp*sin(j*2*pi/n)
                tracer_new(i,j) = tracer(iand((i + k),(n-1)), j)
            end do
        end do
    tracer = tracer_new

        do i = 1, n ! y
            do j = 1, n ! x
                l = amp*sin(i*2*pi/n)
                tracer_new(i,j) = tracer(i, iand((j + l), (n-1)))
            end do
        end do
    tracer = tracer_new
    !end do

    ! Formatting for CSV
101 format(1x, *(g0, ", "))
    ! Open connection (i.e. create file where to write)
    open(unit = 10, access = "sequential", action = "write", &
        status = "replace", file = "data.csv", form = "formatted")
    ! Loop across rows
    do i= 1, n
        write(10, 101) tracer(i,:)
    end do
    ! Close connection
    close(10)

end program sine_flow