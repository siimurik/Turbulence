program sine_flow

    implicit none
    integer(4)                  :: i, j, k, l, nn
    integer(4), parameter       :: n = 2**2, amp = 300     ! n = 2**8
    !integer(4), parameter       :: dim = 2, nx = 100, ny = 100, niter = 1000, amp = 2
    real(8), parameter          :: pi = 4.0*atan(1.0)!, amp = 180.0 !amp = 300.0
    real(8), dimension(n)       :: x, y
    real(8), dimension(n, n)    :: c, c_new, mat, mat_new
    real(8)                     :: Lx, Ly, hx, hy, w1, w2, amo
    logical                     :: a, b, x_dir

    Lx = 1.0
    Ly = 1.0
    hx = Lx/(float(n)-1)
    hy = Ly/(float(n)-1)
!================================================================
!================================================================
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
!================================================================
    amo = 0.7 
    do i = 1, n
        do j = 1, n
            c(i,j) = sin(i * pi * 2.0/n)
        end do
    end do
    c_new = 0.0

    do i = 1, n ! y
        do j = 1, n ! x
            l = amp*cos(i*2*pi/n)
            c_new(i,j) = c(i, iand((j + l), (n-1)))
        end do
    end do
    mat = c_new
!===================================================================
    x_dir = .true.
    a = .true.
    b = .false.
    do i = 1, n
        do j = 1, n
            if (x_dir .eqv. a) then
                w1 = amp*cos(j*pi*2/n) - floor(amp*cos(j*pi*2/n))
                print *, 'eqv w1 =', w1
                w2 = 1 - w1
                mat(i,j) = w1*mat(i,j-1) + w2*mat(i,j)
            else
                w1 = amp*cos(i*pi*2/n) - floor(amp*cos(i*pi*2/n))
                !print *, 'neqv w1 =', w1
                w2 = 1 - w1
                mat(i,j) = w1*mat(i-1,j) + w2*mat(i,j)
            end if
        end do
    end do
    mat_new = 0
    if (x_dir .eqv. a) then
        do i = 1, n
            do j = 1, n
                nn = n+1
                mat_new(i,j) = mat(nn-j,nn-i)
            end do
        end do
        mat = mat_new
        do i = 1, n
            print *, mat(i,:)
        end do
    else
        do i = 1, n
            do j = 1, n
                nn = n+1
                mat_new(i,j) = mat(nn-i,nn-j)               
            end do
        end do
        mat = mat_new
    end if


    do i = 1, n ! y
        do j = 1, n ! x
            if (x_dir .eqv. a) then
                w1 = amp*cos(j*pi*2/n) - floor(amp*cos(j*pi*2/n))
                w2 = 1 - w1
                mat(i,j) = w1*mat(i,j-1) + w2*mat(i,j)
            else
                w1 = amp*cos(j*pi*2/n) - floor(amp*cos(j*pi*2/n))
                w2 = 1 - w1
                mat(i,j) = w1*mat(i-1,j) + w2*mat(i,j)
            end if
        end do
    end do  
!====================================================================
        ! Formatting for CSV
101 format(1x, *(g0, ", "))
    ! Open connection (i.e. create file where to write)
    open(unit = 10, access = "sequential", action = "write", &
        status = "replace", file = "data1.csv", form = "formatted")
    ! Loop across rows
    do i= 1, n
        write(10, 101) mat(i,:)
    end do
    ! Close connection
    close(10)

end program sine_flow