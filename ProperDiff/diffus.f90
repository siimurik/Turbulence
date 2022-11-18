!=================================================================
! Compile with:
!   > gfortran main_mix.f90 -o main
!   > ./main
!=================================================================
! Code for a 'simple' sinusoidal mixing in a hypothetical liquid 
! or gasueous medium. The tracer or pigment inside the medium is 
! mixed by a simple trigonometric function, whose values are then
! bitwise manipulated to store new values into a matrix, which
! will be plotted by the 'plot_mix.py' code. Results are stored
! in the file 'dataMAIN.csv', which contains the final values of
! of the matrixRes. 
!
! Every new iteration after c1 matrix must initialized as in the
! code below. This code is currenty able to calculate up to 4 
! iterations (c5) of mixing, most of which are currently commented out. 
! Of course new matrices can be stored like c6, c7 and so on, but
! that has to be done by the user. For current purposes the four
! different cycles are quite sufficient to show how mixing 
! progresses. 
!
! To avoid too caotic mixing and the plots looking like static,
! a subroutine named 'mixit(...)' has been made to keep the 
! mixing process under control. This means that thin lines that 
! occur during the mixing process will still look like lines, 
! not static.
!=================================================================
program main
    !==========================================================
        implicit none
        integer(8)              :: i, j, l
        integer(8), parameter   :: N = 512
        real(8)                 :: amp, w1, w2
        real(8), parameter      :: pi = 4.D0*atan(1.D0)
        real(8), dimension(N,N) :: c1, c2, c3, c4, c5
        real(8), dimension(N,N) :: matRev, matInter, matResInit
        real(8), dimension(N,N) :: matrixInit, matrixRes
    !==========================================================
    ! Initial conditions
        amp = 500
        c1 = 0.0
        do i = 1, N
            do j = 1, N
                c1(i,j) = -cos(j*pi*2/N)
            end do
        end do
    ! Print out c1
        !do i = 1, N
        !    print *, c1(i,:)
        !end do
    !--------------------------------------------------------------------------
    ! Initialize c2 matrix elements as zero
        c2 = 0.0
    ! Start generating elements for c2
        do i = 1, N ! y
            do j = 1, N ! x
                l = amp*sin(i*2*pi/N)
                !print *, "l =", l
                c2(i,j) = c1(i, iand((j+l-1), (N-1))+1)
                !print *, (j+l-1), "AND", N-1, "  =  ", iand((j+l-1), (N-1))+1
            end do
        end do
    ! Print out c2
        !do i = 1, N
        !    print *, c2(i,:)
        !end do
    !------------------------------------------------------------------------
    ! Initialize matrixInit for subroutine 'mixit'    
        do i = 1, N
            do j = 1, N
                matrixInit(i,j) = c2(i,j)
            end do
        end do
    ! This is where the fun begins
        call mixit(amp, w1, w2, matrixInit, matRev, matInter, matResInit, matrixRes)
    ! END OF FIRST ITERATION
    !========================================================================
    ! BEGINNING OF SECOND ITERATION
        !c3 = 0.0
        !do i = 1, N ! y
        !    do j = 1, N ! x
        !        l = amp*sin(i*2*pi/N)
        !        !print *, "l =", l
        !        c3(i,j) = matrixRes(i, iand((j+l-1), (N-1))+1)
        !        !print *, (j+l-1), "AND", N-1, "  =  ", iand((j+l-1), (N-1))+1
        !    end do
        !end do
        !
        !do i = 1, N
        !    do j = 1, N
        !        matrixInit(i,j) = c3(i,j)
        !    end do
        !end do
    !   !
        !call mixit(amp, w1, w2, matrixInit, matRev, matInter, matResInit, matrixRes)

! END OF SECOND ITERATION
!========================================================================
! BEGINNING OF THIRD ITERATION
        !c4 = 0.0
        !do i = 1, N ! y
        !    do j = 1, N ! x
        !        l = amp*sin(i*2*pi/N)
        !        !print *, "l =", l
        !        c4(i,j) = matrixRes(i, iand((j+l-1), (N-1))+1)
        !        !print *, (j+l-1), "AND", N-1, "  =  ", iand((j+l-1), (N-1))+1
        !    end do
        !end do
    !   !
        !do i = 1, N
        !    do j = 1, N
        !        matrixInit(i,j) = c4(i,j)
        !    end do
        !end do
    !   !
        !call mixit(amp, w1, w2, matrixInit, matRev, matInter, matrixRes)
! END OF THIRD ITERATION
!========================================================================
! BEGINNING OF FOURTH ITERATION
        !c4 = 0.0
        !do i = 1, N ! y
        !    do j = 1, N ! x
        !        l = amp*sin(i*2*pi/N)
        !        !print *, "l =", l
        !        c5(i,j) = matrixRes(i, iand((j+l-1), (N-1))+1)
        !        !print *, (j+l-1), "AND", N-1, "  =  ", iand((j+l-1), (N-1))+1
        !    end do
        !end do
    !   !
        !do i = 1, N
        !    do j = 1, N
        !        matrixInit(i,j) = c5(i,j)
        !    end do
        !end do
    !   !
        !call mixit(amp, w1, w2, matrixInit, matRev, matInter, matrixRes)
! END OF FOURTH ITERATION
!------------------------------------------------------------------------
! Formatting for CSV
        101 format(1x, *(g0, ", "))
        ! Open connection (i.e. create file where to write)
        open(unit = 10, access = "sequential", action = "write", &
            status = "replace", file = "dataMAIN.csv", form = "formatted")
        ! Loop across rows
        do i= 1, n
            write(10, 101) matrixRes(i,:)
        end do
        ! Close connection
        close(10)
!=========================================================================
end program main
    
subroutine mixit(amp, w1, w2, matrixInit, matRev, matInter, matResInit, matrixRes)
!=========================================================================
!   The main heart of the code uses the diffusing 'mixit' subroutine.
!   Uses a simplistic way of diffusing the contcentrated areas.    
!-------------------------------------------------------------------------
!   HOW TO USE:
!-------------------------------------------------------------------------
! INPUTS:
!   amp - amplitude of the wave values
!   matrixInit - the initial matrix that will be mixed
!-------------------------------------------------------------------------
! INTERNAL VARIABLES:
!   w1      - diffusion coefficient; real number
!   w2      - diffusion coefficient; real number
!   matRev  - after first two do loops matrixInit values are rewritten
!             into a new matrix called matRev, which stores the matrxiInit
!             values in reverse order, hence the name. If this ignored, it 
!             will cause a bug, bc old values from matrixInit may be reused.
!   matInter - basically an internal variable that is used for storing
!              temporary matrix values. Used for the if condition where
!              the new first element is the last old element. In python
!              this is just A[i] = A[i-1].
!-------------------------------------------------------------------------
! OUTPUT:
!   matrixRes - final output of the subroutine, which is also written in
!               the dataMAIN.csv file.
!               NB! If there is a need to flip the result image,
!               instead of the 2 do loops at the end of the subroutine
!               just use the transpose command and comment out the do loops.
!-------------------------------------------------------------------------
! EXAMPLE:
!-------------------------------------------------------------------------
!    ! Initialize c2 matrix elements as zero
!    c2 = 0.0
!    ! Start generating elements for c2
!        do i = 1, N ! y
!            do j = 1, N ! x
!                l = amp*sin(i*2*pi/N)
!                !print *, "l =", l
!                c2(i,j) = c1(i, iand((j+l-1), (N-1))+1)
!                !print *, (j+l-1), "AND", N-1, "  =  ", iand((j+l-1), (N-1))+1
!            end do
!        end do
!    !------------------------------------------------------------------------
!    ! Initialize matrixInit for subroutine 'mixit'    
!       do i = 1, N
!           do j = 1, N
!               matrixInit(i,j) = c2(i,j)
!           end do
!       end do
!       call mixit(amp, w1, w2, matrixInit, matRev, matInter, matrixRes)    
!-------------------------------------------------------------------------
!=========================================================================
    implicit none
    integer(8)              :: i, j, k
    integer(8), parameter   :: N = 512
    real(8)                 :: amp, w1, w2
    real(8), parameter      :: pi = 4.D0*atan(1.D0)
    real(8), dimension(N,N) :: matRev, matInter, matResInit
    real(8), dimension(N,N) :: matrixInit, matrixRes
!=========================================================================
! Introduce new matrix matInter as matrix with
! elements set to zero to avoid reusing old values.

    matInter = 0.0
    do i = 1, N
        do j = 1, N
            w1 = amp*cos((j-1)*pi*2/N) - floor(amp*cos((j-1)*pi*2/N))
            !print *, ' eqv w1 = ', w1
            w2 = 1 - w1
            !print *, "w2 =", w2
            if (j == 1) then
                k = 4
                matInter(i,j) = w1*matrixInit(i,k)
            else
                k = j-1
                matInter(i,j) = w1*matrixInit(i,k)
            end if
            matrixInit(i,j) =  matInter(i,j) + w2*matrixInit(i,j)
        end do
    end do

    !do i = 1, N
    !    print *, matrixInit(i,:)
    !end do

    matRev = 0.0
    do i = 1, N
        do j = 1, N
            matRev(i,j) = matrixInit(N+1-i,N+1-j) ! reverse martrix
        end do
    end do

    !do i = 1, N
    !    print *, matRev(i, :)
    !end do
    matInter = 0.0
    matrixRes = 0.0
    do i = 1, N
        do j = 1, N
            w1 = amp*cos((j-1)*pi*2/N) - floor(amp*cos((j-1)*pi*2/N)) 
            !print *, ' eqv w1 new= ', w1
            w2 = w1 - 1
            !print *, "w2 =", w2
            if (j == 1) then
                k = 4
                matInter(i,j) = w1*matRev(i,k)
            else
                k = j-1
                matInter(i,j) = w1*matRev(i,k)
            end if
            matRev(i,j) =  matInter(i,j) + w2*matRev(i,j)
        end do
    end do

    do i = 1, N
        do j = 1, N
            matResInit(i,j) = matRev(i,j)
        end do
    end do
    !matResInit = transpose(matRev)
    matInter = 0.0
    do i = 1, N
        do j = 1, N
            w1 = amp*cos((j-1)*pi*2/N) - floor(amp*cos((j-1)*pi*2/N))
            !print *, ' eqv w1 = ', w1
            w2 = 1 - w1
            !print *, "w2 =", w2
            if (j == 1) then
                k = 4
                matInter(i,j) = w1*matResInit(i,k)
            else
                k = j-1
                matInter(i,j) = w1*matResInit(i,k)
            end if
            matResInit(i,j) =  matInter(i,j) + w2*matResInit(i,j)
        end do
    end do

    do i = 1, N
        do j = 1, N
            matrixRes(i,j) = matResInit(i,j)
        end do
    end do

end subroutine mixit
