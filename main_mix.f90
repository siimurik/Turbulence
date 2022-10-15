program main
!==========================================================
    implicit none
    integer(8)              :: i, j, l, nn, k
    integer(8), parameter   :: N = 512
    real(8)                 :: a, amo, w1, w2
    real(8), parameter      :: pi = 4.D0*atan(1.D0)
    real(8), dimension(N,N) :: c1, c2, c3, c4, c_new, c_miin
    real(8), dimension(N,N) :: matrixInit, matrixRes
!==========================================================
! Initial conditions
    a = 400
    amo = 0.7
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
            l = a*sin(i*2*pi/N)
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
    call mixit_hor(matrixInit, matrixRes)
! END OF FIRST ITERATION
!========================================================================
! BEGINNING OF SECOND ITERATION
    do i = 1, N ! y
        do j = 1, N ! x
            l = a*sin(i*2*pi/N)
            !print *, "l =", l
            c3(i,j) = matrixRes(i, iand((j+l-1), (N-1))+1)
            !print *, (j+l-1), "AND", N-1, "  =  ", iand((j+l-1), (N-1))+1
        end do
    end do

    do i = 1, N
        do j = 1, N
            matrixInit(i,j) = c3(i,j)
        end do
    end do
!
    call mixit_hor(matrixInit, matrixRes)
! END OF SECOND ITERATION
!========================================================================
! BEGINNING OF THIRD ITERATION

    !do i = 1, N ! y
    !    do j = 1, N ! x
    !        l = a*sin(i*2*pi/N)
    !        !print *, "l =", l
    !        c4(i,j) = matrixRes(i, iand((j+l-1), (N-1))+1)
    !        !print *, (j+l-1), "AND", N-1, "  =  ", iand((j+l-1), (N-1))+1
    !    end do
    !end do
!
    !do i = 1, N
    !    do j = 1, N
    !        matrixInit(i,j) = c4(i,j)
    !    end do
    !end do
!!
    !call mixit_hor(matrixInit, matrixRes)
! END OF THIRD ITERATION
!========================================================================

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

subroutine setup_mat(matrixInit, matrixRes)
!===========================================================================
    implicit none
    integer(8)              :: i, j, l, k
    integer(8), parameter   :: N = 512
    real(8)                 :: a
    real(8), parameter      :: pi = 4.D0*atan(1.D0)
    real(8), dimension(N,N) :: c_miin
    real(8), dimension(N,N) :: matrixInit, matrixRes
!===========================================================================
    do i = 1, N ! y
        do j = 1, N ! x
            l = a*sin(i*2*pi/N)
            !print *, "l =", l
            matrixInit(i,j) = matrixRes(i, iand((j+l-1), (N-1))+1)
            !print *, (j+l-1), "AND", N-1, "  =  ", iand((j+l-1), (N-1))+1
        end do
    end do
end subroutine setup_mat


subroutine mixit(matrixInit, matrixRes)
!=========================================================================
    implicit none
    integer(8)              :: i, j, nn, k
    integer(8), parameter   :: N = 512
    real(8)                 :: a, amo, w1, w2
    real(8), parameter      :: pi = 4.D0*atan(1.D0)
    real(8), dimension(N,N) :: c1, c2, c_new, c_miin
    real(8), dimension(N,N) :: matrixInit, matrixRes
!=========================================================================
! Introduce new matrix c_miin as matrix with
! elements set to zero to avoid reusing old values.
    c_miin = 0.0
!------------------------------------------------------------------------
! FIRST LOOP
    do i = 1, N
        do j = 1, N
            w1 = a*cos((j-1)*pi*2/N) - floor(a*cos((j-1)*pi*2/N))
            !print *, ' eqv w1 = ', w1
            w2 = 1 - w1
            !print *, "w2 =", w2
            if (j == 1) then                        ! Whole point of this segment
                k = 4                               ! is to start writing column
                c_miin(i,j) = w1*matrixInit(i,k)    ! elements from the end and
            else                                    ! then continue on from from
                k = j-1                             ! first index. NB! Why this is 
                c_miin(i,j) = w1*matrixInit(i,k)    ! necessary in the fisrt place, 
            end if                                  ! I legitimately have no clue :(
            c2(i,j) =  c_miin(i,j) + w2*matrixInit(i,j)
        end do
    end do
! Print out elements of c2
    !do i = 1, N
    !    print *, c2(i,:)
    !end do
!------------------------------------------------------------------------
! Introduce new matrix c_new as matrix with
! elements set to zero to avoid reusing old values.
    c_new = 0.0
    do i = 1, N
        do j = 1, N
            nn = n+1
            c_new(i,j) = c2(nn-i,nn-j)
        end do
    end do
! Print out the result.
    !do i = 1, N
    !    print *, c_new(i, :)
    !end do
!------------------------------------------------------------------
! Introduce new matrix c_miin as matrix with
! elements set to zero to avoid reusing old values.
    c_miin = 0.0
    matrixRes = 0.0
!------------------------------------------------------------------
! SECOND LOOP
    do i = 1, N
        do j = 1, N
            w1 = a*cos((j-1)*pi*2/N) - floor(a*cos((j-1)*pi*2/N)) 
            !print *, ' eqv w1 new= ', w1
            w2 = w1 - 1
            !print *, "w2 =", w2
            if (j == 1) then                ! same story here 
                k = 4                       ! "It just works" - Todd Howard
                c_miin(i,j) = w1*c2(i,k)
            else
                k = j-1
                c_miin(i,j) = w1*c2(i,k)
            end if
            matrixRes(i,j) =  c_miin(i,j) + w2*c2(i,j)
        end do
    end do
! Print out matrixRes
    !do i = 1, N
    !    print *, matrixRes(i, :)
    !end do

end subroutine mixit

subroutine mixit_hor(matrixInit, matrixRes)
    !=========================================================================
        implicit none
        integer(8)              :: i, j, nn, k
        integer(8), parameter   :: N = 512
        real(8)                 :: a, w1, w2
        real(8), parameter      :: pi = 4.D0*atan(1.D0)
        real(8), dimension(N,N) :: c2, c_new, c_miin, c_new_T
        real(8), dimension(N,N) :: matrixInit, matrixRes
    !=========================================================================
    ! Introduce new matrix c_miin as matrix with
    ! elements set to zero to avoid reusing old values.
        c_miin = 0.0
    !------------------------------------------------------------------------
    ! FIRST LOOP
        do i = 1, N
            do j = 1, N
                w1 = a*cos((j-1)*pi*2/N) - floor(a*cos((j-1)*pi*2/N))
                !print *, ' eqv w1 = ', w1
                w2 = 1 - w1
                !print *, "w2 =", w2
                if (j == 1) then                        ! Whole point of this segment
                    k = 4                               ! is to start writing column
                    c_miin(i,j) = w1*matrixInit(k,j)    ! elements from the end and
                else                                    ! then continue on from from
                    k = j-1                             ! first index. NB! Why this is 
                    c_miin(i,j) = w1*matrixInit(k,j)    ! necessary in the fisrt place, 
                end if                                  ! I legitimately have no clue :(
                c2(i,j) =  c_miin(i,j) + w2*matrixInit(i,j)
            end do
        end do
    ! Print out elements of c2
        !do i = 1, N
        !    print *, c2(i,:)
        !end do
    !------------------------------------------------------------------------
    ! Introduce new matrix c_new as matrix with
    ! elements set to zero to avoid reusing old values.
        c_new = 0.0
        do i = 1, N
            do j = 1, N
                nn = n+1
                c_new(i,j) = c2(nn-i,nn-j)
            end do
        end do

        !c_new_T = transpose(c_new)
    ! Print out the result.
        !do i = 1, N
        !    print *, c_new(i, :)
        !end do
    !------------------------------------------------------------------
    ! Introduce new matrix c_miin as matrix with
    ! elements set to zero to avoid reusing old values.
        c_miin = 0.0
        matrixRes = 0.0
    !------------------------------------------------------------------
    ! SECOND LOOP
        do i = 1, N
            do j = 1, N
                w1 = a*cos((j-1)*pi*2/N) - floor(a*cos((j-1)*pi*2/N)) 
                !print *, ' eqv w1 new= ', w1
                w2 = w1 - 1
                !print *, "w2 =", w2
                if (j == 1) then                ! same story here 
                    k = 4                       ! "It just works" - Todd Howard
                    c_miin(i,j) = w1*c2(k,j)
                else
                    k = j-1
                    c_miin(i,j) = w1*c2(k,j)
                end if
                matrixRes(i,j) =  c_miin(i,j) + w2*c2(i,j)
            end do
        end do
    ! Print out matrixRes
        !do i = 1, N
        !    print *, matrixRes(i, :)
        !end do

        matrixRes = transpose(matrixRes)
    
    end subroutine mixit_hor
