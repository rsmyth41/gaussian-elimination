    program mat_inv_nxn
    implicit none

    real*8, allocatable :: a(:,:), l(:,:), b(:), x(:), id(:,:), ainv(:,:)
    real*8, allocatable :: a_safe(:,:), b_safe(:), id_safe(:,:)
    real*8, allocatable :: check(:,:), checkvec(:)
    real*8 :: z, det, mult, temp, total, totalinv, gaussstart, gaussfinish
    integer :: i, j, k, N, piv, count

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!! DEFINING MATRIX (N X N) !!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    print *, 'Enter dimension of square matrix'
    read *, N

    allocate(a(1:N,1:N),b(1:N),x(1:N),l(1:N,1:N),id(1:N,1:N),ainv(1:N,1:N))
    allocate(a_safe(1:N,1:N),b_safe(1:N),id_safe(1:N,1:N),check(1:N,1:N),checkvec(1:N))

    a=0.0           ! matrix A - to be written over
    a_safe=0.0      ! matrix A - saved version
    l=0.0           ! lower triangular matrix L
    id=0.0          ! identity matrix - to be written over
    id_safe=0.0     ! identity matrix - saved version
    ainv=0.0        ! inverse of matrix A
    b=0.0           ! column vector - to be written over
    b_safe=0.0      ! Saved column vector
    x=0.0           ! Solution vector
    check=0.0       ! matrix to check calculations
    checkvec=0.0    ! vector to check solution

    ! Diagonal elements of L and identity matrix = 1 !
    do i=1, N
        l(i,i)=1.0
        id(i,i)=1.0
        id_safe(i,i)=1.0
    end do

    print *, ' '
    do i=1, N
        do j=1, N
            print 1000, 'Enter value of matrix element A(', i, ',', j, ')'
            read *, z
            a(i,j)= z
            a_safe(i,j)=z
            print *, ' '
        end do
    end do

    print *, ' '
    do i=1, N
        print 1000, 'Enter value of vector element b(', i, ')'
        read *, z
        b(i)= z
        b_safe(i)=z
        print *, ' '
    end do

    print *, 'A='
    do i=1, N
        print *, a_safe(i, 1:N)
    end do
    print *, ' '

    print *, 'b='
    do i=1, N
        print *, b_safe(i)
    end do
    print *, ' '

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!! GAUSSIAN ELIMINATION !!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    call cpu_time(gaussstart)

    do i=1, N-1
        if(a(i,i)==0) then
            ! Finding maximal value on pivot row !
            ! (i-1) is to ensure correct position since we are considering !
            ! an i:N portion of the matrix i.e. 1:N, 3:N, 4:N where i=1,3,4 for example !
            ! Thus, the size of this i:N portion is reduced each time !
            ! and (i-1) ensures we get correct column number / posion in the matrix !
            piv=maxval(maxloc(abs(a(i,i:N))))+(i-1) ! Pivot column number !
            if(abs(a(i,piv))<0.00000000001) then
                print *, 'Matrix is Singular.'
                stop
            end if
            print *, 'Exchanging column', i, 'with column', piv
            do j=1, N
                temp=a(j,piv)
                a(j,piv)=a(j,i)   ! Swaps columns so maximal value is on pivot !
                a(j,i)=temp
            end do
            do j=1, N
                temp=a_safe(j,piv)          ! Bookkeeping to ensure consistency when column swaps occur !
                a_safe(j,piv)=a_safe(j,i)   ! Swaps columns so maximal value is on pivot !
                a_safe(j,i)=temp
            end do
        end if

        do j=i+1, N                  ! From row i+1 to row N i.e down column i!
            mult=a(j,i)/a(i,i)       ! Elimination multiplier !
            l(j,i)=mult              ! Lower triangular matrix elements !
            do k=i, N                ! from column i to N i.e. across row j !
                a(j,k)=a(j,k)-mult*a(i,k)
            end do
            do k=1, N
                id(j,k)=id(j,k)-mult*id(i,k) ! Similar for identity matrix for calculation of a^-1 !
            end do
            b(j)=b(j)-mult*b(i)    ! Same multiplication operations to column vector b !
        end do
    end do

    det=a(1,1) ! Calculates determinant of matrix !
    do i=2, N
        det=det*a(i,i)
    end do
    print *, 'Determinant = ', det
    print *, ' '
    if(det==0.00000000000000) then
        print *, 'Matrix has no inverse'
    go to 100
    end if

    do i=N, 1, -1       ! Back substitution starting from bottom of matrix / vector!
        total=0.0
        do j=i, N
            total=total+x(j)*a(i,j)
        end do
        x(i)=(b(i)-total)/a(i,i)     ! Solution vector !
    end do

    do k=1, N    ! Back substitution for each of the k columns !
        do i=N, 1, -1
            totalinv=0.0
            do j=i, N
                totalinv=totalinv+ainv(j,k)*a(i,j)
            end do
            ainv(i,k)=(id(i,k)-totalinv)/a(i,i) ! The inverse matrix !
        end do
    end do

    call cpu_time(gaussfinish)
    print *, 'Gaussian elimination time = ', gaussfinish-gaussstart, 'seconds'
    print *, ' '

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!! PRINT CHECKS OF MATRICES !!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
100 continue

    print *, 'U='
    do i=1, N
        print *, a(i, 1:N)
    end do
    print *, ' '

    print *, 'L='
    do i=1, N
        print *, l(i, 1:N)
    end do
    print *, ' '

    print *, 'L^-1='
    do i=1, N
        print *, id(i, 1:N)
    end do
    print *, ' '

    print *, 'B='
    do i=1, N
        print *, b(i)
    end do
    print *, ' '

    print *, 'Solution='
    do i=1, N
        print *, x(i)
    end do
    print *, ' '

    print *, 'A^-1='  ! Inverse matrix
    do i=1, N
        print *, ainv(i, 1:N)
    end do
    print *, ' '


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!! SANITY CHECKS !!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Calculations to ensure A=LU, AA^-1=identity,
! LL^-1=identity and Ax=b are satisfied.
! Flags up problems and stops program if not satisfied

!!!!!!!!!!!!!!! A=LU CHECK !!!!!!!!!!!!!!!!!!!

    check=0.0
    do i=1, N
        do j=1, N
            do k=1, N
                check(i,j)=check(i,j)+l(i,k)*a(k,j) ! Matrix multiplication
            end do
        end do
    end do

    print *, 'check='
    do i=1, N
        print *, check(i, 1:N)
    end do
    print *, ' '

    check=check-a_safe

    print *, 'check='
    do i=1, N
        print *, check(i, 1:N)
    end do
    print *, ' '

    count=0
    do i=1, N
        do j=1, N
            if(check(i,j)<0.00000000001) then !10-12
                continue
                count=count+1
            else
                print *, 'Error. A=/=LU.'
                go to 101
            end if
        end do
    end do

    print *, 'A=LU satisfied.'
    print*, count, 'matrix elements checked.'

!!!!!!!!!!!!!!! AA^-1=I CHECK !!!!!!!!!!!!!!!!
101 continue

    check=0.0
    do i=1, N
        do j=1, N
            do k=1, N
                check(i,j)=check(i,j)+ainv(i,k)*a_safe(k,j)
            end do
        end do
    end do

    print *, 'check='
    do i=1, N
        print *, check(i, 1:N)
    end do
    print *, ' '

    check=check-id_safe

    print *, 'check='
    do i=1, N
        print *, check(i, 1:N)
    end do
    print *, ' '

    count=0
    do i=1, N
        do j=1, N
            if(abs(check(i,j))<0.00000000001) then !10-12
                continue
                count=count+1
            else
                print *, 'Error. AA^-1=/=I.'
                go to 102
            end if
        end do
    end do

    print *, 'AA^-1=I satisfied.'
    print *, count, 'matrix elements checked.'

!!!!!!!!!!!!!!! LL^-1=I CHECK !!!!!!!!!!!!!!!!
102 continue

    check=0.0
    do i=1, N
        do j=1, N
            do k=1, N
                check(i,j)=check(i,j)+l(i,k)*id(k,j)
            end do
        end do
    end do

    print *, 'check='
    do i=1, N
        print *, check(i, 1:N)
    end do
    print *, ' '

    check=check-id_safe

    print *, 'check='
    do i=1, N
        print *, check(i, 1:N)
    end do
    print *, ' '

    count=0
    do i=1, N
        do j=1, N
            if(abs(check(i,j))<0.00000000001) then !10-12
                continue
                count=count+1
            else
                print *, 'Error. LL^-1=/=I.'
                go to 103
            end if
        end do
    end do

    print *, 'LL^-1=I satisfied.'
    print *, count, 'matrix elements checked.'
    print *, ' '

!!!!!!!!!!!!!!!! Ax=b CHECK !!!!!!!!!!!!!!!!!!
103 continue

    checkvec=0.0
    do i=1, N
        do j=1, N
            checkvec(i)=checkvec(i)+a_safe(i,j)*x(j)
        end do
    end do

    print *, 'check='
    do i=1, N
        print *, checkvec(i)
    end do
    print *, ' '

    checkvec=checkvec-b_safe
    do i=1, N
        print *, checkvec(i)
    end do

    count=0
    do i=1, N
        if(abs(checkvec(i))<0.00000000001) then
            continue
            count=count+1
        else
            print *, 'Error. Ax=/=b.'
            go to 104
        end if
    end do

    print *, 'Ax=b satisfied.'
    print *, count, 'vector elements checked.'
    print *, ' '

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
104 continue

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!! FORMATTING !!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

1000 format(1x, a, 1x, i2, a, 1x, i2, a)

    end program mat_inv_nxn
