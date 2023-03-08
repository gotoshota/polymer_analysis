!This program calculate Aspherity parameter (AS).  !The algorism is based on following paper.
!K. Alim and E. Frey, "Shapes of Semiflexible Polymer Rings", PRL (2007).
program main
    use omp_lib
    use prm_list
    use read_prm
    use read_dump
    use calc_poscm
    use calc_inertia_tensor
    
    implicit none
    DOUBLE PRECISION tmp, tmptmp
    DOUBLE PRECISION ave(2)
    DOUBLE PRECISION var(2)
    DOUBLE PRECISION mat(3,3)
    double precision t, dummy, tread, tcalc
    integer h, i, j, k, l, m, n
    integer tbegin, tend, CountPerSec, CountMax
    INTEGER normalize
    INTEGER id_atom
    
    DOUBLE PRECISION, ALLOCATABLE :: A(:,:)
    DOUBLE PRECISION, ALLOCATABLE :: P(:,:)
   
    DOUBLE PRECISION eigenval(3)
    integer, parameter :: lwork = 10000
    double precision :: work(lwork)
    integer :: info
!############################################################################

    call readprm()
    call system_clock(tbegin, CountPerSec, CountMax)
    call readdump()
    call system_clock(tend)
    tread = real(tend - tbegin) / CountPerSec

    CALL calcposcm()
    CALL calcinertiatensor()
    DEALLOCATE(pos, pos_cm)
    !! calculate A and P
    ALLOCATE(A(nmol,0:nframe))
    ALLOCATE(P(nmol,0:nframe))

    !$ omp parallel default(share)
    !$ omp do private(j, eigenval, tmp)
    do i = 0, nframe
        do j = 1, nmol
            eigenval(1) = inertia_eigenval(1, j, i)
            eigenval(2) = inertia_eigenval(2, j, i)
            eigenval(3) = inertia_eigenval(3, j, i)
            
            !! delta
            A(j,i) = (eigenval(1) - eigenval(2)) * (eigenval(1) - eigenval(2))
            A(j,i) = A(j,i) + (eigenval(1) - eigenval(3)) * (eigenval(1) - eigenval(3))
            A(j,i) = A(j,i) + (eigenval(2) - eigenval(3)) * (eigenval(2) - eigenval(3))
            A(j,i) = A(j,i) * 0.50d0 / ( (eigenval(1) + eigenval(2) + eigenval(3)) * (eigenval(1) + eigenval(2) + eigenval(3)))
            !! sigma
            P(j,i) = 2.0d0*eigenval(1) - eigenval(2) - eigenval(3)
            P(j,i) = P(j,i) * (2.0d0*eigenval(2) - eigenval(1) - eigenval(3))
            P(j,i) = P(j,i) * (2.0d0*eigenval(3) - eigenval(1) - eigenval(2))
            tmp = eigenval(1)*eigenval(1) + eigenval(2)*eigenval(2) + eigenval(3)*eigenval(3) &
                - eigenval(1)*eigenval(2) - eigenval(1)*eigenval(3) - eigenval(2)*eigenval(3)
            tmp = tmp**(1.50d0)
            P(j,i) = P(j,i) * 0.50d0 / tmp
        enddo
    enddo
    !$ omp end do
    !$ omp end parallel

    !! calculate average
    ave = 0.0d0
    do i = 0, nframe
        do j = 1, nmol
            ave(1) = ave(1) + A(j,i)
            ave(2) = ave(2) + P(j,i)
        enddo
    enddo
    ave = ave / DBLE(nframe + 1) / DBLE(nmol)

    ! -- calculate variance -- !
    var = 0.0d0 ! -- represents mean square -- !
    do i = 0, nframe
        do j = 1, nmol
            var(1) = var(1) + A(j,i)*A(j,i)
            var(2) = var(2) + P(j,i)*P(j,i)
        end do
    end do
    var(1) = var(1) - ave(1)*ave(1)
    var(2) = var(2) - ave(2)*ave(2)


    call system_clock(tend)
    tcalc = real(tend - tbegin) / CountPerSec
    
    open(15, file='aspherity.txt', status='replace')
    write(15,*)'# aspherity'
    write(15,*)'# average, variance'
    write (15,*)ave(1), var(1)
    close(15)

    open(15, file='prolateness.txt', status='replace')
    write(15,*)'# prolateness'
    write(15,*)'# average, variance'
    write (15,*)ave(2), var(2)
    close(15)

    open(17, file='log.as_param', status='replace')
    write (17,*)'read time is ',tread,'sec'
    write (17,*)'calclation time is ',tcalc,'sec'
    close(17)
    print *,'# aspherity'
    print *,'# A is:'
    print *,ave(1)
    print *,'# prolateness'
    print *,'# P is:'
    print *,ave(2)
    


end program main
