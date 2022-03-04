!This program calculate Aspherity parameter (AS).
!The algorism is based on following paper.
!J. Swetlana, et al., "Diffusion- and reaction-limited cluster aggregation revisited"
!PCCP, 2019.
program main
    use prm_list
    use read_prm
    use read_dump
    use calc_poscm
    use calc_inertia_tensor
    
    implicit none
    DOUBLE PRECISION tmp(3), tmptmp
    DOUBLE PRECISION as
    DOUBLE PRECISION ave(3)
    double precision t, dummy, tread, tcalc
    integer h, i, j, k, l, m, n
    integer tbegin, tend, CountPerSec, CountMax
    INTEGER normalize
    
!############################################################################

    call readprm()
    call system_clock(tbegin, CountPerSec, CountMax)
    call readdump()
    call system_clock(tend)
    tread = real(tend - tbegin) / CountPerSec

    CALL calcposcm()
    CALL calcinertiatensor()
    DEALLOCATE(pos, pos_cm)

    ! calculate average eigenval 
    ave=0.0d0
    do i = 0, nframe
        tmp = 0.0d0
        do j = 1, nmol
            do k = 1,3 
                tmp(k) = tmp(k) + eigenval(l,j,i)*eigenval(l,j,k)
            enddo
        enddo
        do k = 1, 3
            ave(k) = ave(k) + tmp(k) / DBLE(nmol)
        enddo
    enddo
    do k = 1, 3
        ave(k) = ave(k) / DBLE(nframe +1)
    enddo
    
    ! calculate AS parameter
    as = ave(3) - 0.50d0*(ave(1) + ave(2))

    call system_clock(tend)
    tcalc = real(tend - tbegin) / CountPerSec
    
    deallocate(inertia_tensor, inertia_eigenval)
    DEALLOCATE(work)

    open(15, file='AS.txt', status='replace')

    write(15,*)'# aspherity parameter is:'
    write (15,*)as
    close(15)
    open(17, file='log.order_param', status='replace')
    write (17,*)'read time is ',tread,'sec'
    write (17,*)'calclation time is ',tcalc,'sec'
    close(17)
    
    print*, '# AS param is:'
    print*, as


end program main
