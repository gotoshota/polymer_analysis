program main
    use prm_list
    use read_prm
    use read_dump
    use calc_poscm
    !use calc_inertia_tensor
    
    implicit none
    DOUBLE PRECISION tmp, tmptmp
    DOUBLE PRECISION, ALLOCATABLE :: order_para(:)
    DOUBLE PRECISION, ALLOCATABLE :: dpos(:,:,:)
    DOUBLE PRECISION ave
    double precision t, dummy, tread, tcalc
    integer h, i, j, k, l, m, n
    integer tbegin, tend, CountPerSec, CountMax
    INTEGER normalize
   
    call readprm()
    call system_clock(tbegin, CountPerSec, CountMax)
    call readdump()
    call system_clock(tend)
    tread = real(tend - tbegin) / CountPerSec

    CALL calcposcm()
    !CALL calcinertiatensor()
    DEALLOCATE(pos, pos_cm)
    ALLOCATE(order_para(nframe))
    ALLOCATE(dpos(3,nmol,nframe)

    !calculate dpos
    do i = 1, nframe
        do j = 1, nmol
            dpos(:, j, i) = pos(:, j, i) - pos(:, j, i-1)
        enddo
    enddo
    

    normalize = nmol * (nmol - 1) / 2 
    do i = 1, nframe
        tmp = 0.0d0
        do j = 1, nmol - 1
            do k = j, nmol
                tmptmp = dpos(1, j, i) * dpos(1, k, i) &
                + dpos(2, j, i) * dpos(2, k, i) &
                + dpos(3, j, i) * dpos(3, k, i)
                tmptmp = tmptmp / SQRT(DOT_PRODUCT(dpos(:,j,i),dpos(:,j,i))) &
                / SQRT(DOT_PRODUCT(dpos(:,k,i),dpos(:,k,i)))
                tmp = tmp + tmptmp
            enddo
        enddo
        tmp = tmp / DBLE(normalize)
        order_para(i) = 0.50d0 * (3.0d0 * tmp - 1)
    enddo
    ave = 0.0d0
    do i = 1, nframe
        ave = ave + order_para(i)
    enddo
    ave = ave / DBLE(nframe)

    call system_clock(tend)
    tcalc = real(tend - tbegin) / CountPerSec
    
    !deallocate(inertia_tensor, inertia_eigenval)

    open(15, file='order_parameter.xvg', status='replace')

    write(15,*)'# average of order parameter for displacement of center of mass =', ave
    write(15,*)'# Time, order param (dispalcement_COM)'
    do i = 1 , nframe 
        t = i* dt * nfreq
        write (15,'(f20.3,2X,4(e40.30,2X))')t,order_para(i)
    enddo
    close(15)
    open(17, file='log.order_param', status='replace')
    write (17,*)'read time is ',tread,'sec'
    write (17,*)'calclation time is ',tcalc,'sec'
    close(17)

end program main
