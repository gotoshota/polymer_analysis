program main
    use prm_list
    use read_prm
    use read_dump
    use calc_poscm
    use calc_inertia_tensor
    
    implicit none
    DOUBLE PRECISION tmp, tmptmp, tmp_tensor(3,3)
    DOUBLE PRECISION order_param(3)
    DOUBLE PRECISION tensor(3,3)
    DOUBLE PRECISION ave
    double precision t, dummy, tread, tcalc
    integer h, i, j, k, l, m, n
    integer tbegin, tend, CountPerSec, CountMax
    INTEGER normalize
    
    !parameter of lapack
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

    ! make order parameter tensor
    do i = 1, 3
        do j = 1,3
            tensor(j,i) = 0.0d0
        enddo
    enddo
    do i = 0, nframe
        do j = 1, nmol
            do k = 1, 3
                do l = 1, 3
                    tmp_tensor(l,k) = 0.0d0
                enddo
            enddo
            do k = 1, 3
                do l = 1, 3
                    if (k .eq. l) then
                        tmp_tensor(l,k) = tmp_tensor(l,k) + inertia_tensor(l,3,j,i) * inertia_tensor(k,3,j,i) * 1.50d0 - 0.50d0
                    else
                        tmp_tensor(l,k) = tmp_tensor(l,k) + inertia_tensor(l,3,j,i) * inertia_tensor(k,3,j,i) * 1.50d0
                    endif
                enddo
            enddo
        enddo
        do k = 1, 3
            do l = 1, 3
                tensor(l,k) = tensor (l,k) + tmp_tensor(l,k) / DBLE(nmol)
            enddo
        enddo
    enddo
    do k = 1, 3
        do l = 1,3
            tensor(l,k) = tensor(l,k) / DBLE(nframe + 1)
        enddo
    enddo

    !digonalization of order parameter tensor
    CALL dsyev('V', 'U', 3, tensor, 3, order_param, work, lwork, info)

    ! #old code
    !normalize = nmol * (nmol -1) / 2
    !do i = 0 , nframe 
    !    tmp = 0.0d0
    !    do j = 1 , nmol
    !        do k = j + 1 , nmol
    !            tmptmp = &
    !            + inertia_tensor(1,3,j,i) * inertia_tensor(1,3,k,i) &
    !            + inertia_tensor(2,3,j,i) * inertia_tensor(2,3,k,i) &
    !            + inertia_tensor(3,3,j,i) * inertia_tensor(3,3,k,i)
    !            tmp = tmp + tmptmp*tmptmp
    !        enddo
    !    enddo
    !    tmp = tmp / DBLE(normalize)
    !    order_para(i) = 0.50d0 * (tmp * 3.0d0 - 1)
    !enddo
    !ave = 0.0d0
    !do i = 0, nframe
    !    ave = ave + order_para(i)
    !enddo
    !ave = ave / DBLE(nframe + 1)

    call system_clock(tend)
    tcalc = real(tend - tbegin) / CountPerSec
    
    deallocate(inertia_tensor, inertia_eigenval)
    DEALLOCATE(work)

    open(15, file='order_parameter.xvg', status='replace')

    write(15,*)'# order param is:'
    do i = 1, 3 
        write (15,*)order_para(i)
    enddo
    close(15)
    open(17, file='log.order_param', status='replace')
    write (17,*)'read time is ',tread,'sec'
    write (17,*)'calclation time is ',tcalc,'sec'
    close(17)
    
    print*, '# order param is:'
    do i = 1, 3 
        print*, order_para(i)
    enddo


end program main
