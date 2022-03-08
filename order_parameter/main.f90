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
    DOUBLE PRECISION order_param_simple
    DOUBLE PRECISION test
    double precision t, dummy, tread, tcalc
    integer h, i, j, k, l, m, n
    integer tbegin, tend, CountPerSec, CountMax
    INTEGER normalize
    
    !parameter of lapack
    integer, parameter :: lwork = 100000
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
        do k = 1, 3
            do l = 1, 3
                tmp_tensor(l,k) = 0.0d0
            enddo
        enddo
        do j = 1, nmol
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
                tensor(l,k) = tensor(l,k) + tmp_tensor(l,k) / DBLE(nmol)
            enddo
        enddo
    enddo
    do k = 1, 3
        do l = 1,3
            tensor(l,k) = tensor(l,k) / DBLE(nframe + 1)
        enddo
    enddo
    print*,'=============================='
    print*,tensor

    !digonalization of order parameter tensor
    CALL dsyev('V', 'U', 3, tensor, 3, order_param, work, lwork, info)

    ! #old code
    order_param_simple = 0.0d0
    do i = 0 , nframe 
        tmp = 0.0d0
        normalize = 0
        do j = 1 , nmol
            do k = j + 1 , nmol
                tmptmp = &
                + inertia_tensor(1,3,j,i) * inertia_tensor(1,3,k,i) &
                + inertia_tensor(2,3,j,i) * inertia_tensor(2,3,k,i) &
                + inertia_tensor(3,3,j,i) * inertia_tensor(3,3,k,i)
                tmp = tmp + 1.50d0*tmptmp*tmptmp-0.50d0
    !            print*,k,j
    !            print*,tmptmp
                normalize = normalize + 1
            enddo
        enddo
        tmp = tmp / DBLE(normalize)
        order_param_simple = order_param_simple + tmp
    enddo
    order_param_simple = order_param_simple / DBLE(nframe + 1)

    ! #test code
    !test = 0.0d0
    !normalize = nmol
    !do i = 0 , nframe 
    !    tmp = 0.0d0
    !    do j = 1 , nmol
    !        tmptmp = &
    !        + inertia_tensor(1,3,j,i) * inertia_tensor(1,3,j,i) &
    !        + inertia_tensor(2,3,j,i) * inertia_tensor(2,3,j,i) &
    !        + inertia_tensor(3,3,j,i) * inertia_tensor(3,3,j,i)
    !        tmp = tmp + tmptmp*tmptmp
    !    enddo
    !    tmp = tmp / DBLE(normalize)
    !    order_param_simple = order_param_simple + 0.50d0 * (tmp * 3.0d0 - 1)
    !enddo
    !order_param_simple = order_param_simple / DBLE(nframe + 1)

    call system_clock(tend)
    tcalc = real(tend - tbegin) / CountPerSec
    
    deallocate(inertia_tensor, inertia_eigenval)

    open(15, file='order_parameter.xvg', status='replace')

    write(15,*)'# order param is:'
    do i = 1, 3 
        write (15,*)order_param(i)
    enddo
    write(15,*)'# order param defined as simple dot product is:'
    write(15,*)order_param_simple
!    write(15,*)'# test order param is:'
!    write(15,*)test
    close(15)
    open(17, file='log.order_param', status='replace')
    write (17,*)'read time is ',tread,'sec'
    write (17,*)'calclation time is ',tcalc,'sec'
    close(17)
    
    print*, '# order param is:'
    do i = 1, 3 
        print*, order_param(i)
    enddo
    print *,'# order param defined as simple dot product is:'
    print *,order_param_simple
!    print *,'# test order param is:'
!    print *,test


end program main
