program msdcalc
    use prm_list
    use read_prm
    use coarse_grain
    use read_dump
    use calc_poscm
    use calc_inertia_tensor
    
    implicit none
    double precision,allocatable :: acf(:,:)
    DOUBLE PRECISION tmp(3)
    double precision t, dummy, tread, tcalc
    integer h, i, j, k, l, m, n
    integer tbegin, tend, CountPerSec, CountMax
   
    call readprm()
    call coarsegrain()
    ALLOCATE(acf(4,0:npoint))
    call system_clock(tbegin, CountPerSec, CountMax)
    call readdump()
    call system_clock(tend)
    tread = real(tend - tbegin) / CountPerSec

    CALL calcposcm()
    CALL calcinertiatensor()
    DEALLOCATE(pos, pos_cm)
    do i = 1, 3
        tmp(i) = 0.0d0
    enddo
    do i = 1, 4
        do j = 0, npoint
            acf(i,j) = 0.0d0
        enddo  
    enddo
    !inertia_tensor = ABS(inertia_tensor)
    do i = 0 , npoint !何個目の点か（時刻差に相当）
        m=0
        do j = 0 , nframe - Targetframe(i) !基準時刻
            tmp = 0.0d0
            do k = 1 , nmol
                do l = 1, 3
                    tmp(l) = tmp(l) &
                    + inertia_tensor(1,l,k,j) * inertia_tensor(1,l,k,j+TargetFrame(i)) &
                    + inertia_tensor(2,l,k,j) * inertia_tensor(2,l,k,j+TargetFrame(i)) &
                    + inertia_tensor(3,l,k,j) * inertia_tensor(3,l,k,j+TargetFrame(i))
                enddo
            enddo
            do l = 1, 3
                tmp(l) = tmp(l) / DBLE(nmol)
                acf(l,i)= acf(l,i) + tmp(l)
                acf(4,i)= acf(4,i) + tmp(l)
            enddo
            m = m + 1
        enddo
        do l = 1, 4 
            acf(l,i) = acf(l,i) / DBLE(m) !基準時刻jに対して平均化
        enddo
        acf(4,i) = acf(4,i) / 3.0d0
    enddo
    do i = 0, npoint
        do j = 1, 4
            acf (j, i) = acf (j, i) / acf (j, 0)
        enddo
    enddo
    call system_clock(tend)
    tcalc = real(tend - tbegin) / CountPerSec
    
    deallocate(inertia_tensor, inertia_eigenval)

    open(15, file='rg_vector_acf.xvg', status='replace')
    write(15,*)'# lambda 1 < lambda 2 < lambda 3'
    write(15,*)'# Coordinate index 3 means longitudinal direction.'
    write(15,*)'# Time, acf of x1, acf of x2, acf of x3, acf total'
    do i = 0 , npoint 
        t = Targetframe(i)* dt * nfreq
        write (15,'(f20.3,2X,4(e40.30,2X))')t,acf(1,i), acf(2,i), acf(3,i), acf(4,i)
    enddo
    close(15)
    open(17, file='log.rg_vecor_acf', status='replace')
    write (17,*)'read time is ',tread,'sec'
    write (17,*)'calclation time is ',tcalc,'sec'
    close(17)

    deallocate(acf)
end program msdcalc
