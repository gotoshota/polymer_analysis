program main
    use prm_list
    use read_prm
    use coarse_grain
    use read_dump
    use calc_poscm
    use calc_inertia_tensor
    
    implicit none
    integer h, i, j, k, l, m, n

    double precision,allocatable :: acf(:), angle(:,:), angle_distri(:), i_angle_distri(:)
    INTEGER n_angle, i_angle
    DOUBLE PRECISION tmp
    DOUBLE PRECISION size_ree, size_major_axis
    double precision t, dummy, tread, tcalc
    integer tbegin, tend, CountPerSec, CountMax
   
    double precision, ALLOCATABLE :: ree(:,:,:)
    INTEGER id_atom

    DOUBLE PRECISION, PARAMETER :: d_angle = 0.10d0
    DOUBLE PRECISION, PARAMETER :: pi = acos(-1.0d0)

!#################################################################

    ! Read and calculate parameters
    call readprm()
    call coarsegrain()
    
    !Read dumpfile
    call system_clock(tbegin, CountPerSec, CountMax)
    call readdump()
    call system_clock(tend)
    tread = real(tend - tbegin) / CountPerSec
    
    ! Calculate C.O.M and Rg Tensor
    CALL calcposcm()
    CALL calcinertiatensor()
    
    ! Calculate R_ee vector
    ALLOCATE(ree(3, nmol, 0:nframe))
    do i = 0, nframe
        do j = 1, nmol
            id_atom = (j - 1) * natom
            ree(:,j,i) = pos(:, id_atom + 1, i) - pos(:, id_atom + natom, i)
        enddo
    enddo
    DEALLOCATE(pos, pos_cm)

    ! Calculate AutoCorrelation Function of R_ee
    ALLOCATE(acf(0:npoint))
    acf = 0.0d0
    do i = 0, npoint
        m = 0
        do j = 0, nframe - TargetFrame(i)
            tmp = 0.0d0
            do k = 1, nmol
                tmp = tmp &
                + ree(1, k, j) * ree (1, k, j + TargetFrame(i)) &
                + ree(2, k, j) * ree (2, k, j + TargetFrame(i)) &
                + ree(3, k, j) * ree (3, k, j + TargetFrame(i))
            enddo
            tmp = tmp / DBLE(nmol)
            acf(i) = acf(i) + tmp
            m = m + 1
        enddo
        acf(i) = acf(i) / DBLE(m)
    enddo
    do i = 0, npoint
        acf(i) = acf(i) / acf(0)
    enddo 

    !Calculate angle of R_ee and Major axis
    ALLOCATE(angle(nmol, 0:nframe))
    do i = 0, nframe
        do j = 1, nmol
            angle(j, i) = inertia_tensor(1,3,j,i) * ree(1,j,i) &
            + inertia_tensor(2,3,j,i) * ree(2,j,i) &
            + inertia_tensor(3,3,j,i) * ree(3,j,i)
            size_ree        = 0.0d0
            size_major_axis = 0.0d0
            do k = 1, 3
                size_ree        = size_ree + ree(k,j,i) ** 2.0d0
                size_major_axis = size_major_axis + inertia_tensor(k,3,j,i)**2.0d0
            enddo
            size_ree        = SQRT(size_ree)
            size_major_axis = SQRT(size_major_axis)
            angle(j,i) = angle(j,i) / size_ree / size_major_axis
            angle(j,i) = ACOS(angle(j,i))
            angle(j,i) = angle(j,i) / pi * 180.0d0
            ! 0 < angle < 90 (上下を区別しない)
            if (angle(j,i) .gt. 90.0d0) angle(j,i) = 180.0d0 - angle(j,i)
        enddo
    enddo
    
    !Calculate distribution function
    n_angle = CEILING(90 / d_angle)
    ALLOCATE(i_angle_distri(0:n_angle))
    do i = 0 , nframe
        do j = 1, nmol
            i_angle = CEILING(angle(j,i) / d_angle)
            i_angle_distri(i_angle) = i_angle_distri(i_angle) + 1
        enddo
    enddo
    !normalization
    ALLOCATE(angle_distri(0:n_angle))
    do i = 0, n_angle
        angle_distri(i) = DBLE(i_angle_distri(i)) / DBLE(nframe + 1) / DBLE(nmol)
    enddo
    
    !check
    tmp = 0.0d0
    do i = 0, n_angle
        tmp = tmp + angle_distri(i)
    enddo
    print *, 'check normalizetion', tmp

    call system_clock(tend)
    tcalc = real(tend - tbegin) / CountPerSec
    
    deallocate(inertia_tensor, inertia_eigenval)
    DEALLOCATE(angle)

    ! output
    open(15, file='ree_acf.xvg', status='replace')
    write(15,*)'# Time, acf'
    do i = 0 , npoint 
        t = Targetframe(i)* dt * nfreq
        write (15,'(f20.3,2X,e20.13)')t,acf(i)
    enddo
    close(15)
    
    open(16, file='angle_ree_majoraxis_distribution.xvg', status='replace')
    WRITE(16,*)'# angle, P(angle)  [0 < angle < 90]'
    do i = 0, n_angle
        t = i * 90.0d0
        WRITE(16,'(f7.4,2X,e20.13)')t,angle_distri(i)
    enddo

    open(17, file='log.rg_vecor_acf', status='replace')
    write (17,*)'read time is ',tread,'sec'
    write (17,*)'calclation time is ',tcalc,'sec'
    close(17)

    deallocate(acf)
end program main
