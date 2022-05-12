program main
    use prm_list
    use read_prm
    use coarse_grain
    use read_dump
    use calc_poscm
    use calc_inertia_tensor
    
    implicit none
    INTEGER(KIND=8) i, j, k
    DOUBLE PRECISION summation
    DOUBLE PRECISION rg2_ave
    DOUBLE PRECISION tmp
    DOUBLE PRECISION, ALLOCATABLE :: rg2(:,:)
    DOUBLE PRECISION, ALLOCATABLE :: tacf(:)

    INTEGER(KIND=4) outfile
    DOUBLE PRECISION time

    CALL readprm()
    CALL readdump()
    CALL coarsegrain()
    CALL calcposcm()
    CALL calcinertiatensor()

    ALLOCATE(rg2(nmol,0:nframe))
    summation = 0.0d0
    do i = 0, nframe
        do j = 1, nmol
            rg2(j,i) = SUM(inertia_eigenval(:,j,i))
            summation = summation + rg2(j,i)
        enddo
    enddo
    
    rg2_ave = summation/ dble(nmol) / dble(nframe - 1)

    ALLOCATE(tacf(0:npoint))
    do i = 0, npoint
        tacf(i) = 0.0d0
    enddo
    do i = 0, npoint 
        do j = 0, nframe - TargetFrame(i) 
            tmp = 0.0d0
            do k = 1, nmol
                tmp = tmp + (rg2(k, j + TargetFrame(i)) - rg2_ave) * (rg2(k, j) - rg2_ave)
            enddo
            tacf(i) = tacf(i) + tmp / DBLE(nmol)
        enddo
        tacf(i) = tacf(i) / DBLE(nframe - TargetFrame(i) + 1)
    enddo
    do i = 1, npoint
        tacf(i) = tacf(i) / tacf(0)
    enddo
    tacf(0) = 1.0d0
  
    outfile = 15
    OPEN(outfile, file='rg2_fluc_tacf.txt', status = 'replace')
    do i = 0, npoint
        time = Targetframe(i)* dt * nfreq
        write(outfile,'(F15.4,2x,E30.23)') time, tacf(i)
    enddo
    CLOSE(outfile)

end program
