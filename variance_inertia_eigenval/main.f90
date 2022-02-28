program main
    use prm_list
    use read_prm
    use read_dump
    use calc_poscm
    use calc_inertia_tensor
    
    implicit none
    INTEGER(KIND=8) i, j, k
    INTEGER outfile
    DOUBLE PRECISION, ALLOCATABLE :: average(:,:)
    DOUBLE PRECISION, ALLOCATABLE :: variance(:,:)
    

    CALL readprm()
    CALL readdump()
    CALL calcposcm()
    CALL calcinertiatensor()
    DEALLOCATE(pos,pos_cm,inertia_tensor)
    
    ! Calculate average of eigenvalue
    ALLOCATE(average(3,nmol))
    average = 0.0d0
    do i = 0, nframe
        do j = 1, nmol
            do k = 1, 3
                average(k,j) = average(k,j) + inertia_eigenval(k,j,i)
            enddo
        enddo
    enddo
    do i = 1, nmol
        do j = 1, 3
        average(j,i) = average(j,i) / DBLE(nframe + 1)
        enddo
    enddo

    ! Calculate variance of eigenval 
    ALLOCATE(variance(3,nmol))
    variance = 0.0d0
    do i = 0, nframe
        do j = 1, nmol
            do k = 1, 3
                variance(k,j) = variance(k,j) + (inertia_eigenval(k,j,i) - average(k,j)) ** 2.0d0
            enddo
        enddo
    enddo
    do i = 1, nmol
        do j = 1,3
        variance(j,i) = variance(j,i) / DBLE(nframe + 1)
        enddo
    enddo
    DEALLOCATE(inertia_eigenval)

    outfile = 11
    open(outfile, file="sd_eigenval.xvg")
    WRITE(outfile,*) "# ID (mol), standard deviates of lambda1, lambda2, lambda3"
    do i = 1, nmol
        WRITE(outfile, '(I5, 1x, 3(e10.3, 1x))') i, SQRT(variance(1,i))/average(1,i), SQRT(variance(2,i))/average(2,i) &
        ,SQRT(variance(3,i))/average(3,i)
    enddo
    close(outfile)
    DEALLOCATE(variance,average)

    !print *,"The average of:"
    !print *,"lambda 1 =",average(1)
    !print *,"lambda 2 =",average(2)
    !print *,"lambda 3 =",average(3)
    !print *,"The variance of:"
    !print *,"lambda 1 =",variance(1)
    !print *,"lambda 2 =",variance(2)
    !print *,"lambda 3 =",variance(3)
    !print *,"The standard deviation of:"
    !print *,"lambda 1 =",SQRT(variance(1))
    !print *,"lambda 2 =",SQRT(variance(2))
    !print *,"lambda 3 =",SQRT(variance(3))


end program
