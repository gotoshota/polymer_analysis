program main
    use prm_list
    use read_prm
    use read_dump
    use calc_poscm
    use calc_inertia_tensor
    
    implicit none
    INTEGER(KIND=8) i, j,k
    DOUBLE PRECISION summation(3)
    

    CALL readprm()
    CALL readdump()
    CALL calcposcm()
    CALL calcinertiatensor()
    
    summation = 0.0d0
    do i = 0, nframe
        do j = 1, nmol
            do k = 1, 3
                summation(k) = summation(k) + inertia_eigenval(k,j,i)
            enddo
        enddo
    enddo
    print *, "# Rg is"
    print *, sum(summation)/ dble(nmol) / dble(nframe + 1)
    print *, summation(1)/ dble(nmol) / dble(nframe + 1)
    print *, summation(2)/ dble(nmol) / dble(nframe + 1)
    print *, summation(3)/ dble(nmol) / dble(nframe + 1)

    

end program
