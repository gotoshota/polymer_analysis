program main
    use prm_list
    use read_prm
    use read_dump
    use calc_poscm
    use calc_inertia_tensor
    
    implicit none
    INTEGER(KIND=8) i, j
    DOUBLE PRECISION summation
    

    CALL readprm()
    CALL readdump()
    CALL calcposcm()
    CALL calcinertiatensor()
    
    summation = 0.0d0
    do i = 0, nframe
        do j = 1, nmol
            summation = summation + SUM(inertia_eigenval(:,j,i))
        enddo
    enddo
    print *, "# Rg is"
    print *, summation/ dble(nmol) / dble(nframe - 1)
    

end program
