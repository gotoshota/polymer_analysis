Program main
    use prm_list
    use read_prm
    use coarse_grain
    use read_dump
    implicit none

    double precision tmp
    DOUBLE PRECISION tmptmp
    double precision dummy
    DOUBLE PRECISION time
    DOUBLE PRECISION size_ree

    double precision,allocatable :: ree(:,:,:), ree_save(:), p_ree(:)
    DOUBLE PRECISION, ALLOCATABLE :: tacf(:)
    integer i, j, k
    integer id_atom
    INTEGER(KIND=8) normalize
    INTEGER(KIND=4) id_atom_read
    INTEGER(KIND=4) outfile

    CALL  readprm()
    CALL  coarsegrain()
    CALL  readdump()

    ALLOCATE(ree(3, nparticle/2, 0:nframe))
    ALLOCATE(tacf(0:npoint))
    
    do i = 0, npoint
        tacf(i) = 0.0d0
    enddo

    ! calculate Ree vector
    do k = 0, nframe
        do j = 1, nmol
            do i = 1, natom/2
                id_atom_read = (j - 1) * natom + i 
                id_atom      = id_atom_read - natom/2 * (j-1) 
                ree(:, id_atom, k) = pos(:, id_atom_read , k) - pos(:, id_atom_read + CEILING(natom/2.0d0), k)
                size_ree = ree(1, id_atom, k)**2.0d0 + ree(2, id_atom, k)**2.0d0 + ree(3, id_atom, k)**2.0d0
                ree(:, id_atom, k) = ree(:, id_atom, k) / SQRT(size_ree)
            enddo
        enddo
    enddo

    ! calculate TACF of Ree
    do k = 0, npoint
        normalize = 0
        do j = 0, nframe - TargetFrame(k)
            tmp = 0.0d0
            do i = 1, nparticle/2
                tmptmp = DOT_PRODUCT(ree(:,i,TargetFrame(k) + j), ree(:,i,j))
                tmp = tmp + tmptmp 
            enddo
            tacf(k) = tacf(k) + tmp / DBLE(nparticle / 2)
            normalize = normalize + 1
        enddo
        tacf(k) = tacf(k) / DBLE(normalize)
    enddo

    do i = 1, npoint
        tacf(i) = tacf(i) / tacf(0)
    enddo
    tacf(0) = 1.0d0


    outfile = 15
    OPEN(outfile, file='ree_tacf.txt', status = 'replace')
    do i = 0, npoint
        time = Targetframe(i)* dt * nfreq
        write(outfile,'(F15.4,2x,E30.23)') time, tacf(i)
    enddo
    CLOSE(outfile)
    
    deallocate(ree_save)
end program main
