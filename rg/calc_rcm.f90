module calc_rcm
    implicit none
    double PRECISION, ALLOCATABLE :: rcm(:,:,:)

    contains
    subroutine calcrcm
    use prm_list

    implicit none
    INTEGER(KIND=8) i, j, k, id_atom
    DOUBLE PRECISION tmp(3)

    ALLOCATE(rcm(3, nmol, 0:nframe))
    
    do i = 0, nframe
        do j = 1, nmol
            tmp(:) = 0.0d0
            do k = 1, natom
                id_atom = (j - 1) * natom + k
                tmp(:) = tmp(:) + pos(:, id_atom, i)
            enddo
            rcm(:, j, i) = tmp(:) / DBLE(natom)
        enddo
    enddo

    end subroutine
end module calc_rcm
