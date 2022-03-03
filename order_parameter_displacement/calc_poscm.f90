module calc_poscm
    implicit none
    double PRECISION, ALLOCATABLE :: pos_cm(:,:,:)

    contains
    subroutine calcposcm
    use prm_list

    implicit none
    INTEGER(KIND=8) i, j, k, id_atom
    DOUBLE PRECISION tmp(3)

    ALLOCATE(pos_cm(3, nmol, 0:nframe))
    
    do i = 0, nframe
        do j = 1, nmol
            tmp(:) = 0.0d0
            do k = 1, natom
                id_atom = (j - 1) * natom + k
                tmp(:) = tmp(:) + pos(:, id_atom, i)
            enddo
            pos_cm(:, j, i) = tmp(:) / DBLE(natom)
        enddo
    enddo

    end subroutine
end module calc_poscm
