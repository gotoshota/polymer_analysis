!######################################
!  This module needs LAPACK
!######################################
!  Add option "-llapack -lblas"
!  when complile
!######################################
module calc_inertia_tensor
    implicit none
    DOUBLE PRECISION, ALLOCATABLE :: inertia_tensor(:,:,:,:) !! output as eigenvector 
    DOUBLE PRECISION, ALLOCATABLE :: inertia_eigenval(:,:,:) !! output as eigenvaue 


    contains
    subroutine calcinertiatensor
    use prm_list
    use calc_poscm 

    implicit none 
    DOUBLE PRECISION eigenval(3), eigenvec(3,3)
    INTEGER(KIND=8) i, j, k, l, m
    INTEGER(KIND=8) id_atom
    integer, parameter :: lwork = 10000
    double precision :: work(lwork)
    integer :: info
    INTEGER(KIND=4) target_frame

    ALLOCATE(inertia_tensor(3, 3, nmol, nframe)) !! 1=x, 2=y, 3=z
    ALLOCATE(inertia_eigenval(3,nmol, nframe)) !! 1=x, 2=y, 3=z
    inertia_tensor(:,:,:,:) = 0.0d0

    !! calculate inertia ternsor
    do i = 1, nframe
        do j = 1, nmol 
            do k = 1, natom 
                id_atom = (j - 1) * natom + k
                do l = 1, 3 !xyz
                    do m = 1, 3 !xyz
                        inertia_tensor(m, l, j, i) = &
                        inertia_tensor(m, l, j, i) + &
                        (pos(m, id_atom, i) - pos_cm(m, j, i)) * &
                        (pos(l, id_atom, i) - pos_cm(l, j, i))
                    enddo
                enddo
            enddo 
            inertia_tensor(:, :, j, i) = inertia_tensor(:, :, j, i) / DBLE(natom)
        enddo
    enddo

    !! call LAPACK to diagonalize
    do i = 1, nframe 
        do j = 1, nmol
            do k = 1, 3
                do l = 1,3
                    eigenvec(l, k) = inertia_tensor(l, k, j, i)
                enddo
            enddo
            CALL dsyev('V', 'U', 3, eigenvec, 3, eigenval, work, lwork, info)
            do k = 1, 3
                do l = 1,3
                    inertia_tensor(l, k, j, i) = eigenvec(l, k)
                enddo
            enddo
        enddo
    enddo

    end subroutine
end module
