! This program calculate bond correlation function of linear polymer
program main
    use prm_list
    use read_prm
    use read_dump

    implicit none

    INTEGER(KIND=4)                 :: i, j, k
    DOUBLE PRECISION, ALLOCATABLE   :: bond(:,:)
    INTEGER(KIND=4)                 :: id_atom

    CALL readprm()
    CALL readdump()

    ALLOCATE (bond(3, natom))
    do i = 0, nframe
        do j = 1, nmol
            do k = 1, natom-1
                id_atom = (j-1) * natom + k
                bond(:, id_atom) = pos(:, id_atom + 1, i) - pos(:, id_atom, i)
                !! bond(i) = pos(i+1) - pos(i)
            enddo
        enddo
    enddo




end program main
