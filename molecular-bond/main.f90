program main
    !$ use omp_lib
    use prm_list
    use read_prm
    use coarse_grain
    use read_dump
    use calc_poscm

    implicit none

    INTEGER(KIND=4)                 :: i, j, k, l
    INTEGER(KIND=4)                 :: ii
    INTEGER(KIND=4)                 :: id_particle
    INTEGER(KIND=4)                 :: total_num_bond

    INTEGER(KIND=4), ALLOCATABLE    :: num_bond(:,:)
    INTEGER(KIND=2), ALLOCATABLE    :: flag_bond(:,:,:)

    DOUBLE PRECISION                :: normalize
    DOUBLE PRECISION                :: rg2_ave
    DOUBLE PRECISION                :: tmp
    DOUBLE PRECISION                :: bond_length
    DOUBLE PRECISION                :: dpos(3)
    DOUBLE PRECISION                :: distance
    DOUBLE PRECISION                :: prefactor
    DOUBLE PRECISION                :: summation

    DOUBLE PRECISION, ALLOCATABLE   :: rg2(:,:)
    DOUBLE PRECISION, ALLOCATABLE   :: p_num_bond(:)
    DOUBLE PRECISION, ALLOCATABLE   :: cf_bond(:)

    INTEGER(KIND=4), PARAMETER      :: outfile = 17
    CHARACTER(LEN=256)              :: outfilename
    CHARACTER(LEN=256)              :: command
    CHARACTER(LEN=256)              :: chara

    CALL readprm()
    CALL coarsegrain()
    CALL readdump()
    CALL calcposcm()

    ALLOCATE(rg2(nmol, 0:nframe))

    ! -- calculate Rg^2 -- !
    rg2 = 0.0d0
    normalize = 1/DBLE(natom)
    !$ omp palallel
    !$ omp do private (j,k,id_particle)
    do i = 0, nframe
        do j = 1, nmol
            do k = 1, natom
                id_particle = (j-1)*natom + k
                rg2(j,i) = rg2(j,i) + (pos(1, id_particle, i) - pos_cm(1, j, i))*(pos(1, id_particle, i) - pos_cm(1, j, i)) &
                + (pos(2, id_particle, i) - pos_cm(2, j, i))*(pos(2, id_particle, i) - pos_cm(2, j, i)) &
                + (pos(3, id_particle, i) - pos_cm(3, j, i))*(pos(3, id_particle, i) - pos_cm(3, j, i)) 
            end do
            rg2(j,i) = rg2(j,i) * normalize
        end do
    end do
    !$ omp end do
    DEALLOCATE(pos)

    ! -- calculate the average of Rg^2 -- !
    !$ omp single 
    normalize = 1 / DBLE(nmol)
    rg2_ave = 0.0d0
    do i = 0, nframe
        tmp = 0.0d0
        do j = 1, nmol
            tmp = tmp + rg2(j,i)
        end do
        tmp = tmp * normalize
        rg2_ave = rg2_ave + tmp
    end do
    rg2_ave = rg2_ave / DBLE(nframe + 1)
    outfilename = "rg2.txt"
    OPEN(outfile, file=outfilename, status="replace", form="formatted")
        WRITE(outfile,*) rg2_ave
    CLOSE(outfile)
    !$ omp end single

    ! -- calculate the number of bonded molecules-- !
    ALLOCATE(num_bond(nmol, 0:nframe))
    ALLOCATE(flag_bond(nmol, nmol, 0:nframe))
    ALLOCATE(cf_bond(1:npoint))
    ALLOCATE(p_num_bond(0:nmol))
    command = "if [ ! -d CF ] ; then mkdir -p CF; fi"
    CALL system(command)
    command = "if [ ! -d PF ] ; then mkdir -p PF; fi"
    CALL system(command)
    do ii = 1, 20
        prefactor = ii / 10.0d0
        bond_length = SQRT(rg2_ave) * prefactor
        num_bond = 0
        flag_bond = 0

        !$ omp do private(i,j,k,dpos,distance)
        do i = 0, nframe
            do j = 1, nmol
                do k = 1, nmol
                    if (j .ne. k) then
                        dpos(:) = pos_cm(:,k,i) - pos_cm(:,j,i)
                        dpos(:) = ABS(dpos(:) - box_l * NINT(dpos(:) / box_l))
                        distance = SQRT(dpos(1)*dpos(1) + dpos(2)*dpos(2) + dpos(3)*dpos(3))
                    end if 
                    if (distance .le. bond_length) then
                        num_bond(j,i) = num_bond(j,i) + 1
                        flag_bond(k,j,i) = 1
                    end if 
                end do
            end do 
        end do
        !$ omp end do

        !$ omp single
        total_num_bond = 0
        do i = 0, nframe
            do j = 1, nmol
                total_num_bond = total_num_bond + num_bond(j,i)
            end do
        end do
        !$ omp end single


        ! -- calculate the probability of the number of neighbor polymer -- !
        p_num_bond = 0.0d0
        !$ omp do private(j)
        do i = 0, nframe
            do j = 1, nmol
                p_num_bond(num_bond(j,i)) = p_num_bond(num_bond(j,i)) + 1
            end do
        end do
        summation = 0.0d0
        do i = 0, nmol
            p_num_bond(i) = p_num_bond(i) / (nmol * (nframe + 1))
            summation = summation + p_num_bond(i)
        end do
        print *, summation
        !$ omp end do

        !$ omp single
        WRITE(chara, "(I3.3)") ii
        outfilename = "PF/" // TRIM(ADJUSTL(chara))// ".txt"
        OPEN(outfile, file=outfilename, status="replace", form="formatted")
            do i = 0, nmol
                WRITE(outfile, "(I4, 1X, F15.7)") i, p_num_bond(i)
            end do
        CLOSE(outfile)
        !$ omp end single

        ! -- calculate the life time of molecular bond -- !
        cf_bond = 0.0d0
        !$ omp do private(j,k,l)
        do i = 1, npoint
            do j = 0, nframe - TargetFrame(i)
                do k = 1, nmol 
                    do l = 1, nmol
                        cf_bond(i) = cf_bond(i) + flag_bond(l,k,j) * flag_bond(l,k,j+TargetFrame(i))
                    end do
                end do    
            end do
            cf_bond(i) = cf_bond(i) / DBLE(nframe - TargetFrame(i) + 1)
        end do
        !$ omp end do

        !$ omp single
        do i = 1, npoint
            cf_bond(i) = cf_bond(i) / DBLE(total_num_bond) * DBLE(nframe + 1)
        end do

        WRITE(chara, "(I3.3)") ii
        outfilename = "CF/" // TRIM(ADJUSTL(chara))// ".txt"
        OPEN(outfile, file=outfilename, status="replace", form="formatted")
            WRITE(outfile, *) "# The difinision of bond length = ", bond_length
            WRITE(outfile, *) "# The radius of gyration sqrt(Rg^2) = ", SQRT(rg2_ave)
            WRITE(outfile, *) "# The ratio of bondlength and Rg = ", bond_length/SQRT(rg2_ave)
            do i = 1, npoint
                WRITE(outfile, "(F15.4, 1X, F15.7)") TargetFrame(i)*nfreq*dt, cf_bond(i)
            end do
        CLOSE(outfile)
        !$ omp end single
    end do
    !$ omp end parallel
end program main
