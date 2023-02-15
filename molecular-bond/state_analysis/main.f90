! This programs calculates the time how long the state of N_o maintain.

program main
    !$ use omp_lib
    use prm_list
    use read_prm
    use read_dump
    use calc_poscm

    implicit none

    INTEGER(KIND=4)                 :: i, j, k, l
    INTEGER(KIND=4)                 :: ii
    INTEGER(KIND=4)                 :: id_particle
    INTEGER(KIND=4)                 :: max_n_o

    INTEGER(KIND=4), ALLOCATABLE    :: i_normalize(:)
    INTEGER(KIND=4), ALLOCATABLE    :: total_num_bond(:)
    INTEGER(KIND=4), ALLOCATABLE    :: num_bond(:,:)
    INTEGER(KIND=4), ALLOCATABLE    :: int_tmp_ave_time_of_state(:)
    INTEGER(KIND=4), ALLOCATABLE    :: time_of_state(:,:)
    INTEGER(KIND=1), ALLOCATABLE    :: flag_bond(:,:,:)

    DOUBLE PRECISION                :: normalize
    DOUBLE PRECISION                :: rg2_ave
    DOUBLE PRECISION                :: tmp
    DOUBLE PRECISION                :: bond_length
    DOUBLE PRECISION                :: dpos(3)
    DOUBLE PRECISION                :: distance
    DOUBLE PRECISION                :: prefactor
    DOUBLE PRECISION                :: summation
    DOUBLE PRECISION                :: n_ave

    DOUBLE PRECISION, ALLOCATABLE   :: rg2(:,:)
    DOUBLE PRECISION, ALLOCATABLE   :: p_num_bond(:)
    DOUBLE PRECISION, ALLOCATABLE   :: p_time_of_state(:,:)
    DOUBLE PRECISION, ALLOCATABLE   :: ave_time_of_state(:)

    INTEGER(KIND=4), PARAMETER      :: outfile = 17
    CHARACTER(LEN=256)              :: outfilename
    CHARACTER(LEN=256)              :: command
    CHARACTER(LEN=256)              :: chara
    CHARACTER(LEN=256)              :: chara_num_ol

    CALL readprm()
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
    ALLOCATE(total_num_bond(0:nmol+1))
    ALLOCATE(flag_bond(nmol, nmol, 0:nframe))
    ALLOCATE(p_num_bond(0:nmol))
    ALLOCATE(i_normalize(0:nmol+1))
    ALLOCATE(time_of_state(nmol, 0:nframe))
    ALLOCATE(p_time_of_state(nframe, 0:nmol+1))
    ALLOCATE(ave_time_of_state(0:nmol+1))
    ALLOCATE(int_tmp_ave_time_of_state(0:nmol+1))

    command = "if [ ! -d maintain_time ] ; then mkdir -p maintain_time; fi"
    CALL system(command)
    command = "if [ ! -d pdf_maintain_time ] ; then mkdir -p pdf_maintain_time; fi"
    CALL system(command)
    command = "if [ ! -d PF ] ; then mkdir -p PF; fi"
    CALL system(command)

    ! -- the variable "ii" is the prefactor of r_overlap (i.e. bond length) -- !
    do ii = 1, 20
        prefactor = ii / 10.0d0
        bond_length = SQRT(rg2_ave) * prefactor
        num_bond = 0
        flag_bond = 0

        ! -- Determin overlapped (bonded)  or not --  !
        !$ omp do private(i,j,k,dpos,distance)
        do i = 0, nframe
            do j = 1, nmol
                do k = 1, nmol
                    if (j .ne. k) then
                        dpos(:) = pos_cm(:,k,i) - pos_cm(:,j,i)
                        dpos(:) = ABS(dpos(:) - box_l * NINT(dpos(:) / box_l))
                        distance = SQRT(dpos(1)*dpos(1) + dpos(2)*dpos(2) + dpos(3)*dpos(3))
                    else 
                        flag_bond(k,j,i) = 0
                    end if 
                    if (distance .le. bond_length) then
                        num_bond(j,i) = num_bond(j,i) + 1
                        flag_bond(k,j,i) = 1
                    else 
                        flag_bond(k,j,i) = 0
                    end if 
                end do
            end do 
        end do
        !$ omp end do

        !$ omp single
        ! -- calculate the number of polymer whose N_o = n -- !
        ! -- nmol + 1 represents the summation of N_o,i(t) over i and t
        total_num_bond = 0
        do i = 0, nframe
            do j = 1, nmol ! polymer index
                total_num_bond(num_bond(j,i)) = total_num_bond(num_bond(j,i)) + 1
                total_num_bond(nmol+1) = total_num_bond(nmol+1) + num_bond(j,i)
            end do
        end do
        n_ave = total_num_bond(nmol+1) / DBLE(nmol * (nframe+1))
        !$ omp end single


        ! -- calculate the probability density function of the number of overlapped polymers -- !
        p_num_bond = 0.0d0
        !$ omp do private(j) reduction(+: p_num_bond, summation)
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
        !print *, summation
        !$ omp end do

        !$ omp single
        WRITE(chara, "(I3.3)") ii
        outfilename = "PF/" // TRIM(ADJUSTL(chara))// ".txt"
        OPEN(outfile, file=outfilename, status="replace", form="formatted")
            do i = 0, nmol
                WRITE(outfile, "(I4, 1X, F15.7)") i, p_num_bond(i)
            end do
        CLOSE(outfile)

        command = "if [ ! -d PF/norm ] ; then mkdir -p PF/norm; fi"
        CALL system(command)
        outfilename = "PF/norm/" // TRIM(ADJUSTL(chara))// ".txt"
        OPEN(outfile, file=outfilename, status="replace", form="formatted")
            WRITE(outfile, *) "# N / N_ave, p(N) * N_ave"
            WRITE(outfile, *) "# N_ave = ", n_ave
            do i = 0, nmol
                WRITE(outfile, "(F15.7, 1X, F15.7)") DBLE(i)/DBLE(n_ave), p_num_bond(i)*n_ave
            end do 
        CLOSE(outfile)
        !$ omp end single


        ! -- Test -- !
        !$ omp do private (j,k)
        do i = 0, nframe
            do j = 1, nmol
                do k = 1, nmol
                    if (flag_bond(k,j,i) .ne. 1 .and. flag_bond(k,j,i) .ne. 0)then
                        print*,"Error: flag_bond =", flag_bond(k,j,i)
                        stop
                    end if
                end do
            enddo
        enddo
        !$ omp end do
        
        ! -- calculate the time of maintaining the state as N_o = n -- !
        !$ omp do private(j,k)
        time_of_state = -1
        do i = 0, nframe - 1  ! -- "i" represents the starting time -- !
            do j = 1, nmol    ! -- "j" represents the index of the polymer of interest -- !
                do k = 1, nframe - i     ! -- "k" represents lag time -- !
                    if (num_bond(j,i+k) .ne. num_bond(j,i)) then 
                        time_of_state(j,i) = k
                        exit  
                    end if
                end do
            end do
        end do
        !$ omp end do

        ! -- calculate probability distribution function of time of state -- !
        p_time_of_state= 0
        !$ omp do private(j) reduction(+: p_time_of_state)
        do i = 0, nframe
            do j = 1, nmol
                p_time_of_state(time_of_state(j,i), num_bond(j,i)) = p_time_of_state(time_of_state(j,i), num_bond(j,i)) + 1
                p_time_of_state(time_of_state(j,i), nmol+1) = p_time_of_state(time_of_state(j,i), nmol+1) + 1
            end do
        enddo
        !$ omp end do
        !$ omp single
        do j = 0, nmol
            p_time_of_state(:,j) = p_time_of_state(:,j) / total_num_bond(j)
        enddo
        p_time_of_state(:,nmol+1) = p_time_of_state(:,nmol+1) / ((nframe+1)*nmol)

        ! -- output -- !
        WRITE(chara, "(I3.3)") ii
        command = "if [ ! -d pdf_maintain_time/total ] ; then mkdir -p pdf_maintain_time/total; fi"
        CALL system(command)
        outfilename = "pdf_maintain_time/total/" // TRIM(ADJUSTL(chara))// ".txt"
        OPEN(outfile, file=outfilename, status="replace", form="formatted")
            WRITE(outfile, *) "# The difinision of overlap length = ", bond_length
            WRITE(outfile, *) "# The radius of gyration sqrt(Rg^2) = ", SQRT(rg2_ave)
            WRITE(outfile, *) "# The ratio of bondlength and Rg = ", bond_length/SQRT(rg2_ave)
            do i = 1, nframe
                WRITE(outfile, "(F15.4, 1X, F15.7)") i*dt*nfreq, p_time_of_state(i,nmol+1)
            end do
        CLOSE(outfile)

        do j = 1, nmol
            WRITE(chara_num_ol, "(I3.3)") j
            command = "if [ ! -d pdf_maintain_time/" // TRIM(ADJUSTL(chara_num_ol)) //  &
            " ] ; then mkdir -p pdf_maintain_time/" // TRIM(ADJUSTL(chara_num_ol)) // "; fi"
            CALL system(command)
            outfilename = "pdf_maintain_time/"// TRIM(ADJUSTL(chara_num_ol)) // "/" // TRIM(ADJUSTL(chara))// ".txt"
            OPEN(outfile, file=outfilename, status="replace", form="formatted")
                WRITE(outfile, *) "# The difinision of bond length = ", bond_length
                WRITE(outfile, *) "# The radius of gyration sqrt(Rg^2) = ", SQRT(rg2_ave)
                WRITE(outfile, *) "# The ratio of bondlength and Rg = ", bond_length/SQRT(rg2_ave)
                WRITE(outfile, "(A92, I3)") "# This is probability distribution &
                function of maintain time the polymers whose N_overlap = ",j
                do i = 1, nframe
                    WRITE(outfile, "(F15.4, 1X, F15.7)") i*dt*nfreq, p_time_of_state(i,j)
                end do
            CLOSE(outfile)
        enddo
        !$ omp end single


        ! -- calculate the maintaining time of N_o = n -- !
        ! -- "nmol + 1" represents totally averaged over N_o-- !
        ave_time_of_state = 0.0d0
        !$ omp do private(j) reduction(+: ave_time_of_state)
        do i = 0, nframe
            int_tmp_ave_time_of_state = 0
            do j = 1, nmol
                if (time_of_state(j,i) .ne. -1) then
                    int_tmp_ave_time_of_state(num_bond(j,i)) = int_tmp_ave_time_of_state(num_bond(j,i)) + time_of_state(j,i)
                end if
            end do
            do j = 0, nmol
                ave_time_of_state(j) = ave_time_of_state(j) + DBLE(int_tmp_ave_time_of_state(j) * nfreq) * dt 
                ave_time_of_state(nmol+1) = ave_time_of_state(nmol+1) + DBLE(int_tmp_ave_time_of_state(j) * nfreq) * dt 
            enddo
        enddo
        !$ omp end do
        !$ omp single
        do j = 0, nmol
            ave_time_of_state(j) = ave_time_of_state(j) / DBLE(total_num_bond(j))
        enddo
        ave_time_of_state(nmol+1) = ave_time_of_state(nmol+1) / DBLE((nframe + 1) * nmol)


        ! -- output -- !
        !WRITE(chara, "(I3.3)") ii
        !command = "if [ ! -d maintain_time/total ] ; then mkdir -p maintain_time/total; fi"
        !CALL system(command)
        !outfilename = "maintain_time/total/" // TRIM(ADJUSTL(chara))// ".txt"
        !OPEN(outfile, file=outfilename, status="replace", form="formatted")
        !    WRITE(outfile, *) "# The difinision of overlap length = ", bond_length
        !    WRITE(outfile, *) "# The radius of gyration sqrt(Rg^2) = ", SQRT(rg2_ave)
        !    WRITE(outfile, *) "# The ratio of bondlength and Rg = ", bond_length/SQRT(rg2_ave)
        !    do i = 1, npoint
        !        WRITE(outfile, "(F15.4, 1X, F15.7)") times(i), cf_bond(i,nmol+1)
        !    end do
        !CLOSE(outfile)

        outfilename = "maintain_time/"// TRIM(ADJUSTL(chara))// ".txt"
        OPEN(outfile, file=outfilename, status="replace", form="formatted")
            WRITE(outfile, *) "# The difinision of bond length = ", bond_length
            WRITE(outfile, *) "# The radius of gyration sqrt(Rg^2) = ", SQRT(rg2_ave)
            WRITE(outfile, *) "# The ratio of bondlength and Rg = ", bond_length/SQRT(rg2_ave)
            WRITE(outfile, *) "# N_o, averaged maintain time"
            WRITE(outfile, *) "# N_o = nmol+1 represents averaged over N_o"
            do i = 0, nmol+1
                if (ave_time_of_state(i) .ge. 0.0d0) then
                    WRITE(outfile, "(I4, 1X, E15.7)") i, ave_time_of_state(i)
                endif
            end do
        CLOSE(outfile)
        !$ omp end single
    end do
    !$ omp end parallel
end program main
