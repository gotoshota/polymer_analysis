! -- This programs calculates Bond-brakage correlations.                  -- !
! -- Reffer the following paper,                                          -- !
! -- H. Shiba, T. Kawasaki, and A. Onuki, Phys. Rev. E 86, 041504 (2012). -- !

program main
    !$ use omp_lib
    use prm_list
    use read_prm
    use coarse_grain
    use read_dump
    use calc_poscm

    implicit none

    INTEGER(KIND=4)                 :: i, j, k, l
    INTEGER(KIND=4)                 :: ii, iii
    INTEGER(KIND=4)                 :: id_particle
    INTEGER(KIND=4)                 :: max_n_o
    INTEGER(KIND=4)                 :: delta_t
    INTEGER(KIND=4)                 :: delta_frame

    INTEGER(KIND=4), ALLOCATABLE    :: i_normalize(:)
    INTEGER(KIND=4), ALLOCATABLE    :: total_num_bond(:)
    INTEGER(KIND=4), ALLOCATABLE    :: num_bond(:,:)
    INTEGER(KIND=4), ALLOCATABLE    :: int_tmp_ave_time_of_state(:)
    INTEGER(KIND=4), ALLOCATABLE    :: time_of_state(:,:)
    INTEGER(KIND=1), ALLOCATABLE    :: flag_bond(:,:,:)
    INTEGER(KIND=1), ALLOCATABLE    :: flag_breakage(:,:,:)
    INTEGER(KIND=4), ALLOCATABLE    :: broken_bond(:,:,:)
    INTEGER(KIND=4), ALLOCATABLE    :: num_broken_bond(:,:)

    DOUBLE PRECISION                :: normalize
    DOUBLE PRECISION                :: rg2_ave
    DOUBLE PRECISION                :: tmp
    DOUBLE PRECISION                :: bond_length
    DOUBLE PRECISION                :: breakage_length
    DOUBLE PRECISION                :: bond_prefactor
    DOUBLE PRECISION                :: breakage_prefactor
    DOUBLE PRECISION                :: dpos(3)
    DOUBLE PRECISION                :: distance
    DOUBLE PRECISION                :: summation
    DOUBLE PRECISION                :: n_ave
    DOUBLE PRECISION                :: density
    DOUBLE PRECISION                :: t

    DOUBLE PRECISION, ALLOCATABLE   :: rg2(:,:)
    DOUBLE PRECISION, ALLOCATABLE   :: p_num_bond(:)
    DOUBLE PRECISION, ALLOCATABLE   :: p_time_of_state(:,:)
    DOUBLE PRECISION, ALLOCATABLE   :: ave_time_of_state(:)
    DOUBLE PRECISION, ALLOCATABLE   :: frac_of_num_of_broken_bond(:,:)
    DOUBLE PRECISION, ALLOCATABLE   :: cf_broken_bond(:)
    DOUBLE PRECISION, ALLOCATABLE   :: susceptibility(:)

    INTEGER(KIND=4), PARAMETER      :: outfile = 17
    CHARACTER(LEN=256)              :: outfilename
    CHARACTER(LEN=256)              :: dir_name, dir_name_header
    CHARACTER(LEN=256)              :: command
    CHARACTER(LEN=256)              :: chara1, chara2, chara3
    INTEGER(KIND=4)                 :: access

    CHARACTER(LEN=256)              :: arg

! ------------------------------------------------------------------ !
! ------------------------------------------------------------------ !

    CALL readprm()
    CALL coarsegrain()
    CALL readdump()
    CALL calcposcm()

    ! -- number of density of polymer chain -- !
    density = DBLE(nmol) / (box_l**3.0d0)

    ! -- Define Rg^2 -- !
    if (iargc() > 0) then
        print *, "Get the value of Rg^2 from the file whose path is of cmmand line argument"
        CALL GET_COMMAND_ARGUMENT(NUMBER=1, VALUE=arg)
        OPEN(10, file=arg, status="old")
            READ(10, *) rg2_ave
        CLOSE(10)

    else if ( access( "rg2.txt", " " ) .eq. 0 ) then
        print *, "Get the value of Rg^2 from output file rg2.txt"
        OPEN(10, file="rg2.txt", status="old")
            READ(10, *) rg2_ave
        CLOSE(10)

    else if ( access( "rg2_long.txt", " " ) .eq. 0 ) then
        print *, "Get the value of Rg^2 from output file rg2_long.txt"
        OPEN(10, file="rg2_long.txt", status="old")
            READ(10, *) rg2_ave
        CLOSE(10)

    else
        print *, "Calculate the value of Rg^2"
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
        DEALLOCATE(rg2)
        !$ omp end single
    end if

    ! -- calculate the number of bonded molecules-- !
    ALLOCATE(num_bond(nmol, 0:nframe))
    ALLOCATE(total_num_bond(0:nmol+1))
    ALLOCATE(flag_bond(nmol, nmol, 0:nframe))
    ALLOCATE(flag_breakage(nmol, nmol, 0:nframe))
    ALLOCATE(p_num_bond(0:nmol))
    ALLOCATE(i_normalize(0:nmol+1))
    ALLOCATE(time_of_state(nmol, 0:nframe))
    ALLOCATE(p_time_of_state(nframe, 0:nmol+1))
    ALLOCATE(ave_time_of_state(0:nmol+1))
    ALLOCATE(int_tmp_ave_time_of_state(0:nmol+1))
    ALLOCATE(broken_bond(nmol, 0:nframe-TargetFrame(1), npoint))
    ALLOCATE(num_broken_bond(0:nmol, npoint))
    ALLOCATE(frac_of_num_of_broken_bond(0:nmol, npoint))
    ALLOCATE(cf_broken_bond(npoint))
    ALLOCATE(susceptibility(npoint))

    command = "touch BOND_LENGTH"
    CALL system(command)
    !$ omp do default(private) 
    do iii = 1, 20
        WRITE(chara1,"(I3.3)") iii
        dir_name = chara1
        CALL MKDIR(dir_name)

        ! -- Define bond length -- !
        bond_prefactor = DBLE(iii) * 0.10d0
        bond_length = SQRT(rg2_ave) * bond_prefactor

        ! -- Calculate bond flag --  !
        CALL CALC_BOND_FLAG(nframe, box_l, nmol, pos_cm, bond_length, flag_bond)
        print *,"Calculated flag of bond." 

        ! -- Count up the number of bond of j-th particle -- !
        num_bond = 0
        do i = 0, nframe
            do j = 1, nmol
                do k = 1, nmol
                    num_bond(j,i) = num_bond(j,i) + flag_bond(k,j,i)
                enddo
            enddo
        enddo

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

        ! -- calculate the probability density function of the number of overlapped polymers -- !
        !dir_name = TRIM(ADJUSTL(chara1)) // "/PF"
        !CALL MKDIR(dir_name)
        !dir_name = TRIM(ADJUSTL(chara1)) // "/PF/norm"
        !CALL MKDIR(dir_name)
        CALL CALC_PDF_NUM_OF_BOND(nframe, nmol, ii, n_ave, num_bond,  p_num_bond, chara1)

        command = "touch " // TRIM(ADJUSTL(chara1)) // "/BREAKAGE_LENGTH"
        CALL system(command)
        do ii = iii, 20
            WRITE(chara2,"(I3.3)") ii
            dir_name = TRIM(ADJUSTL(chara1)) // "/"  //  TRIM(ADJUSTL(chara2))
            CALL MKDIR(dir_name)
            dir_name_header = dir_name

            ! -- Define bond-breakage length -- !
            breakage_prefactor = DBLE(ii) * 0.10d0
            breakage_length = SQRT(rg2_ave) * breakage_prefactor
          
            ! -- Calculate bond flag --  !
            CALL CALC_BOND_FLAG(nframe, box_l, nmol, pos_cm, breakage_length, flag_breakage)
            print *,"Calculated flag of breakage." 

            ! -- calculate bond-breakage number B_i and the number of the particle whose B_i = k -- !
            CALL CALC_BROKEN_BOND(nmol, nframe, npoint, TargetFrame, flag_bond, flag_breakage, broken_bond)
            print *,"Calculated broken bond." 

            CALL CALC_NUM_OF_BROKEN_BOND(nmol, nframe, npoint, TargetFrame, broken_bond, num_broken_bond)
            print *,"Calculated the number of broken bond." 

            ! -- fraction of the number of particle whose bond are broken in time window t -- !
            do i = 1, npoint ! -- lag time (second arg) -- !
                do j = 0, nmol ! -- the number of broken bond (first arg) -- !
                    frac_of_num_of_broken_bond(j,i) = num_broken_bond(j,i) / DBLE(nmol * (nframe + 1 - TargetFrame(i))) 
                enddo
            enddo
            do i = 1, npoint
                summation = 0.0d0
                do j = 0, nmol
                    summation = summation + frac_of_num_of_broken_bond(j,i)
                enddo
            enddo
            print *, "Calculated the fraction of broken bond in time lag t"
            ! -- output -- !
            dir_name = dir_name_header // "/frac_broken_bond"
            CALL MKDIR(dir_name)
            do j = 0, nmol
                WRITE(chara3, "(I3.3)") j
                outfilename = TRIM(ADJUSTL(dir_name)) // "/" // TRIM(ADJUSTL(chara3)) // ".txt"
                OPEN (outfile, file=outfilename, status="replace", form="formatted")
                WRITE(outfile, *) "# time , fraction of the particle with broken bond = ", j
                do i = 1, npoint
                    WRITE(outfile, 100) times(i), frac_of_num_of_broken_bond(j,i)
                enddo
                CLOSE(outfile)
            enddo
            print *, "Output the fraction of broken bond was completed."
            
            ! -- time correlation function of bond breakage -- !
            cf_broken_bond = 0.0d0
            do i = 1, npoint
                do j = 1, nmol
                    cf_broken_bond(i) = cf_broken_bond(i) + j * frac_of_num_of_broken_bond(j,i)
                enddo
                cf_broken_bond(i) = cf_broken_bond(i) * 0.50d0 * density
            enddo
            print *, "Calculated the total fraction of broken bond in time lag t"

            ! -- output -- !
            outfilename = dir_name_header // "/cf_broken_bond.txt"
            OPEN (outfile, file=outfilename, status="replace", form="formatted")
            WRITE(outfile, *) "# time , TCF "
            do i = 1, npoint
                WRITE(outfile, 100) times(i), cf_broken_bond(i)
            enddo
            CLOSE(outfile)
            print *, "Output the total fraction of broken bond was completed."

            CALL CALC_SUSCEPTIBILITY(box_l, nmol, nframe, npoint, TargetFrame, broken_bond, cf_broken_bond, susceptibility)
            print *, "Calculated the susceptibility"

            ! -- output -- !
            outfilename = dir_name_header // "/susceptibility.txt"
            OPEN (outfile, file=outfilename, status="replace", form="formatted")
            WRITE(outfile, *) "# time , susceptibility "
            do i = 1, npoint
                WRITE(outfile, 100) times(i), susceptibility(i)
            enddo
            CLOSE(outfile)
            print *, "Output the susceptibility was completed."

            t = 10.0d0**5
            CALL VIS_BOND_BREAKAGE(box_l, nmol, t, dt, nfreq, nframe, flag_bond, pos_cm)
        enddo
    enddo    
    !$ omp end do
    !$ omp end parallel
    ! -- format -- !
    100 format(E15.7, 1X, E15.7)

! ------------------------------------------------------------------------------------------------------- !
! ------------------------------------------------------------------------------------------------------- !
! ------------------------------------------------------------------------------------------------------- !

    contains
        ! -- Calculate the flag whether pair of particles is bonded or not -- !
        subroutine CALC_BOND_FLAG(nframe, box_l, nmol, pos, threshold, flag)
            implicit none
            INTEGER(KIND=4), INTENT(IN)                 :: nmol, nframe
            DOUBLE PRECISION, INTENT(IN)                :: box_l
            DOUBLE PRECISION, INTENT(IN)                :: pos(1:, 1:, 0:)
            DOUBLE PRECISION, INTENT(IN)                :: threshold 
            
            INTEGER(KIND=1), INTENT(OUT)                :: flag(1:, 1:, 0:)

            DOUBLE PRECISION                            :: dpos(3)
            DOUBLE PRECISION                            :: distance


            flag = 0
            do i = 0, nframe
                do j = 1, nmol
                    do k = 1, nmol
                        if (j .ne. k) then
                            dpos(:) = pos(:,k,i) - pos(:,j,i)
                            dpos(:) = ABS(dpos(:) - box_l * NINT(dpos(:) / box_l))
                            distance = SQRT(dpos(1)*dpos(1) + dpos(2)*dpos(2) + dpos(3)*dpos(3))
                        else 
                            flag(k,j,i) = 0
                        end if 
                        ! -- if bonded flag_bond = 1, if not flag_bond = 0 -- !
                        if (distance .le. threshold) then
                            flag(k,j,i) = 1
                        else 
                            flag(k,j,i) = 0
                        end if 
                    end do
                end do 
            end do

        end subroutine

        subroutine CALC_PDF_NUM_OF_BOND(nframe, nmol, int_bond_prefactor, n_ave, num_bond, pdf, chara1)
            implicit none    

            INTEGER(KIND=4), INTENT(IN)                 :: nmol, nframe
            INTEGER(KIND=4), INTENT(IN)                 :: int_bond_prefactor
            DOUBLE PRECISION, INTENT(IN)                :: n_ave
            INTEGER, INTENT(IN)                         :: num_bond(1:, 0:)
            CHARACTER(LEN=256), INTENT(IN)              :: chara1

            DOUBLE PRECISION, INTENT(OUT)               :: pdf(0:)

            INTEGER(KIND=4)                             :: i, j
            INTEGER(KIND=4)                             :: outfile = 16
            CHARACTER(LEN=256)                          :: chara, outfilename, command

            p_num_bond = 0.0d0
            do i = 0, nframe
                do j = 1, nmol
                    pdf(num_bond(j,i)) = pdf(num_bond(j,i)) + 1
                end do
            end do
            do i = 0, nmol
                pdf(i) = pdf(i) / (nmol * (nframe + 1))
            end do

            outfilename = TRIM(ADJUSTL(chara1)) // "/probability_nb.txt"
            OPEN(outfile, file=outfilename, status="replace", form="formatted")
            WRITE(outfile, *) "# N_b, pdf(N_b)"
                do i = 0, nmol
                    WRITE(outfile, "(I4, 1X, F15.7)") i, pdf(i)
                end do
            CLOSE(outfile)

            outfilename = "PF/norm/" // TRIM(ADJUSTL(chara))// ".txt"
            outfilename = TRIM(ADJUSTL(chara1)) // "/probability_nb_norm.txt"
            OPEN(outfile, file=outfilename, status="replace", form="formatted")
            WRITE(outfile, *) "# N / N_ave, p(N) * N_ave"
            WRITE(outfile, *) "# N_ave = ", n_ave
                do i = 0, nmol
                    WRITE(outfile, "(F15.7, 1X, F15.7)") DBLE(i)/DBLE(n_ave), pdf(i)*n_ave
                end do
            CLOSE(outfile)

        end subroutine

        ! -- Calculate the function B_i which represents how many bond are broken whose i-th polymer has. -- !
        subroutine CALC_BROKEN_BOND(nmol, nframe, npoint, TargetFrame, flag1, flag2, broken_bond)
            implicit none

            INTEGER(KIND=4), INTENT(IN)                 :: nmol, nframe, npoint, TargetFrame(0:)
            INTEGER(KIND=1), INTENT(IN)                 :: flag1(1:, 1:, 0:)
            INTEGER(KIND=1), INTENT(IN)                 :: flag2(1:, 1:, 0:)
            INTEGER(KIND=4), INTENT(OUT)                :: broken_bond(1:,0:,1:) ! -- args = (i-th polymer, t0, t) -- !
                                                                                 ! -- t0 represents initial time   -- !
                                                                                 ! -- t represents lag time        -- !
            INTEGER(KIND=4)                             :: i, j, k, l

            broken_bond = 0
            do i = 1, npoint ! -- lag time (third arg) -- !
                do j = 0, nframe-TargetFrame(i) ! -- initial time (second arg) -- !
                    do k = 1, nmol ! -- polymer index (first arg) -- !
                        do l = 1, nmol
                            broken_bond(k, j, i) = broken_bond(k, j, i) + flag1(l,k,j)*(1-flag2(l,k,j+TargetFrame(i)))
                        end do
                    end do
                end do
            end do

        end subroutine

        ! -- calculate the number of the particle with broken_bond = k -- !
        subroutine CALC_NUM_OF_BROKEN_BOND(nmol, nframe, npoint, TargetFrame, broken_bond, num_broken_bond)
            implicit none

            INTEGER(KIND=4), INTENT(IN)                 :: nmol, nframe, npoint, TargetFrame(0:), broken_bond(1:,0:,1:)
            INTEGER(KIND=4), INTENT(OUT)                :: num_broken_bond(0:,1:) ! -- 1st arg is the number of broken bond -- !
                                                                                  ! -- 2nd arg is lag time                  -- !

            INTEGER(KIND=4)                             :: i, j, k, l

            
            num_broken_bond = 0
            do i = 1, npoint
                do j = 0, nframe - TargetFrame(i)
                    do k = 1, nmol
                        num_broken_bond(broken_bond(k,j,i) ,i) = num_broken_bond(broken_bond(k,j,i), i) + 1
                    enddo
                enddo
            enddo
        end subroutine

        ! -- calculate the susceptibility of overall degree of the bond-breakage -- !
        subroutine CALC_SUSCEPTIBILITY(box_l, nmol, nframe, npoint, TargetFrame, broken_bond, cf_broken_bond, susceptibility)
            implicit none   
            INTEGER(KIND=4), INTENT(IN)                 :: nmol, nframe, npoint, TargetFrame(0:), broken_bond(1:,0:,1:)
            DOUBLE PRECISION, INTENT(IN)                :: cf_broken_bond(1:), box_l
            DOUBLE PRECISION, INTENT(OUT)               :: susceptibility(1:)
            
            INTEGER(KIND=4)                             :: i, j, k, l
            DOUBLE PRECISION                            :: tmp, density, vol

            vol = box_l ** 3.0d0
            density = nmol / vol
            susceptibility = 0.0d0
            do i = 1, npoint
                tmp = 2.0d0 * cf_broken_bond(i) / density
                do j = 0, nframe - TargetFrame(i)
                    do k = 1, nmol
                        do l = 1, nmol
                            susceptibility(i) = susceptibility(i) + (broken_bond(k,j,i) - tmp) * (broken_bond(l,j,i) - tmp)
                        enddo
                    enddo
                enddo
                susceptibility(i) = susceptibility(i) * 0.250d0 / vol
            enddo

        end subroutine

        ! -- Visualize the bond breakage -- !
        ! -- Note: pos_cm is unwrapped coordinates -- !
        subroutine VIS_BOND_BREAKAGE(box_l, nmol, t, dt, nfreq, nframe, flag_bond, pos_cm)
            implicit none   
            INTEGER(KIND=4), INTENT(IN)                 :: nmol, nfreq,nframe
            INTEGER(KIND=1), INTENT(IN)                 :: flag_bond(1:, 1:, 0:)
            DOUBLE PRECISION, INTENT(IN)                :: dt, t, pos_cm(1:, 1:, 0:), box_l

            INTEGER(KIND=4)                             :: frame0, frame1, delta_frame
            INTEGER(KIND=4)                             :: i, j, k, ii
            INTEGER(KIND=4)                             :: num_breakage
            INTEGER(KIND=4)                             :: num_frame
            INTEGER(KIND=4)                             :: outfile = 17
            DOUBLE PRECISION                            :: break_point(3)
            DOUBLE PRECISION                            :: dpos(3)
            DOUBLE PRECISION                            :: pos_tmp(3)
            DOUBLE PRECISION,ALLOCATABLE                :: arr_break_point(:,:)
            CHARACTER(LEN=256)                          :: dumpfilename = "bond_breakage.lammpstrj"
            CHARACTER(LEN=10)                           :: ii_chara
                
            delta_frame = INT(t / dt / nfreq)
            command = "if [ ! -d bond_breakage_traj ] ; then mkdir -p bond_breakage_traj ; fi"
            CALL system(command)
            do ii = 0, nframe-delta_frame
                frame0 = INT(ii / dt / nfreq)
                frame1 = frame0 + delta_frame
                ! -- Count up the number of breakage point -- !
                num_breakage = 0
                do i = 1, nmol
                    do j = 1, nmol
                        if (flag_bond(j,i,frame0) .eq. 1 .and. flag_bond(j,i,frame1) .eq. 0) then
                            num_breakage = num_breakage + 1
                        endif
                    enddo
                enddo

                ALLOCATE(arr_break_point(3,num_breakage))
                k = 0
                do i = 1, nmol
                    do j = 1, nmol
                        if (flag_bond(j,i,frame0) .eq. 1 .and. flag_bond(j,i,frame1) .eq. 0) then
                            k = k + 1
                            dpos(:) = pos_cm(:,j,frame1) - pos_cm(:,i,frame1)
                            pos_tmp(:) = pos_cm(:,j,frame1) - box_l * INT(dpos(:) / box_l)
                            arr_break_point(1,k) = 0.50d0 * (pos_tmp(1) + pos_cm(1,i,frame1))
                            arr_break_point(2,k) = 0.50d0 * (pos_tmp(2) + pos_cm(2,i,frame1))
                            arr_break_point(3,k) = 0.50d0 * (pos_tmp(3) + pos_cm(3,i,frame1))
                        endif
                    enddo
                enddo

                ! -- Write out break point -- !
                WRITE(ii_chara, "(I5.5)") ii
                dumpfilename = "bond_breakage_traj/"// TRIM(ADJUSTL(ii_chara)) // ".lammpstrj"
                CALL mktraj_lmp_format(dumpfilename, box_l, arr_break_point, num_breakage, frame1)
                DEALLOCATE(arr_break_point)
            enddo
        end subroutine

        subroutine mktraj_lmp_format(dumpfilename, box_l, pos, nmol, frame)
            
            implicit none

            CHARACTER(LEN=256), INTENT(IN)  :: dumpfilename
            INTEGER(KIND=4),    INTENT(IN)  :: nmol, frame
            DOUBLE PRECISION,   INTENT(IN)  :: box_l
            DOUBLE PRECISION,   INTENT(IN)  :: pos(1:,1:)
            
            INTEGER(KIND=4)                 :: i,j
            INTEGER(KIND=4), PARAMETER      :: dump = 16

            CHARACTER(LEN=50)              :: header(4)

            ! -- define headers -- !
            header(1) = "ITEM: TIMESTEP"
            header(2) = "ITEM: NUMBER OF ATOMS"
            header(3) = "ITEM: BOX BOUNDS pp pp pp"
            header(4) = "ITEM: ATOMS id xu yu zu"

            ! -- write out -- !
            OPEN(dump, file=dumpfilename, status="replace")
                WRITE(dump, 101) header(1)
                WRITE(dump, *) frame
                WRITE(dump, 102) header(2)
                WRITE(dump, "(I10)") nmol
                WRITE(dump, 103) header(3)
                WRITE(dump, "(2(E15.7,1X))") 0.0d0, box_l
                WRITE(dump, "(2(E15.7,1X))") 0.0d0, box_l
                WRITE(dump, "(2(E15.7,1X))") 0.0d0, box_l
                WRITE(dump, 104) header(4)
                do j = 1, nmol
                    WRITE(dump, "(I4, 1X, F15.7, 1X, F15.7, 1X, F15.7)") j, pos(1, j), pos(2, j), pos(3, j)
                enddo
            CLOSE(dump)
            
            101 format(A14)
            102 format(A21)
            103 format(A25)
            104 format(A23)

        end subroutine

        subroutine MKDIR(path)
            implicit none

            CHARACTER(LEN=256), INTENT(IN)      :: path

            CHARACTER(LEN=256)                  :: command

            command = "if [ ! -d " // TRIM(ADJUSTL(path)) // " ] ; then mkdir -p " // TRIM(ADJUSTL(path)) // "; fi"
            CALL system(command)
        
        end subroutine

end program main

