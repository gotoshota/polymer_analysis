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
    INTEGER(KIND=4)                 :: ii
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
    INTEGER(KIND=4), ALLOCATABLE    :: broken_bond(:,:,:)
    INTEGER(KIND=4), ALLOCATABLE    :: num_broken_bond(:,:)

    DOUBLE PRECISION                :: normalize
    DOUBLE PRECISION                :: rg2_ave
    DOUBLE PRECISION                :: tmp
    DOUBLE PRECISION                :: bond_length
    DOUBLE PRECISION                :: dpos(3)
    DOUBLE PRECISION                :: distance
    DOUBLE PRECISION                :: prefactor
    DOUBLE PRECISION                :: summation
    DOUBLE PRECISION                :: n_ave
    DOUBLE PRECISION                :: density

    DOUBLE PRECISION, ALLOCATABLE   :: rg2(:,:)
    DOUBLE PRECISION, ALLOCATABLE   :: p_num_bond(:)
    DOUBLE PRECISION, ALLOCATABLE   :: p_time_of_state(:,:)
    DOUBLE PRECISION, ALLOCATABLE   :: ave_time_of_state(:)
    DOUBLE PRECISION, ALLOCATABLE   :: frac_of_num_of_broken_bond(:,:)
    DOUBLE PRECISION, ALLOCATABLE   :: cf_broken_bond(:)
    DOUBLE PRECISION, ALLOCATABLE   :: susceptibility(:)

    INTEGER(KIND=4), PARAMETER      :: outfile = 17
    CHARACTER(LEN=256)              :: outfilename
    CHARACTER(LEN=256)              :: command
    CHARACTER(LEN=256)              :: chara
    CHARACTER(LEN=256)              :: chara_num_ol

    INTEGER(KIND=4)                 :: access

! ------------------------------------------------------------------ !
! ------------------------------------------------------------------ !

    CALL readprm()
    CALL coarsegrain()
    CALL readdump()
    CALL calcposcm()

    ! -- number of density of polymer chain -- !
    density = DBLE(nmol) / (box_l**3.0d0)

    if ( access( "rg2.txt", " " ) .eq. 0 ) then
        OPEN(10, file="rg2.txt", status="old")
            READ(10, *) rg2_ave
        CLOSE(10)
    else
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

    prefactor = 1.0d0
    bond_length = SQRT(rg2_ave) * prefactor
    num_bond = 0
    flag_bond = 0

    ! -- Determin overlapped (bonded) or not --  !
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
                ! -- if bonded flag_bond = 1, if not flag_bond = 0 -- !
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

    ! -- calculate the number of polymer whose N_o = n -- !
    ! -- nmol + 1 represents the summation of N_o,i(t) over i and t
    !$ omp single
    total_num_bond = 0
    do i = 0, nframe
        do j = 1, nmol ! polymer index
            total_num_bond(num_bond(j,i)) = total_num_bond(num_bond(j,i)) + 1
            total_num_bond(nmol+1) = total_num_bond(nmol+1) + num_bond(j,i)
        end do
    end do
    n_ave = total_num_bond(nmol+1) / DBLE(nmol * (nframe+1))
    !$ omp end single
    !$ omp end parallel

    ! -- calculate bond-breakage number B_i and the number of the particle whose B_i = k -- !
    CALL CALC_BROKEN_BOND(nmol, nframe, npoint, TargetFrame, flag_bond, broken_bond)
    CALL CALC_NUM_OF_BROKEN_BOND(nmol, nframe, npoint, TargetFrame, broken_bond, num_broken_bond)

    ! -- fraction of the number of particle whose bond are broken in time window t -- !
    do i = 1, npoint ! -- lag time (second arg) -- !
        do j = 0, nmol ! -- the number of broken bond (first arg) -- !
            frac_of_num_of_broken_bond(j,i) = num_broken_bond(j,i) / DBLE(nmol * (nframe + 1 -TargetFrame(i))) 
        enddo
    enddo
    do i = 1, npoint
        summation = 0.0d0
        do j = 0, nmol
            summation = summation + frac_of_num_of_broken_bond(j,i)
        enddo
    enddo
    
    ! -- time correlation function of bond breakage -- !
    cf_broken_bond = 0.0d0
    do i = 1, npoint
        do j = 1, nmol
            cf_broken_bond(i) = cf_broken_bond(i) + j * frac_of_num_of_broken_bond(j,i)
        enddo
        cf_broken_bond(i) = cf_broken_bond(i) * 0.50d0 * density
    enddo

    CALL CALC_SUSCEPTIBILITY(box_l, nmol, nframe, npoint, TargetFrame, broken_bond, cf_broken_bond, susceptibility)
    ! -- output -- !
    outfilename = "susceptibility.txt"
    OPEN (outfile, file=outfilename, status="replace", form="formatted")
    WRITE(outfile, *) "# time , susceptibility "
    do i = 1, npoint
        WRITE(outfile, 100) times(i), susceptibility(i)
    enddo
    
    ! -- format -- !
    100 format(E15.7, 1X, E15.7)
! ------------------------------------------------------------------------------------------------------- !
! ------------------------------------------------------------------------------------------------------- !
! ------------------------------------------------------------------------------------------------------- !

    contains

    ! -- Calculate the function B_i which represents how many bond are broken whose i-th polymer has. -- !
    subroutine CALC_BROKEN_BOND(nmol, nframe, npoint, TargetFrame, flag_bond, broken_bond)
        !$ use omp_lib
        implicit none

        INTEGER(KIND=4), INTENT(IN)                 :: nmol, nframe, npoint, TargetFrame(0:)
        INTEGER(KIND=1), INTENT(IN)                 :: flag_bond(1:, 1:, 0:)
        INTEGER(KIND=4), INTENT(OUT)                :: broken_bond(1:,0:,1:) ! -- args = (i-th polymer, t0, t) -- !
                                                                             ! -- t0 represents initial time   -- !
                                                                             ! -- t represents lag time        -- !

        INTEGER(KIND=4)                             :: i, j, k, l

        broken_bond = 0
        !$ omp parallel do private(j,k,l)
        do i = 1, npoint ! -- lag time (third arg) -- !
            do j = 0, nframe-TargetFrame(i) ! -- initial time (second arg) -- !
                do k = 1, nmol ! -- polymer index (first arg) -- !
                    do l = 1, nmol
                        broken_bond(k, j, i) = broken_bond(k, j, i) + flag_bond(l,k,j)*(1-flag_bond(l,k,j+TargetFrame(i)))
                    end do
                end do
            end do
        end do
        !$ omp end do
        !$ omp end parallel
    end subroutine

    ! -- calculate the number of the particle with broken_bond = k -- !
    subroutine CALC_NUM_OF_BROKEN_BOND(nmol, nframe, npoint, TargetFrame, broken_bond, num_broken_bond)
        !$ use omp_lib
        implicit none

        INTEGER(KIND=4), INTENT(IN)                 :: nmol, nframe, npoint, TargetFrame(0:), broken_bond(1:,0:,1:)
        INTEGER(KIND=4), INTENT(OUT)                :: num_broken_bond(0:,1:) ! -- 1st arg is the number of broken bond -- !
                                                                              ! -- 2nd arg is lag time                  -- !

        INTEGER(KIND=4)                             :: i, j, k, l

        
        num_broken_bond = 0
        !$ omp parallel do private(j,k)
        do i = 1, npoint
            do j = 0, nframe - TargetFrame(i)
                do k = 1, nmol
                    num_broken_bond(broken_bond(k,j,i) ,i) = num_broken_bond(broken_bond(k,j,i), i) + 1
                enddo
            enddo
        enddo
        !$ omp end do parallel
    end subroutine

    ! -- calculate the susceptibility of overall degree of the bond-breakage -- !
    subroutine CALC_SUSCEPTIBILITY(box_l, nmol, nframe, npoint, TargetFrame, broken_bond, cf_broken_bond, susceptibility)
        !$ use omp_lib
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
end program main

