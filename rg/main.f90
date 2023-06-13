Program main
    use prm_list
    use read_prm
    use read_dump
    use calc_poscm
    implicit none

    double precision s(3), s2, ave, squared_ave, var
    double precision drg, normalize, dummy, tmp, tmp2
    double precision,allocatable :: p_rg(:), rg2(:)
    integer i, j, k
    integer drg_max, i_rg, i_peak
    integer id_atom
    integer, allocatable :: i_p_rg(:)

    CALL  readprm()
    CALL  readdump()

    CALL calcposcm()
    ALLOCATE(rg2(nmol))

    ave = 0.0d0
    squared_ave = 0.0d0
    drg = 0.01d0
    drg_max = 0
    rg2 = 0.0d0

    open(18,file='rg_time.xvg',status='replace')
    do k = 0, nframe
        tmp = 0.0d0
        tmp2 = 0.0d0
        ! -- calculate Rg2 for each polymer chain -- !
        do j = 1, nmol
            s2= 0.0d0
            do i = 1, natom
                id_atom = (j-1)*natom + i
                s(:) = pos(:, id_atom, k) - pos_cm(:, j, k)
                s2 = s2 + dot_product(s,s)
            enddo
            rg2(j)   = s2 / dble(natom)
            tmp      = tmp + rg2(j) !: summation of rg2 on polymer ID
            tmp2     = tmp2 + rg2(j)*rg2(j) !: summation of squared Rg2 on polymer ID
        enddo
        tmp = tmp / nmol !: averaged Rg2 of one snapshot
        ave = ave + tmp  !: summation of rg2 on snapshot ID
        write(18,'(F15.7, 1X, F15.7)')dt * k * nfreq, tmp
        tmp2 = tmp2 / nmol !: averaged squared Rg2 of one snapshot
        squared_ave = squared_ave + tmp2  !: summation of rg2 on snapshot ID

        ! -- initialization variables for distribution function -- !
        if(k .eq. 0)then
            drg_max = ceiling(2.0d0 * maxval(rg2) / drg)
            allocate(p_rg(drg_max))  !: distribution function for Rg2
            allocate(i_p_rg(drg_max))
            do i = 1, drg_max
                p_rg(i) = 0.0d0
                i_p_rg(i) = 0
            enddo
        endif

        ! -- calculate distribution function -- !
        do i = 1, nmol
            i_rg = ceiling(rg2(i) / drg)
            if (i_rg .eq. 0)then
                i_rg = 1
            endif
            !i_p_rg(i_rg) = i_p_rg(i_rg) + 1
        enddo
    enddo
    close(18)
    ave = ave / (nframe + 1)
    squared_ave = squared_ave / (nframe + 1)

    ! -- calculate variance of Rg2 -- !
    var = squared_ave - ave * ave

    OPEN(15, file='rg2.txt', status = 'replace')
    write(15,'(A20,1X,F18.9)')'#average of Rg, variance of Rg '
    WRITE(15, "(F15.7, 1x, F15.7)")ave, var
    CLOSE(15)
    
    do i = 1, drg_max
        p_rg(i) = dble(i_p_rg(i)) / dble((nmol * drg * (nframe+1)))
    enddo
    open(19,file='distribution_of_rg2.xvg',status='replace')
    do i = 1, drg_max
        normalize = normalize + p_rg(i)*drg
    enddo
    write(19,'(A20,1X,F17.15)')'#integral p_rg ==',normalize
    do i = 1, drg_max
        write(19,'(F30.4,2X,F30.15,I20)')(drg * i), p_rg(i)
    enddo
    close(19)
    end program main
