Program distribution_of_r_of_gyration
    use prm_list
    use read_prm
    use read_dump

    implicit none

    double precision r(3), s(3), rcm(3), s2, rg2, rg_ave
    double precision drg, normalize, dummy
    double precision,allocatable :: rg_save(:), p_rg(:)
    integer i, j, k
    integer drg_max, i_rg, i_peak
    integer, allocatable :: i_p_rg(:)

    CALL  readprm()
    CALL  readdump()

    open(18,file='rg_time.xvg',status='replace')

    allocate(rg_save(nmol))

    rg_ave = 0.0d0
    drg = 0.01d0
    drg_max = 0

    do k = 0, nframe
        rg2 = 0.0d0
        do i = 1, nmol
            rg_save(i) = 0.0d0 
        enddo

        do i = 1, 9
            read(17,'()')
        enddo

        do j = 1, nmol
            s2= 0.0d0
            rcm = 0.0d0
            do i = 1, natom
                read(17,*)dummy,r(:)
                rcm(:) = rcm(:) + r(:)
            enddo
            rcm(:) = rcm(:) / natom
            do i = 1, natom
                backspace(17)
            enddo
            do i = 1, natom
                read(17,*)dummy,r(:)
                s(:) = r(:) - rcm(:)
                s2 = s2 + dot_product(s,s)
            enddo
            s2 = s2 /dble(natom)
            rg2 = rg2 + s2
            rg_save(j) = s2**0.50d0
        enddo
        rg2 = rg2 / dble(nmol)
        write(18,'(F25.3,2X,F30.8,2X,F30.8)')dt * k * nfreq, rg2
        rg_ave = rg_ave + rg2
        if(k .eq. 0)then
            drg_max = ceiling(maxval(rg_save) / drg) + int(10000 / drg) !yoyuu motasu
            allocate(p_rg(drg_max))
            allocate(i_p_rg(drg_max))
            do i = 1, drg_max
                p_rg(i) = 0.0d0
                i_p_rg(i) = 0
            enddo
        endif
        do i = 1, nmol
            i_rg = ceiling(rg_save(i) / drg)
            if (i_rg .eq. 0)then
                i_rg = 1
            endif
            i_p_rg(i_rg) = i_p_rg(i_rg) + 1
        enddo
        print*,k,'th frame'
    enddo

    rg_ave = rg_ave / (nframe + 1)

    OPEN(15, file='rg.txt', status = 'replace')
    write(15,'(A20,1X,F18.9)')'#average of rg =',rg_ave
    CLOSE(15)
    
    do i = 1, drg_max
        p_rg(i) = dble(i_p_rg(i)) / dble((nmol * drg * (nframe+1)))
    enddo
    open(19,file='distri_rg.xvg',status='replace')
    do i = 1, drg_max
        normalize = normalize + p_rg(i)*drg
    enddo
    write(19,'(A20,1X,F17.15)')'#integral p_rg ==',normalize
    open(20,file='distri_rg_scaled.xvg')
    do i = 1, drg_max
        write(19,'(F30.4,2X,F30.15,I20)')(drg * i), p_rg(i),i_p_rg(i)
        write(20,'(F30.4,2X,F30.15,I20)')(drg * i / rg_ave), p_rg(i) *rg_ave
    enddo
    print*,rg_ave
    deallocate(rg_save,p_rg,i_p_rg)
    close(18)
    close(19)
    close(20)
    end program distribution_of_r_of_gyration
