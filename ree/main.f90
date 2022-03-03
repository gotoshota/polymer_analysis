Program main
    use prm_list
    use read_prm
    use read_dump
    implicit none

    double precision ree(3), tmp, ree2, ree_ave
    double precision dree, normalize, dummy
    double precision,allocatable :: ree_save(:), p_ree(:)
    integer i, j, k
    integer dree_max, i_ree, i_peak
    integer id_atom
    integer, allocatable :: i_p_ree(:)

    CALL  readprm()
    CALL  readdump()

    open(18,file='ree_time.xvg',status='replace')
    allocate(ree_save(nmol))

    ree_ave = 0.0d0
    dree = 0.01d0
    dree_max = 0

    do k = 0, nframe
        ree2 = 0.0d0
        do j = 1, nmol
            id_atom = (j - 1) * natom
            ree(:) = pos(:, id_atom + 1, k) - pos(:, id_atom + natom, k)
            tmp = dot_product(ree,ree)
            ree2 = ree2 + tmp
            !ree_save(j) = SQRT(tmp)
        enddo
        ree2 = ree2 / dble(nmol)
        write(18,'(F25.3,2X,F30.8,2X,F30.8)')dt * k * nfreq, ree2
        ree_ave = ree_ave + ree2
        !if(k .eq. 0)then
        !    dree_max = ceiling(maxval(rg_save) / drg) + int(10000 / drg) !yoyuu motasu
        !    allocate(p_rg(drg_max))
        !    allocate(i_p_rg(drg_max))
        !    do i = 1, drg_max
        !        p_rg(i) = 0.0d0
        !        i_p_rg(i) = 0
        !    enddo
        !endif
        !do i = 1, nmol
        !    i_rg = ceiling(rg_save(i) / drg)
        !    if (i_rg .eq. 0)then
        !        i_rg = 1
        !    endif
        !    i_p_rg(i_rg) = i_p_rg(i_rg) + 1
        !enddo
        !print*,k,'th frame'
    enddo
    close(18)

    ree_ave = ree_ave / dble(nframe + 1)

    OPEN(15, file='ree.txt', status = 'replace')
    write(15,'(A20,1X,F18.9)')'#average of rg =',ree_ave
    CLOSE(15)
    
    !do i = 1, drg_max
    !    p_rg(i) = dble(i_p_rg(i)) / dble((nmol * drg * (nframe+1)))
    !enddo
    !open(19,file='distri_rg.xvg',status='replace')
    !do i = 1, drg_max
    !    normalize = normalize + p_rg(i)*drg
    !enddo
    !write(19,'(A20,1X,F17.15)')'#integral p_rg ==',normalize
    !open(20,file='distri_rg_scaled.xvg')
    !do i = 1, drg_max
    !    write(19,'(F30.4,2X,F30.15,I20)')(drg * i), p_rg(i),i_p_rg(i)
    !    write(20,'(F30.4,2X,F30.15,I20)')(drg * i / rg_ave), p_rg(i) *rg_ave
    !enddo
    !print*,rg_ave
    !close(19)
    !close(20)
    deallocate(ree_save)
end program main
