!This program calculate Aspherity parameter (AS).
!The algorism is based on following paper.
!J. Swetlana, et al., "Diffusion- and reaction-limited cluster aggregation revisited"
!PCCP, 2019.
program main
    use prm_list
    use read_prm
    use read_dump
    use calc_poscm
    use calc_inertia_tensor
    
    implicit none
    DOUBLE PRECISION, ALLOCATABLE :: inertia_tensor(:,:,:,:) !! output as eigenvector 
    DOUBLE PRECISION tmp(3), tmptmp
    DOUBLE PRECISION as
    DOUBLE PRECISION ave(3)
    double precision t, dummy, tread, tcalc
    integer h, i, j, k, l, m, n
    integer tbegin, tend, CountPerSec, CountMax
    INTEGER normalize
    
    
!############################################################################

    call readprm()
    call system_clock(tbegin, CountPerSec, CountMax)
    call readdump()
    call system_clock(tend)
    tread = real(tend - tbegin) / CountPerSec

    CALL calcposcm()

    !! calculate inertia ternsor
    ALLOCATE(inertia_tensor(3, 3, nmol, 0:nframe)) !! 1=x, 2=y, 3=z
    do i = 0, nframe
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
    DEALLOCATE(pos, pos_cm)

    call system_clock(tend)
    tcalc = real(tend - tbegin) / CountPerSec
    
    deallocate(inertia_tensor, inertia_eigenval)

    open(15, file='aspherity_parameter.txt', status='replace')

    write(15,*)'# aspherity parameter is:'
    write (15,*)as
    close(15)
    open(17, file='log.as_param', status='replace')
    write (17,*)'read time is ',tread,'sec'
    write (17,*)'calclation time is ',tcalc,'sec'
    close(17)
    
    print*, '# AS param is:'
    print*, as


end program main
