!This program calculate Aspherity parameter (AS).  !The algorism is based on following paper.
!K. Alim and E. Frey, "Shapes of Semiflexible Polymer Rings", PRL (2007).
program main
    use prm_list
    use read_prm
    use read_dump
    use calc_poscm
    use calc_inertia_tensor
    
    implicit none
    DOUBLE PRECISION, ALLOCATABLE :: inertia_tensor(:,:,:,:) !! output as eigenvector 
    DOUBLE PRECISION, ALLOCATABLE :: inertia_tensor_hat(:,:,:,:) !! output as eigenvector 
    DOUBLE PRECISION, ALLOCATABLE :: inertia_tensor_hat_square(:,:,:,:) !! output as eigenvector 
    DOUBLE PRECISION tmp(3), tmptmp, tmp_tensor
    DOUBLE PRECISION ave(2)
    DOUBLE PRECISION mat(3,3)
    double precision t, dummy, tread, tcalc
    integer h, i, j, k, l, m, n
    integer tbegin, tend, CountPerSec, CountMax
    INTEGER normalize
    
    DOUBLE PRECISION, ALLOCATABLE :: delta(:,:), sigma(:,:)
    DOUBLE PRECISION trace, tracetrace
    DOUBLE PRECISION det

    DOUBLE PRECISION eigenval(3)
    integer, parameter :: lwork = 10000
    double precision :: work(lwork)
    integer :: info
!############################################################################

    call readprm()
    call system_clock(tbegin, CountPerSec, CountMax)
    call readdump()
    call system_clock(tend)
    tread = real(tend - tbegin) / CountPerSec

    CALL calcposcm()

    !! calculate inertia ternsor
    ALLOCATE(inertia_tensor(3, 3, nmol, 0:nframe)) !! 1=x, 2=y, 3=z
    ALLOCATE(inertia_tensor_hat(3, 3, nmol, 0:nframe)) !! 1=x, 2=y, 3=z
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
    !calculate inertia tensor hat
    do i = 0, nframe
        do j = 1, nmol
            trace = inertia_tensor(1,1,j,i) + inertia_tensor(2,2,j,i) + inertia_tensor(3,3,j,i)
            trace = trace / 3.0d0
            do k = 1, 3
                do l = 1, 3
                    inertia_tensor_hat(l,k,j,i) = inertia_tensor(l,k,j,i)
                enddo
                inertia_tensor_hat(k,k,j,i) = inertia_tensor_hat(k,k,j,i) - trace
            enddo
        enddo
    enddo
    !! calculate square inertia tensor hat
    do i = 0, nframe
        do j = 1, nmol
            mat(:,:) = inertia_tensor_hat(:,:,j,i)
            inertia_tensor_hat_square(:,:,j,i) = MATMUL(mat,mat)
        enddo
    enddo
    !! calculate delta
    ALLOCATE(delta(nmol,0:nframe))
    do i = 0, nframe
        do j = 1, nmol
            !! delta
            trace = inertia_tensor(1,1,j,i) + inertia_tensor(2,2,j,i) + inertia_tensor(3,3,j,i)
            tracetrace = inertia_tensor_hat_square(1,1,j,i) + inertia_tensor_hat_square(2,2,j,i) + inertia_tensor_hat_square(3,3,j,i) 
            delta(j,i) = 1.50d0*tracetrace/trace/trace
            !! sigma
            mat = inertia_tensor_hat(:,:,j,i)
            CALL dsyev('V', 'U', 3, mat, 3, eigenval, work, lwork, info)
            det = eigenval(1) + eigenval(2) + eigencal(3)
            tmp = 2.0d0 / 3.0d0 *tracetrace
            tmp = tmp ** 1.50d0
            sigma(j,i) = 4.0d0*det / tmp
        enddo
    enddo

    !! calculate average
    ave = 0.0d0
    do i = 0, nframe
        do j = 1, nmol
            ave(1) = ave(1) + delta(j,i)
            ave(2) = ave(2) + sigma(j,i)
        enddo
    enddo
    ave = ave / DBLE(nframe + 1) / DBLE(nmol)

    call system_clock(tend)
    tcalc = real(tend - tbegin) / CountPerSec
    
    deallocate(inertia_tensor, inertia_eigenval)

    open(15, file='aspherity_parameter.txt', status='replace')

    write(15,*)'# aspherity parameter'
    write(15,*)'# delta is:'
    write (15,*)ave(1)
    write(15,*)'# sigma is:'
    write (15,*)ave(2)
    close(15)
    open(17, file='log.as_param', status='replace')
    write (17,*)'read time is ',tread,'sec'
    write (17,*)'calclation time is ',tcalc,'sec'
    close(17)
    print *,'# aspherity parameter'
    print *,'# delta is:'
    print *,ave(1)
    print *,'# sigma is:'
    print *,ave(2)
    


end program main
