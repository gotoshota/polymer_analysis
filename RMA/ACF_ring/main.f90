Program main
    use prm_list
    use read_prm
    use coarse_grain
    use read_dump
    implicit none

    character output*50
    CHARACTER mode_number*5
    double precision dummy, Xp2, t
    double precision,allocatable :: Xp(:,:,:), func(:), ampli(:)
    double precision,parameter :: pi = acos(-1.0d0)
    DOUBLE PRECISION tread, tcalc
    integer n, i, j, k, l, p, a, frag
    integer tbegin, tend, CountPerSec, CountMax

    call readprm()
    call coarsegrain()

    
    call system_clock(tbegin, CountPerSec, CountMax)
    call readdump()
    call system_clock(tend)
    tread = real(tend - tbegin) / CountPerSec

    allocate(Xp(3,nparticle,0:nframe)) 
    allocate(func(0:npoint))
    allocate(ampli(natom - 1))

    call system_clock(tbegin, CountPerSec, CountMax)
    a=0
    Xp=0.0d0
    ampli=0.0d0
    !Caluculate Xp
    do i = 0, nframe !frame
        do j = 0,nmol-1 !polychain
            do p = 0,natom-1 !mode
                do k = 1,natom !beads in a chain
                    Xp(:,j*natom+p+1,i)= Xp(:,j*natom+p+1,i) + pos(:,j*natom+k,i)* &
                    (cos(2.0d0 *dble(k) * DBLE(p) * pi / dble(natom)) + &
                    sin(2.0d0 *dble(k) *DBLE(p) *pi / dble(natom)))
                enddo
                Xp(:,j*natom+p+1,i) = Xp(:,j*natom+p+1,i) / SQRT(dble(natom))
            enddo
        enddo
    enddo
    DEALLOCATE(pos)
    !Calculate correlate function
    do p = 1, natom - 1
        func=0.0d0
        do i = 0,npoint !jikokusa
            a=0
            do j=0,nframe - Targetframe(i) !initial time
                a=a+1
                Xp2=0.0d0
                do k = 0,nmol-1 !which chain
                    Xp2 = Xp2 + dot_product(Xp(:,k*natom+p+1,j),Xp(:,k*natom+p+1,j+TargetFrame(i))) 
                enddo
                Xp2 = Xp2/DBLE(nmol)
                func(i) = func(i) + Xp2
            enddo
            func(i)=func(i)/DBLE(a)
        enddo
        
        WRITE(mode_number, '(i4)') p
        output = 'rmf_'//TRIM(ADJUSTL(mode_number))//'.xvg'
        open(17,file=output,status='replace')
        do i = 1,npoint
            t = Targetframe(i)* dt * nfreq
            write(17,'(F15.3,1X,e30.23e2)')t,func(i)/func(0)
        enddo
        close(17)
        ampli(p) = func(0)
    enddo
    
    output = 'rouse_ampli.xvg'
    open(17, file=output, status='replace')
    write(17,*)'#mode p, Rouse Amplitude Xp**2, sin(pi*p/2N)**2'
    do p = 1, natom - 1
       WRITE(17,'(I4, 1X, e30.23e2, 1X, F30.15)') p, ampli(p), sin(DBLE(p)*pi/DBLE(natom)/2.0d0)**2.0d0
    enddo
end Program main
