Program corre
    implicit none

    character filename*50
    double precision dummy, dt, Xp2, taup, func_m, func_l, func_h, func_03
    double precision rangemax, trange
    double precision,allocatable :: pos(:,:,:), Xp(:,:,:), func(:)
    double precision,parameter :: pi = acos(-1.0d0)
    integer n, atoms, i, mol, j, k, l, p, dump, step, frame, a, frag
    integer iframe, dframe

    print*,'Input trajectory filename.'
    read*,filename
    
    print*,'Input number of steps (Integer)'
    read*,step

    print*,'Input timestep (dt : float)'
    read*,dt

    print*,'Input dump frequency (integer)'
    read*,dump

    ! print*,'Input the number of particles in this system (integer)'
    ! read*,nparticle

    print*,'Input the number of atoms in a molecule (integer)'
    read*,n

    !nmol = nparticle / natom
    frame = step / dump

    rangeMax = (frame - 14) ** (1.0d0 / 15.0d0)
    print*,'Input the time range of point (float)'
    print*,'Reccomended value is (total points are 15)',rangemax
    read*,trange
    iframe = int(log(dble(frame)) / log(trange))
    if (trange .eq. 1.0d0) then
        iframe = frame
    endif

    print*,'#################################################'
    open(17,file=filename,status='old')

    print*,'Reading...'

    do i = 1,3
        read(17,'()')
    enddo
    read(17,*)atoms !Number of all atoms in simulation
    do i =1,5
        read(17,'()')
    enddo
    mol = atoms / N
    allocate(pos(3,atoms,0:frame))
    allocate(Xp(3,atoms,0:frame)) 
    allocate(func(0:iframe))
    

    do i = 1 , atoms
        read(17,*)dummy,pos(:,i,0)
    enddo
    do l = 1,frame
        do j=1,9
            read(17,'()')
        enddo
        do i=1,atoms
            read(17,*)dummy,pos(:,i,l)
        enddo
    enddo

    a=0
    !Caluculate Xp
    do i = 0, frame !frame
        do j = 0,mol-1 !polychain
            do p = 0,N-1 !mode
                do k = 1,N !beads in a chain
                    Xp(:,j*N+p+1,i)= Xp(:,j*N+p+1,i) + pos(:,j*N+k,i)* &
                    (cos(2.0d0 *dble(k) * DBLE(p) * pi / dble(n)) + &
                    sin(2.0d0 *dble(k) *DBLE(p) *pi / dble(n)))
                enddo
                Xp(:,j*n+p+1,i) = Xp(:,j*n+p+1,i) / sqrt(DBLE(N))
            enddo
        enddo
        !Xp(:,:,i) = Xp (:,:,i) / mol
    enddo
    close(17)

    deallocate(pos)

    !Calculate correlate function
    func=0.0d0
    print*,'Which mode? (p=0,1,2,...)'
    read*,p
    print*,'Input filename'
    read*,filename
    
    p=p+1
    do i = 0,iframe !jikokusa
        a=0
        dframe = int(trange ** i)
        do j=0,frame - dframe !initial time
            a=a+1
            Xp2=0.0d0
            do k = 0,mol-1 !which chain
                Xp2 = Xp2 + dot_product(Xp(:,k*N+p,j),Xp(:,k*N+p,j+i)) 
            enddo
            Xp2 = Xp2/mol
            func(i) = func(i) + Xp2
        enddo
        func(i)=func(i)/a
    enddo

    open(17,file=filename,status='replace')
    frag = 0
    do i = 0,iframe
        write(17,'(F30.3,1X,F30.15)')i*dt*dump,func(i)/func(0)
    enddo

    close(17)
end Program corre
