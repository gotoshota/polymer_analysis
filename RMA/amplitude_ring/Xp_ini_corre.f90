Program ini_corre
    implicit none

    character filename*50
    double precision dummy, dt
    double precision,allocatable :: pos(:,:,:), Xp(:,:), Xp2(:), Xp2_temp(:)
    double precision,parameter :: pi = acos(-1.0d0)
    integer nparticle, natom, i, nmol, j, k, l, dump, step, nframe, delta
    integer skip

    print*,'Input trajectory filename.'
    read*,filename
    
    print*,'Input number of steps (Integer)'
    read*,step

    print*,'Input dump frequency (integer)'
    read*,dump
    
    print*,'Input timestep (dt : float)'
    read*,dt

    print*,'Input the number of particles in this system (integer)'
    read*,nparticle

    print*,'Input the number of atoms in a molecule (integer)'
    read*,natom

    nmol = nparticle / natom
    nframe = step / dump

    print*,'#################################################'

    open(17,file=filename,status='old')

    allocate(pos(3,nmol,natom))
    allocate(Xp(3,0:natom-1)) !(xyz , modes)
    allocate(Xp2(0:natom-1))
    allocate(Xp2_temp(0:natom-1))
    Xp2=0.0d0
    Xp2_temp=0.0d0
    
    !Input coordinates
    do l = 0, nframe
        do j=1,9
            read(17,'()')
        enddo
        do j=1,nmol
            do k = 1 , natom
                read(17,*)dummy,pos(:,j,k)
            enddo
        enddo
        do i = 1,nmol !each molecule
            do j = 0,natom-1 !Rouse mode
                Xp = 0.0d0
                do k = 1,natom  !atoms in a molecule
                    Xp(:,j) = Xp(:,j) + pos(:,i,k)* (cos(2.0d0 *dble(k) * DBLE(j) * pi / dble(natom)) + &
                    sin(2.0d0 *dble(k) *DBLE(j) *pi / dble(natom)))
                enddo
                !if(j == 0) then
                !    delta = 1
                !elseif (mod(natom,2) == 0 .and. j == natom/2) then
                !    delta = 1
                !else
                !    delta = 0
                !endif
                delta=1
                Xp(:,j) = Xp(:,j)*sqrt(dble(2-delta)/dble(natom))
                Xp2_temp(j) = Xp2_temp(j) + dot_product(Xp(:,j),Xp(:,j))
            enddo
        enddo
        do j = 0,natom-1
            Xp2_temp(j) = Xp2_temp(j)/dble(nmol)
            Xp2(j) = Xp2(j) + Xp2_temp(j)
            Xp2_temp(j)=0.0d0
        enddo
    enddo
    Xp2(:)=Xp2(:)/DBLE(nframe+1)

    close(17)
    open(17,file='Xp_inicorre.xvg',status='replace')
    write(17,*)'#mode p, Rouse Amplitude Xp**2, sin(pi*p/N)**2'
    do i = 0,natom-1
        write(17,'(I6,1X,F30.15,1X,F30.15)')i,Xp2(i),(sin(DBLE(i)*pi/DBLE(natom)))**2.0
    enddo
    close(17)

end Program ini_corre
