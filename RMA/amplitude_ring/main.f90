Program ini_corre
    use prm_list
    use read_prm
    use read_dump
    implicit none

    double precision,allocatable :: Xp(:,:), Xp2(:), Xp2_temp(:)
    double precision,parameter :: pi = acos(-1.0d0)
    integer i, j, k, l, dump, delta, id_atom

    CALL readprm()
    CALL readdump()

    allocate(Xp(3,0:natom-1)) !(xyz , modes)
    allocate(Xp2(0:natom-1))
    allocate(Xp2_temp(0:natom-1))
    Xp2=0.0d0
    Xp2_temp=0.0d0
    
    !Input coordinates
    do l = 0, nframe
        do i = 1,nmol !each molecule
            Xp = 0.0d0
            do j = 0,natom-1 !Rouse mode
                do k = 1,natom  !atoms in a molecule
                    id_atom = (i - 1) * natom + k
                    Xp(:,j) = Xp(:,j) + pos(:,id_atom, l)* &
                    (cos(2.0d0 *dble(k) * DBLE(j) * pi / dble(natom)) + &
                    sin(2.0d0 *dble(k) *DBLE(j) *pi / dble(natom)))
                enddo
                if(j == 0) then
                    delta = 1
                elseif (mod(natom,2) == 0 .and. j == natom/2) then
                    delta = 1
                else
                    delta = 0
                endif
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
    write(17,*)'#mode p, Rouse Amplitude Xp**2, sin(pi*p/Natom)**2'
    do i = 0,natom-1
        write(17,'(I6,1X,F30.15,1X,F30.15)')i,Xp2(i),(sin(DBLE(i)*pi/DBLE(natom)))**2.0
    enddo
    close(17)

end Program ini_corre
