program pdr1

! module
use xdr , only : xtcfile
use info
implicit none 
type(xtcfile) :: xtc

!======================================================================!
!============================== variable ==============================!
!======================================================================!
! parameter 
integer :: M2, N, OP
real*8, parameter :: logn = 1.3d0

! allocatable
integer*1, allocatable :: n_j1(:), n_j2(:)
integer*1, allocatable :: l_j1(:), l_j2(:)
integer*4, allocatable :: h(:,:,:,:)
integer*8, allocatable :: nt(:,:)
integer*1, allocatable :: otype(:,:,:,:)
integer*4, allocatable :: flag(:,:,:,:)
real*8, allocatable :: pos(:,:)
integer*4, allocatable :: list_call(:), list(:,:)

! almost parameter from module
character(50) :: trjfile
integer :: npol, nsol, nstep
real*8 :: dt
real*8 :: box_dim
integer :: logstep
integer :: na_p, tot_p,tot_o, natoms

! others
integer :: i, j, k, l,m
integer :: kk, j2, j3
integer :: step, kstep, interval
integer :: posow, posh1, posh2, poso
real*8 :: v_oo(3), v_oh1(3), v_oh2(3)
real*8 :: r_oo, r_oh1, r_oh2
real*8 :: beta1, beta2
real*8 :: n0(3)
real*8 :: pdr, t
character(50) :: num, OUT1
!======================================================================!


print *, '*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*'
print*, ''
print *, '         == TCF of H_bond broken but near=='
print *, ''

! MDinfo
call mdinfo(trjfile,npol,nsol,nstep,dt)

! open xtc and read 1st flame
call xtc % init(trjfile)
call xtc % read

! parameter 
call pol_info(M2,N,OP, nsol, npol, na_p, tot_p, tot_o, natoms)
box_dim = xtc % box(1,1)
logstep = log_info(logn, nstep)

! allocate
allocate(n_j1(nsol))
allocate(n_j2(nsol))
allocate(l_j1(nsol))
allocate(l_j2(nsol))
allocate(h(10,nsol,nstep + 1,2))
allocate(nt(3, 0:logstep))
allocate(otype(10, nsol, nstep + 1, 2))
allocate(flag(20,nsol, nstep + 1, 2))
allocate(pos(3,natoms))
allocate(list_call(npol*N))
allocate(list(npol*N,3))

! make index of polymer
do i = 1, 3
  call make_ndx(M2, N, OP, npol,i, list_call)
  list(:,i) = list_call(:)
end do

! XXXX = 0

n_j1(:)=0
n_j2(:)=0
l_j1(:)=0
l_j2(:)=0
h(:,:,:,:)=0
nt(:,:) = 0
otype(:,:,:,:) = 0
flag(:,:,:,:) = 0
n0(:)=0
print *, 'Start Checking Hbond ...'

do step = 1, nstep + 1
  pos = xtc % pos(:,:)
  print * , step
!$omp parallel shared(nsol, pos, otype, h, npol, na_p, tot_p, box_dim, flag, natoms,step) &
!$omp & private(v_oh1, v_oh2,v_oo,poso, posow, posh1, posh2,r_oh1,r_oh2,r_oo, beta1, beta2, kk) &
!$omp & reduction(+:n_j1, n_j2, l_j1, l_j2)
!$omp do
  do i = 1, nsol
    call pos_water(i, tot_p, posow, posh1, posh2)

    v_oh1 =  pos (:, posh1 ) - pos(:, posow )
    v_oh2 =  pos (:, posh2 ) - pos(:, posow )

    call length(box_dim, v_oh1, r_oh1)
    call length(box_dim, v_oh2, r_oh2)

    kk = 3
    do j = 1, nsol
      
      poso = tot_p + (j-1)*4 + 1
      v_oo = pos(:,poso) - pos(:,posow)
      call length(box_dim, v_oo, r_oo)

      call degree(v_oh1, v_oo, r_oh1, r_oo, beta1)
      call degree(v_oh2, v_oo, r_oh2, r_oo, beta2)

      if (r_oo > 0.0 .and. r_oo <= 0.35d0) then
        if( beta1 <= 30.0d0) then
          n_j1(i) = n_j1(i) + 1
          h(n_j1(i), i, step, 1) = poso
          otype(n_j1(i),i,step,1) = kk
        else
          l_j1(i) = l_j1(i) + 1
          flag(l_j1(i), i, step, 1) = poso
        end if

        if (beta2 <= 30.0d0 ) then
          n_j2(i) = n_j2(i) + 1
          h(n_j2(i), i, step, 2) = poso
          otype(n_j2(i),i,step,2) = kk
        else
          l_j2(i) = l_j2(i) + 1
          flag(l_j2(i), i, step, 2) = poso
        end if

      end if
    end do

    do k = 1,2
      kk = 1 + (k - 1)*2
      do j = 1, npol*N


          poso = list(j,kk)
          v_oo = pos(:, poso ) - pos(:, posow)
          call length(box_dim, v_oo, r_oo)

          call degree(v_oh1, v_oo, r_oh1, r_oo, beta1)
          call degree(v_oh2, v_oo, r_oh2, r_oo, beta2)

          if (r_oo <= 0.35d0) then
    
            if( beta1 <= 30.0d0) then
              n_j1(i) = n_j1(i) + 1
              h(n_j1(i), i, step, 1) = poso
              otype(n_j1(i),i,step,1) = k
            else
              l_j1(i) = l_j1(i) + 1
              flag(l_j1(i), i, step, 1) = poso
            end if
    
            if (beta2 <= 30.0d0 ) then
              n_j2(i) = n_j2(i) + 1
              h(n_j2(i), i, step, 2) = poso
              otype(n_j2(i),i,step,2) = k
            else
              l_j2(i) = l_j2(i) + 1
              flag(l_j2(i), i, step, 2) = poso
            end if
    
          end if


      end do
    end do
  end do
!$omp end do
!$omp end parallel
  n_j1(:) = 0
  n_j2(:) = 0
  l_j1(:) = 0
  l_j2(:) = 0
  call xtc % read
end do

print *, 'END Checking Hbond !!!'
print *, 'Start Counting H-Bond Time ...'


!$omp parallel shared(nstep, nsol, otype, h) &
!$omp & private(i,j,l, kk) &
!$omp & reduction(+: n0)
!$omp do
do k = 1, nstep + 1
  do i = 1, nsol
    do l = 1, 2
      do j2 = 1, 10
        j = h(j2,i,k,l)
        kk = otype(j2,i,k,l)
        if (j /= 0 ) then
          if ( kk == 3) then
            n0(3) = n0(3) + 1d0
          elseif ( kk == 1 ) then
            n0(1) = n0(1) + 1d0
          elseif ( kk == 2 ) then
            n0(2) = n0(2) + 1d0
          end if
        else
          exit
        end if
      end do
    end do
  end do
end do
!$omp end do
!$omp end parallel

print *, n0
!$omp parallel shared(nstep, nsol, otype, flag, h) &
!$omp & private(j, kk, kstep, interval) &
!$omp & reduction(+: nt)
!$omp do
do step = 0, logstep
  interval = int(logn**step)
  kstep = nstep + 1 - interval
  do k = 1, kstep
    do i = 1, nsol
      do l = 1, 2
        do j2 = 1, 10
          j = h(j2,i,k,l)
          kk = otype(j2,i,k,l)
          if(j /= 0) then
            do j3 = 1, 20
              if (flag(j3,i,k+interval,l) == 0 ) then 
                exit
              end if

              if (flag(j3,i,k+interval,l) == j ) then
                if ( kk == 1) then
                  nt(1,step) = nt(1,step) + 1
                elseif ( kk == 2 ) then
                  nt(2,step) = nt(2,step) + 1
                elseif ( kk == 3 ) then
                  nt(3,step) = nt(3,step) + 1
                end if
                exit
              end if
            end do
          end if
        end do
      end do
    end do
  end do
print*, step, nt(:,step)
end do
!$omp end do
!$omp end parallel

do kk = 1,3
  t=0d0
  print *, n0(kk)
  n0(kk) = (n0(kk))/dble(nstep+1)
  print *, n0(kk)
  num=''
  write(num,'(I0)') kk
  OUT1 = 'pdr'//trim(num)//'.xvg'
  open(17,file=OUT1)

  write(17,'(2f10.3)') t, 0d0

  do i = 0, logstep
    interval = int(logn**i)
    t = interval * dt
    pdr = dble(nt(kk,i))/dble(nstep +1 - interval)
    pdr = pdr/n0(kk)
    write(17,'(2f10.3)') t, pdr
  end do

  close(17)
end do

print *,''
print *, '             == MISSION COMPLETED !!! =='
print *, ''
print *, '*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*'
end program pdr1
