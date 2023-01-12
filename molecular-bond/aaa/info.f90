module info
implicit none
contains

subroutine mdinfo(trjfile,npol,nsol,nstep,dt)

integer, intent(out) :: nsol, nstep, npol
real*8, intent(out) :: dt
character(50), intent(out) :: trjfile

namelist /system_info/trjfile, npol, nsol, nstep, dt

open (17,file='MDinfo', status = 'old')
read(17,nml=system_info) 
close(17)

print *, ''
print *, '= MDinfo ='
print*,  'polymer : ', npol, 'molecules'
print*,  'water   : ', nsol, 'molecules'
write(*,'(a11, i8)') ' flames   : ', nstep
write(*,'(a11, f9.3)') ' dt [ps] : ', dt
write(*,'(a16, f10.2)') ' Total t [ns] : ', nstep*dt/1000

end subroutine mdinfo


subroutine pos_water(i, tot_p, posow, posh1, posh2)
  integer, intent(in) :: tot_p, i
  integer, intent(out) :: posow, posh1, posh2

  posow = tot_p + (i-1)*4 + 1
  posh1 = tot_p + (i-1)*4 + 2
  posh2 = tot_p + (i-1)*4 + 3
end subroutine pos_water

integer function func_poso(j, tot_p)
  implicit none
  integer, intent(in) :: tot_p, j
  func_poso = tot_p + (j-1)*4 + 1
end function func_poso

subroutine pol_info(M2,N,OP, nsol, npol, na_p, tot_p, tot_o, natoms )
  integer , intent(in) :: npol, nsol
  integer, intent(out) :: M2, N,OP
  integer, intent(out) :: na_p, tot_p, natoms, tot_o

  namelist /polymer_info/N,M2,OP

  open(17,file='POLinfo', status='old')
  read(17,nml=polymer_info)
  close(17)

  na_p = M2*N +2
  tot_p = na_p * npol
  tot_o = npol*N*2 + nsol
  natoms = tot_p + 4* nsol
end subroutine pol_info

integer function log_info(logn, nstep)
  implicit none
  real*8, intent(in) :: logn
  integer, intent(in) :: nstep

  log_info = int(log(dble(nstep))/log(logn))
end function log_info

subroutine length(box_dim, v, r)
  real*8, intent(in) :: box_dim
  real*8, intent(inout) :: v(3)
  real*8, intent(out) :: r

   v = v - box_dim * nint(v / box_dim)
   r = sqrt(sum(v**2))
end subroutine length

subroutine degree(v_oh, v_oo, r_oh, r_oo, beta)
real*8, intent(in) :: v_oh(3), v_oo(3), r_oh, r_oo
real*8, intent(out) :: beta

real*8, parameter :: pi = dacos(-1d0)
real*8 :: cosb

  cosb = (sum(v_oh*v_oo))/(r_oo*r_oh)
  if ( cosb > 1d0 ) then
    cosb = 1d0
  elseif ( cosb < -1d0 ) then
    cosb = -1d0
  end if

  beta = dacos(cosb)
  beta = 180*beta/pi
end subroutine degree

subroutine make_ndx(M2, N, OP, npol, k, list)
  integer, intent(in) :: M2, N, npol, OP, k
  integer, allocatable, intent(out) :: list(:)
  integer :: cn, i, j, na_p

  allocate(list(npol*N))
  na_p = M2*N +2
  
  cn = 0
  do i = 1, npol
    cn = cn + 1
    list(cn) =  na_p*(i-1) + OP + k
    do j = 2, N
      cn = cn + 1
      list(cn) = na_p*(i-1) + 1 + M2*(j-1) + OP + k 
    end do
  end do

end subroutine make_ndx

integer function func_bw(card, bwl, bwr)
  implicit none
  integer, intent(inout) :: bwl, bwr 
  integer, intent(in) :: card
  integer :: dumy, j, bw

  bw = 0
  if ( bwl > bwr ) then
    dumy = bwl
    bwl = bwr
    bwr = dumy
  end if

  if ( bwl /= 0 ) then
    do j = 0, bwl - 1
        bw = bw + (card - j)
    end do ! j
  end if

  bw = bw + bwr - bwl
  func_bw = bw
end function func_bw


end module info
