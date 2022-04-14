Program main
    use prm_list
    use read_prm
    use read_dump
    use calc_poscm
    use calc_inertia_tensor

    implicit none
    integer tbegin, tend, CountPerSec, CountMax
    integer i, j, k, m, r_max, i_r, start
    integer delta_frame
    integer, allocatable :: i_vhf(:,:), i_vhf_self(:,:)
    double precision calc_t, t, box_l, x_1, x_2, dummy, box_half, dpos(3), distance, dv, t_vanhove
    DOUBLE PRECISION theta !! angule of molecular axis  at time t
    INTEGER theta_max, i_theta
    DOUBLE PRECISION limit_ri, r, integ, rg_estimate
    double precision, parameter :: dr = 0.20d0
    double precision, parameter :: pi = acos(-1.0d0)
    double precision, parameter :: b2 = 3.0d0  !! bead-spring model
    DOUBLE PRECISION, PARAMETER :: dtheta = 0.50d0
    double precision, allocatable :: vhf(:,:), vhf_self(:), vhf_nonself(:)
    
    call readprm()
    call system_clock(tbegin, CountPerSec, CountMax)
    call readdump()
    CALL calcposcm()
    CALL calcinertiatensor()

    !! Determin t
    Print*,'This program output self part of van hove function Gs(r,theta,t) as 2D by fixing time. (parameter is only distance r)'
    print*,'Time t of self van hove function Gs(r,theta,t)'
    print*,'should be multiple of i*nfreq*dt'
    print*,'Input the prefactor i (integer). imax=',nframe 
    read(*,*)delta_frame
    t_vanhove = nfreq *delta_frame *dt
    if (delta_frame .gt. nframe) then
        print*,'Error: prefactor i is larger then nframe'
        stop
    endif

    !!Read box_l
    open(16, file = filename, status = 'old')
    do i = 1, 5
        read(16,'()')
    enddo
    read(16,'(F18.16,1X,I3,1X,F18.16,1X,I3)')x_1, i, x_2, j
    close(16)
    box_l = x_2 * 10.0d0 ** j - x_1 * 10.0d0 ** i
    box_half = box_l / 2.0d0
    rg_estimate = SQRT(b2 *natom/6)  
    r_max = ceiling(rg_estimate*4 / dr)
    !! end reading box_l

    
    theta_max = 90 / dtheta
    allocate(vhf(r_max, theta_max))
    allocate(i_vhf(r_max, theta_max))
    vhf             = 0.0d0
    i_vhf           = 0 
    !! Calculate distance and count particles
    if (nframe .eq. 1)then
        start = 1
    else
        start = 0
    endif
    do i = start, nframe - delta_frame
        do j = 1, nmol !!set j th particle on origin
            !! Calculate displacement 
            dpos(:) = pos_cm(:,j,i + delta_frame) - pos_cm(:,j,i)
            distance = sqrt(dpos(1)**2.0d0 + dpos(2)**2.0d0 + dpos(3)**2.0d0)
            i_r = ceiling(distance / dr)
            if (i_r .eq. 0)then
                i_r = 1
            endif 

            !! Calculate angle 
            theta = inertia_tensor(1,1,j,i+delta_frame) * inertia_tensor(1,1,j,i) &
                  + inertia_tensor(2,1,j,i+delta_frame) * inertia_tensor(2,1,j,i) &
                  + inertia_tensor(3,1,j,i+delta_frame) * inertia_tensor(3,1,j,i)
            theta = ACOS(theta) / pi * 180.0d0
            if (theta .lt. 0 .or. theta .gt. 90) print *, 'Error: the value of theta'
            i_theta = theta / dtheta

            if (i_r <= r_max .and. i_theta <= theta_max) then
                i_vhf(i_r,i_theta) = i_vhf(i_r,i_theta) + 1
            else
                print*,'Error: Missing particle (line 87) '
            endif
        enddo
    enddo
    


    vhf(:,:)  = dble(i_vhf(:,:)) / dble(nmol) /dble(nframe - delta_frame)

    !print*,sum(vhf)
    vhf(:,:) = vhf(:,:) /dr /dtheta
    !! Caluculate vhf
    !do i = 1, r_max
    !    dv                = (4.0d0 / 3.0d0)*pi* (dble(i) ** 3.0d0 - dble(i-1) ** 3.0d0) * dr **3.0d0
    !    vhf(i)            = vhf(i) / dv
    !enddo
    call system_clock(tend)
    calc_t = real(tend - tbegin) / CountPerSec

    open (16, file ='vanhove_func.xvg',status='replace')
    write(16,'(A17,1X,F10.3,1X,A3)')'##elapsed time is',calc_t,'sec'
    write(16,*)'## t = ',t_vanhove
    write(16,*)'## r, theta,  van hove function(VHF) self'
    do i = 1, r_max
        do j = 1, theta_max
            write(16,'(2(F15.7,1X),F30.20)')dr*(i-0.5), dtheta*(j-0.5), vhf(i,j)
        enddo
        WRITE(16,'()')
    enddo
    close(16)
end Program main
