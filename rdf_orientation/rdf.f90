Program rdf_calc
    use prm_list
    use read_prm
    use read_dump
    use calc_poscm
    use calc_inertia_tensor


    implicit none

    integer tbegin, tend, CountPerSec, CountMax
    integer i, j, k, m, r_max, i_r, start, skip_frame
    INTEGER id_atom, theta_max, i_theta

    double precision t, box_l, x_1, x_2, dummy, rho, box_half, dpos(3), distance, dv
    DOUBLE PRECISION pos_cm_tmp(3), rho_cm_cm
    DOUBLE PRECISION dot_pro, theta, size_dot_pro

    integer, allocatable :: i_rdf(:), i_rdf_dotpro(:), i_rdf_dotpro2(:)
    double precision, allocatable :: rdf_cm(:), rdf_cm_dotpro(:), rdf_cm_dotpro2(:)

    INTEGER, parameter :: num_calc_frame = 2000
    double precision, parameter :: dr=0.020d0
    double precision, parameter :: pi = acos(-1.0d0)

! #######################################################################

    call readprm()
    call system_clock(tbegin, CountPerSec, CountMax)
    call readdump()
    call calcposcm()
    CALL calcinertiatensor()

    !!Read box_l
    open(16, file = filename, status = 'old')
    do i = 1, 5
        read(16,'()')
    enddo
    read(16,'(F18.16,1X,I3,1X,F18.16,1X,I3)')x_1, i, x_2, j
    close(16)
    box_l = x_2 * 10.0d0 ** j - x_1 * 10.0d0 ** i
    rho        = dble(nparticle) / box_l**3.0d0 !!Density
    rho_cm_cm  = dble(nmol) / box_l**3.0d0 
    box_half = box_l / 2.0d0
    r_max = ceiling(box_half / dr)
    !! end reading box_l
    
    allocate(rdf_cm(r_max), rdf_cm_dotpro(r_max), rdf_cm_dotpro2(r_max))
    allocate(i_rdf(r_max), i_rdf_dotpro(:), i_rdf_dotpro2(:))

    rdf_cm_         = 0.0d0
    rdf_cm_2        = 0.0d0
    i_rdf           = 0 
    i_rdf_dotpro    = 0 
    i_rdf_dotpro2   = 0 
    !! Calculate distance and count particles
    if (nframe .eq. 1)then
        start = 1
    else
        start = 0
    endif
    if (nframe .le. num_calc_frame) then
        skip_frame = 1
    else
        skip_frame = CEILING(DBLE(nframe) / DBLE(num_calc_frame))
    endif

    
    !! cm_cm
    do i = start, nframe, skip_frame
        do j = 1, nmol !!set j th particle on origin
            do k = 1, nmol
                if (j .ne. k) then
                    !! calculate distance
                    dpos(:) = pos_cm(:,j,i) - pos_cm(:,k,i)
                    dpos(:) = abs(dpos(:) - box_l * nint(dpos(:) / box_l))
                    distance = dsqrt(dpos(1)**2.0d0 + dpos(2)**2.0d0 + dpos(3)**2.0d0)
                    i_r = ceiling(distance / dr)

                    !! calculate theta
                    dot_pro = inertia_tensor(1,3,j,i)* inertia_tensor(1,3,k,i) &
                            + inertia_tensor(2,3,j,i)* inertia_tensor(2,3,k,i) &
                            + inertia_tensor(3,3,j,i)* inertia_tensor(3,3,k,i)

                    if (i_r .eq. 0)then
                        i_r = 1
                    endif 
                    if (i_r <= r_max) then
                        i_rdf(i_r) = i_rdf(i_r) + 1
                        i_rdf_dotpro(i_r) = i_rdf_dotpro(i_r) + dot_pro
                        i_rdf_dotpro2(i_r) = i_rdf_dotpro2(i_r) + 0.50d0*(3.0d0*dot_pro*dot_pro - 1.0d0)
                    endif
                endif
            enddo
        enddo
    enddo
    rdf_cm            = dble(i_rdf) / (nmol * num_calc_frame)
    rdf_cm_dotpro     = dble(i_rdf_dotpro) / (nmol * num_calc_frame)
    rdf_cm_dotpro2    = dble(i_rdf_dotpro2) / (nmol * num_calc_frame)
    
    !! Caluculate RDF
    do i = 1, r_max
        dv                  = (4.0d0 / 3.0d0)*pi* (dble(i) ** 3.0d0 - dble(i-1) ** 3.0d0) * dr **3.0d0
        rdf_cm_cm(i)              = rdf_cm_cm(i) / (dv * rho_cm_cm)
        rdf_orient(i)             = rdf_orient(i) / (dv * rho_cm_cm)
    enddo
    call system_clock(tend)
    t = real(tend - tbegin) / CountPerSec
    
    open (16, file ='rdf_cm.xvg',status='replace')
    write(16,'(A17,1X,F10.3,1X,A3)')'##elapsed time is',t,'sec'
    write(16,*)'## distance, rdf cm-cm, rdf cm-mon'
    do i = 1, r_max
        write(16,'(F15.7,1X,F15.7,1X,F15.7,1X,F15.7)')(i-0.5d0)*dr,rdf_cm(i),rdf_cm_dotpro(i),rdf_cm_dotpro2(i)
    enddo
    close(16)
end Program rdf_calc
