Program rdf_calc
    use prm_list
    use read_prm
    use read_dump

    implicit none
    integer tbegin, tend, CountPerSec, CountMax
    integer i, j, k, m, r_max, i_r, start, skip_frame
    integer, allocatable :: i_rdf(:), i_rdf_bonded(:), i_rdf_nonbonded(:)
    INTEGER, parameter :: num_calc_frame = 2000
    double precision t, box_l, x_1, x_2, dummy, rho, box_half, dpos(3), distance, dv
    double precision, parameter :: dr=0.020d0
    double precision, parameter :: pi = acos(-1.0d0)
    double precision, allocatable :: rdf(:), rdf_bonded(:), rdf_nonbonded(:)
    
    call readprm()
    call system_clock(tbegin, CountPerSec, CountMax)
    call readdump()

    !!Read box_l
    open(16, file = filename, status = 'old')
    do i = 1, 5
        read(16,'()')
    enddo
    read(16,'(F18.16,1X,I3,1X,F18.16,1X,I3)')x_1, i, x_2, j
    close(16)
    box_l = x_2 * 10.0d0 ** j - x_1 * 10.0d0 ** i
    rho = dble(nparticle) / box_l**3.0d0 !!Density
    box_half = box_l / 2.0d0
    r_max = ceiling(box_half / dr)
    !! end reading box_l
    
    allocate(rdf(r_max), rdf_bonded(r_max), rdf_nonbonded(r_max))
    allocate(i_rdf(r_max), i_rdf_bonded(r_max), i_rdf_nonbonded(r_max))
    rdf             = 0.0d0
    rdf_bonded      = 0.0d0
    rdf_nonbonded   = 0.0d0
    i_rdf           = 0 
    i_rdf_bonded    = 0
    i_rdf_nonbonded = 0
    !! Calculate distance and count particles
    if (nframe .eq. 1)then
        start = 1
    else
        start = 0
    endif
    if (nframe .le. num_calc_frame) then
        skip_frame = 1
    else
        skip_frame = CEILING(nframe / num_calc_frame)
    endif
    do i = start, nframe, skip_frame
        do j = 1, nparticle !!set j th particle on origin
            m = ceiling(real(j) / natom)
            do k = 1, nparticle
                if (j .ne. k) then
                    dpos(:) = pos(:,j,i) - pos(:,k,i)
                    dpos(:) = abs(dpos(:) - box_l * nint(dpos(:) / box_l))
                    distance = sqrt(dpos(1)**2.0d0 + dpos(2)**2.0d0 + dpos(3)**2.0d0)
                    i_r = ceiling(distance / dr)
                    if (i_r .eq. 0)then
                        i_r = 1
                    endif 
                    if (i_r <= r_max) then
                        i_rdf(i_r) = i_rdf(i_r) + 1
                        if (ceiling(real(k) / natom) .eq. m)then
                            i_rdf_bonded(i_r) = i_rdf_bonded(i_r) + 1
                            !if (i==1 .and. j <= 10)then
                            !    print*,j,m,k
                            !endif
                        elseif (ceiling(real(k) / natom) .ne. m)then 
                        !if(distance <0.98) print *, j,k
                            i_rdf_nonbonded(i_r) = i_rdf_nonbonded(i_r) + 1
                        else
                        print*,'Error: molecule number'
                        endif
                    endif
                endif
            enddo
        enddo
    enddo
    
    rdf               = dble(i_rdf) / (nparticle * num_calc_frame)
    rdf_bonded        = dble(i_rdf_bonded) / (nparticle * num_calc_frame)
    rdf_nonbonded     = dble(i_rdf_nonbonded) / (nparticle * num_calc_frame)

    !! Caluculate RDF
    do i = 1, r_max
        dv                  = (4.0d0 / 3.0d0)*pi* (dble(i) ** 3.0d0 - dble(i-1) ** 3.0d0) * dr **3.0d0
        rdf(i)              = rdf(i) / (dv * rho)
        rdf_bonded(i)       = rdf_bonded(i) / (dv * rho)
        rdf_nonbonded(i)    = rdf_nonbonded(i) / (dv * rho)
    enddo
    call system_clock(tend)
    t = real(tend - tbegin) / CountPerSec
    
    open (16, file ='rdf.xvg',status='replace')
    write(16,'(A17,1X,F10.3,1X,A3)')'##elapsed time is',t,'sec'
    write(16,*)'## distance, rdf(ALL), rdf(Bonded), rdf(nonbonded)'
    do i = 1, r_max
        write(16,'(F15.7,1X,F15.7,1X,F15.7,1X,F15.7)')(i-0.5d0)*dr,rdf(i),rdf_bonded(i),rdf_nonbonded(i)
    enddo
    close(16)
end Program rdf_calc
