module prm_list
    implicit none 
    integer nfreq, nparticle, natom, nmol, nstep, nframe
    double precision dt
    character(50) filename
    double precision, allocatable :: pos(:,:,:) 
    integer, parameter :: skip = 9
end module prm_list 
