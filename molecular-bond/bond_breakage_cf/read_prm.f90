module read_prm 
    implicit none
    contains
    subroutine readprm()
        use prm_list
        implicit none 
        integer :: prmfile = 11
        integer*4 access
        
        NAMELIST /param/ filename, nstep, nfreq, dt, nparticle, natom
        if (access( './param.nml',' ') .eq. 0) then
            OPEN(prmfile, file="param.nml", status="old")
            READ(prmfile, nml=param)
            CLOSE(prmfile)
        endif

        !! for old parameter file format
        if (nparticle .eq. 0) then
            if (access( './param.txt',' ') .eq. 0) then
                open (prmfile,file='param.txt',status='old')

                read(prmfile,*)filename
                
                read(prmfile,*)nstep
                
                read(prmfile,*)nfreq
                
                read(prmfile,*)dt
                
                read(prmfile,*)nparticle
                
                read(prmfile,*)natom

            else

                print*,''

                print*,'Input trajectory filename'
                read*,filename
                print*,''
                
                print*,'Input the number of Step'
                read*,nstep
                print*,''
                
                print*,'Input the frequency of the dump file'
                read*,nfreq
                print*,''
                
                print*,'Input dt'
                read*,dt
                print*,''
                
                print*,'Input the number of particles'
                read*,nparticle
                print*,''
                
                print*,'Input the number of atoms in a molecule'
                read*,natom
                print*,''

            endif
        endif

        nmol    =   nparticle / natom
        nframe  =   nstep / nfreq
        allocate (pos(3,nparticle,0:nframe))
    end subroutine
end module read_prm
