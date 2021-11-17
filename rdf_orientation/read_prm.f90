module read_prm 
    implicit none
    contains
    subroutine readprm()
        use prm_list
        implicit none 
        integer prmfile
        integer*4 access
        
        if (access( './param.txt',' ') .eq. 0) then
            prmfile=11
            open (prmfile,file='param.txt',status='old')

            read(prmfile,*)filename
            
            read(prmfile,*)nstep
            
            read(prmfile,*)nfreq
            
            read(prmfile,*)dt
            
            read(prmfile,*)nparticle
            
            read(prmfile,*)natom

            nmol    =   nparticle / natom
            nframe  =   nstep / nfreq
            allocate (pos(3,nparticle,0:nframe))
        
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

            nmol    =   nparticle / natom
            nframe  =   nstep / nfreq
            allocate (pos(3,nparticle,0:nframe))
        endif
    end subroutine
end module read_prm
