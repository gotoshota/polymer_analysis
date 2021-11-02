module read_prm 
    implicit none
    contains
    subroutine readprm()
        use prm_list
        implicit none 
        
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
    end subroutine
end module read_prm
