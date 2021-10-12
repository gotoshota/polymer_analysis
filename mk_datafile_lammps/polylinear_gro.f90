Program ini_str
    implicit none
    integer     n, i
    double precision b_r

    print*,'Input Number of atoms'
    read*,n

    open(17,file='ini_str.gro',status='replace')
    write(17,'(A16,/)')'LAMMPS data file'
    write(17,'(I5,1X,A5)')n,'atoms'
    write(17,'(I5,1X,A5)')n-1,'bonds'
    write(17,'(I5,1X,A6,/)')n-2,'angles'    
    write(17,'(I1,1X,A10)')1,'atom types'
    write(17,'(I1,1X,A10)')1,'bond types'
    write(17,'(I1,1X,A11,/)')1,'angle types'
    b_r = n*0.96 + 1
    
    write(17,'(F4.2,1X,F6.2,1X,A3,1X,A3)')0.0,b_r,'xlo','xhi'
    write(17,'(F4.2,1X,F6.2,1X,A3,1X,A3)')0.0,b_r,'ylo','yhi'
    write(17,'(F4.2,1X,F6.2,1X,A3,1X,A3,/)')0.0,b_r,'zlo','zhi'
   
    write(17,'(A6,/)')'Masses'
    write(17,'(I1,1X,F3.1,/)')1,1.0

    write(17,'(A5,/)')'Atoms'
    do i = 1,n
        write(17,'(I3,1X,I1,1X,I1,1X,F3.1,1X,F3.1,1X,F6.2)')i,1,1,0.0,0.0,0.96*(i-1)
    enddo
    
    write(17,'(/,A5,/)')'Bonds'
    do i = 1,n-1
        write(17,'(I3,1X,I1,1X,I3,1X,I3)')i,1,i,i+1
    enddo

    write(17,'(/,A6,/)')'Angles'
    do i = 2, n-1
        write(17,'(I3,1X,I1,1X,I3,1X,I3,1X,I3)')i,1,i-1,i,i+1
    enddo

    close(17)
end Program ini_str
