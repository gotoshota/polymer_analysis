Program ini_str
    implicit none
    integer     n, i
    double precision b_r,r
    double precision,parameter :: pi = acos(-1.0d0)

    print*,'Input Number of atoms'
    read*,n

    open(17,file='ini_str.gro',status='replace')
    write(17,'(A16,/)')'LAMMPS data file'
    write(17,'(I6,1X,A5)')n,'atoms'
    write(17,'(I6,1X,A5)')n,'bonds'
    write(17,'(I6,1X,A6,/)')n,'angles'
    write(17,'(I3,1X,A10)')1,'atom types'
    write(17,'(I3,1X,A10)')1,'bond types'
    write(17,'(I3,1X,A11,/)')1,'angle types'
    r = n*0.96/pi/2
    b_r = r+1
    
    write(17,'(F9.3,1X,F9.3,1X,A3,1X,A3)')-b_r,b_r,'xlo','xhi'
    write(17,'(F9.3,1X,F9.3,1X,A3,1X,A3)')-b_r,b_r,'ylo','yhi'
    write(17,'(F9.3,1X,F9.3,1X,A3,1X,A3,/)')-b_r,b_r,'zlo','zhi'
   
    write(17,'(A6,/)')'Masses'
    write(17,'(I1,1X,F3.1,/)')1,1.0

    write(17,'(A5,/)')'Atoms'
    do i = 1,n
        write(17,'(I6,1X,I1,1X,I1,1X,F9.3,1X,F9.3,1X,F9.3)')i,1,1,r*cos(2*pi*i/n),r*sin(2*pi*i/n),0.0
    enddo
    
    write(17,'(/,A5,/)')'Bonds'
    do i = 1,n-1
        write(17,'(I6,1X,I1,1X,I6,1X,I6)')i,1,i,i+1
    enddo
    write(17,'(I6,1X,I1,1X,I6,1X,I6)')n,1,n,1
    
    write(17,'(/,A6,/)')'Angles'
    write(17,'(I1,1X,I1,1X,I6,1X,I3,1X,I3)')1,1,n,1,2
    do i = 2, n-1
        write(17,'(I6,1X,I1,1X,I6,1X,I6,1X,I6)')i,1,i-1,i,i+1
    enddo
    write(17,'(I6,1X,I1,1X,I6,1X,I6,1X,I3)')n,1,n-1,n,1
    close(17)
end Program ini_str
