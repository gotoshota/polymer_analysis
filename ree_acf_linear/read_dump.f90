module read_dump
    implicit none
    contains
    
    subroutine readdump()
    use prm_list
    implicit none
    double precision dummy
    integer i, j

    open(16,file=filename,status='old')

    print*,'######################################'
    print*,'Reading Trajectory File...'
    !! Read Box Size
    do i = 1, 5
        READ(16,'()')
    enddo
    read(16, '(E22.15,1x,e22.15)') box_l
    REWIND(16)
    do i = 0, nframe
        do j = 1, skip
            read(16,'()')
        enddo
        do j = 1, nparticle
            read(16,*)dummy, pos(:,j,i)
        enddo
    enddo
    print*,'Finished Reading file!'
    print*,'######################################'

    end subroutine
end module read_dump
