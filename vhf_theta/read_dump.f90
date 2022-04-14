module read_dump
    implicit none
    !double precision box_l
    !double precision, allocatable :: pos_wrap(:,:,:)
    contains
    
    subroutine readdump()
    use prm_list
    implicit none
    double precision dummy
    integer i, j

    open(16,file=filename,status='old')

    print*,'######################################'
    print*,'Reading Trajectory File...'
    
    pos(:,:,:)=0.0d0
    if (nframe .eq. 1)then
        do j = 1, skip
            read(16,'()')
        enddo
        do j = 1, nparticle
            read(16,*)dummy,pos(:,j,1)
        enddo
    else
        do i = 0, nframe
            do j = 1, skip
                read(16,'()')
            enddo
            do j = 1, nparticle
                read(16,*)dummy, pos(:,j,i)
            enddo
        enddo
    endif
    close(16)
    print*,'Finished Reading file!'
    print*,'######################################'
    end subroutine

!    subroutine readdumpwrap()
!    use prm_list
!    implicit none
!    double precision dummy
!    integer i, j
!
!    open(16,file=filename,status='old')
!
!    print*,'######################################'
!    print*,'Reading Trajectory File...'
!
!    pos(:,:,:)=0.0d0
!    !! read box size using first snapshot
!    do j = 1, 5
!        read(16,'()')
!    enddo
!    read(16,'(F18.16,5X,F18.16,1X,I3)')dummy,box_l,j
!    box_l = box_l * 10.0d0 ** j
!    do j = 1, 3
!        read(16,'()')
!    enddo
!    do i = 1, nparticle
!        read(16,*)dummy, pos(:,i,0)
!    enddo
!    !! end reading box size and first snapshot
!    !! read the rest of dump file
!    do i = 1, nframe
!        do j = 1, skip
!            read(16,'()')
!        enddo
!        do j = 1, nparticle
!            read(16,*)dummy, pos(:,j,i)
!        enddo
!    enddo
!    !! end reading dump file
!    !! Calculate wrap trajectory
!    do i = 0, nframe
!        do j = 1, nparticle
!            do k = 1, nparticle
!                if (j .ne. k)then
!                    dpos(:) = pos(:,j) -
!    !! end Calculating wrap trajectory
!    close(16)
!    print*,'Finished Reading file!'
!    print*,'######################################'
!
!    end subroutine  

end module read_dump
