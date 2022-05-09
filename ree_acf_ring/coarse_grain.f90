module coarse_grain
    implicit none
    public
    integer npoint
    integer, allocatable ::TargetFrame(:)


    contains
    subroutine coarsegrain
    use prm_list
    implicit none
    double precision a, fn, fi
    integer DoubleNumber, dframe, i, j
    integer, allocatable :: DoubleFlag(:), AllTargetFrame(:)

    print*,''
    print*,'########################################################################'
    print*,'Determine coarse-grained parameter.'
    print*,'How many points do you want?'
    read*,npoint
    print*,'########################################################################'
    
    allocate(AllTargetFrame(npoint), DoubleFlag(npoint))
    fn =dble(nframe) ** (1/dble(npoint - 1))
    
    AllTargetFrame(:) = 1
    fi = 1
    !!Calculate All Target Frame (This can take the same value)
    do i = 1, npoint - 1
        dframe = int((fn - 1)*fi)
        AllTargetFrame(i+1) = AllTargetFrame(i) + dframe
        fi = fi * fn
    enddo
    
    !!Check doubled frame   
    DoubleFlag(:) = 0
    DoubleNumber = 0
    do i = 1, npoint - 1
        if (AllTargetFrame(i) == AllTargetFrame(i+1)) then
            DoubleFlag(i+1) = 1
            DoubleNumber = DoubleNumber + 1
        endif
    enddo

    !!Remove doubled frame
    allocate(TargetFrame(npoint - DoubleNumber))
    j = 1 
    do i = 1, npoint
        if (DoubleFlag(i) == 0) then
            TargetFrame(j) = AllTargetFrame(i)
            j = j + 1
        endif
    enddo
    if (j - 1 .ne. npoint - DoubleNumber) then
        print*,'Error: Net number of target frame is wrong.'
        print*,'j',j
        print*,'nopoint',npoint - DoubleNumber
        stop
    elseif (TargetFrame(npoint - DoubleNumber) .gt. nframe) then
        print*,'Error: Last target frame is greater than nframe.'
        print*,'TargetFrame',TargetFrame(npoint - DoubleNumber)
        print*,'nframe',nframe
        stop
    endif

    npoint = npoint - DoubleNumber !!Net point number
    print*,''
    print*,'Coarse-grained parameter has defined.'
    print*,'Net number of points is ',npoint
    print*,'########################################################################'
    print*,''
    print*,'check'
    
    deallocate(AllTargetFrame, DoubleFlag)
    print*,'check'
    end subroutine 
end module coarse_grain
