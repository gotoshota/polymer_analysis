module coarse_grain
    implicit none
    public
    integer npoint
    integer, allocatable ::Targetpoint(:)


    contains
    subroutine coarsegrain(max_val)
    use prm_list
    implicit none
    INTEGER, INTENT(IN) :: max_val
    double precision a, fn, fi
    integer DoubleNumber, delta, i, j
    integer, allocatable :: DoubleFlag(:), AllTargetpoint(:)

    print*,''
    print*,'########################################################################'
    print*,'Determine coarse-grained parameter.'
    print*,'How many points do you want?'
    read*,npoint
    print*,'########################################################################'
    
    allocate(AllTargetpoint(npoint), DoubleFlag(npoint))
    fn =dble(max_val) ** (1/dble(npoint - 1))
    
    AllTargetpoint(:) = 1
    fi = 1
    !!Calculate All Targetpoint (This can take the same value)
    do i = 1, npoint - 1
        delta = int((fn - 1)*fi)
        AllTargetpoint(i+1) = AllTargetpoint(i) + delta
        fi = fi * fn
    enddo
    
    !!Check doubled frame   
    DoubleFlag(:) = 0
    DoubleNumber = 0
    do i = 1, npoint - 1
        if (AllTargetpoint(i) == AllTargetpoint(i+1)) then
            DoubleFlag(i+1) = 1
            DoubleNumber = DoubleNumber + 1
        endif
    enddo

    !!Remove doubled frame
    allocate(Targetpoint(npoint - DoubleNumber))
    j = 1 
    do i = 1, npoint
        if (DoubleFlag(i) == 0) then
            Targetpoint(j) = AllTargetpoint(i)
            j = j + 1
        endif
    enddo
    if (j - 1 .ne. npoint - DoubleNumber) then
        print*,'Error: Net number of target point is wrong.'
        print*,'j',j
        print*,'nopoint',npoint - DoubleNumber
        stop
    elseif (Targetpoint(npoint - DoubleNumber) .gt. max_val) then
        print*,'Error: Last target frame is greater than max_val.'
        print*,'Targetpoint',Targetpoint(npoint - DoubleNumber)
        print*,'max_val',max_val
        stop
    endif

    npoint = npoint - DoubleNumber !!Net point number
    print*,''
    print*,'Coarse-grained parameter has defined.'
    print*,'Net number of points is ',npoint
    print*,'########################################################################'
    print*,''
    
    deallocate(AllTargetpoint, DoubleFlag)
    end subroutine 
end module coarse_grain
