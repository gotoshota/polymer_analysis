module coarse_grain
    implicit none
    public
    integer :: npoint = 0
    integer, allocatable ::TargetFrame(:)
    DOUBLE PRECISION, ALLOCATABLE :: times(:)


    contains
    subroutine coarsegrain
    use prm_list
    implicit none
    double precision a, fn, fi
    integer DoubleNumber, dframe, i, j
    integer, allocatable :: DoubleFlag(:), AllTargetFrame(:)

    CHARACTER(LEN=256), PARAMETER :: nmlfile = "./param.nml"
    integer   access
    INTEGER, PARAMETER :: id_nml = 11

    namelist /cg_para/ npoint

    OPEN(id_nml, file=nmlfile, status="old")
    READ(id_nml, nml=cg_para)
    CLOSE(id_nml)
    
    if (npoint .eq. 0) then
        print*,''
        print*,'########################################################################'
        print*,'Determine coarse-grained parameter.'
        print*,'How many points do you want?'
        read*,npoint
        print*,'########################################################################'
    endif
    
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
    allocate(TargetFrame(0:npoint - DoubleNumber))
    TargetFrame(0) = 0
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
    ALLOCATE(times(0:npoint))
    do i = 0, npoint
        times(i) = TargetFrame(i) *dt * nfreq
    enddo
    
    print*,''
    print*,'Coarse-grained parameter has defined.'
    print*,'Net number of points is ',npoint
    print*,'########################################################################'
    print*,''
    
    deallocate(AllTargetFrame, DoubleFlag)
    end subroutine 
end module coarse_grain
