module read_dump
    use prm_list
    implicit none

    INTEGER(KIND=4) :: dumpfile = 666
    
    contains
    subroutine readdump()
        CHARACTER(LEN=256) :: ext
        INTEGER(KIND=4) period

        print*,'######################################'
        print*,'Reading Trajectory File...'
        ! get ext
        period = INDEX(filename, ".", BACK=.true.)
        ext = filename(period+1 : LEN(TRIM(filename))) 

        if ( ext .eq. "lammpstrj" ) then
            print *, "The extension of dumpfile is 'lammpstrj'."
            CALL readdump_txt()
        else if ( ext .eq. "dcd" ) then
            print *, "The extension of dumpfile is 'dcd'."
            CALL readdump_dcd()
        else if ( ext .eq. "xtc" ) then
            print *, "The extension of dumpfile is 'xtc'."
            CALL readdump_xtc()
        else if ( ext .eq. "mydumpbin" ) then
            print *, "The extension of dumpfile is 'mydumpbin'."
            CALL readdump_mydumpbin()
        else
            print *,"dumpfile is NOT found."
            stop
        end if
        print*,'Finished Reading file!'
        print*,'######################################'
    end subroutine
   
    subroutine readdump_txt()
        implicit none
        double precision :: dummy
        DOUBLE PRECISION :: a, b
        integer i, j

        open(dumpfile,file=filename,status='old')

        !! Read Box Size
        do i = 1, 5
            READ(dumpfile,'()')
        enddo
        read(dumpfile, '(E22.15,1x,e22.15)') a, b
        box_l = b - a
        REWIND(dumpfile)
        do i = 0, nframe
            do j = 1, skip
                read(dumpfile,'()')
            enddo
            do j = 1, nparticle
                read(dumpfile,*)dummy, pos(:,j,i)
            enddo
        enddo
    end subroutine

    subroutine readdump_dcd()
        implicit none
        integer              :: i, natom
        integer(4)           :: dcdinfo(20), ntitle
        real(4), allocatable :: coord(:,:)
        real(8)              :: boxsize(1:6)
        character(len=4)     :: header
        character(len=80)    :: title(10)

        open(dumpfile, file=filename, form='unformatted', status='old')
        read(dumpfile) header, dcdinfo(1:20)
        read(dumpfile) ntitle, (title(i),i=1,ntitle)
        read(dumpfile) natom

        allocate(coord(1:3,1:natom))

        do i = 0, dcdinfo(1)
        read(dumpfile) boxsize(1:6)
        read(dumpfile) pos(1,1:natom,i)
        read(dumpfile) pos(2,1:natom,i)
        read(dumpfile) pos(3,1:natom,i)
        end do
        close(dumpfile)
    end subroutine

    subroutine readdump_xtc()
        use xdr, only : xtcfile

        implicit none

        INTEGER(KIND=4) i, j

        TYPE(xtcfile) :: xtcf

        call xtcf % init(filename)

        call xtcf % read
        do i = 0, nframe
            do j = 1, nparticle
                !pos(1,j,i) = xtcf % pos(1,j)
                !pos(2,j,i) = xtcf % pos(2,j)
                !pos(3,j,i) = xtcf % pos(3,j)
                pos(:,j,i) = xtcf % pos(:,j)
            enddo
            call xtcf % read
        end do
        call xtcf % close

    end subroutine

    subroutine readdump_mydumpbin()
        implicit none
        INTEGER(KIND=4) i, j
        REAL(KIND=4) x, y, z

        open(dumpfile, file=filename, form='unformatted', status='old')

        do i = 0, nframe
            READ(dumpfile) box_l
            do j = 1, nparticle
                READ(dumpfile) x, y, z
                pos(1,j,i) = x
                pos(2,j,i) = y 
                pos(3,j,i) = z
            enddo
        enddo
    end subroutine
end module read_dump
