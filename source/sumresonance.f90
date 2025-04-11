subroutine sumresonance(type)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose: Summarize resonance data
!
! Revision    Date      Author      Quality  Description
! ======================================================
!    1     2025-02-08   A.J. Koning    A     Original code
!-----------------------------------------------------------------------------------------------------------------------------------
!
! *** Use data from other modules
!
  use A0_resonancetables_mod
  use A1_error_handling_mod
!
! *** Declaration of local data
!
  implicit none
  integer, parameter :: Nsource=10
  character(len=132) :: line
  character(len=132) :: resfile   ! nuclide file
  character(len=20)  :: msource(Nsource)
  character(len=132) :: sourcefile(Nsource)
  logical            :: lexist
  integer            :: istat
  integer            :: isource
  integer            :: ifile
  integer            :: ix
  real               :: xs
  real               :: dxs
  real               :: ratio
  character(len=2)   :: iso        ! extension
  character(len=3)   :: Astring    ! mass string
  character(len=6)   :: nuclide
  character(len=132) :: nucfile    ! nuclide file
  character(len=18)  :: rfile
  character(len=18)  :: react      ! reaction
  character(len=15)  :: col(8)     ! header
  character(len=15)  :: un(8)      ! units
  character(len=80)  :: quantity   ! quantity
  character(len=132) :: topline    ! topline
  integer            :: k          ! counter
  integer            :: N          ! counter
  integer            :: Ncol       ! number of columns
  integer            :: type       ! reaction type
!
! **************** Write databases for thermal cross sections *****
!
  react=restype(type)
  msource(1) = 'RIPL-2'
  msource(2) = 'RIPL-3'
  msource(3) = 'JUKO'
  msource(4) = 'Mughabghab_2006'
  msource(5) = 'Mughabghab_2016'
  msource(6) = 'cendl3.2'
  msource(7) = 'jendl5.0'
  msource(8) = 'tendl.2023'
  msource(9) = 'endfb8.1'
  msource(10) = 'jeff4.0'
  do isource = 1, Nsource
    sourcefile(isource)=trim(respath)//trim(react)//'/'//trim(msource(isource))//'.res'
    ifile = 10 + isource
    open (unit = ifile, status = 'unknown', file = trim(sourcefile(isource)))
    write(ifile,'("##       Z              A            Liso           Value         dValue          Ratio")')
    write(ifile,'("##      []             []             []             [b]            [b]            []")')
  enddo
  quantity='resonance data'
  topline=trim(react)//' '//trim(quantity)
  rfile=trim(reac(type))//iso
  nucfile=trim(respath)//trim(react)//'/'//trim(react)//'.res'
  write(*,*) " Writing to ", trim(nucfile) 
  open (unit = 1, status = 'unknown', file = trim(nucfile))
  call write_header(topline,source,user,date,oformat)
  call write_reaction(react,0.D0,0.D0,0,0)
  col(1) = 'Z'
  col(2) = 'A'
  col(3) = 'Liso'
  col(4) = 'Value'
  col(5) = 'dValue'
  col(6) = 'Reference'
  col(7) = '#Experiments'
  col(8) = 'Nuclide'
  Ncol = 8
  un = ''
  un(4) = 'eV'
  un(5) = 'eV'
  N = Nsave
  call write_quantity(quantity)
  call write_datablock(Ncol,N,col,un)
  do k = 1, N
    Astring='   '
    write(Astring(1:3),'(i3.3)') Asave(k)
    nuclide=trim(nuc(Zsave(k)))//Astring
    if (Lisosave(k) == 1) nuclide = trim(nuclide)//'m'
    if (Lisosave(k) == 2) nuclide = trim(nuclide)//'n'
    write(1, '(3(6x,i4,5x),2es15.6,2x,a15,4x,i4,11x,a6)') Zsave(k), Asave(k), Lisosave(k), xssave(k), dxssave(k), refsave(k), &
 &     Nexpsave(k), nuclide
    resfile=trim(respath)//trim(react)//'/'//trim(nuclide)//'.res'
    inquire (file = resfile, exist = lexist)
    if (lexist) then
      do isource = 1, Nsource
        ifile = 10 + isource
        open (unit = 2, status = 'unknown', file = trim(resfile))
        do
          read(2, '(a)', iostat = istat) line
          if (istat == -1) exit
          if (istat > 0) call read_error(resfile, istat)
          ix = index(line(1:25),trim(msource(isource)))
          if (ix > 0) then
            read(line(61:75), *) xs
            read(line(76:90), *) dxs
            if (xssave(k) > 0.) then
              ratio = xs / xssave(k)
            else
              ratio = 0.
            endif
            write(ifile, '(3(6x,i4,5x),3es15.6)') Zsave(k), Asave(k), Lisosave(k), xs, dxs, ratio
            exit
          endif
        enddo
        close(2)
      enddo
    endif
  enddo
  close(1)
  do isource = 1, Nsource
    close (unit = 10 + isource)
  enddo
  return
end subroutine sumresonance
! Copyright A.J. Koning 2025
