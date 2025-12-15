subroutine sumresonance(type)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose: Summarize resonance data
!
! Revision    Date      Author      Quality  Description
! ======================================================
!    1     2025-08-09   A.J. Koning    A     Original code
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
  integer, parameter :: Nsource=12
  character(len=200) :: line
  character(len=200) :: fline(10000)
  character(len=132) :: resfile   ! nuclide file
  character(len=20)  :: msource(Nsource)
  character(len=132) :: sourcefile(Nsource)
  logical            :: lexist
  integer            :: istat
  integer            :: isource
  integer            :: ifile
  integer            :: ix
  integer            :: indent
  real               :: xs
  real               :: dxs
  real               :: ratio
  character(len=2)   :: iso        ! extension
  character(len=3)   :: Astring    ! mass string
  character(len=6)   :: nuclide
  character(len=132) :: nucfile    ! nuclide file
  character(len=18)  :: rfile
  character(len=15)  :: ref
  character(len=30)  :: auth
  character(len=18)  :: react      ! reaction
  character(len=15)  :: col(12)     ! header
  character(len=15)  :: un(12)      ! units
  character(len=80)  :: quantity   ! quantity
  character(len=132) :: topline    ! topline
  integer            :: k          ! counter
  integer            :: N          ! counter
  integer            :: Ncol       ! number of columns
  integer            :: type       ! reaction type
!
! **************** Write databases for thermal cross sections *****
!
  indent = 0
  react=restype(type)
  msource(1) = 'RIPL-2'
  msource(2) = 'RIPL-3'
  msource(3) = 'JUKO'
  msource(4) = 'Mughabghab-2006'
  msource(5) = 'Mughabghab-2018'
  msource(6) = 'cendl3.2'
  msource(7) = 'jendl5.0'
  msource(8) = 'tendl.2025'
  msource(9) = 'endfb8.1'
  msource(10) = 'jeff4.0'
  msource(11) = 'EXFOR'
  msource(12) = 'TARES'
  do isource = 1, Nsource
    sourcefile(isource)=trim(respath)//trim(react)//'/all/'//trim(msource(isource))//'_'//trim(react)//'.txt'
    ifile = 10 + isource
    open (unit = ifile, status = 'unknown', file = trim(sourcefile(isource)))
  enddo
  quantity='resonance data'
  topline=trim(react)//' '//trim(quantity)
  rfile=trim(reac(type))//iso
  nucfile=trim(respath)//trim(react)//'/all/selected_'//trim(react)//'.txt'
  write(*,*) " Writing to ", trim(nucfile) 
  open (unit = 1, status = 'unknown', file = trim(nucfile))
  call write_header(indent,topline,source,user,date,oformat)
  call write_reaction(indent,react,0.D0,0.D0,0,0)
  col(1) = 'Z'
  col(2) = 'A'
  col(3) = 'Liso'
  col(4) = 'Value'
  col(5) = 'dValue'
  col(6) = 'Reference'
  col(7) = 'Rel. dev. comp.'
  col(8) = 'Rel. dev. NDL'
  col(9) = 'Rel. dev. EXFOR'
  col(10) = 'Rel. dev. all'
  col(11) = '#Experiments'
  col(12) = 'Nuclide'
  Ncol = 12
  un = ''
  if (type == 2 .or. type == 5) then
    un(4) = '*e-4'
    un(5) = '*e-4'
  else
    un(4) = 'eV'
    un(5) = 'eV'
  endif
  if (type == 8) then
    un(5) = 'fm'
    un(6) = 'fm'
  endif
  un(7) = '%'
  un(8) = '%'
  un(9) = '%'
  un(10) = '%'
  N = Nsave
  call write_quantity(indent,quantity)
  call write_datablock(indent,Ncol,N,col,un)
  do k = 1, N
    Astring='   '
    write(Astring(1:3),'(i3.3)') Asave(k)
    nuclide=trim(nuc(Zsave(k)))//Astring
    if (Lisosave(k) == 1) nuclide = trim(nuclide)//'m'
    if (Lisosave(k) == 2) nuclide = trim(nuclide)//'n'
    write(1, '(3(6x,i4,5x),2es15.6,2x,a15,4(4x,f7.2,4x),4x,i4,11x,a6)') Zsave(k), Asave(k), Lisosave(k), &
 &     xssave(k), dxssave(k), refsave(k), compsave(k), NDLsave(k), expsave(k), varsave(k), Nexpsave(k), nuclide
    resfile=trim(respath)//trim(react)//'/nuc/'//trim(nuclide)//'_'//trim(react)//'.txt'
    inquire (file = resfile, exist = lexist)
    if (lexist) then
      do isource = 1, Nsource
        ifile = 10 + isource
        open (unit = 2, status = 'unknown', file = trim(resfile))
        do
          read(2, '(a)', iostat = istat) line
          if (istat == -1) exit
          if (istat > 0) call read_error(resfile, istat)
          ix = index(line(1:45),trim(msource(isource)))
          if (ix > 0 .and. line(1:1) /= '#') then
            read(line(61:75), *) xs
            read(line(76:90), *) dxs
            if (xssave(k) > 0.) then
              ratio = xs / xssave(k)
            else
              ratio = 0.
            endif
            if (isource == 11) then
              read(line(1:30), '(a)') auth
              read(line(91:105), '(a)') ref
              write(ifile, '(3(6x,i4,5x),2es15.6,f15.4,a15,a30,6x,a6)') Zsave(k), Asave(k), Lisosave(k), xs, dxs, ratio, &
 &              ref, auth, nuclide
            else
              write(ifile, '(3(6x,i4,5x),2es15.6,f15.4,6x,a6)') Zsave(k), Asave(k), Lisosave(k), xs, dxs, ratio, nuclide
            endif
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
  col(6) = 'Ratio'
  fline = ''
  do isource = 1, Nsource
    sourcefile(isource)=trim(respath)//trim(react)//'/all/'//trim(msource(isource))//'_'//trim(react)//'.txt'
    open (unit = 10, status = 'unknown', file = trim(sourcefile(isource)))
    k = 1
    do
      read(10, '(a)', iostat = istat) fline(k)
      if (istat == -1) exit
      k = k + 1
    enddo
    N = k - 1
    close(10)
    if (isource == 11) then
      col(7) = 'Reference'
      col(8) = 'Author'
      col(9) = ''
      col(10) = 'Nuclide'
      Ncol = 10
    else
      col(7) = 'Nuclide'
      Ncol = 7
    endif
    open (unit = 1, status = 'unknown', file = trim(sourcefile(isource)))
    if (N > 0) then
      call write_header(indent,trim(topline)//' for '//trim(msource(isource)),source,user,date,oformat)
      call write_reaction(indent,react,0.D0,0.D0,0,0)
      call write_quantity(indent,quantity)
      call write_datablock(indent,Ncol,N,col,un)
      do k = 1, N
        write(1, '(a)') trim(fline(k))
      enddo
      close(1)
    else
      close(1, status = 'delete')
    endif
  enddo
  return
end subroutine sumresonance
! Copyright A.J. Koning 2025
