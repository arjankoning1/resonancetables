subroutine sumthermal(Riso, type)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose: Summarize thermal cross section data
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
  character(len=2)   :: iso        ! extension
  character(len=132) :: thermfile   ! nuclide file
  character(len=20)  :: msource(Nsource)
  character(len=200) :: line
  character(len=200) :: fline(10000)
  character(len=132) :: sourcefile(Nsource)
  character(len=3)   :: Astring    ! mass string
  character(len=6)   :: nuclide
  character(len=132) :: nucfile    ! nuclide file
  character(len=18)  :: rfile
  character(len=15)  :: ref
  character(len=30)  :: auth
  character(len=18)  :: react      ! reaction
  character(len=15)  :: col(13)     ! header
  character(len=15)  :: un(13)      ! units
  character(len=80)  :: quantity   ! quantity
  character(len=132) :: topline    ! topline
  logical            :: lexist
  integer            :: istat
  integer            :: indent
  integer            :: k          ! counter
  integer            :: N          ! counter
  integer            :: Ncol       ! number of columns
  integer            :: Riso       ! residual isomer
  integer            :: type       ! reaction type
  integer            :: isource
  integer            :: ifile
  integer            :: ix
  real               :: xs
  real               :: dxs
  real               :: ratio
!
! **************** Write databases for thermal cross sections *****
!
  indent = 0
  iso = ''
  if (Riso == 0) iso='-g'
  if (Riso == 1) iso='-m'
  if (Riso == 2) iso='-n'
  if (type <= 6) then
    quantity='thermal cross section'
  else
    quantity='thermal neutron multiplicity'
  endif
  react=trim(reaction(type))//iso
  topline=trim(react)//' '//trim(quantity)
  rfile=trim(reac(type))//trim(iso)
  msource(1) = 'RIPL-3'
  msource(2) = 'Kayzero'
  msource(3) = 'Sukhoruchkin'
  msource(4) = 'Mughabghab-2006'
  msource(5) = 'Mughabghab-2018'
  msource(6) = 'Firestone-2022'
  msource(7) = 'cendl3.2'
  msource(8) = 'jendl5.0'
  msource(9) = 'tendl.2025'
  msource(10) = 'endfb8.1'
  msource(11) = 'jeff4.0'
  msource(12) = 'EXFOR'
  do isource = 1, Nsource
    sourcefile(isource)=trim(thermalpath)//trim(rfile)//'/all/'//trim(msource(isource))//'_'//trim(rfile)//'.txt'
    ifile = 10 + isource
    open (unit = ifile, status = 'unknown', file = trim(sourcefile(isource)))
  enddo
  nucfile=trim(thermalpath)//trim(rfile)//'/all/selected_'//trim(rfile)//'.txt'
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
  col(12) = 'Spectrum'
  col(13) = 'Nuclide'
  Ncol = 13
  un = ''
  if (type <= 6) then
    un(4) = 'b'
    un(5) = 'b'
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
    write(1, '(3(6x,i4,5x),2es15.6,2x,a15,4(4x,f7.2,4x),4x,i4,11x,a9,10x,a6)') Zsave(k), Asave(k), Lisosave(k), &
 &     xssave(k), dxssave(k), refsave(k), compsave(k), NDLsave(k), expsave(k), varsave(k), Nexpsave(k), avsave(k), nuclide 
    thermfile=trim(thermalpath)//trim(rfile)//'/nuc/'//trim(nuclide)//'_'//trim(rfile)//'.txt'
    inquire (file = thermfile, exist = lexist)
    if (lexist) then
      do isource = 1, Nsource
        ifile = 10 + isource
        open (unit = 2, status = 'unknown', file = trim(thermfile))
        do
          read(2, '(a)', iostat = istat) line
          if (istat == -1) exit
          if (istat > 0) call read_error(thermfile, istat)
          ix = index(line(1:45),trim(msource(isource)))
          if (ix > 0 .and. line(1:1) /= '#') then
            read(line(61:75), *) xs
            read(line(76:90), *) dxs
            if (xssave(k) > 0.) then
              ratio = xs / xssave(k)
            else
              ratio = 0.
            endif
            if (isource == Nsource) then
              read(line(1:30), '(a)') auth
              read(line(91:105), '(a)') ref
              write(ifile, '(3(6x,i4,5x),2es15.6,f15.4,a15,a30,6x,a6)') Zsave(k), Asave(k), Lisosave(k), xs, dxs, ratio, &
 &              ref, auth, nuclide
            else
              write(ifile, '(3(6x,i4,5x),2es15.6,f16.4,6x,a6)') Zsave(k), Asave(k), Lisosave(k), xs, dxs, ratio, nuclide
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
    sourcefile(isource)=trim(thermalpath)//trim(rfile)//'/all/'//trim(msource(isource))//'_'//trim(rfile)//'.txt'
    open (unit = 10, status = 'unknown', file = trim(sourcefile(isource)))
    k = 1
    do
      read(10, '(a)', iostat = istat) fline(k)
      if (istat == -1) exit
      k = k + 1
    enddo
    N = k - 1
    close(10)
    if (isource == Nsource) then
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
end subroutine sumthermal
! Copyright A.J. Koning 2025
