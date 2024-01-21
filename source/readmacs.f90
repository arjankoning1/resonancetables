subroutine readmacs
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose: Read data from MACS cross section databases
!
! Revision    Date      Author      Quality  Description
! ======================================================
!    1     17-08-2023   A.J. Koning    A     Original code
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
  character(len=1)   :: isosym     ! isomeric symbol
  character(len=4)   :: year       ! year
  character(len=9)   :: sub        ! subentry
  character(len=24)  :: author     ! author
  character(len=132) :: macsfile  ! macsfile
  character(len=1200):: line       ! input line
  integer            :: Z          ! charge number
  integer            :: A          ! mass number
  integer            :: N          ! counter
  integer            :: ix
  integer            :: isoT       ! target isomer
  integer            :: isoR       ! residual isomer
  integer            :: lib        ! library number
  integer            :: type       ! reaction type
  integer            :: istat      ! error code
  real(sgl)          :: Einc       ! incident energy
  real(sgl)          :: xs         ! cross section
  real(sgl)          :: dxs        ! cross section uncertainty
  real(sgl)          :: rochread
  real(sgl)          :: CE
  real(sgl)          :: chi2
!
! **************** Read RIPL-3 database for MACS cross sections *****
!
  macs_xs = 0.
  macs_dxs = 0.
  macs_ref = ''
  macs_exist = .false.
  ratiomacs_xs = 0.
  Nmacs_xs = 0
  Emacs_xs = 0.
  Emacs_dxs = 0.
  Emacs_author = ''
  Emacs_year = ''
  Emacs_subentry = ''
  Emacs_ref = ''
  Emacs_xs = 0.
  Emacs_exist = .false.
  Lmacs_xs = 0.
  Lratiomacs_xs = 0.
!
! KADONIS database: experimental values only
!
  write(*, *) "Reading data from MACS databases....."
  lib = 1
  type = 4
  macsfile = trim(filespath)//'macs_kadonis.ng'
  open (unit = 1, status = 'old', file = macsfile)
  read(1,'()')
  read(1,'()')
  do
    read(1, '(a)', iostat = istat) line
    if (istat == -1) exit
    if (istat > 0) call read_error(macsfile, istat)
    isosym = ''
    read(line(1:9), * ) Z, A
    if (line(10:10) == 'm') isoT = 1
    read(line(78:89), * ) xs
    read(line(150:161), * ) dxs
    isoT = 0
    isoR = -1
    if (isosym == 'm') isoT = 1
    macs_xs(lib, type, Z, A, isoT, isoR) = 0.001 * xs
    macs_dxs(lib, type, Z, A, isoT, isoR) = 0.001 * dxs
    macs_ref(lib, type, Z, A, isoT, isoR) = 'Kadonis'
    macs_exist(lib, type, isoR) = .true.
  enddo
  close (1)
!
! ASTRAL database: experimental values only
!
  write(*, *) "Reading data from MACS databases....."
  lib = 1
  type = 4
  macsfile = trim(filespath)//'macs_astral.ng'
  open (unit = 1, status = 'old', file = macsfile)
  do
    read(1, '(a)', iostat = istat) line
    if (istat == -1) exit
    if (istat > 0) call read_error(macsfile, istat)
    if (line(1:1) == '#') cycle
    isosym = ''
    read(line, * ) Z, A
    ix=index(line,'_0')
    if (ix > 0) isosym = 'g'
    ix=index(line,'_1')
    if (ix > 0) isosym = 'm'
    read(line(40:80), * ) xs,dxs
    isoT = 0
    isoR = -1
    if (isosym == 'g') isoR = 0
    if (isosym == 'm') isoR = 1
    macs_xs(lib, type, Z, A, isoT, isoR) = 0.001 * xs
    macs_dxs(lib, type, Z, A, isoT, isoR) = 0.001 * dxs
    macs_ref(lib, type, Z, A, isoT, isoR) = 'Astral'
    macs_exist(lib, type, isoR) = .true.
  enddo
  close (1)
!
! Mughabghab 2016 MACS database
!
  lib = 3
  macsfile = trim(filespath)//'thermal.mugh18'
  open (unit = 1, status = 'old', file = macsfile)
  do
    read(1, '(a)', iostat = istat) line
    if (istat == -1) exit
    if (istat > 0) call read_error(macsfile, istat)
    if (line(4:4) < 'A' .or. line(4:4) > 'Z') cycle
    read(line(10:24), * ) A, Z, isoT
    if (A == 0) cycle
    isoR = -1
    type = 4
    xs = rochread(line(839:849))
    dxs = rochread(line(853:860))
    if (xs > 0.) then
      macs_xs(lib, type, Z, A, isoT, isoR) = xs
      macs_dxs(lib, type, Z, A, isoT, isoR) = dxs
      macs_ref(lib, type, Z, A, isoT, isoR) = 'Mugh18'
      macs_exist(lib, type, isoR) = .true.
    endif
  enddo
  close (1)
!
! Sukhoruchkin 2015 MACS database
!
  lib = 5
  macsfile = trim(filespath)//'sukhoruchkin.txt'
  open (unit = 1, status = 'old', file = macsfile)
  do
    read(1, '(a)', iostat = istat) line
    if (istat == -1) exit
    if (istat > 0) call read_error(macsfile, istat)
    if (line(4:4) < 'A' .or. line(4:4) > 'Z') cycle
    read(line(10:24), * ) A, Z, isoT
    if (A == 0) cycle
    isoR = -1
    type = 4
    xs = rochread(line(187:199))
    dxs = rochread(line(201:210))
    if (xs > 0.) then
      macs_xs(lib, type, Z, A, isoT, isoR) = xs
      macs_dxs(lib, type, Z, A, isoT, isoR) = dxs
      macs_ref(lib, type, Z, A, isoT, isoR) = 'Sukhoruchkin'
      macs_exist(lib, type, isoR) = .true.
    endif
  enddo
  close (1)
!
! EXFOR MACS database
!
  do type = 4,4
    macsfile = trim(filespath)//'macs_exfor.'//reac(type)
    open (unit = 1, status = 'old', file = macsfile)
    do
      read(1, '(a)', iostat = istat) line
      if (istat == -1) exit
      if (istat > 0) call read_error(macsfile, istat)
      read(line, * ) Z, A, isoT, isoR, xs ,dxs, Einc, author, year, sub
      Nmacs_xs(type, Z, A, isoT, isoR) = Nmacs_xs(type, Z, A, isoT, isoR) + 1
      N = Nmacs_xs(type, Z, A, isoT, isoR)
      Emacs_xs(type, Z, A, isoT, isoR, N) = 0.001*xs
      Emacs_dxs(type, Z, A, isoT, isoR, N) = 0.001*dxs
      Emacs_author(type, Z, A, isoT, isoR, N) = author
      Emacs_year(type, Z, A, isoT, isoR, N) = year
      Emacs_subentry(type, Z, A, isoT, isoR, N) = sub
      Emacs_ref(type, Z, A, isoT, isoR, N) = trim(author)//'_'//trim(year)//'_'//trim(sub)
      Emacs_exist(type, isoR) = .true.
    enddo
    close (1)
  enddo
!
! MACS from Nuclear data libraries
!
  do type = 4,4
    do lib = 1, numndlib
      macsfile = trim(libspath)//trim(ndlib(lib))//'.MACS'
      open (unit = 1, status = 'old', file = macsfile)
      do
        read(1, '(a)', iostat = istat) line
        if (istat == -1) exit
        if (istat > 0) call read_error(macsfile, istat)
        if (line(1:1) == '#') cycle
        read(line, * ) Z, A, isoT, CE, chi2, xs
        Lmacs_xs(lib, type, Z, A, isoT, isoR) = xs
      enddo
      close (1)
    enddo
  enddo
  return
end subroutine readmacs
! Copyright A.J. Koning 2019
