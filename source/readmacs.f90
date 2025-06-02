subroutine readmacs(Z, A, Liso, Riso)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose: Read data from MACS cross section databases
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
  character(len=1)   :: isosym     ! isomeric symbol
  character(len=2)   :: Unuc
  character(len=3)   :: exten
  character(len=9)   :: sub        ! subentry
  character(len=9)   :: ref        ! reference
  character(len=24)  :: author     ! author
  character(len=132) :: macsfile  ! macsfile
  character(len=1200):: line       ! input line
  integer            :: iz         ! charge number
  integer            :: ia         ! mass number
  integer            :: Z          ! charge number
  integer            :: year       ! year
  integer            :: A          ! mass number
  integer            :: N          ! counter
  integer            :: ix
  integer            :: iav
  integer            :: i
  integer            :: k
  integer            :: isoT       ! target isomer
  integer            :: Liso       ! target isomer
  integer            :: Riso       ! residual isomer
  integer            :: isoR       ! residual isomer
  integer            :: lib        ! library number
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
  res_author = ''
  res_type = ''
  res_year = 0
  res_ref = ''
  res_xs = 0.
  res_dxs = 0.
  res_E = 0.03
  res_av = ''
  ref = ''
  Nres_exp = 0
  Nres = 0
  k = 0
!
! KADONIS database: experimental values only
!
  xs = 0.
  dxs = 0.
  macsfile = trim(filespath)//'macs_kadonis.ng'
  open (unit = 2, status = 'old', file = macsfile)
  read(2,'()')
  read(2,'()')
  do
    read(2, '(a)', iostat = istat) line
    if (istat == -1) exit
    if (istat > 0) call read_error(macsfile, istat)
    isoT = 0
    read(line(1:9), * ) iz, ia
    if (line(10:10) == 'm') isoT = 1
    if (iz == Z .and. ia == A .and. Liso == isoT .and. Riso == -1) then
      read(line(78:89), * ) xs
      read(line(150:161), * ) dxs
      if (xs > 0.) then
        k = k + 1
        res_author(k) = 'Kadonis'
        res_type(k) = 'Compilation'
        res_year(k) = 2000 
        res_ref(k) = ref
        res_av(k) = 'MXW'
        res_xs(k) = 0.001 * xs
        res_dxs(k) = 0.001 * dxs
        res_exist = .true.
      endif
      exit
    endif
  enddo
  close (2)
!
! ASTRAL database: experimental values only
!
  xs = 0.
  dxs = 0.
  macsfile = trim(filespath)//'macs_astral.ng'
  open (unit = 2, status = 'old', file = macsfile)
  do
    read(2, '(a)', iostat = istat) line
    if (istat == -1) exit
    if (istat > 0) call read_error(macsfile, istat)
    if (line(1:1) == '#') cycle
    isosym = ''
    read(line, * ) iz, ia
    ix=index(line,'_0')
    if (ix > 0) isosym = 'g'
    ix=index(line,'_1')
    if (ix > 0) isosym = 'm'
    if (iz == Z .and. ia == A .and. Liso == 0) then
      read(line(40:80), * ) xs,dxs
      isoR = -1
      if (isosym == 'g') isoR = 0
      if (isosym == 'm') isoR = 1
      if (xs > 0.) then
        k = k + 1
        res_author(k) = 'Astral'
        res_type(k) = 'Compilation'
        res_year(k) = 2020
        res_ref(k) = ref
        res_av(k) = 'MXW'
        res_xs(k) = 0.001 * xs
        res_dxs(k) = 0.001 * dxs
        res_exist = .true.
      endif
      exit
    endif
  enddo
  close (2)
!
! Bao database: experimental values only
!
  xs = 0.
  dxs = 0.
  macsfile = trim(filespath)//'macs_bao.ng'
  open (unit = 2, status = 'old', file = macsfile)
  read(2,'(////////////////)')
  do
    read(2, '(a)', iostat = istat) line
    if (istat == -1) exit
    if (istat > 0) call read_error(macsfile, istat)
    isoT = 0
    read(line(1:4), * ) ia
    read(line(5:6), '(a2)') Unuc
    if (Unuc(2:2) /= ' ') Unuc(2:2)=achar(iachar(Unuc(2:2)) + 32)
    do i=1,numZ
      if (trim(Unuc) == trim(nuc(i))) then
        iz = i
        exit
      endif
    enddo
    if (iz == Z .and. ia == A .and. Liso == isoT .and. Riso == -1) then
      read(line(62:68), * ) xs
      read(line(120:126), * ) dxs
      if (xs > 0.) then
        k = k + 1
        res_author(k) = 'Bao'
        res_type(k) = 'Compilation'
        res_year(k) = 2000 
        res_ref(k) = ref
        res_av(k) = 'MXW'
        res_xs(k) = 0.001 * xs
        res_dxs(k) = 0.001 * dxs
        res_exist = .true.
      endif
      exit
    endif
  enddo
  close (2)
!
! Mughabghab 2018 MACS database
!
  xs = 0.
  dxs = 0.
  macsfile = trim(filespath)//'global.2023.txt'
  open (unit = 2, status = 'old', file = macsfile)
  do
    read(2, '(a)', iostat = istat) line
    if (istat == -1) exit
    if (istat > 0) call read_error(macsfile, istat)
    if (line(4:4) < 'A' .or. line(4:4) > 'Z') cycle
    read(line(10:24), * ) ia, iz, isoT
    if (ia == 0) cycle
    isoR = -1
    if (iz == Z .and. ia == A .and. Liso == isoT .and. Riso == isoR) then
      xs = rochread(line(839:849))
      dxs = rochread(line(853:860))
      if (xs > 0.) then
        k = k + 1
        res_author(k) = 'Mughabghab-2018'
        res_type(k) = 'Compilation'
        res_year(k) = 2018
        res_ref(k) = ref
        res_xs(k) = xs
        res_dxs(k) = dxs
        res_av(k) = 'MXW'
        res_exist = .true.
      endif
      exit
    endif
  enddo
  close (2)
!
! Sukhoruchkin 2015 MACS database
!
  macsfile = trim(filespath)//'sukhoruchkin.txt'
  open (unit = 2, status = 'old', file = macsfile)
  do
    read(2, '(a)', iostat = istat) line
    if (istat == -1) exit
    if (istat > 0) call read_error(macsfile, istat)
    if (line(4:4) < 'A' .or. line(4:4) > 'Z') cycle
    read(line(10:24), * ) ia, iz, isoT
    if (ia == 0) cycle
    isoR = -1
    if (iz == Z .and. ia == A .and. Liso == isoT .and. Riso == isoR) then
      xs = rochread(line(187:199))
      dxs = rochread(line(201:210))
      if (xs > 0.) then
        k = k + 1
        res_author(k) = 'Sukhoruchkin'
        res_type(k) = 'Compilation'
        res_year(k) = 2015
        res_ref(k) = ref
        res_xs(k) =  xs
        res_dxs(k) = dxs
        res_av(k) = 'MXW'
        res_exist = .true.
      endif
      exit
    endif
  enddo
  close (2)
!
! MACS from Nuclear data libraries
!
  xs = 0.
  dxs = 0.
  do lib = 1, numndlib
    macsfile = trim(libspath)//trim(ndlib(lib))//'.MACS'
    open (unit = 2, status = 'old', file = macsfile)
    do
      read(2, '(a)', iostat = istat) line
      if (istat == -1) exit
      if (istat > 0) call read_error(macsfile, istat)
      if (line(1:1) == '#') cycle
      read(line, * ) iz, ia, isoT, CE, chi2, xs
      if (iz == Z .and. ia == A .and. Liso == isoT .and. Riso == -1) then
        if (xs > 0.) then
          k = k + 1
          res_author(k) = ndlib(lib)
          res_type(k) = 'NDL'
          res_year(k) = ndyear(lib)
          res_ref(k) = ref
          res_xs(k) = xs
          res_dxs(k) = dxs
          res_av(k) = 'MXW'
          res_exist = .true.
        endif
        exit
      endif
    enddo
    close (2)
  enddo
!
! EXFOR MACS database
!
  do iav= 0, 1
    if (iav == 1) then   
      exten='_av'      
    else
      exten=''
    endif
    xs = 0.
    dxs = 0.
    macsfile = trim(exforpath)//'exfor_30keV'//trim(exten)//'.txt'
    open (unit = 2, status = 'old', file = macsfile)
    do
      read(2, '(a)', iostat = istat) line
      if (istat == -1) exit
      if (istat > 0) call read_error(macsfile, istat)
      read(line, * ) iz, ia, isoT,isoR, xs ,dxs, Einc, author, year, sub
      if (iz == Z .and. ia == A .and. Liso == isoT .and. Riso == isoR) then
        Nres_exp = Nres_exp + 1
        N = Nres_exp
        if (N > numex) then
          Nres_exp = numex
          exit
        endif
        if (xs > 0.) then
          k = k + 1
          res_author(k) = author
          res_type(k) = 'EXFOR'
          res_year(k) = year
          res_ref(k) = sub 
          res_E(k) = Einc 
          res_xs(k) = 0.001 * xs
          res_dxs(k) = 0.001 * dxs
          if (iav == 1) then   
            res_av(k)='MXW'      
          else
            res_av(k)=''
          endif
          res_exist = .true.
        endif
      endif
    enddo
    close (2)
  enddo
  Nres = k
  return
end subroutine readmacs
! Copyright A.J. Koning 2025
