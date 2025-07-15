subroutine readthermal
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose: Read data from thermal cross section databases
!
! Revision    Date      Author      Quality  Description
! ======================================================
!    1     01-01-2019   A.J. Koning    A     Original code
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
  character(len=4)   :: year       ! year
  character(len=9)   :: ref        ! reference
  character(len=9)   :: sub        ! subentry
  character(len=24)  :: author     ! author
  character(len=132) :: thermfile  ! thermfile
  character(len=1200):: line       ! input line
  integer            :: Z          ! charge number
  integer            :: A          ! mass number
  integer            :: N          ! counter
  integer            :: IZA        ! ZA identifier
  integer            :: isoT       ! target isomer
  integer            :: isoR       ! residual isomer
  integer            :: MT         ! MT number
  integer            :: lib        ! library number
  integer            :: type       ! reaction type
  integer            :: istat      ! error code
  real(sgl)          :: Einc       ! incident energy
  real(sgl)          :: xs0        ! cross section
  real(sgl)          :: xs         ! cross section
  real(sgl)          :: dxs        ! cross section uncertainty
  real(sgl)          :: rochread
  real(sgl)          :: CE
  real(sgl)          :: chi2
!
! **************** Read database for thermal cross sections *****
! type = 1: (n,tot)
! type = 2: (n,el)
! type = 3: (n,f)
! type = 4: (n,g)
! type = 5: (n,p)
! type = 6: (n,a)
! type = 7: nubar total
! type = 8: nubar delayed
! type = 9: nubar prompt
!
  therm_xs = 0.
  therm_dxs = 0.
  therm_ref = ''
  therm_exist = .false.
  ratio_xs = 0.
  Ntherm_xs = 0
  Etherm_xs = 0.
  Etherm_dxs = 0.
  Etherm_author = ''
  Etherm_year = ''
  Etherm_subentry = ''
  Etherm_ref = ''
  Eratio_xs = 0.
  Etherm_exist = .false.
  capratio_xs = 0.
  Ltherm_xs = 0.
  Lratio_xs = 0.
!
! RIPL-3 thermal database
!
  write(*, *) "Reading data from databases....."
  lib = 1
  thermfile = trim(filespath)//'thermal.ripl'
  open (unit = 1, status = 'old', file = thermfile)
  do
    read(1, '(a)', iostat = istat) line
    if (istat == -1) exit
    if (istat > 0) call read_error(thermfile, istat)
    if (line(1:2) == 'ZA') cycle
    read(line, * ) IZA, isoT, MT, isoR, xs ,dxs, ref
    Z=IZA/1000
    A=IZA-1000*Z
    if (isoR == 99) isoR = -1
    type = 0
    if (MT == 18) type = 3
    if (MT == 102) type = 4
    if (MT == 103) type = 5
    if (MT == 107) type = 6
    if (type == 0) cycle
    therm_xs(lib, type, Z, A, isoT, isoR) = xs
    therm_dxs(lib, type, Z, A, isoT, isoR) = dxs
    therm_ref(lib, type, Z, A, isoT, isoR) = ref
    therm_exist(lib, type, isoR) = .true.
  enddo
  close (1)
!
! Correct or remove wrong/suspicious values from RIPL database
!
  therm_xs(lib, 3, 86, 223, 1, -1) = 0.
  therm_xs(lib, 4, 99, 254, 1, -1) = 0.
  therm_xs(lib, 4, 24, 51, 0, -1) = 0.
  therm_xs(lib, 5, 35, 76, 0, -1) = 1000. * therm_xs(lib, 5, 35, 76, 0, -1) 
  therm_dxs(lib, 5, 35, 76, 0, -1) = 1000. * therm_dxs(lib, 5, 35, 76, 0, -1) 
  therm_xs(lib, 6, 48, 113, 0, -1) = 0.
  therm_xs(lib, 6, 24, 51, 0, -1) = 0.
  therm_xs(lib, 6, 26, 55, 0, -1) = 0.
  therm_xs(lib, 6, 26, 59, 0, -1) = 0.
  therm_xs(lib, 6, 54, 127, 0, -1) = 0.
  therm_xs(lib, 6, 71, 176, 0, -1) = 0.
  therm_xs(lib, 6, 76, 192, 0, -1) = 0.
  therm_xs(lib, 6, 81, 203, 0, -1) = 0.
!
! Mughabghab 2006 thermal database
!
  lib = 2
  thermfile = trim(filespath)//'thermal.mugh06'
  open (unit = 1, status = 'old', file = thermfile)
  do
    read(1, '(a)', iostat = istat) line
    if (istat == -1) exit
    if (istat > 0) call read_error(thermfile, istat)
    if (line(1:1) == '#') cycle
    read(line(1:16), * ) A, Z, isoT
    isoR = -1
    type = 4
    xs = rochread(line(29:37))
    dxs = rochread(line(41:49))
    if (xs > 0.) then
      therm_xs(lib, type, Z, A, isoT, isoR) = xs
      therm_dxs(lib, type, Z, A, isoT, isoR) = dxs
      therm_ref(lib, type, Z, A, isoT, isoR) = 'Mugh06'
      therm_exist(lib, type, isoR) = .true.
    endif
    xs = rochread(line(51:59))
    dxs = rochread(line(63:69))
    if (xs > 0.) then
      therm_xs(lib, type, Z, A, isoT, isoR) = xs
      therm_dxs(lib, type, Z, A, isoT, isoR) = dxs
      therm_ref(lib, type, Z, A, isoT, isoR) = 'Mugh06'
      therm_exist(lib, type, isoR) = .true.
    endif
    type = 2
    xs = rochread(line(74:79))
    dxs = rochread(line(85:90))
    if (xs > 0.) then
      therm_xs(lib, type, Z, A, isoT, isoR) = xs
      therm_dxs(lib, type, Z, A, isoT, isoR) = dxs
      therm_ref(lib, type, Z, A, isoT, isoR) = 'Mugh06'
      therm_exist(lib, type, isoR) = .true.
    endif
    type = 6
    xs = rochread(line(97:103))
    dxs = rochread(line(109:114))
    if (xs > 0.) then
      therm_xs(lib, type, Z, A, isoT, isoR) = xs
      therm_dxs(lib, type, Z, A, isoT, isoR) = dxs
      therm_ref(lib, type, Z, A, isoT, isoR) = 'Mugh06'
      therm_exist(lib, type, isoR) = .true.
    endif
    type = 3
    xs = rochread(line(117:123))
    dxs = rochread(line(129:134))
    if (xs > 0.) then
      therm_xs(lib, type, Z, A, isoT, isoR) = xs
      therm_dxs(lib, type, Z, A, isoT, isoR) = dxs
      therm_ref(lib, type, Z, A, isoT, isoR) = 'Mugh06'
      therm_exist(lib, type, isoR) = .true.
    endif
  enddo
  close (1)
!
! Mughabghab 2016 thermal database
!
  lib = 3
  thermfile = trim(filespath)//'thermal.mugh18'
  open (unit = 1, status = 'old', file = thermfile)
  do
    read(1, '(a)', iostat = istat) line
    if (istat == -1) exit
    if (istat > 0) call read_error(thermfile, istat)
    if (line(4:4) < 'A' .or. line(4:4) > 'Z') cycle
    read(line(10:24), * ) A, Z, isoT
!
! (n,g)
!
    type = 4
    isoR = -1
    xs = rochread(line(61:69))
    dxs = rochread(line(73:82))
    if (xs > 0.) then
      therm_xs(lib, type, Z, A, isoT, isoR) = xs
      therm_dxs(lib, type, Z, A, isoT, isoR) = dxs
      therm_ref(lib, type, Z, A, isoT, isoR) = 'Mugh18'
      therm_exist(lib, type, isoR) = .true.
    endif
    xs0 = xs
    xs = rochread(line(133:143))
    dxs = rochread(line(147:157))
    if (xs > 0.) then
      therm_xs(lib, type, Z, A, isoT, isoR) = xs
      therm_dxs(lib, type, Z, A, isoT, isoR) = dxs
      therm_ref(lib, type, Z, A, isoT, isoR) = 'Mugh18'
      therm_exist(lib, type, isoR) = .true.
    endif
    if (xs0 > 0. .and. xs >0.) write(*,*) "Warning: perhaps wrong column in Mugh18 ", Z, A
    isoR = 1
    xs = rochread(line(85:93))
    dxs = rochread(line(97:106))
    if (xs > 0.) then
      therm_xs(lib, type, Z, A, isoT, isoR) = xs
      therm_dxs(lib, type, Z, A, isoT, isoR) = dxs
      therm_ref(lib, type, Z, A, isoT, isoR) = 'Mugh18'
      therm_exist(lib, type, isoR) = .true.
    endif
    xs0 = xs
    xs = rochread(line(160:170))
    dxs = rochread(line(174:184))
    if (xs > 0.) then
      therm_xs(lib, type, Z, A, isoT, isoR) = xs
      therm_dxs(lib, type, Z, A, isoT, isoR) = dxs
      therm_ref(lib, type, Z, A, isoT, isoR) = 'Mugh18'
      therm_exist(lib, type, isoR) = .true.
    endif
    if (xs0 > 0. .and. xs >0.) write(*,*) "Warning: perhaps wrong column in Mugh18 ", Z, A
    isoR = 0
    xs = rochread(line(109:117))
    dxs = rochread(line(121:130))
    if (xs > 0.) then
      therm_xs(lib, type, Z, A, isoT, isoR) = xs
      therm_dxs(lib, type, Z, A, isoT, isoR) = dxs
      therm_ref(lib, type, Z, A, isoT, isoR) = 'Mugh18'
      therm_exist(lib, type, isoR) = .true.
    endif
    xs0 = xs
    xs = rochread(line(187:197))
    dxs = rochread(line(201:211))
    if (xs > 0.) then
      therm_xs(lib, type, Z, A, isoT, isoR) = xs
      therm_dxs(lib, type, Z, A, isoT, isoR) = dxs
      therm_ref(lib, type, Z, A, isoT, isoR) = 'Mugh18'
      therm_exist(lib, type, isoR) = .true.
    endif
    if (xs0 > 0. .and. xs >0.) write(*,*) "Warning: perhaps wrong column in Mugh18 ", Z, A
    isoR = -1
!
! (n,el)
!
    type = 2
    xs = rochread(line(214:225))
    dxs = rochread(line(229:236))
    if (xs > 0.) then
      therm_xs(lib, type, Z, A, isoT, isoR) = xs
      therm_dxs(lib, type, Z, A, isoT, isoR) = dxs
      therm_ref(lib, type, Z, A, isoT, isoR) = 'Mugh18'
      therm_exist(lib, type, isoR) = .true.
    endif
    xs0 = xs
    xs = rochread(line(239:250))
    dxs = rochread(line(254:261))
    if (xs > 0.) then
      therm_xs(lib, type, Z, A, isoT, isoR) = xs
      therm_dxs(lib, type, Z, A, isoT, isoR) = dxs
      therm_ref(lib, type, Z, A, isoT, isoR) = 'Mugh18'
      therm_exist(lib, type, isoR) = .true.
    endif
!
! (n,a)
!
    type = 6
    xs = 0.
    xs = rochread(line(264:272))
    dxs = rochread(line(276:285))
    if (xs > 0.) then
      therm_xs(lib, type, Z, A, isoT, isoR) = xs
      therm_dxs(lib, type, Z, A, isoT, isoR) = dxs
      therm_ref(lib, type, Z, A, isoT, isoR) = 'Mugh18'
      therm_exist(lib, type, isoR) = .true.
    endif
!
! (n,f)
!
    type = 3
    xs = rochread(line(288:295))
    dxs = rochread(line(299:308))
    if (xs > 0.) then
      therm_xs(lib, type, Z, A, isoT, isoR) = xs
      therm_dxs(lib, type, Z, A, isoT, isoR) = dxs
      therm_ref(lib, type, Z, A, isoT, isoR) = 'Mugh18'
      therm_exist(lib, type, isoR) = .true.
    endif
!
! (n,p)
!
    type = 5
    xs = rochread(line(743:753))
    dxs = rochread(line(757:764))
    if (xs > 0.) then
      therm_xs(lib, type, Z, A, isoT, isoR) = xs
      therm_dxs(lib, type, Z, A, isoT, isoR) = dxs
      therm_ref(lib, type, Z, A, isoT, isoR) = 'Mugh18'
      therm_exist(lib, type, isoR) = .true.
    endif
!
! nubar
!
    type = 7
    xs = rochread(line(311:318))
    dxs = rochread(line(322:331))
    if (xs > 0.) then
      therm_xs(lib, type, Z, A, isoT, isoR) = xs
      therm_dxs(lib, type, Z, A, isoT, isoR) = dxs
      therm_ref(lib, type, Z, A, isoT, isoR) = 'Mugh18'
      therm_exist(lib, type, isoR) = .true.
    endif
!
! prompt nubar
!
    type = 9
    xs = rochread(line(334:341))
    dxs = rochread(line(354:354))
    if (xs > 0.) then
      therm_xs(lib, type, Z, A, isoT, isoR) = xs
      therm_dxs(lib, type, Z, A, isoT, isoR) = dxs
      therm_ref(lib, type, Z, A, isoT, isoR) = 'Mugh18'
      therm_exist(lib, type, isoR) = .true.
    endif
  enddo
  close (1)
!
! Kayzero thermal database
!
  lib = 4
  thermfile = trim(filespath)//'kayzero.txt'
  open (unit = 1, status = 'old', file = thermfile)
  do
    read(1, '(a)', iostat = istat) line
    if (istat == -1) exit
    if (istat > 0) call read_error(thermfile, istat)
    if (line(4:4) < 'A' .or. line(4:4) > 'Z') cycle
    read(line(10:24), * ) A, Z, isoT
!
! (n,g)
!
    type = 4
    isoR = -1
    xs = rochread(line(89:98))
    dxs = rochread(line(101:110))
    if (xs > 0.) then
      therm_xs(lib, type, Z, A, isoT, isoR) = xs
      therm_dxs(lib, type, Z, A, isoT, isoR) = dxs
      therm_ref(lib, type, Z, A, isoT, isoR) = 'Kayzero'
      therm_exist(lib, type, isoR) = .true.
    endif
  enddo
  close (1)
!
! Sukhoruchkin 2015 thermal database
!
  lib = 5
  thermfile = trim(filespath)//'sukhoruchkin.txt'
  open (unit = 1, status = 'old', file = thermfile)
  do
    read(1, '(a)', iostat = istat) line
    if (istat == -1) exit
    if (istat > 0) call read_error(thermfile, istat)
    if (line(4:4) < 'A' .or. line(4:4) > 'Z') cycle
    read(line(10:24), * ) A, Z, isoT
!
! (n,g)
!
    type = 4
    isoR = -1
    xs = rochread(line(87:99))
    dxs = rochread(line(102:111))
    if (xs > 0.) then
      therm_xs(lib, type, Z, A, isoT, isoR) = xs
      therm_dxs(lib, type, Z, A, isoT, isoR) = dxs
      therm_ref(lib, type, Z, A, isoT, isoR) = 'Sukhoruchkin'
      therm_exist(lib, type, isoR) = .true.
    endif
!
! (n,a)
!
    type = 6
    xs = 0.
    xs = rochread(line(137:147))
    dxs = rochread(line(152:161))
    if (xs > 0.) then
      therm_xs(lib, type, Z, A, isoT, isoR) = xs
      therm_dxs(lib, type, Z, A, isoT, isoR) = dxs
      therm_ref(lib, type, Z, A, isoT, isoR) = 'Sukhoruchkin'
      therm_exist(lib, type, isoR) = .true.
    endif
!
! (n,p)
!
    type = 5
    xs = 0.
    xs = rochread(line(163:172))
    dxs = rochread(line(177:185))
    if (xs > 0.) then
      therm_xs(lib, type, Z, A, isoT, isoR) = xs
      therm_dxs(lib, type, Z, A, isoT, isoR) = dxs
      therm_ref(lib, type, Z, A, isoT, isoR) = 'Sukhoruchkin'
      therm_exist(lib, type, isoR) = .true.
    endif
!
! (n,f)
!
    type = 3
    xs = 0.
    xs = rochread(line(38:48))
    dxs = rochread(line(49:60))
    if (xs > 0.) then
      therm_xs(lib, type, Z, A, isoT, isoR) = xs
      therm_dxs(lib, type, Z, A, isoT, isoR) = dxs
      therm_ref(lib, type, Z, A, isoT, isoR) = 'Sukhoruchkin'
      therm_exist(lib, type, isoR) = .true.
    endif
  enddo
  close (1)
!
! EXFOR thermal database
!
  do type = 1, numtype
    thermfile = trim(filespath)//'thermal_exfor.'//reac(type)
    open (unit = 1, status = 'old', file = thermfile)
    do
      read(1, '(a)', iostat = istat) line
      if (istat == -1) exit
      if (istat > 0) call read_error(thermfile, istat)
      read(line, * ) Z, A, isoT, isoR, xs ,dxs, Einc, author, year, sub
      Ntherm_xs(type, Z, A, isoT, isoR) = Ntherm_xs(type, Z, A, isoT, isoR) + 1
      N = Ntherm_xs(type, Z, A, isoT, isoR)
      if (N > numex) then
        Ntherm_xs(type, Z, A, isoT, isoR) = numex
        exit
      endif
      if (type <= 6) then
        Etherm_xs(type, Z, A, isoT, isoR, N) = 0.001*xs
        Etherm_dxs(type, Z, A, isoT, isoR, N) = 0.001*dxs
      else
        Etherm_xs(type, Z, A, isoT, isoR, N) = xs
        Etherm_dxs(type, Z, A, isoT, isoR, N) = dxs
      endif
      Etherm_author(type, Z, A, isoT, isoR, N) = author
      Etherm_year(type, Z, A, isoT, isoR, N) = year
      Etherm_subentry(type, Z, A, isoT, isoR, N) = sub
      Etherm_ref(type, Z, A, isoT, isoR, N) = trim(author)//'_'//trim(year)//'_'//trim(sub)
      Etherm_exist(type, isoR) = .true.
    enddo
    close (1)
  enddo
!
! Nuclear data libraries
!
  isoR = -1
! do type = 1, numtype
  do type = 4, 4
    do lib = 1, numndlib
      thermfile = trim(libspath)//trim(ndlib(lib))//'.therm'
      open (unit = 1, status = 'old', file = thermfile)
      do
        read(1, '(a)', iostat = istat) line
        if (istat == -1) exit
        if (istat > 0) call read_error(thermfile, istat)
        if (line(1:1) == '#') cycle
        read(line, * ) Z, A, isoT, CE, chi2, xs
        Ltherm_xs(lib, type, Z, A, isoT, isoR) = xs
      enddo
      close (1)
    enddo
  enddo
  return
end subroutine readthermal
! Copyright A.J. Koning 2019
