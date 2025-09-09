subroutine readthermal(Z, A, Liso, Riso, type)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose: Read data from thermal cross section databases
!
! Revision    Date      Author      Quality  Description
! ======================================================
!    1     2025-08-19   A.J. Koning    A     Original code
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
  logical            :: lexist
  logical            :: shin
  integer            :: year       ! year
  character(len=2)   :: Unuc
  character(len=3)   :: exten      ! extension
  character(len=3)   :: Zstring
  character(len=3)   :: Astring
  character(len=4)   :: ext
  character(len=9)   :: reacdir(6) 
  character(len=9)   :: ref        ! reference
  character(len=12)  :: sub        ! subentry
  character(len=10)   :: nucstring
  character(len=24)  :: author     ! author
  character(len=132) :: thermfile  ! thermfile
  character(len=1200):: line       ! input line
  integer            :: iz         ! charge number
  integer            :: ia         ! mass number
  integer            :: Z          ! charge number
  integer            :: A          ! mass number
  integer            :: k        
  integer            :: iav
  integer            :: L
  integer            :: mxw1
  integer            :: mxw2
  integer            :: mxw3
  integer            :: recom
  integer            :: go       
  integer            :: IZA        ! ZA identifier
  integer            :: Liso       ! target isomer
  integer            :: Riso       ! residual isomer
  integer            :: isoT       ! target isomer
  integer            :: isoR       ! residual isomer
  integer            :: MT         ! MT number
  integer            :: lib        ! library number
  integer            :: type       ! reaction type
  integer            :: istat      ! error code
  real(sgl)          :: Einc       ! incident energy
  real(sgl)          :: dEinc 
  real(sgl)          :: G
  real(sgl)          :: dum
  real(sgl)          :: xs0        ! cross section
  real(sgl)          :: xs1        ! cross section
  real(sgl)          :: xs         ! cross section
  real(sgl)          :: dxs        ! cross section uncertainty
  real(sgl)          :: dxs0       ! cross section uncertainty
  real(sgl)          :: dxs1       ! cross section uncertainty
  real(sgl)          :: rochread
  real(sgl)          :: CE
  real(sgl)          :: chi2
!
! **************** Read database for thermal cross sections *****
!
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
  res_author = ''
  res_type = ''
  res_year = 0
  res_ref = ''
  res_xs = 0.
  res_dxs = 0.
  res_G = 0.
  res_av = ''
  res_E = 2.53e-8
  ref = ''
  Nres_exp = 0
  Nres = 0
  k = 0
!
! RIPL-3 thermal database
!
  xs = 0.
  dxs = 0.
  if (type >= 3 .and. type <= 6) then
    thermfile = trim(filespath)//'thermal.ripl'
    open (unit = 2, status = 'old', file = thermfile)
    read(2, '()')
    do
      read(2, '(a)', iostat = istat) line
      if (istat == -1) exit
      if (istat > 0) call read_error(thermfile, istat)
      read(line, * ) IZA, isoT, MT, isoR, xs ,dxs
      iz=IZA/1000
      ia=IZA-1000*Z
      if (isoR == 99) isoR = -1
      if (iz == Z .and. ia == A .and. Liso == isoT .and. Riso == isoR) then
        go = 0
        if (type == 3 .and. MT == 18) go = 1
        if (type == 4 .and. MT == 102) go = 1
        if (type == 5 .and. MT == 103) go = 1
        if (type == 6 .and. MT == 107) go = 1
        if (go == 1) then
          k = k + 1
          res_author(k) = 'RIPL-3'
          res_type(k) = 'Compilation'
          res_year(k) = 2004
          res_ref(k) = ref
          res_av(k) = ''
          res_xs(k) = xs
          res_dxs(k) = dxs
          res_exist = .true.
        endif
        exit
      endif
    enddo
    close (2)
  endif
!
! Mughabghab 2006 thermal database
!
  do iav = 0, 1
    xs = 0.
    dxs = 0.
    ref = ''
    if (iav == 0 .and. (type == 2 .or. type == 3 .or. type == 6) .or. type == 4) then
      thermfile = trim(filespath)//'atlas_2006.txt'
      open (unit = 2, status = 'old', file = thermfile)
      do
        read(2, '(a)', iostat = istat) line
        if (istat == -1) exit
        if (istat > 0) call read_error(thermfile, istat)
        if (line(1:1) == '#') cycle
        read(line(1:16), * ) ia, iz, isoT
        if (iz == Z .and. ia == A .and. Liso == isoT .and. Riso == -1) then
          xs = 0.
          dxs = 0.
          go = 0
          if (type == 4) then
            if (iav == 1) then
              xs = rochread(line(30:36))
              dxs = rochread(line(42:48))
            else
              xs = rochread(line(51:59))
              dxs = rochread(line(63:69))
            endif
            if (xs > 0.) go = 1
          endif
          if (type == 2) then
            xs = rochread(line(74:79))
            dxs = rochread(line(85:90))
            if (xs > 0.) go = 1
          endif
          if (type == 6) then
            xs = rochread(line(97:103))
            dxs = rochread(line(109:114))
            if (xs > 0.) go = 1
          endif
          if (type == 3) then
            xs = rochread(line(117:123))
            dxs = rochread(line(129:134))
            if (xs > 0.) go = 1
          endif
          if (go == 1) then
            k = k + 1
            res_author(k) = 'Mughabghab-2006'
            res_type(k) = 'Compilation'
            res_year(k) = 2006
            res_ref(k) = ref
            res_xs(k) = xs
            res_dxs(k) = dxs
            if (iav == 1) then
              res_av(k) = 'MXW'
            else
              res_av(k) = ''
            endif
            res_exist = .true.
          endif
          exit
        endif
      enddo
      close (2)
    endif
  enddo
!
! Mughabghab 2018 thermal database
!
  xs = 0.
  dxs = 0.
  thermfile = trim(filespath)//'global.2023.txt'
  open (unit = 2, status = 'old', file = thermfile)
  do
    read(2, '(a)', iostat = istat) line
    if (istat == -1) exit
    if (istat > 0) call read_error(thermfile, istat)
    if (line(4:4) < 'A' .or. line(4:4) > 'Z') cycle
    read(line(10:24), * ) ia, iz, isoT
    if (iz == Z .and. ia == A .and. Liso == isoT) then
      xs0 = 0.
      dxs0 = 0.
      xs1 = 0.
      dxs1 = 0.
      xs = 0.
      dxs = 0.
      go = 0
!
! (n,g)
!
      if (type == 4) then
        if (Riso == -1) then
          xs0 = rochread(line(61:69))
          dxs0 = rochread(line(73:82))
          xs1 = rochread(line(133:143))
          dxs1 = rochread(line(147:157))
          xs = max(xs1, xs0)
          dxs = max(dxs1, dxs0)
          if (xs > 0.) go = 1
        endif
        if (Riso == 1) then
          xs0 = rochread(line(85:93))
          dxs0 = rochread(line(97:106))
          xs = max(rochread(line(160:170)),xs0)
          dxs = max(rochread(line(174:184)),dxs0)
          xs = max(xs1, xs0)
          dxs = max(dxs1, dxs0)
          if (xs > 0.) go = 1
        endif
        if (Riso == 0) then
          xs0 = rochread(line(109:117))
          dxs0 = rochread(line(121:130))
          xs1 = rochread(line(187:197))
          dxs1 = rochread(line(201:211))
          xs = max(xs1, xs0)
          dxs = max(dxs1, dxs0)
          if (xs > 0.) go = 1
        endif
      else
        if (Riso /= -1) cycle
      endif
!
! (n,el)
!
      if (type == 2) then
        xs0 = rochread(line(214:225))
        dxs0 = rochread(line(229:236))
        xs1 = rochread(line(239:250))
        dxs1 = rochread(line(254:261))
        xs = max(xs1, xs0)
        dxs = max(dxs1, dxs0)
        if (xs > 0.) go = 1
      endif
!
! (n,a)
!
      if (type == 6) then
        xs = rochread(line(264:272))
        dxs = rochread(line(276:285))
        if (xs > 0.) go = 1
      endif
!
! (n,f)
!
      if (type == 3) then
        xs = rochread(line(288:295))
        dxs = rochread(line(299:308))
        if (xs > 0.) go = 1
      endif
!
! (n,p)
!
      if (type == 5) then
        xs = rochread(line(743:753))
        dxs = rochread(line(757:764))
        if (xs > 0.) go = 1
      endif
!
! nubar
!
      if (type == 7) then
        xs = rochread(line(311:318))
        dxs = rochread(line(322:331))
        if (xs > 0.) go = 1
      endif
!
! prompt nubar
!
      if (type == 9) then
        xs = rochread(line(334:341))
        dxs = rochread(line(354:354))
        if (xs > 0.) go = 1
      endif
      if (xs <= 0.) go = 0
!
! Add to database
!
      if (go == 1) then
        k = k + 1
        res_author(k) = 'Mughabghab-2018'
        res_type(k) = 'Compilation'
        res_year(k) = 2018
        res_ref(k) = ref
        res_xs(k) = xs
        res_dxs(k) = dxs
        res_av(k) = ''
        res_exist = .true.
      endif
      if (xs0 > 0. .and. xs1 > 0.) write(*,*) "Warning: perhaps double data in Mughabghab-2018 ", Z, A, " xs0 ",xs0," xs1 ",xs1
      exit
    endif
  enddo
  close (2)
!
! Kayzero thermal database
!
  xs = 0.
  dxs = 0.
  if (type == 4) then
    thermfile = trim(filespath)//'kayzero.txt'
    open (unit = 2, status = 'old', file = thermfile)
    do
      read(2, '(a)', iostat = istat) line
      if (istat == -1) exit
      if (istat > 0) call read_error(thermfile, istat)
      if (line(4:4) < 'A' .or. line(4:4) > 'Z') cycle
      read(line(10:24), * ) ia, iz, isoT
!
! (n,g)
!
      if (iz == Z .and. ia == A .and. Liso == isoT .and. Riso == -1) then
        xs = rochread(line(89:98))
        dxs = rochread(line(101:110))
        if (xs > 0.) then
          k = k + 1
          res_author(k) = 'Kayzero'
          res_type(k) = 'Compilation'
          res_year(k) = 2018
          res_ref(k) = ref
          res_xs(k) = xs
          res_dxs(k) = dxs
          res_av(k) = ''
          res_exist = .true.
        endif
        exit
      endif
    enddo
    close (2)
  endif
!
! Firestone thermal database
!
  xs = 0.
  dxs = 0.
  if (type == 4 .and. Liso == 0 .and. Riso == -1) then
    thermfile = trim(filespath)//'firestone2022.txt'
    open (unit = 2, status = 'old', file = thermfile)
    do
      read(2, '(a)', iostat = istat) line
      if (istat == -1) exit
      if (istat > 0) call read_error(thermfile, istat)
      if (line(1:1) == '#') cycle
      read(line(12:14), * ) ia
      read(line(17:20), * ) iz
!
! (n,g)
!
      if (iz == Z .and. ia == A) then
        xs = rochread(line(63:72))
        dxs = rochread(line(77:86))
        if (xs > 0.) then
          k = k + 1
          res_author(k) = 'Firestone-2022'
          res_type(k) = 'Compilation'
          res_year(k) = 2022
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
  endif
!
! Sukhoruchkin 2015 thermal database
!
  do iav= 0, 1
    xs = 0.
    dxs = 0.
    if (type == 4 .or. (iav == 0 .and. (type == 3 .or. type == 5 .or. type == 6))) then
      thermfile = trim(filespath)//'sukhoruchkin.txt'
      open (unit = 2, status = 'old', file = thermfile)
      do
        read(2, '(a)', iostat = istat) line
        if (istat == -1) exit
        if (istat > 0) call read_error(thermfile, istat)
        if (line(4:4) < 'A' .or. line(4:4) > 'Z') cycle
        read(line(10:24), * ) ia, iz, isoT
!
! (n,g)
!
        if (iz == Z .and. ia == A .and. Liso == isoT .and. Riso == -1) then
          if (type == 4) then
            if (iav == 1) then
              if (line(87:87) == ' ') exit
            else
              if (line(87:87) == 'W') exit
            endif
            xs = rochread(line(88:99))
            dxs = rochread(line(102:111))
            if (xs > 0.) go = 1
          endif
!
! (n,a)
!
          if (type == 6) then
            xs = rochread(line(137:147))
            dxs = rochread(line(152:161))
            if (xs > 0.) go = 1
          endif
!
! (n,p)
!
          if (type == 5) then
            xs = rochread(line(163:172))
            dxs = rochread(line(177:185))
            if (xs > 0.) go = 1
          endif
!
! (n,f)
!
          if (type == 3) then
            xs = rochread(line(38:48))
            dxs = rochread(line(49:60))
            if (xs > 0.) go = 1
          endif
          if (go == 1) then
            k = k + 1
            res_author(k) = 'Sukhoruchkin'
            res_type(k) = 'Compilation'
            res_year(k) = 2015
            res_ref(k) = ref
            res_xs(k) = xs
            res_dxs(k) = dxs
            if (iav == 1) then
              res_av(k) = 'MXW'
            else
              res_av(k) = ''
            endif
            res_exist = .true.
          endif
          exit
        endif
      enddo
      close (2)
    endif
  enddo
!
! Nuclear data libraries
!
  xs = 0.
  dxs = 0.
  if (type <= 6) then
    if (type == 1) ext = 'tot'
    if (type == 2) ext = 'el'
    if (type == 3) ext = 'nf'
    if (type == 4) ext = 'ng'
    if (type == 5) ext = 'np'
    if (type == 6) ext = 'na'
    if (type >= 4) then
      if (Riso == 0) ext=trim(ext)//'_g'
      if (Riso == 1) ext=trim(ext)//'_m'
      if (Riso == 2) ext=trim(ext)//'_n'
    endif
    do lib = 1, numndlib
      if (type <= 3 .and. Riso /= -1) cycle
      thermfile = trim(libspath)//trim(ndlib(lib))//'.therm_'//ext
      inquire (file = thermfile, exist = lexist)
      if (lexist) then
        open (unit = 2, status = 'old', file = thermfile)
        do
          read(2, '(a)', iostat = istat) line
          if (istat == -1) exit
          if (istat > 0) call read_error(thermfile, istat)
          if (line(1:1) == '#') cycle
          read(line, * ) iz, ia, isoT, CE, chi2, xs, dum, dum, dum, G
          if (iz == Z .and. ia == A .and. Liso == isoT) then
            k = k + 1
            res_author(k) = ndlib(lib)
            res_type(k) = 'NDL'
            res_year(k) = ndyear(lib)
            res_ref(k) = ref
            res_xs(k) = xs
            res_dxs(k) = dxs
            res_G(k) = G
            res_av(k) = ''
            res_exist = .true.
            exit
          endif
        enddo
        close (2)
      endif
    enddo
  endif
!
! Shin Okumura EXFOR mining
!
! Get thermal c.s. from Shin Okumra repo, and nubar from EXFORtables
!
  reacdir(1) = 'n-tot'
  reacdir(2) = 'n-el'
  reacdir(3) = 'n-f'
  reacdir(4) = 'n-g'
  reacdir(5) = 'n-p'
  reacdir(6) = 'n-a'
  if (type <= 6) then
    shin = .true.
  else
    shin = .false.
  endif
  if (shin) then
    isoR = -1
    xs = 0.
    dxs = 0.
    nucstring='          '
    Zstring='   '
    write(Zstring,'(i3)') Z
    Astring='   '
    write(Astring,'(i3)') A
    Unuc=nuc(Z)
    if (Unuc(2:2) /= ' ') Unuc(2:2) = achar(iachar(Unuc(2:2)) - 32) 
    write(nucstring,'(a,"-",a,"-",a)') trim(adjustl(Zstring)),trim(Unuc),trim(adjustl(Astring))
    if (Liso == 1) nucstring=trim(nucstring)//'-M'
    thermfile = trim(filespath)//'exforfiles/'//trim(reacdir(type))//'/'//trim(nucstring)//'.txt'
    inquire (file = thermfile, exist = lexist)
    if (lexist .and. Liso == 0) then
      open (unit = 2, status = 'old', file = thermfile)
      do
        read(2, '(a)', iostat = istat) line
        if (istat == -1) exit
        if (istat > 0) call read_error(thermfile, istat)
        if (line(1:10) == '# Residual') then
          read(2, '(//,a)', iostat = istat) line
          L=len_trim(line)
          if (line(L-1:L) == '-G') isoR = 0
          if (line(L-1:L) == '-M') isoR = 1
          if (line(L-1:L) == '-N') isoR = 2
        endif
        if (line(1:1) == '#') cycle
        if (line(1:1) == ' ') cycle
        if (Riso /= isoR) cycle
        recom=index(line,'RECOM') 
        recom=max(recom,index(line,'DERIV') )
        recom=max(recom,index(line,'RAW') )
        if (recom > 0) cycle
        mxw1=index(line,'MXW') 
        mxw2=index(line,'AV') 
        mxw3=index(line,'SPA') 
        sub=line(1:12)
        author=line(21:44)
        read(line(45:120), *, iostat=istat ) year, Einc, dEinc, xs, dxs
        if (istat > 0) then
          write(*,*) "EXFOR problem: ",trim(line)
          cycle
        endif
        if (xs <= 0.) cycle
        Nres_exp = Nres_exp + 1
        if (Nres_exp  > numex - 20) then
          Nres_exp = numex - 20
          exit
        endif
        k = k + 1
        res_author(k) = author
        res_type(k) = 'EXFOR'
        res_year(k) = year
        res_ref(k) = sub
        res_xs(k) = xs
        res_dxs(k) = dxs
        res_E(k) = Einc
        res_av(k) = ''
        if (mxw1 > 0) res_av(k) = 'MXW'
        if (mxw2 > 0) res_av(k) = 'AV'
        if (mxw3 > 0) res_av(k) = 'SPA'
        res_exist = .true.
      enddo
      close (2)
    endif
  else
!
! EXFOR thermal database
!
    do iav = 0, 1
      if (iav == 0) then
        exten=''
      else
        exten='_av'
      endif
      xs = 0.
      dxs = 0.
      thermfile = trim(exforpath)//'exfor_thermal'//trim(exten)//'_'//trim(reac(type))//'.txt'
      inquire (file = thermfile, exist = lexist)
      if (lexist) then
        open (unit = 2, status = 'old', file = thermfile)
        do
          read(2, '(a)', iostat = istat) line
          if (istat == -1) exit
          if (istat > 0) call read_error(thermfile, istat)
          read(line, * ) iz, ia, isoT, isoR, xs ,dxs, Einc, author, year, sub
          if (iz == Z .and. ia == A .and. Liso == isoT .and. Riso == isoR) then
            Nres_exp = Nres_exp + 1
            if (Nres_exp  > numex) then
              Nres_exp = numex
              exit
            endif
            if (type <= 6) then
              xs = 0.001 * xs
              dxs = 0.001 * dxs
            endif
            k = k + 1
            res_author(k) = author
            res_type(k) = 'EXFOR'
            res_year(k) = year
            res_ref(k) = sub
            res_xs(k) = xs
            res_dxs(k) = dxs
            if (iav == 1) then
              res_av(k) = 'MXW'
            else
              res_av(k) = ''
            endif
            res_exist = .true.
          endif
        enddo
        close (2)
      endif
    enddo
  endif
  Nres = k
  return
end subroutine readthermal
! Copyright A.J. Koning 2025
