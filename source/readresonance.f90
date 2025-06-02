subroutine readresonance(Z, A, Liso, type)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose: Read data from resonance databases
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
  character(len=2)   :: cdum       !
  character(len=2)   :: Unuc
  character(len=3)   :: Zstring
  character(len=3)   :: Astring
  character(len=10)   :: nucstring
  character(len=9)   :: ref        ! reference
  character(len=9)   :: sub        ! subentry
  character(len=24)  :: author     ! author
  character(len=132) :: Nrfile
  character(len=132) :: resfile    ! resfile
  character(len=2000):: line       ! input line
  logical            :: lexist
  integer            :: year       ! year
  integer            :: Z          ! charge number
  integer            :: iz         ! charge number
  integer            :: A          ! mass number
  integer            :: k
  integer            :: ix
  integer            :: lib
  integer            :: iripl
  integer            :: isoR
  integer            :: ia         ! mass number
  integer            :: Nrr
  integer            :: isoT       ! target isomer
  integer            :: Liso       ! target isomer
  integer            :: type       ! reaction type
  integer            :: istat      ! error code
  real(sgl)          :: xs        
  real(sgl)          :: dxs        
  real(sgl)          :: CE        
  real(sgl)          :: Emin
  real(sgl)          :: Emax
  real(sgl)          :: rL
  real(sgl)          :: rJ
  real(sgl)          :: rP
  real(sgl)          :: chi2        
  real(sgl)          :: rochread
  real(sgl)          :: Einc       ! incident energy
!
! **************** Read RIPL-3 database for average resonance parameters *****
!
! type = 1: D0
! type = 2: S0
! type = 3: gamgam0
! type = 4: D1
! type = 5: S1
! type = 6: gamgam1
! type = 7: D2
! type = 8: Ig
! type = 9: If
!
  res_author = ''
  res_type = ''
  res_year = 0
  res_ref = ''
  res_xs = 0.
  res_dxs = 0.
  res_Nrr = ''
  res_Emin = ''
  res_Emax = ''
  res_L = ''
  res_J = ''
  res_P = ''
  ref = ''
  Nres_exp = 0
  Nres = 0
  k = 0
  if (type <= 6) then
!
! RIPL resonance database
!
    do iripl= 2, 3
      xs = 0.
      dxs = 0.
      if (iripl == 2) then
        if (type <= 3) then
          resfile = trim(filespath)//'resonances_0.ripl2'
        else
          resfile = trim(filespath)//'resonances_1.ripl2'
        endif
      else
        if (type <= 3) then
          resfile = trim(filespath)//'resonances_0.ripl3'
        else
          resfile = trim(filespath)//'resonances_1.ripl3'
        endif
      endif
      inquire (file = resfile, exist = lexist)
      if (.not.lexist) then
        write(*,*) trim(resfile), ' does not exist'
        cycle
      endif
      open (unit = 2, status = 'old', file = resfile)
      do
        read(2, '(a)', iostat = istat) line
        if (istat == -1) exit
        if (istat > 0) call read_error(resfile, istat)
        if (line(1:1) == '#') cycle
        if (iripl == 2) then
          read(line, * ) iz, ia
        else
          read(line, * ) iz, cdum, ia
        endif
        if (iz == Z .and. ia == A .and. Liso == 0) then
          if (type == 1 .or. type == 4) then
            xs = 1000. * rochread(line(25:34))
            dxs = 1000. * rochread(line(35:44))
          endif
          if (type == 2 .or. type == 5) then
            xs = rochread(line(46:51))
            dxs = rochread(line(52:57))
          endif
          if (type == 3 .or. type == 6) then
            xs = 0.001 * rochread(line(58:63))
            dxs = 0.001 * rochread(line(64:69))
          endif
          if (xs > 0.) then
            k = k + 1
            ref = line(70:73)
            if (iripl == 2) then
              res_author(k) = 'RIPL-2'
              res_year(k) = 2000
              Nrfile = trim(filespath)//'obninsk.dat'
              open (unit = 3, status = 'old', file = Nrfile)
              read(3,'(//)')
              do
                read(3, '(i4,4x,i3,11x,i3)', iostat = istat) iz, ia, Nrr
                if (istat == -1) exit
                if (istat > 0) call read_error(resfile, istat)
                if (iz == Z .and. ia == A .and. Liso == 0) then
                  write(res_Nrr(k)(1:15),'(6x,i3)') Nrr
                  exit
                endif
              enddo
              close (3)
            else
              res_author(k) = 'RIPL-3'
              res_year(k) = 2009
            endif
            res_type(k) = 'Compilation'
            res_ref(k) = ref
            res_xs(k) = xs
            res_dxs(k) = dxs
            res_exist = .true.
          endif
          exit
        endif
      enddo
      close (2)
    enddo
  endif
!
! JUKO database for resonance integrals
!
  xs = 0.
  dxs = 0.
  ref = ''
  if (type >= 8) then
    resfile = trim(filespath)//'resint.juko'
    open (unit = 2, status = 'old', file = resfile)
    do
      read(2, '(a)', iostat = istat) line
      if (istat == -1) exit
      if (istat > 0) call read_error(resfile, istat)
      if (line(1:1) == '#') cycle
      read(line, * ) ia, iz
      if (iz == Z .and. ia == A .and. Liso == 0) then
        if (type == 8 .or. type == 9 .and. line(37:41) == '(n,f)') then
          xs = rochread(line(49:56))
          dxs = rochread(line(66:73))
          if (xs > 0.) then
            k = k + 1
            res_author(k) = 'JUKO'
            res_type(k) = 'Compilation'
            res_year(k) = 2000
            res_ref(k) = ref
            res_xs(k) = xs
            res_dxs(k) = dxs
            res_exist = .true.
          endif
          exit
        endif
      endif
    enddo
    close (2)
  endif
!
! Mughabghab 2018 database
!
  xs = 0.
  dxs = 0.
  ref = ''
  resfile = trim(filespath)//'global.2023.txt'
  open (unit = 2, status = 'old', file = resfile)
  do
    read(2, '(a)', iostat = istat) line
    if (istat == -1) exit
    if (istat > 0) call read_error(resfile, istat)
    if (line(4:4) < 'A' .or. line(4:4) > 'Z') cycle
    read(line(10:24), * ) ia, iz, isoT
    if (ia == 0) cycle
    if (iz == Z .and. ia == A .and. Liso == isoT) then
      if (type == 1) then
        xs = rochread(line(1031:1041))
        dxs = rochread(line(1045:1052))
      endif
      if (type == 2) then
        xs = rochread(line(1055:1065))
        dxs = rochread(line(1069:1076))
      endif
      if (type == 3) then
        xs = rochread(line(1175:1185))
        dxs = rochread(line(1189:1196))
      endif
      if (type == 4) then
        xs = rochread(line(935:945))
        dxs = rochread(line(949:956))
      endif
      if (type == 5) then
        xs = rochread(line(959:969))
        dxs = rochread(line(973:980))
      endif
      if (type == 6) then
        xs = rochread(line(1199:1209))
        dxs = rochread(line(1213:1220))
      endif
      if (type == 7) then
        xs = rochread(line(1103:1112))
        dxs = rochread(line(1117:1124))
      endif
      if (type == 8) then
        xs = rochread(line(607:615))
        dxs = rochread(line(619:626))
      endif
      if (type == 9) then
        xs = rochread(line(357:364))
        dxs = rochread(line(368:377))
      endif
      if (xs > 0.) then
        k = k + 1
        res_author(k) = 'Mughabghab-2018'
        res_type(k) = 'Compilation'
        res_year(k) = 2018
        res_ref(k) = ref
        res_xs(k) = xs
        res_dxs(k) = dxs
        res_exist = .true.
      endif
      exit
    endif
  enddo
  close (2)
!
! Kayzero database
!
  xs = 0.
  dxs = 0.
  if (type == 8) then
    resfile = trim(filespath)//'kayzero.txt'
    open (unit = 2, status = 'old', file = resfile)
    do
      read(2, '(a)', iostat = istat) line
      if (istat == -1) exit
      if (istat > 0) call read_error(resfile, istat)
      if (line(4:4) < 'A' .or. line(4:4) > 'Z') cycle
      read(line(10:24), * ) ia, iz, isoT
      if (ia == 0) cycle
      if (iz == Z .and. ia == A .and. Liso == isoT) then
        xs = rochread(line(112:122))
        dxs = rochread(line(127:136))
        if (xs > 0.) then
          k = k + 1
          res_author(k) = 'Kayzero'
          res_type(k) = 'Compilation'
          res_year(k) = 2018
          res_ref(k) = ref
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
! Sukhoruchkin database
!
  xs = 0.
  dxs = 0.
  if (type >= 8) then
    resfile = trim(filespath)//'sukhoruchkin.txt'
    open (unit = 2, status = 'old', file = resfile)
    do
      read(2, '(a)', iostat = istat) line
      if (istat == -1) exit
      if (istat > 0) call read_error(resfile, istat)
      if (line(4:4) < 'A' .or. line(4:4) > 'Z') cycle
      read(line(10:24), * ) ia, iz, isoT
      if (ia == 0) cycle
      if (iz == Z .and. ia == A .and. Liso == isoT) then
        if (type == 8) then
          xs = rochread(line(112:123))
          dxs = rochread(line(127:136))
        endif
        if (type == 9) then
          xs = rochread(line(61:72))
          dxs = rochread(line(76:84))
        endif
        if (xs > 0.) then
          k = k + 1
          res_author(k) = 'Sukhoruchkin'
          res_type(k) = 'Compilation'
          res_year(k) = 2015
          res_ref(k) = ref
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
! Nuclear data libraries
!
  xs = 0.
  dxs = 0.
  if (type == 8) then
    do lib = 1, numndlib
      resfile = trim(libspath)//trim(ndlib(lib))//'.RI'
      open (unit = 2, status = 'old', file = resfile)
      do
        read(2, '(a)', iostat = istat) line
        if (istat == -1) exit
        if (istat > 0) call read_error(resfile, istat)
        if (line(1:1) == '#') cycle
        read(line, * ) iz, ia, isoT, CE, chi2, xs
        if (iz == Z .and. ia == A .and. Liso == isoT) then
          k = k + 1
          res_author(k) = ndlib(lib)
          res_type(k) = 'NDL'
          res_year(k) = ndyear(lib)
          res_ref(k) = ref
          res_xs(k) = xs
          res_dxs(k) = dxs
          res_av(k) = ''
          res_exist = .true.
          exit
        endif
      enddo
      close (2)
    enddo
  endif
!
! Number of cases
!
  Nres = k
!
! EXFOR database
!
  if (type /= 2) then
    xs = 0.
    dxs = 0.
    if (type == 1 .or. type == 4 .or. type == 7) then
      nucstring='          '
      Zstring='   '
      write(Zstring,'(i3)') Z
      Astring='   '
      write(Astring,'(i3)') A        
      Unuc=nuc(Z)      
      if (Unuc(2:2) /= ' ') Unuc(2:2) = achar(iachar(Unuc(2:2)) - 32)
      write(nucstring,'(a,"-",a,"-",a)') trim(adjustl(Zstring)),trim(Unuc),trim(adjustl(Astring))
      if (Liso == 1) nucstring=trim(nucstring)//'-M'
      resfile = trim(filespath)//'exforfiles/n-0/'//trim(nucstring)//'.txt'
      inquire (file = resfile, exist = lexist)
      if (lexist .and. Liso <= 1) then
        open (unit = 2, status = 'old', file = resfile)
        do
          read(2, '(a)', iostat = istat) line
          if (istat == -1) exit
          if (istat > 0) call read_error(resfile, istat)
          if (line(1:1) == '#') cycle
          if (line(1:1) == ' ') cycle
          if (line(1:4) == 'RIPL') exit
          sub=line(1:12)
          author=line(21:44)
          read(line(45:160), *, iostat=istat ) year, Emin, Emax, xs, dxs, rL, rJ, rP
          if (isnan(xs)) xs = 0.
          if (isnan(dxs)) dxs = 0.
          if (isnan(rL)) rL = 0.
          if (type == 1 .and. rL /= 0.) cycle
          if (type == 4 .and. rL /= 1.) cycle
          if (type == 7 .and. rL /= 2.) cycle
          if (istat > 0) then
            write(*,*) "EXFOR problem: ",trim(line)
            cycle
          endif
          Nres_exp = Nres_exp + 1 
          if (Nres_exp  > numex) then
            Nres_exp = numex
            exit
          endif
          k = k + 1
          res_author(k) = author
          res_type(k) = 'EXFOR'
          res_year(k) = year
          res_ref(k) = sub
          res_xs(k) = xs
          res_dxs(k) = dxs
          if (Emax > 1.) write(*,*) "Warning: for Z=", Z," A=",A," resonance range given as ", Emin, " - ", Emax," MeV"
          if (.not.isnan(Emin)) write(res_Emin(k)(1:15),'(es15.6)') Emin
          if (.not.isnan(Emax)) write(res_Emax(k)(1:15),'(es15.6)') Emax
          if (.not.isnan(rL)) write(res_L(k)(7:10),'(f4.1)') rL
          if (.not.isnan(rJ)) write(res_J(k)(7:10),'(f4.1)') rJ
          if (.not.isnan(rP)) write(res_P(k)(7:10),'(f4.1)') rP
          res_av(k) = ''
          res_exist = .true.
        enddo
        close (2)
      endif
    else
      resfile = trim(exforpath)//'exfor_'//trim(restype(type))//'.txt'
      inquire (file = resfile, exist = lexist)
      if (lexist) then
        open (unit = 2, status = 'old', file = resfile)
        do
          read(2, '(a)', iostat = istat) line
          if (istat == -1) exit
          if (istat > 0) call read_error(resfile, istat)
          read(line, * ) iz, ia, isoT, isoR, xs ,dxs, Einc, author, year, sub
          if (iz == Z .and. ia == A .and. Liso == isoT) then
            Nres_exp = Nres_exp + 1
            if (Nres_exp  > numex) then
              Nres_exp = numex
              exit         
            endif
            k = k + 1      
            res_author(k) = author       
            res_type(k) = 'EXFOR'        
            res_year(k) = year
            res_ref(k) = sub
            res_xs(k) = xs 
            res_dxs(k) = dxs
            res_exist = .true.
          endif            
        enddo
        close (2)          
      endif
    endif
  endif
!
! resbase (TARES database)
!
  xs = 0.
  dxs = 0.
  resfile = trim(resbasepath)//trim(targetnuclide)//'/files/n-'//trim(targetnuclide)//'.txt'
  inquire (file = resfile, exist = lexist)
  if (lexist .and. type <= 6 .and. Nres > 0) then
    open (unit = 2, status = 'old', file = resfile)
    ix = 0
    do
      read(2, '(a)', iostat = istat) line
      if (istat == -1) exit
      if (istat > 0) call read_error(resfile, istat)
      if (type == 1) ix=index(line,'Ave. D        l=0')
      if (type == 2) ix=index(line,'Ave. S        l=0')
      if (type == 3) ix=index(line,'Ave. Gg       l=0')
      if (type == 4) ix=index(line,'Ave. D        l=1')
      if (type == 5) ix=index(line,'Ave. S        l=1')
      if (type == 6) ix=index(line,'Ave. Gg       l=1')
      if (ix > 0) then
        read(line(22:35),*) xs
        exit
      endif
    enddo
    close (2)
    if (xs > 0.) then
      k = k + 1
      res_author(k) = 'TARES'
      res_type(k) = 'NDL'
      res_year(k) = 2025
      res_ref(k) = ref
      if (type == 2 .or. type == 5) then
        res_xs(k) = 1.e4 * xs
        res_dxs(k) = 1.e4 * dxs
      else
        res_xs(k) = xs
        res_dxs(k) = dxs
      endif
      res_exist = .true.
    endif
  endif
  Nres = k
  return
end subroutine readresonance
! Copyright A.J. Koning 2025
