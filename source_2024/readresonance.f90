subroutine readresonance
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose: Read data from resonance databases
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
  character(len=2)   :: cdum       !
  character(len=9)   :: ref        ! reference
  character(len=132) :: resfile    ! thermfile
  character(len=1200):: line       ! input line
  integer            :: Z          ! charge number
  integer            :: A          ! mass number
  integer            :: isoT       ! target isomer
  integer            :: isoR       ! residual isomer
  integer            :: lib        ! library number
  integer            :: type       ! reaction type
  integer            :: istat      ! error code
  real(sgl)          :: D0        
  real(sgl)          :: dD0        
  real(sgl)          :: S0        
  real(sgl)          :: dS0        
  real(sgl)          :: gamgam        
  real(sgl)          :: dgamgam        
  real(sgl)          :: RI        
  real(sgl)          :: dRI        
  real(sgl)          :: rochread
!
! **************** Read RIPL-3 database for average resonance parameters *****
!
! type = 1: D0
! type = 2: S0
! type = 3: gamgam
! type = 4: Ig
! type = 5: If
!
  res_R = 0.
  res_dR = 0.
  res_ref = ''
  res_exist = .false.
!
! RIPL-3 resonance database
!
  write(*, *) "Reading data from resonance databases....."
  lib = 1
  resfile = trim(filespath)//'resonance.ripl'
  open (unit = 1, status = 'old', file = resfile)
  do
    read(1, '(a)', iostat = istat) line
    if (istat == -1) exit
    if (istat > 0) call read_error(resfile, istat)
    if (line(1:1) == '#') cycle
    read(line, * ) Z, cdum, A
    D0 = rochread(line(25:34))
    dD0 = rochread(line(36:45))
    res_R(lib, 1, Z, A) = D0
    res_dR(lib, 1, Z, A) = dD0
    S0 = rochread(line(46:51))
    dS0 = rochread(line(53:57))
    res_R(lib, 2, Z, A) = S0
    res_dR(lib, 2, Z, A) = dS0
    gamgam = rochread(line(59:63))
    dgamgam = rochread(line(65:69))
    res_R(lib, 3, Z, A) = gamgam
    res_dR(lib, 3, Z, A) = dgamgam
    ref = line(70:73)
    res_ref(lib, 1, Z, A) = ref
    res_ref(lib, 2, Z, A) = ref
    res_ref(lib, 3, Z, A) = ref
    res_exist(lib, 1) = .true.
    res_exist(lib, 2) = .true.
    res_exist(lib, 3) = .true.
  enddo
  close (1)
!
! JUKO database for resonance integrals
!
  resfile = trim(filespath)//'resint.juko'
  open (unit = 1, status = 'old', file = resfile)
  do
    read(1, '(a)', iostat = istat) line
    if (istat == -1) exit
    if (istat > 0) call read_error(resfile, istat)
    if (line(1:1) == '#') cycle
    read(line, * ) A, Z
    RI = rochread(line(49:56))
    dRI = rochread(line(66:73))
    type = 4
    if (line(37:41) == '(n,f)') type = 5
    res_R(lib, type, Z, A) = RI
    res_dR(lib, type, Z, A) = dRI
    res_ref(lib, type, Z, A) = 'JUKO'
    res_exist(lib, type) = .true.
  enddo
!
! Mughabghab 2016 database
!
  lib = 3
  resfile = trim(filespath)//'thermal.mugh18'
  open (unit = 1, status = 'old', file = resfile)
  do
    read(1, '(a)', iostat = istat) line
    if (istat == -1) exit
    if (istat > 0) call read_error(resfile, istat)
    if (line(4:4) < 'A' .or. line(4:4) > 'Z') cycle
    read(line(10:24), * ) A, Z, isoT
    if (A == 0) cycle
    isoR = -1
    type = 1
    D0 = rochread(line(1031:1041))
    dD0 = rochread(line(1045:1052))
    if (D0 > 0.) then
      res_R(lib, type, Z, A) = 0.001 * D0
      res_dR(lib, type, Z, A) = 0.001 * dD0
      res_ref(lib, type, Z, A) = 'Mugh18'
      res_exist(lib, type) = .true.
    endif
    type = 2
    S0 = rochread(line(1055:1065))
    dS0 = rochread(line(1069:1076))
    if (S0 > 0.) then
      res_R(lib, type, Z, A) = S0
      res_dR(lib, type, Z, A) = dS0
      res_ref(lib, type, Z, A) = 'Mugh18'
      res_exist(lib, type) = .true.
    endif
    type = 3
    gamgam = rochread(line(1175:1185))
    dgamgam = rochread(line(1189:1196))
    if (gamgam > 0.) then
      res_R(lib, type, Z, A) = 1000. * gamgam
      res_dR(lib, type, Z, A) = 1000. * dgamgam
      res_ref(lib, type, Z, A) = 'Mugh18'
      res_exist(lib, type) = .true.
    endif
    type = 4
    RI = rochread(line(607:615))
    dRI = rochread(line(619:626))
    if (RI > 0.) then
      res_R(lib, type, Z, A) = RI
      res_dR(lib, type, Z, A) = dRI
      res_ref(lib, type, Z, A) = 'Mugh18'
      res_exist(lib, type) = .true.
    endif
    type = 5
    RI = rochread(line(357:364))
    dRI = rochread(line(368:377))
    if (RI > 0.) then
      res_R(lib, type, Z, A) = RI
      res_dR(lib, type, Z, A) = dRI
      res_ref(lib, type, Z, A) = 'Mugh18'
      res_exist(lib, type) = .true.
    endif
  enddo
  close (1)
!
! Kayzero database
!
  lib = 4
  resfile = trim(filespath)//'kayzero.txt'
  open (unit = 1, status = 'old', file = resfile)
  do
    read(1, '(a)', iostat = istat) line
    if (istat == -1) exit
    if (istat > 0) call read_error(resfile, istat)
    if (line(4:4) < 'A' .or. line(4:4) > 'Z') cycle
    read(line(10:24), * ) A, Z, isoT
    if (A == 0) cycle
    isoR = -1
    type = 4
    RI = rochread(line(112:122))
    dRI = rochread(line(127:136))
    if (RI > 0.) then
      res_R(lib, type, Z, A) = RI
      res_dR(lib, type, Z, A) = dRI
      res_ref(lib, type, Z, A) = 'Kayzero'
      res_exist(lib, type) = .true.
    endif
  enddo
  close (1)
!
! Sukhoruchkin database
!
  lib = 5
  resfile = trim(filespath)//'sukhoruchkin.txt'
  open (unit = 1, status = 'old', file = resfile)
  do
    read(1, '(a)', iostat = istat) line
    if (istat == -1) exit
    if (istat > 0) call read_error(resfile, istat)
    if (line(4:4) < 'A' .or. line(4:4) > 'Z') cycle
    read(line(10:24), * ) A, Z, isoT
    if (A == 0) cycle
    isoR = -1
    type = 4
    RI = rochread(line(112:123))
    dRI = rochread(line(127:136))
    if (RI > 0.) then
      res_R(lib, type, Z, A) = RI
      res_dR(lib, type, Z, A) = dRI
      res_ref(lib, type, Z, A) = 'Sukhoruchkin'
      res_exist(lib, type) = .true.
    endif
    type = 5
    RI = rochread(line(61:72))
    dRI = rochread(line(76:84))
    if (RI > 0.) then
      res_R(lib, type, Z, A) = RI
      res_dR(lib, type, Z, A) = dRI
      res_ref(lib, type, Z, A) = 'Sukhoruchkin'
      res_exist(lib, type) = .true.
    endif
  enddo
  close (1)
  return
end subroutine readresonance
! Copyright A.J. Koning 2019
