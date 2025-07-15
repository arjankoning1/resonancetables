subroutine resonancetablesinitial
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose: Initialization
!
! Revision    Date      Author      Quality  Description
! ======================================================
!    1     01-08-2020   A.J. Koning    A     Original code
!-----------------------------------------------------------------------------------------------------------------------------------
!
! *** Use data from other modules
!
  use A0_resonancetables_mod
!
! *** Declaration of local data
!
  implicit none
  character(len=80) :: cmd       ! command
  integer           :: isys      ! counter
  integer           :: i         ! counter
  integer           :: system    ! system command
!
! **************** Initialization of constants and arrays **************
!
  nuc =  (/ 'H ', 'He', 'Li', 'Be', 'B ', 'C ', 'N ', 'O ', 'F ', 'Ne', &
    'Na', 'Mg', 'Al', 'Si', 'P ', 'S ', 'Cl', 'Ar', 'K ', 'Ca', 'Sc', 'Ti', 'V ', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', &
    'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y ', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', &
    'Sb', 'Te', 'I ', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', &
    'Lu', 'Hf', 'Ta', 'W ', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th', &
    'Pa', 'U ', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm'/)
  reac(1) = 'tot'
  reac(2) = 'el'
  reac(3) = 'nf'
  reac(4) = 'ng'
  reac(5) = 'np'
  reac(6) = 'na'
  reac(7) = 'nu'
  reac(8) = 'nud'
  reac(9) = 'nup'
  restype(1) = 'D0'
  restype(2) = 'S0'
  restype(3) = 'gamgam'
  restype(4) = 'Ig'
  restype(5) = 'If'
  ext(0) = 'final'
  ext(1) = 'ripl'
  ext(2) = 'mugh06'
  ext(3) = 'mugh18'
  ext(4) = 'kayzero'
  ext(5) = 'sukhoruchkin'
  ext(6) = 'exfor '
  ndlib(1) = 'cendl3.2  '
  ndlib(2) = 'endfb8.1  '
  ndlib(3) = 'jeff3.3   '
  ndlib(4) = 'jendl5.0  '
  ndlib(5) = 'tendl.2023'
!
! Cleanup of previous results
!
  write(*, *) "Removing directories from previous run....."
  cmd = 'rm -r '//trim(thermalpath)
  isys = system(cmd)
  cmd = 'mkdir -p '//trim(thermalpath)//'nuc'
  isys = system(cmd)
  cmd = 'rm -r '//trim(macspath)
  isys = system(cmd)
  cmd = 'mkdir -p '//trim(macspath)//'nuc'
  isys = system(cmd)
  cmd = 'rm -r '//trim(respath)
  isys = system(cmd)
  cmd = 'mkdir '//trim(respath)
  isys = system(cmd)
  do i = 1, numtype
    cmd = 'mkdir '//trim(thermalpath)//reac(i)
    isys = system(cmd)
    if (i == 4) then
      cmd = 'mkdir '//trim(macspath)//reac(i)
      isys = system(cmd)
    endif
    if (i <= 5) then
      cmd = 'mkdir '//trim(respath)//restype(i)
      isys = system(cmd)
    endif
  enddo
  return
end subroutine resonancetablesinitial
! Copyright A.J. Koning 2020
