subroutine resonancetablesinitial
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose: Initialization
!
! Revision    Date      Author      Quality  Description
! ======================================================
!    1     2025-09-14   A.J. Koning    A     Original code
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
  light = (/ 1,   3,   6,   9,  10,  12,  14,  16,  19,  20, &
     23,   24,   27,   28,   31,   32,   35,   36,   39,   40,   45,   46,   50,   50,   55,   54,   59,   58,   63,   64,  &
     69,   70,   75,   74,   79,   78,   85,   84,   89,   90,   93,   92,   99,   96,  103,  102,  107,  106,  113,  112,  &
    121,  120,  127,  124,  133,  130,  138,  136,  141,  142,  147,  144,  151,  152,  159,  156,  165,  162,  169,  168, &
    175,  174,  180,  180,  185,  184,  191,  190,  197,  196,  203,  204,  209,  208,  208,  211,  212,  223,  225,  228, &
    229,  232,  235,  236,  239,  240,  242,  245,  248,  252 /)
  heavy = (/ 2,    4,    7,    9,   11,   13,   15,   18,   19,   22, &
     23,   26,   27,   30,   31,   36,   37,   40,   41,   48,   45,   50,   51,   54,   55,   58,   59,   64,   65,   70, &
     71,   76,   75,   82,   81,   86,   87,   88,   89,   96,   93,  100,   99,  104,  103,  110,  109,  116,  115,  124, &
    123,  130,  127,  136,  133,  138,  139,  142,  141,  150,  149,  154,  153,  160,  159,  164,  165,  170,  169,  176, &
    176,  180,  181,  186,  187,  192,  193,  198,  197,  204,  205,  208,  209,  210,  211,  222,  222,  226,  227,  233, &
    233,  238,  238,  242,  244,  250,  250,  253,  256,  258 /)
  reaction(1) = '(n,tot)'
  reaction(2) = '(n,el)'
  reaction(3) = '(n,f)'
  reaction(4) = '(n,g)'
  reaction(5) = '(n,p)'
  reaction(6) = '(n,a)'
  reaction(7) = '(n,nubar)'
  reaction(8) = '(n,nudel)'
  reaction(9) = '(n,nuprompt)'
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
  restype(4) = 'D1'
  restype(5) = 'S1'
  restype(6) = 'gamgam1'
  restype(7) = 'D2'
  restype(8) = 'R'
  restype(9) = 'Ig'
  restype(10) = 'If'
  ndlib(1) = 'cendl3.2  '
  ndyear(1) = 2019
  ndlib(2) = 'endfb8.1  '
  ndyear(2) = 2024
  ndlib(3) = 'jeff4.0   '
  ndyear(3) = 2025
  ndlib(4) = 'jendl5.0  '
  ndyear(4) = 2021
  ndlib(5) = 'tendl.2025'
  ndyear(5) = 2025
  user = ''
! user = 'Arjan Koning'
  source = 'Resonancetables'
! oformat = 'YANDF-0.4'
  oformat = ''
!
! Cleanup of previous results
!
  write(*, *) "Removing directories from previous run....."
  cmd = 'rm -r '//trim(thermalpath)
  isys = system(cmd)
  cmd = 'rm -r '//trim(macspath)
  isys = system(cmd)
  cmd = 'rm -r '//trim(respath)
  isys = system(cmd)
  do i = 1, numtype
    cmd = 'mkdir -p '//trim(thermalpath)//trim(reac(i))//'/nuc'
    isys = system(cmd)
    cmd = 'mkdir -p '//trim(thermalpath)//trim(reac(i))//'/all'
    isys = system(cmd)
    if (i == 4 .or. i == 6) then
      cmd = 'mkdir -p '//trim(thermalpath)//trim(reac(i))//'-g/nuc'
      isys = system(cmd)
      cmd = 'mkdir -p '//trim(thermalpath)//trim(reac(i))//'-g/all'
      isys = system(cmd)
      cmd = 'mkdir -p '//trim(thermalpath)//trim(reac(i))//'-m/nuc'
      isys = system(cmd)
      cmd = 'mkdir -p '//trim(thermalpath)//trim(reac(i))//'-m/all'
      isys = system(cmd)
      if (i == 4) then
        cmd = 'mkdir -p '//trim(thermalpath)//trim(reac(i))//'-n/nuc'
        isys = system(cmd)
        cmd = 'mkdir -p '//trim(thermalpath)//trim(reac(i))//'-n/all'
        isys = system(cmd)
      endif
    endif
    if (i == 4) then
      cmd = 'mkdir -p '//trim(macspath)//trim(reac(i))//'/nuc'
      isys = system(cmd)
      cmd = 'mkdir -p '//trim(macspath)//trim(reac(i))//'/all'
      isys = system(cmd)
    endif
    cmd = 'mkdir -p '//trim(respath)//trim(restype(i))//'/nuc'
    isys = system(cmd)
    cmd = 'mkdir -p '//trim(respath)//trim(restype(i))//'/all'
    isys = system(cmd)
  enddo
  return
end subroutine resonancetablesinitial
! Copyright A.J. Koning 2025
