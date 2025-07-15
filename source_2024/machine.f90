subroutine machine
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose: Machine dependent statements
!
! Revision    Date      Author      Quality  Description
! ======================================================
!    1     20-02-2023   A.J. Koning    A     Original code
!-----------------------------------------------------------------------------------------------------------------------------------
!
! *** Use data from other modules
!
  use A0_resonancetables_mod
!
! *** Declaration of local data
!
  implicit none
  character(len=132) :: codedir     ! code directory
!
! ************************ Set directories *****************************
!
  codedir = '/Users/koning/resonancetables/'
  filespath = trim(codedir)//'files/'
  libspath = trim(codedir)//'libs/'
  thermalpath = 'thermal/'
  macspath = 'macs/'
  respath = 'resonance/'
  return
end subroutine machine
! Copyright A.J. Koning 2023
