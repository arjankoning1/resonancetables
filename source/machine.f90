subroutine machine
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose: Machine dependent statements
!
! Revision    Date      Author      Quality  Description
! ======================================================
!    1     2025-02-25   A.J. Koning    A     Original code
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
  character(len=132) :: basedir     ! base directory
  integer :: ix
  integer           :: year    ! year
  integer           :: month   ! month
  integer           :: day     ! day
  integer           :: values(8) ! date and time values
!
! ************************ Set directories *****************************
!
  codedir = '/Users/koning/resonancetables/'
  ix = index(codedir,'/resonancetables/')
  basedir = codedir(1:ix)
  exforpath = trim(basedir)//'exfortables/special/'
  resbasepath = trim(basedir)//'libraries/resbase/'
  filespath = trim(codedir)//'files/'
  libspath = trim(codedir)//'libs/'
  thermalpath = 'thermal/'
  macspath = 'macs/'
  respath = 'resonance/'
  resparpath = trim(basedir)//'tendl/'
!
! Set date
!
  call date_and_time(VALUES=values)
  year=values(1)
  month=values(2)
  day=values(3)
  date='xxxx-xx-xx'
  write(date(1:4),'(i4.4)') year
  write(date(6:7),'(i2.2)') month
  write(date(9:10),'(i2.2)') day
  return
end subroutine machine
! Copyright A.J. Koning 2025
