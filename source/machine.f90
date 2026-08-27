subroutine machine
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose: Machine dependent statements
!
! Revision    Date      Author      Quality  Description
! =====================================================
!    1     2025-02-25   A.J. Koning    A     Original code
!    2     2026-08-27   A.J. Koning    A     Runtime definition of RESONANCETABLES directory
!-----------------------------------------------------------------------------------------------------------------------------------
!
! *** Use data from other modules
!
  use A0_resonancetables_mod
!
! *** Declaration of local data
!
  implicit none
  logical             :: lexist
  character(len=1024) :: code_dir
  character(len=1024) :: base_dir
  character(len=1024) :: resonancetables_dir
  integer             :: envstat
  integer             :: i
  integer             :: n
  integer             :: year
  integer             :: month
  integer             :: day
  integer             :: values(8)
!
! ************************ Set directories *****************************
!
! The preferred option is to set the environment variable
! RESONANCETABLES_DIR, for example in ~/.profile or ~/.zshrc:
!
! export RESONANCETABLES_DIR=/path/to/resonancetables
!
  call get_environment_variable('RESONANCETABLES_DIR', resonancetables_dir, length=n, status=envstat)
  if (envstat == 0 .and. n > 0) then
    code_dir = trim(resonancetables_dir)
  else
!
! If the environment variable cannot be used, the code directory can be
! changed here manually.
!
    code_dir = '/path/to/resonancetables/'
  endif
!
! Remove a trailing slash, if present, and determine the parent directory.
!
  i = len_trim(code_dir)
  if (i > 1) then
    if (code_dir(i:i) == '/') code_dir = code_dir(:i - 1)
  endif
  i = scan(trim(code_dir), '/', back=.true.)
  if (i > 0) then
    base_dir = code_dir(:i)
  else
    base_dir = './'
  endif
!
! Internal RESONANCETABLES data.
!
  filespath = trim(code_dir)//'/files/'
  libspath = trim(code_dir)//'/libs/'
!
! External data are expected as sibling directories.
!
  exforpath = trim(base_dir)//'exfortables/special/'
  resbasepath = trim(base_dir)//'libraries/resbase/'
  resparpath = trim(base_dir)//'tendl/'
!
! Output directories are intentionally relative to the current working
! directory. RESONANCETABLES removes and recreates these directories when
! it starts.
!
  thermalpath = 'thermal/'
  macspath = 'macs/'
  respath = 'resonance/'
!
! Test accessibility of the internal RESONANCETABLES data files.
!
  inquire (file=trim(filespath)//'atlas_2006.txt', exist=lexist)
  if (.not. lexist) then
    write(*, '(a)') 'RESONANCETABLES error: data files not found.'
    write(*, '(2a)') 'Expected file: ', trim(filespath)//'atlas_2006.txt'
    write(*, '(a)') 'Set the RESONANCETABLES_DIR environment variable:'
    write(*, '(a)') '  export RESONANCETABLES_DIR=/path/to/resonancetables'
    write(*, '(a)') 'Alternatively, edit code_dir in source/machine.f90'
    write(*, '(a)') 'and rebuild RESONANCETABLES.'
    error stop 77
  endif
!
! ************************ Set date ***********************************
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
! Copyright A.J. Koning 2026
