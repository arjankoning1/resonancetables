function rochread(string)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose: Read a number from Rochman's file
!
! Revision    Date      Author      Quality  Description
! ======================================================
!    1     01-01-2017   A.J. Koning    A     Original code
!-----------------------------------------------------------------------------------------------------------------------------------
!
  use A0_kinds_mod
!
! *** Declaration of local data
!
  implicit none
  character(len=*)  :: string    ! line with parameter value
  integer           :: k         ! help variable
  integer           :: i         ! help variable
  integer           :: istat     ! help variable
  integer           :: L         ! help variable
  real(sgl)         :: x         ! help variable
  real(sgl)         :: rochread
!
  L = len_trim(string)
  if (L == 0) then
    x = 0.
  else
    k = index(string,'<')
    if (k > 0) string(k:k) = ' '
    i = index(string,'.')
    if (i == 0) i = index(string,'e')
    if (i > 0) then
      read(string, *, iostat = istat) x
      if (istat > 0) x = 0.
    else
      read(string, *, iostat = istat) k
      x = real(k)
      if (istat > 0) x = 0.
    endif
  endif
  rochread = x
  return
end function rochread
! Copyright A.J. Koning 2020
