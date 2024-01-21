module A0_kinds_mod
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose: Definition of single and double precision variables
!
! Revision    Date      Author           Description
! ====================================================
!    1     01-01-2019   A.J. Koning      Original code
!-----------------------------------------------------------------------------------------------------------------------------------
!
  integer, parameter :: sgl = selected_real_kind(6,37)   ! single precision kind
  integer, parameter :: dbl = selected_real_kind(15,307) ! double precision kind
end module A0_kinds_mod
! Copyright (C) 2019 A.J. Koning
