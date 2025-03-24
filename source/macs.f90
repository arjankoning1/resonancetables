subroutine macs
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose: MACS databases
!
! Revision    Date      Author      Quality  Description
! ======================================================
!    1     2025-03-05   A.J. Koning    A     Original code
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
  character(len=3)   :: Astring    ! mass string
  integer            :: Z          ! charge number
  integer            :: A          ! mass number
  integer            :: Liso       ! isomer
  integer            :: Riso       ! isomer
!
!
!
  write(*, *) "MACS databases"
  Riso = -1
  Nsave = 0
  do Z = 1, numZ
    do A = 0, heavy(Z) + 5
      if (A == 0 .or. A >= light(Z) - 5) then
        do Liso = 0, numisom
          res_exist = .false.
          Astring='   '
          write(Astring(1:3),'(i3.3)') A
          targetnuclide=trim(nuc(Z))//Astring
          if (Liso == 1) targetnuclide = trim(targetnuclide)//'m'
          if (Liso == 2) targetnuclide = trim(targetnuclide)//'n'
          call readmacs(Z,A,Liso,Riso)
          if (res_exist) call procmacs(Z,A,Liso)
          if (res_exist) call writemacs(Z,A,Liso,Riso)
        enddo
      endif
    enddo
  enddo
  if (Nsave > 0) call summacs
  return
end subroutine macs
! Copyright A.J. Koning 2025
