subroutine resonance
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose: Resonance databases
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
  character(len=3)   :: Astring    ! mass string
  integer            :: Z          ! charge number
  integer            :: A          ! mass number
  integer            :: Liso       ! isomer
  integer            :: type
!   
! type = 1: D0
! type = 2: S0
! type = 3: gamgam 
! type = 4: D1
! type = 5: S1
! type = 6: gamgam1
! type = 7: Ig
! type = 8: If
!     
  write(*, *) "Resonance databases"
  do type = 1, 8
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
            call readresonance(Z,A,Liso,type)
            if (res_exist) call procresonance(Z,A,Liso,type)
            if (res_exist) call writeresonance(Z,A,Liso,type)
          enddo
        endif
      enddo
    enddo
    if (Nsave > 0) call sumresonance(type)
  enddo
  return
end subroutine resonance
! Copyright A.J. Koning 2025
