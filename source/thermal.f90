subroutine thermal
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose: Thermal cross section databases
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
  logical            :: flagav     ! flag for spectrum average
  character(len=3)   :: Astring    ! mass string
  integer            :: Z          ! charge number
  integer            :: A          ! mass number
  integer            :: Liso       ! target isomer
  integer            :: Riso       ! residual isomer
  integer            :: type       ! reaction type
  integer            :: iav
!
! type = 1: (n,tot)
! type = 2: (n,el)
! type = 3: (n,f)
! type = 4: (n,g)
! type = 5: (n,p)
! type = 6: (n,a)
! type = 7: nubar total
! type = 8: nubar delayed
! type = 9: nubar prompt
!
  write(*, *) "Thermal databases"
  do Riso = -1, numisom
    do iav = 1, 2
      if (iav == 1) then
        flagav = .false.
      else
        flagav = .true.
      endif
      do type = 1, 9
        if ((type <= 3 .or. (type >= 7 .and. type <= 9)) .and. Riso > -1) cycle
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
                call readthermal(Z, A, Liso, Riso, type, flagav)
                if (res_exist) call procthermal(Z, A, Liso, Riso, type, flagav)
                if (res_exist) call writethermal(Z, A, Liso, Riso, type, flagav)
              enddo
            endif
          enddo
        enddo
        call sumthermal(Riso, type, flagav)
      enddo
    enddo
  enddo
  return
end subroutine thermal
! Copyright A.J. Koning 2025
