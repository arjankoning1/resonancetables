subroutine thermal
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose: Thermal cross section databases
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
  integer            :: Liso       ! target isomer
  integer            :: Riso       ! residual isomer
  integer            :: type       ! reaction type
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
    do type = 1, 9
      if ((type <= 3 .or. (type >= 7 .and. type <= 9)) .and. Riso > -1) cycle
      Nsave = 0
      do Z = 1, numZ
        do A = 0, heavy(Z) + 3
          if (A == 0 .or. A >= light(Z) - 3) then
            do Liso = 0, numisom
              res_exist = .false.
              Astring='   '
              write(Astring(1:3),'(i3.3)') A
              targetnuclide=trim(nuc(Z))//Astring
              if (Liso == 1) targetnuclide = trim(targetnuclide)//'m'
              if (Liso == 2) targetnuclide = trim(targetnuclide)//'n'
              call readthermal(Z, A, Liso, Riso, type)
              if (res_exist) call procthermal(Z, A, Liso)
              if (res_exist) call writethermal(Z, A, Liso, Riso, type)
            enddo
          endif
        enddo
      enddo
      if (Nsave > 0) call sumthermal(Riso, type)
    enddo
  enddo
  return
end subroutine thermal
! Copyright A.J. Koning 2025
