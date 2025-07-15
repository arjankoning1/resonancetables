subroutine procresonance
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose: Process average resonance data
!
! Revision    Date      Author      Quality  Description
! ======================================================
!    1     01-08-2020   A.J. Koning    A     Original code
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
  integer            :: Z             ! charge number
  integer            :: A             ! mass number
  integer            :: lib           ! library
  integer            :: type          ! reaction type
!
! **************** Process databases for average resonance data *****
!
  write(*, *) "Processing data from resonance databases....."
  do type = 1, 5
    do Z = 1, numZ
      do A = 0, numA
        if (Z == 1 .and. A > 10) exit
!
! Final dataset
! Rule D0: RIPL > Mughabghab 2016 > EXFOR
! Rule S0, gamgam: Mughabghab 2016 >   RIPL > EXFOR
! Rule Ig, If: Kayzero > Mughabghab 2016 > Sukhoruchkin 2015 >  JUKO > EXFOR
!
        if (type == 1) then
!
! Mughabghab 2016
!
          if (res_R(3, type, Z, A) > 0.) then
            res_R(0, type, Z, A) = res_R(3, type, Z, A)
            res_dR(0, type, Z, A) = res_dR(3, type, Z, A)
            res_ref(0, type, Z, A) = res_ref(3, type, Z, A)
          endif
!
! RIPL
!
          if (res_R(1, type, Z, A) > 0.) then
            res_R(0, type, Z, A) = res_R(1, type, Z, A)
            res_dR(0, type, Z, A) = res_dR(1, type, Z, A)
            res_ref(0, type, Z, A) = res_ref(1, type, Z, A)
          endif
        endif
        if (type == 2 .or. type == 3) then
!
! Mughabghab 2016
!
          if (res_R(3, type, Z, A) > 0.) then
            res_R(0, type, Z, A) = res_R(3, type, Z, A)
            res_dR(0, type, Z, A) = res_dR(3, type, Z, A)
            res_ref(0, type, Z, A) = res_ref(3, type, Z, A)
          endif
!
! RIPL
!
          if (res_R(1, type, Z, A) > 0.) then
            res_R(0, type, Z, A) = res_R(1, type, Z, A)
            res_dR(0, type, Z, A) = res_dR(1, type, Z, A)
            res_ref(0, type, Z, A) = res_ref(1, type, Z, A)
          endif
        endif
        if (type == 4 .or. type == 5) then
!
! JUKO
!
          if (res_R(1, type, Z, A) > 0.) then
            res_R(0, type, Z, A) = res_R(1, type, Z, A)
            res_dR(0, type, Z, A) = res_dR(1, type, Z, A)
            res_ref(0, type, Z, A) = res_ref(1, type, Z, A)
          endif
!
! Sukhoruchkin 2015
!
          if (res_R(5, type, Z, A) > 0.) then
            res_R(0, type, Z, A) = res_R(5, type, Z, A)
            res_dR(0, type, Z, A) = res_dR(5, type, Z, A)
            res_ref(0, type, Z, A) = res_ref(5, type, Z, A)
          endif
!
! Mughabghab 2016
!
          if (res_R(3, type, Z, A) > 0.) then
            res_R(0, type, Z, A) = res_R(3, type, Z, A)
            res_dR(0, type, Z, A) = res_dR(3, type, Z, A)
            res_ref(0, type, Z, A) = res_ref(3, type, Z, A)
          endif
!
! Kayzero
!
          if (res_R(4, type, Z, A) > 0.) then
            res_R(0, type, Z, A) = res_R(4, type, Z, A)
            res_dR(0, type, Z, A) = res_dR(4, type, Z, A)
            res_ref(0, type, Z, A) = res_ref(4, type, Z, A)
          endif
        endif
!
! Ratios between databases
!
        if (res_R(0, type, Z, A) > 0.) then
          do lib = 1, 5
            ratio_R(lib, type, Z, A) = res_R(lib, type, Z, A) / res_R(0, type, Z, A)
          enddo
        endif
      enddo
    enddo
  enddo
  return
end subroutine procresonance
! Copyright A.J. Koning 2019
