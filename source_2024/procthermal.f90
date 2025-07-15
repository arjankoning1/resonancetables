subroutine procthermal
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose: Process thermal cross section data
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
  character(len=4)   :: ytmp          ! year
  character(len=9)   :: stmp          ! subentry
  character(len=24)  :: atmp          ! author
  character(len=40)  :: rtmp          ! reference
  integer            :: Z             ! charge number
  integer            :: A             ! mass number
  integer            :: N             ! counter
  integer            :: i             ! counter
  integer            :: j             ! counter
  integer            :: lib           ! library
  integer            :: isoT          ! target isomer
  integer            :: isoR          ! residual isomer
  integer            :: type          ! reaction type
  real(sgl)          :: xstmp         ! cross section
  real(sgl)          :: dxstmp        ! cross section uncertainty
!
! **************** Process databases for thermal cross sections *****
!
  write(*, *) "Processing data from databases....."
  do type = 1, numtype
    do Z = 1, numZ
      do A = 0, numA
        if (Z == 1 .and. A > 10) exit
        do isoT = 0, numisom
          do isoR = -1, numisom
!
! Sort EXFOR data according to year
!
            N = Ntherm_xs(type, Z, A, isoT, isoR)
            if (N > 0) then
              do i = 1, N
                do j = 1, i
                  if (Etherm_year(type, Z, A, isoT, isoR, i) > Etherm_year(type, Z, A, isoT, isoR, j)) cycle
                  xstmp = Etherm_xs(type, Z, A, isoT, isoR, i)
                  dxstmp = Etherm_dxs(type, Z, A, isoT, isoR, i)
                  ytmp = Etherm_year(type, Z, A, isoT, isoR, i)
                  atmp = Etherm_author(type, Z, A, isoT, isoR, i)
                  stmp = Etherm_subentry(type, Z, A, isoT, isoR, i)
                  rtmp = Etherm_ref(type, Z, A, isoT, isoR, i)
                  Etherm_xs(type, Z, A, isoT, isoR, i) = Etherm_xs(type, Z, A, isoT, isoR, j)
                  Etherm_dxs(type, Z, A, isoT, isoR, i) = Etherm_dxs(type, Z, A, isoT, isoR, j)
                  Etherm_year(type, Z, A, isoT, isoR, i) = Etherm_year(type, Z, A, isoT, isoR, j)
                  Etherm_author(type, Z, A, isoT, isoR, i) = Etherm_author(type, Z, A, isoT, isoR, j)
                  Etherm_subentry(type, Z, A, isoT, isoR, i) = Etherm_subentry(type, Z, A, isoT, isoR, j)
                  Etherm_ref(type, Z, A, isoT, isoR, i) = Etherm_ref(type, Z, A, isoT, isoR, j)
                  Etherm_xs(type, Z, A, isoT, isoR, j) = xstmp 
                  Etherm_dxs(type, Z, A, isoT, isoR, j) = dxstmp 
                  Etherm_year(type, Z, A, isoT, isoR, j) = ytmp 
                  Etherm_author(type, Z, A, isoT, isoR, j) = atmp 
                  Etherm_subentry(type, Z, A, isoT, isoR, j) = stmp 
                  Etherm_ref(type, Z, A, isoT, isoR, j) = rtmp
                enddo
              enddo
            endif
!
! Final dataset
! Rule: Kayzero > Mughabghab 2016 > Sukhoruchkin 2015 > Mughabghab 2006 > RIPL > EXFOR
!
! EXFOR
!
            if (N > 0) then
              therm_xs(0, type, Z, A, isoT, isoR) = Etherm_xs(type, Z, A, isoT, isoR, N)
              therm_dxs(0, type, Z, A, isoT, isoR) = Etherm_dxs(type, Z, A, isoT, isoR, N)
              therm_ref(0, type, Z, A, isoT, isoR) = trim(Etherm_ref(type, Z, A, isoT, isoR, N))
            endif
!
! RIPL
!
            if (therm_xs(1, type, Z, A, isoT, isoR) > 0.) then
              therm_xs(0, type, Z, A, isoT, isoR) = therm_xs(1, type, Z, A, isoT, isoR)
              therm_dxs(0, type, Z, A, isoT, isoR) = therm_dxs(1, type, Z, A, isoT, isoR)
              therm_ref(0, type, Z, A, isoT, isoR) = therm_ref(1, type, Z, A, isoT, isoR)
            endif
!
! Mughabghab 2006
!
            if (therm_xs(2, type, Z, A, isoT, isoR) > 0.) then
              therm_xs(0, type, Z, A, isoT, isoR) = therm_xs(2, type, Z, A, isoT, isoR)
              therm_dxs(0, type, Z, A, isoT, isoR) = therm_dxs(2, type, Z, A, isoT, isoR)
              therm_ref(0, type, Z, A, isoT, isoR) = therm_ref(2, type, Z, A, isoT, isoR)
            endif
!
! Sukhoruchkin 2015
!
            if (therm_xs(5, type, Z, A, isoT, isoR) > 0.) then
              therm_xs(0, type, Z, A, isoT, isoR) = therm_xs(5, type, Z, A, isoT, isoR)
              therm_dxs(0, type, Z, A, isoT, isoR) = therm_dxs(5, type, Z, A, isoT, isoR)
              therm_ref(0, type, Z, A, isoT, isoR) = therm_ref(5, type, Z, A, isoT, isoR)
            endif
!
! Mughabghab 2016
!
            if (therm_xs(3, type, Z, A, isoT, isoR) > 0.) then
              therm_xs(0, type, Z, A, isoT, isoR) = therm_xs(3, type, Z, A, isoT, isoR)
              therm_dxs(0, type, Z, A, isoT, isoR) = therm_dxs(3, type, Z, A, isoT, isoR)
              therm_ref(0, type, Z, A, isoT, isoR) = therm_ref(3, type, Z, A, isoT, isoR)
            endif
!
! Kayzero
!
            if (therm_xs(4, type, Z, A, isoT, isoR) > 0.) then
              therm_xs(0, type, Z, A, isoT, isoR) = therm_xs(4, type, Z, A, isoT, isoR)
              therm_dxs(0, type, Z, A, isoT, isoR) = therm_dxs(4, type, Z, A, isoT, isoR)
              therm_ref(0, type, Z, A, isoT, isoR) = therm_ref(4, type, Z, A, isoT, isoR)
            endif
!
! Ratios between databases
!
            if (therm_xs(0, type, Z, A, isoT, isoR) > 0.) then
              do lib = 1, 5
                ratio_xs(lib, type, Z, A, isoT, isoR) = therm_xs(lib, type, Z, A, isoT, isoR) / therm_xs(0, type, Z, A, isoT, isoR)
              enddo
              do lib = 1, numndlib
                Lratio_xs(lib, type, Z, A, isoT, isoR) = Ltherm_xs(lib, type, Z, A, isoT, isoR)/ therm_xs(0, type, Z, A, isoT, isoR)
              enddo
              if (N > 0) then
                do i = 1, N
                  Eratio_xs(type, Z, A, isoT, isoR, i) = Etherm_xs(type, Z, A, isoT, isoR, i) / therm_xs(0, type, Z, A, isoT, isoR) 
                enddo
              endif
            endif
          enddo
        enddo
      enddo
    enddo
  enddo
!
! Ratios of other thermal cross section to thermal (n,g)
!
  do type = 1, 6
    if (type == 4) cycle
    do Z = 1, numZ
      do A = 0, numA
        if (Z == 1 .and. A > 10) exit
        do isoT = 0, numisom
          do isoR = -1, numisom
            if (therm_xs(0, 4, Z, A, isoT, isoR) > 0.) &
 &            capratio_xs(type, Z, A, isoT, isoR) = therm_xs(0, type, Z, A, isoT, isoR) / therm_xs(0, 4, Z, A, isoT, isoR)
          enddo
        enddo
      enddo
    enddo
  enddo
  return
end subroutine procthermal
! Copyright A.J. Koning 2020
