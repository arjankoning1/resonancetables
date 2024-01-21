subroutine procmacs
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose: Process MACS cross section data
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
! **************** Process databases for MACS cross sections *****
!
  write(*, *) "Processing data from MACS databases....."
  do type = 1, numtype
    if (type /= 4) cycle
    do Z = 1, numZ
      do A = 0, numA
        if (Z == 1 .and. A > 10) exit
        do isoT = 0, numisom
          do isoR = -1, numisom
!
! Sort EXFOR data according to year
!
            N = Nmacs_xs(type, Z, A, isoT, isoR)
            if (N > 0) then
              do i = 1, N
                do j = 1, i
                  if (Emacs_year(type, Z, A, isoT, isoR, i) > Emacs_year(type, Z, A, isoT, isoR, j)) cycle
                  xstmp = Emacs_xs(type, Z, A, isoT, isoR, i)
                  dxstmp = Emacs_dxs(type, Z, A, isoT, isoR, i)
                  ytmp = Emacs_year(type, Z, A, isoT, isoR, i)
                  atmp = Emacs_author(type, Z, A, isoT, isoR, i)
                  stmp = Emacs_subentry(type, Z, A, isoT, isoR, i)
                  rtmp = Emacs_ref(type, Z, A, isoT, isoR, i)
                  Emacs_xs(type, Z, A, isoT, isoR, i) = Emacs_xs(type, Z, A, isoT, isoR, j)
                  Emacs_dxs(type, Z, A, isoT, isoR, i) = Emacs_dxs(type, Z, A, isoT, isoR, j)
                  Emacs_year(type, Z, A, isoT, isoR, i) = Emacs_year(type, Z, A, isoT, isoR, j)
                  Emacs_author(type, Z, A, isoT, isoR, i) = Emacs_author(type, Z, A, isoT, isoR, j)
                  Emacs_subentry(type, Z, A, isoT, isoR, i) = Emacs_subentry(type, Z, A, isoT, isoR, j)
                  Emacs_ref(type, Z, A, isoT, isoR, i) = Emacs_ref(type, Z, A, isoT, isoR, j)
                  Emacs_xs(type, Z, A, isoT, isoR, j) = xstmp 
                  Emacs_dxs(type, Z, A, isoT, isoR, j) = dxstmp 
                  Emacs_year(type, Z, A, isoT, isoR, j) = ytmp 
                  Emacs_author(type, Z, A, isoT, isoR, j) = atmp 
                  Emacs_subentry(type, Z, A, isoT, isoR, j) = stmp 
                  Emacs_ref(type, Z, A, isoT, isoR, j) = rtmp
                enddo
              enddo
            endif
!
! Final dataset
! Rule: ASTRAL > KADONIS> Mughabghab 2016 > Sukhoruchkin 2015 >  EXFOR
!
! EXFOR
!
            if (N > 0) then
              macs_xs(0, type, Z, A, isoT, isoR) = Emacs_xs(type, Z, A, isoT, isoR, N)
              macs_dxs(0, type, Z, A, isoT, isoR) = Emacs_dxs(type, Z, A, isoT, isoR, N)
              macs_ref(0, type, Z, A, isoT, isoR) = trim(Emacs_ref(type, Z, A, isoT, isoR, N))
            endif
!
! Sukhoruchkin 2015
!
            if (macs_xs(5, type, Z, A, isoT, isoR) > 0.) then
              macs_xs(0, type, Z, A, isoT, isoR) = macs_xs(5, type, Z, A, isoT, isoR)
              macs_dxs(0, type, Z, A, isoT, isoR) = macs_dxs(5, type, Z, A, isoT, isoR)
              macs_ref(0, type, Z, A, isoT, isoR) = macs_ref(5, type, Z, A, isoT, isoR)
            endif
!
! Mughabghab 2016
!
            if (macs_xs(3, type, Z, A, isoT, isoR) > 0.) then
              macs_xs(0, type, Z, A, isoT, isoR) = macs_xs(3, type, Z, A, isoT, isoR)
              macs_dxs(0, type, Z, A, isoT, isoR) = macs_dxs(3, type, Z, A, isoT, isoR)
              macs_ref(0, type, Z, A, isoT, isoR) = macs_ref(3, type, Z, A, isoT, isoR)
            endif
!
! KADONIS
!
            if (macs_xs(1, type, Z, A, isoT, isoR) > 0.) then
              macs_xs(0, type, Z, A, isoT, isoR) = macs_xs(1, type, Z, A, isoT, isoR)
              macs_dxs(0, type, Z, A, isoT, isoR) = macs_dxs(1, type, Z, A, isoT, isoR)
              macs_ref(0, type, Z, A, isoT, isoR) = macs_ref(1, type, Z, A, isoT, isoR)
            endif
!
! ASTRAL
!
            if (macs_xs(2, type, Z, A, isoT, isoR) > 0.) then
              macs_xs(0, type, Z, A, isoT, isoR) = macs_xs(2, type, Z, A, isoT, isoR)
              macs_dxs(0, type, Z, A, isoT, isoR) = macs_dxs(2, type, Z, A, isoT, isoR)
              macs_ref(0, type, Z, A, isoT, isoR) = macs_ref(2, type, Z, A, isoT, isoR)
            endif
!
! Ratios between databases
!
            if (macs_xs(0, type, Z, A, isoT, isoR) > 0.) then
              do lib = 1,5
               ratiomacs_xs(lib, type, Z, A, isoT, isoR) = macs_xs(lib, type, Z, A, isoT, isoR) / macs_xs(0, type, Z, A, isoT, isoR)
              enddo
              do lib = 1, numndlib
                Lratiomacs_xs(lib, type, Z, A, isoT, isoR) =Lmacs_xs(lib, type, Z, A, isoT, isoR)/macs_xs(0, type, Z, A, isoT, isoR)
              enddo
              if (N > 0) then
                do i = 1, N
                  Eratiomacs_xs(type, Z, A, isoT, isoR, i) = Emacs_xs(type, Z, A, isoT, isoR, i)/macs_xs(0, type, Z, A, isoT, isoR) 
                enddo
              endif
            endif
          enddo
        enddo
      enddo
    enddo
  enddo
  return
end subroutine procmacs
! Copyright A.J. Koning 2020
