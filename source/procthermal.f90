subroutine procthermal(Z, A, Liso, Riso, type, flagav)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose: Process thermal cross section data
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
  character(len=40)  :: ttmp          ! subentry
  character(len=24)  :: atmp          ! author
  character(len=40)  :: rtmp          ! reference
  integer            :: Z             ! charge number
  integer            :: A             ! mass number
  integer            :: N             ! counter
  integer            :: i             ! counter
  integer            :: j             ! counter
  integer            :: Nsel          ! counter
  integer            :: ytmp          ! year
  integer            :: Liso          ! target isomer
  integer            :: Riso          ! residual isomer
  integer            :: type          ! reaction type
  real(sgl)          :: xstmp         ! cross section
  real(sgl)          :: dxstmp        ! cross section uncertainty
!
! **************** Process databases for thermal cross sections *****
!
! Sort data according to year
!
  N = Nres
  if (N > 0) then
    do i = 1, N
      do j = 1, i
        if (res_year(i) > res_year(j)) cycle
        xstmp = res_xs(i)
        dxstmp = res_dxs(i)
        ytmp = res_year(i)
        atmp = res_author(i)
        ttmp = res_type(i)
        rtmp = res_ref(i)
        res_xs(i) = res_xs(j)
        res_dxs(i) = res_dxs(j)
        res_year(i) = res_year(j)
        res_author(i) = res_author(j)
        res_type(i) = res_type(j)
        res_ref(i) = res_ref(j)
        res_xs(j) = xstmp 
        res_dxs(j) = dxstmp 
        res_year(j) = ytmp 
        res_author(j) = atmp 
        res_type(j) = ttmp 
        res_ref(j) = rtmp
      enddo
    enddo
  endif
!
! Final dataset
! Rule for MXW average: Mughabghab 2016 > Mughabghab 2006 > Firestone 2022 > EXFOR
!
  Nsel = 0
  N = Nres
  if (flagav) then
    Loop1:  do
      do i = 1, N
        if (res_author(i) == 'Mughabghab_2016') then
          Nsel = i
          exit Loop1
        endif
      enddo
      do i = 1, N
        if (res_author(i) == 'Mughabghab_2006') then
          Nsel = i
          exit Loop1
        endif
      enddo
      do i = 1, N
        if (res_author(i) == 'Firestone') then
          Nsel = i
          exit Loop1
        endif
      enddo
      Nsel = Nres_exp
      exit Loop1
    enddo Loop1
  else
!
! Final dataset
! Rule for xs: Kayzero > Mughabghab 2016 > Sukhoruchkin 2015 > Mughabghab 2006 > RIPL-3 > EXFOR
!
    Loop2:  do
      do i = 1, N
        if (res_author(i) == 'Kayzero') then
          Nsel = i
          exit Loop2
        endif
      enddo
      do i = 1, N
        if (res_author(i) == 'Mughabghab_2016') then
          Nsel = i
          exit Loop2
        endif
      enddo
      do i = 1, N
        if (res_author(i) == 'Sukhoruchkin') then
          Nsel = i
          exit Loop2
        endif
      enddo
      do i = 1, N
        if (res_author(i) == 'Mughabghab_2006') then
          Nsel = i
          exit Loop2
        endif
      enddo
      do i = 1, N
        if (res_author(i) == 'RIPL-3') then
          Nsel = i
          exit Loop2
        endif
      enddo
      Nsel = Nres_exp
      exit Loop2
    enddo Loop2
  endif
  if (Nsel > 0) then
    res_xs_sel = res_xs(Nsel)
    res_dxs_sel = res_dxs(Nsel)
    res_author_sel = res_author(Nsel)
  else
    res_xs_sel = 0.
    res_dxs_sel = 0.
    res_author_sel = ''
  endif
  Nsave = Nsave +1
  N = Nsave
  Zsave(N) = Z
  Asave(N) = A
  Lisosave(N) = Liso
  xssave(N) = res_xs_sel 
  dxssave(N) = res_dxs_sel 
  refsave(N) = res_author_sel 
  Nexpsave(N) = Nres_exp 
  return
end subroutine procthermal
! Copyright A.J. Koning 2025
