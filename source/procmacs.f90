subroutine procmacs(Z, A, Liso)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose: Process MACS data
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
  character(len=40)  :: ttmp          ! subentry
  character(len=24)  :: atmp          ! author
  character(len=40)  :: avtmp 
  character(len=40)  :: macschoice 
  character(len=132) :: macsfile 
  character(len=40)  :: rtmp          ! reference
  integer            :: Z             ! charge number
  integer            :: A             ! mass number
  integer            :: N             ! counter
  integer            :: i             ! counter
  integer            :: j             ! counter
  integer            :: Nsel          ! counter
  integer            :: ytmp          ! year
  integer            :: Liso          ! target isomer
  integer            :: iz
  integer            :: ia
  integer            :: istat
  real(sgl)          :: xstmp         ! cross section
  real(sgl)          :: Etmp 
  real(sgl)          :: Gtmp 
  real(sgl)          :: dxstmp        ! cross section uncertainty
  real(sgl)          :: sumxs
  real(sgl)          :: sumxs_comp
  real(sgl)          :: sumxs_av_comp
  real(sgl)          :: sumxs_exfor 
  real(sgl)          :: sumxs_av_exfor
  real(sgl)          :: sumxs_NDL
  real(sgl)          :: varsumxs
  real(sgl)          :: varsumxs_comp
  real(sgl)          :: varsumxs_av_comp
  real(sgl)          :: varsumxs_exfor
  real(sgl)          :: varsumxs_av_exfor
  real(sgl)          :: varsumxs_NDL
  real(sgl)          :: F
!
! **************** Process databases for MACS *****
!
! Sort data according to year
!
  Nexp = 0
  Nexp_av = 0
  Ncomp = 0
  Ncomp_av = 0
  Nlib = 0
  av_xs = 0.
  av_xs_comp = 0.
  av_xs_av_comp = 0.
  av_xs_exfor = 0.
  av_xs_av_exfor = 0.
  av_xs_NDL = 0.
  var_xs = 0.
  var_xs_comp = 0.
  var_xs_av_comp = 0.
  var_xs_exfor = 0.
  var_xs_av_exfor = 0.
  var_xs_NDL = 0.
  N = Nres
  if (N > 0) then
    do i = 1, N
      do j = 1, i
        if (res_year(i) > res_year(j)) cycle
        xstmp = res_xs(i)
        dxstmp = res_dxs(i)
        Gtmp = res_G(i)
        ytmp = res_year(i)
        atmp = res_author(i)
        ttmp = res_type(i)
        rtmp = res_ref(i)
        avtmp = res_av(i)
        Etmp = res_E(i)
        res_xs(i) = res_xs(j)
        res_dxs(i) = res_dxs(j)
        res_year(i) = res_year(j)
        res_author(i) = res_author(j)
        res_type(i) = res_type(j)
        res_ref(i) = res_ref(j)
        res_av(i) = res_av(j)
        res_E(i) = res_E(j)
        res_G(i) = res_G(j)
        res_xs(j) = xstmp 
        res_dxs(j) = dxstmp 
        res_year(j) = ytmp 
        res_author(j) = atmp 
        res_type(j) = ttmp 
        res_ref(j) = rtmp
        res_av(j) = avtmp
        res_E(j) = Etmp
        res_G(j) = Gtmp
      enddo
    enddo
!
! Determine number of each type, average and coefficient of variation
!
    sumxs = 0.
    sumxs_comp = 0.
    sumxs_av_comp = 0.
    sumxs_exfor = 0.
    sumxs_av_exfor = 0.
    sumxs_NDL = 0.
    varsumxs = 0.
    varsumxs_comp = 0.
    varsumxs_av_comp = 0.
    varsumxs_exfor = 0.
    varsumxs_av_exfor = 0.
    varsumxs_NDL = 0.
    do i = 1, N
      if (res_type(i) == 'Compilation') then
        if (res_av(i) == '') then
          Ncomp = Ncomp + 1
          sumxs_comp = sumxs_comp + res_xs(i)
        else
          Ncomp_av = Ncomp_av + 1
          sumxs_av_comp = sumxs_av_comp + res_xs(i)
        endif
      endif
      if (res_type(i) == 'EXFOR') then
        if (res_av(i) == '') then
          Nexp = Nexp + 1
          sumxs_exfor = sumxs_exfor + res_xs(i)
        else
          Nexp_av = Nexp_av + 1
          sumxs_av_exfor = sumxs_av_exfor + res_xs(i)
        endif
      endif
      if (res_type(i) == 'NDL') then
        Nlib = Nlib + 1
        sumxs_NDL = sumxs_NDL + res_xs(i)
      endif
      sumxs = sumxs + res_xs(i)
    enddo
    if (Ncomp > 0) av_xs_comp = sumxs_comp/Ncomp
    if (Ncomp_av > 0) av_xs_av_comp = sumxs_av_comp/Ncomp_av
    if (Nexp > 0) av_xs_exfor = sumxs_exfor/Nexp
    if (Nexp_av > 0) av_xs_av_exfor = sumxs_av_exfor/Nexp_av
    if (Nlib > 0) av_xs_NDL = sumxs_NDL/Nlib
    av_xs = sumxs/N
    do i = 1, N
      if (res_type(i) == 'Compilation') then
        if (res_av(i) == '') then
          varsumxs_comp = varsumxs_comp + (res_xs(i)-av_xs_comp)**2
        else
          varsumxs_av_comp = varsumxs_av_comp + (res_xs(i)-av_xs_av_comp)**2
        endif
      endif
      if (res_type(i) == 'EXFOR') then
        if (res_av(i) == '') then
          varsumxs_exfor = varsumxs_exfor + (res_xs(i)-av_xs_exfor)**2
        else
          varsumxs_av_exfor = varsumxs_av_exfor + (res_xs(i)-av_xs_av_exfor)**2
        endif
      endif
      if (res_type(i) == 'NDL') then
        varsumxs_NDL = varsumxs_NDL + (res_xs(i)-av_xs_NDL)**2
      endif
      varsumxs = varsumxs + (res_xs(i)-av_xs)**2
    enddo
    if (Ncomp > 0) var_xs_comp = 100. * sqrt(varsumxs_comp/Ncomp)/av_xs_comp
    if (Ncomp_av > 0) var_xs_av_comp = 100. * sqrt(varsumxs_av_comp/Ncomp_av)/av_xs_av_comp
    if (Nexp > 0) var_xs_exfor = 100. * sqrt(varsumxs_exfor/Nexp)/av_xs_exfor
    if (Nexp_av > 0) var_xs_av_exfor = 100. * sqrt(varsumxs_av_exfor/Nexp_av)/av_xs_av_exfor
    if (Nlib > 0) var_xs_NDL = 100. * sqrt(varsumxs_NDL/Nlib)/av_xs_NDL
    if (Nres > 0) var_xs = 100. * sqrt(varsumxs/Nres)/av_xs
    if (var_xs_comp < 1.e-10) var_xs_comp = 0.
    if (var_xs_av_comp < 1.e-10) var_xs_av_comp = 0.
    if (var_xs_exfor < 1.e-10) var_xs_exfor = 0.
    if (var_xs_av_exfor < 1.e-10) var_xs_av_exfor = 0.
    if (var_xs_NDL < 1.e-10) var_xs_NDL = 0.
    if (var_xs < 1.e-10) var_xs = 0.
  endif
!
! Final dataset
! Rule: ASTRAL > KADONIS > EXFOR
!
  Nsel = 0
  N = Nres
Loop1:  do
    do i = 1, N
      if (res_author(i) == 'Astral') then
        Nsel = i
        exit Loop1
      endif
    enddo
    do i = 1, N
      if (res_author(i) == 'Kadonis') then
        Nsel = i
        exit Loop1
      endif
    enddo
    do i = 1, N
      if (res_author(i) == 'Bao') then
        Nsel = i
        exit Loop1
      endif
    enddo
    if (A > 0) Nsel = Nres_exp
    exit Loop1
  enddo Loop1
!
! Overrule rule by specific cases
!
  macsfile = trim(filespath)//'macs_koning.ng'
  open (unit = 2, status = 'old', file = macsfile)
  do
    read(2,'(2i4, 1x, a)', iostat = istat) iz, ia, macschoice
    if (istat == -1) exit
    if (istat > 0) call read_error(macsfile, istat)
    if (iz == Z .and. ia == A) then
      do i = 1, N
        if (trim(res_author(i)) == trim(macschoice)) then
          Nsel = i
          exit 
        endif
      enddo
    endif
  enddo
  close(2)
  if (Nsel > 0) then
     res_xs_sel = res_xs(Nsel)
     res_dxs_sel = res_dxs(Nsel)
     res_G_sel = res_G(Nsel)
     res_author_sel = res_author(Nsel)
     res_av_sel = res_av(Nsel)
    Nsave = Nsave +1
    N = Nsave
    Zsave(N) = Z
    Asave(N) = A
    Lisosave(N) = Liso
    xssave(N) = res_xs_sel
    dxssave(N) = res_dxs_sel
    Gsave(N) = res_G_sel
    refsave(N) = res_author_sel
    avsave(N) = res_av_sel
    Nexpsave(N) = Nres_exp
    varsave(N) = var_xs
    compsave(N) = var_xs_av_comp
    NDLsave(N) = var_xs_NDL
    expsave(N) = var_xs_exfor
    do i = 1, Nres
      F = res_xs(i) / res_xs(Nsel)    
      if (F < 0.2 .or. F > 5.) then
        write(*,*) "Warning: Ratio= ", F, "for Z= ", Z , " A= ", A, "author ", res_author(i)
      endif          
    enddo
  else
     res_xs_sel = 0.
     res_dxs_sel = 0.
     res_G_sel = 0.
     res_author_sel = ''
     res_av_sel = ''
  endif
  return
end subroutine procmacs
! Copyright A.J. Koning 2025
