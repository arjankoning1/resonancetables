subroutine writeresonance(Z, A, Liso, type)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose: Write resonance data
!
! Revision    Date      Author      Quality  Description
! ======================================================
!    1     2025-05-25   A.J. Koning    A     Original code
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
  character(len=132) :: nucfile    ! nuclide file
  character(len=20)  :: react      ! reaction
  character(len=15)  :: col(13)    ! header
  character(len=15)  :: un(13)     ! units
  character(len=80)  :: quantity   ! quantity
  character(len=132) :: topline    ! topline
  integer            :: Z          ! charge number
  integer            :: A          ! mass number
  integer            :: k          ! counter
  integer            :: Ncol       ! number of columns
  integer            :: Liso       ! target isomer
  integer            :: type       ! reaction type
  real               :: F
!
! **************** Write databases for resonance data *****
!
  if (.not.res_exist) return
  Ztarget = Z
  Atarget = A
  quantity='resonance data'
  react=restype(type)
  topline=trim(targetnuclide)//' '//trim(react)//' '//trim(quantity)
  nucfile=trim(respath)//trim(react)//'/nuc/'//trim(targetnuclide)//'.res'
  write(*,*) Z, A, Liso, trim(nucfile), " ", Nres
  open (unit = 1, status = 'unknown', file = trim(nucfile))
  call write_header(topline,source,user,date,oformat)
  call write_target
  call write_reaction(react,0.D0,0.D0,0,0)
  un = ''
  col(1) = 'Author'
  col(2) = ''
  col(3) = 'Type'
  col(4) = 'Year'
  col(5) = 'Value'
  col(6) = 'dValue'
  if (type == 2 .or. type == 5) then
    un(5) = '*e-4'
    un(6) = '*e-4'
  else
    un(5) = 'eV'
    un(6) = 'eV'
  endif
  col(7) = 'Reference'
  col(8) = 'Ratio'
  Ncol = 8
  write(1,'("# parameters:")')
  call write_real(2,'selected value ['//trim(un(5))//']',res_xs_sel)
  call write_real(2,'selected value uncertainty ['//trim(un(5))//']',res_dxs_sel)
  call write_char(2,'selected value source',res_author_sel)
  call write_integer(2,'number of values',Nres)
  call write_real(2,'average value [b]',av_xs)
  write(1,'("#   relative standard deviation [%]:",f15.6)') var_xs
  if (Ncomp > 0) then
    if (type == 1) then
      col(9) = '#Resonances'
      Ncol = 9
    endif
    quantity='Compilation'
    call write_quantity(quantity)
    call write_real(2,'average value',av_xs_comp)
    write(1,'("#   relative standard deviation [%]:",f15.6)') var_xs_comp
    call write_datablock(Ncol,Ncomp,col,un)
    do k = 1, Nres
      if (res_type(k) == 'Compilation') then
        F = res_xs(k) / res_xs_sel
        if (type == 1) then
          write(1, '(a30,a15,6x,i4,5x,2es15.6,3x,a12,f15.6,a15)') res_author(k), res_type(k), res_year(k), res_xs(k), &
 &          res_dxs(k), trim(adjustl(res_ref(k))), F, res_Nrr(k)
        else
          write(1, '(a30,a15,6x,i4,5x,2es15.6,3x,a12,f15.6,a15)') res_author(k), res_type(k), res_year(k), res_xs(k), &
 &          res_dxs(k), trim(adjustl(res_ref(k))), F
        endif
      endif
    enddo
  endif
  Ncol = 8
  if (Nexp > 0) then
    if (type == 1 .or. type == 4 .or. type == 7) then
      Ncol = 13
      col(9) = 'E-min'
      un(9) = 'MeV'
      col(10) = 'E-max'
      un(10) = 'MeV'
      col(11) = 'Momentum L'
      un(11) = ''
      col(12) = 'Spin J'
      un(12) = ''
      col(13) = 'Parity'
      un(13) = ''
    endif
    quantity='EXFOR'
    call write_quantity(quantity)
    call write_real(2,'average value',av_xs_exfor)
    write(1,'("#   relative standard deviation [%]:",f15.6)') var_xs_exfor
    call write_datablock(Ncol,Nexp,col,un)
    F = 1.
    do k = 1, Nres
      if (res_type(k) == 'EXFOR') then
        if (res_xs_sel > 0.) F = res_xs(k) / res_xs_sel
        if (type == 1 .or. type == 4 .or. type == 7) then
          write(1, '(a30,a15,6x,i4,5x,2es15.6,3x,a12,f15.6,5a15)') res_author(k), res_type(k), res_year(k), res_xs(k), &
 &          res_dxs(k), trim(adjustl(res_ref(k))), F, res_Emin(k), res_Emax(k), res_L(k), res_J(k), res_P(k)
       else
          write(1, '(a30,a15,6x,i4,5x,2es15.6,3x,a12,f15.6)') res_author(k), res_type(k), res_year(k), res_xs(k), &
 &          res_dxs(k), trim(adjustl(res_ref(k))), F
       endif
      endif
    enddo
  endif
  Ncol = 8
  if (Nlib > 0) then
    quantity='Nuclear data library'
    call write_quantity(quantity)
    call write_real(2,'average value',av_xs_NDL)
    write(1,'("#   relative standard deviation [%]:",f15.6)') var_xs_NDL
    call write_datablock(Ncol,Nlib,col,un)
    F = 1.
    do k = 1, Nres
      if (res_type(k) == 'NDL') then
        if (res_xs_sel > 0.) F = res_xs(k) / res_xs_sel
        write(1, '(a30,a15,6x,i4,5x,2es15.6,3x,a12,f15.6)') res_author(k), res_type(k), res_year(k), res_xs(k), &
 &        res_dxs(k), trim(adjustl(res_ref(k))), F
      endif
    enddo
  endif
  close(1)
  return
end subroutine writeresonance
! Copyright A.J. Koning 2025
