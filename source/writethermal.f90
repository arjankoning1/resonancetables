subroutine writethermal(Z, A, Liso, Riso, type)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose: Write thermal cross section data
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
  character(len=2)   :: iso        ! extension
  character(len=132) :: nucfile    ! nuclide file
  character(len=18)  :: rfile
  character(len=20)  :: react      ! reaction
  character(len=15)  :: col(10)     ! header
  character(len=15)  :: un(10)      ! units
  character(len=80)  :: quantity   ! quantity
  character(len=132) :: topline    ! topline
  integer            :: Z          ! charge number
  integer            :: A          ! mass number
  integer            :: k          ! counter
  integer            :: Ncol       ! number of columns
  integer            :: Liso       ! target isomer
  integer            :: Riso       ! residual isomer
  integer            :: type       ! reaction type
  real               :: F
  real               :: tolerance
!
! **************** Write databases for thermal cross sections *****
!
  if (.not.res_exist) return
  Ztarget = Z
  Atarget = A
  iso = ''
  if (Riso == 0) iso='_g'
  if (Riso == 1) iso='_m'
  if (Riso == 2) iso='_n'
  if (type <= 6) then
    quantity='thermal cross section'
  else
    quantity='thermal neutron multiplicity'
  endif
  react=trim(reaction(type))//iso
  topline=trim(targetnuclide)//trim(react)//' '//trim(quantity)
  rfile=trim(reac(type))//trim(iso)
  nucfile=trim(thermalpath)//trim(rfile)//'/'//trim(targetnuclide)//'.'//trim(rfile)
  write(*,*) Z, A, Liso, Riso, trim(nucfile), " ", Nres
  open (unit = 1, status = 'unknown', file = trim(nucfile))
  call write_header(topline,source,user,date,oformat)
  call write_target
  call write_reaction(react,0.D0,0.D0,0,0)
  write(1,'("# observables:")')
  if (type <= 6) then
    call write_real(2,'selected value [b]',res_xs_sel)
    call write_real(2,'selected value uncertainty [b]',res_dxs_sel)
  else
    call write_real(2,'selected value]',res_xs_sel)
    call write_real(2,'selected value uncertainty]',res_dxs_sel)
  endif
  call write_char(2,'selected value source',res_author_sel)
  call write_integer(2,'number of values',Nres)
  if (type <= 6) then
    call write_real(2,'average value [b]',av_xs)
  else
    call write_real(2,'average value',av_xs)
  endif
  write(1,'("#   relative standard deviation [%]:",f15.6)') var_xs
  col(1) = 'Author'
  col(2) = ''
  col(3) = 'Type'
  col(4) = 'Year'
  col(5) = 'Value'
  col(6) = 'dValue'
  col(7) = 'Reference'
  col(8) = 'Ratio'
  col(9) = 'Spectrum'
  col(10) = 'Energy'
  Ncol = 10
  un = ''
  un(10) = 'MeV'
  if (type <= 6) then
    un(5) = 'b'
    un(6) = 'b'
  endif
  if (Ncomp > 0) then
    quantity='Compilation'
    call write_quantity(quantity)
    call write_real(2,'average value',av_xs_comp)
    write(1,'("#   relative standard deviation [%]:",f15.6)') var_xs_comp
    call write_datablock(Ncol,Ncomp,col,un)
    F = 1.
    do k = 1, Nres
      if (res_type(k) == 'Compilation' .and. res_av(k) == '') then
        if (res_xs_sel > 0.) F = res_xs(k) / res_xs_sel
        write(1, '(a30,a15,6x,i4,5x,2es15.6,3x,a12,f15.6,3x,a12,es15.6)') res_author(k), res_type(k), res_year(k), res_xs(k), & 
 &        res_dxs(k), res_ref(k), F, res_av(k), res_E(k)
      endif
    enddo
  endif
  if (Ncomp_av > 0) then
    quantity='Compilation spectrum-averaged'
    call write_quantity(quantity)
    call write_real(2,'average value',av_xs_av_comp)
    write(1,'("#   relative standard deviation [%]:",f15.6)') var_xs_av_comp
    call write_datablock(Ncol,Ncomp_av,col,un)
    F = 1.
    do k = 1, Nres
      if (res_type(k) == 'Compilation' .and. res_av(k) /= '') then
        if (res_xs_sel > 0.) F = res_xs(k) / res_xs_sel
        write(1, '(a30,a15,6x,i4,5x,2es15.6,3x,a12,f15.6,3x,a12,es15.6)') res_author(k), res_type(k), res_year(k), res_xs(k), & 
 &        res_dxs(k), res_ref(k), F, res_av(k), res_E(k)
      endif
    enddo
  endif
  if (Nexp > 0) then
    quantity='EXFOR'
    call write_quantity(quantity)
    call write_real(2,'average value',av_xs_exfor)
    write(1,'("#   relative standard deviation [%]:",f15.6)') var_xs_exfor
    call write_datablock(Ncol,Nexp,col,un)
    F = 1.
    do k = 1, Nres
      if (res_type(k) == 'EXFOR' .and. res_av(k) == '') then
        if (res_xs_sel > 0.) F = res_xs(k) / res_xs_sel
        write(1, '(a30,a15,6x,i4,5x,2es15.6,3x,a12,f15.6,3x,a12,es15.6)') res_author(k), res_type(k), res_year(k), res_xs(k), & 
 &        res_dxs(k), res_ref(k), F, res_av(k), res_E(k)
      endif
    enddo
  endif
  if (Nexp_av > 0) then
    quantity='EXFOR spectrum-averaged'
    call write_quantity(quantity)
    call write_real(2,'average value',av_xs_av_exfor)
    write(1,'("#   relative standard deviation [%]:",f15.6)') var_xs_av_exfor
    call write_datablock(Ncol,Nexp_av,col,un)
    F = 1.
    do k = 1, Nres
      if (res_type(k) == 'EXFOR' .and. res_av(k) /= '') then
        if (res_xs_sel > 0.) F = res_xs(k) / res_xs_sel
        write(1, '(a30,a15,6x,i4,5x,2es15.6,3x,a12,f15.6,3x,a12,es15.6)') res_author(k), res_type(k), res_year(k), res_xs(k), & 
 &        res_dxs(k), res_ref(k), F, res_av(k), res_E(k)
      endif
    enddo
  endif
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
        write(1, '(a30,a15,6x,i4,5x,2es15.6,3x,a12,f15.6,6x,a9,es15.6)') res_author(k), res_type(k), res_year(k), res_xs(k), & 
 &        res_dxs(k), res_ref(k), F, res_av(k), res_E(k)
      endif
    enddo
  endif
  close(1)
  return
end subroutine writethermal
! Copyright A.J. Koning 2025
