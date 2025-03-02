subroutine writethermal(Z, A, Liso, Riso, type, flagav)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose: Write thermal cross section data
!
! Revision    Date      Author      Quality  Description
! ======================================================
!    1     2025-02-15   A.J. Koning    A     Original code
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
  character(len=2)   :: iso        ! extension
  character(len=3)   :: exten
  character(len=132) :: nucfile    ! nuclide file
  character(len=18)  :: rfile
  character(len=20)  :: react      ! reaction
  character(len=15)  :: col(8)     ! header
  character(len=15)  :: un(8)      ! units
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
!
! **************** Write databases for thermal cross sections *****
!
  if (.not.res_exist) return
  if (flagav) then
    exten='_av'
  else
    exten=''
  endif
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
  rfile=trim(reac(type))//trim(iso)//trim(exten)
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
  col(1) = 'Author'
  col(2) = ''
  col(3) = 'Type'
  col(4) = 'Year'
  col(5) = 'Value'
  col(6) = 'dValue'
  col(7) = 'Reference'
  col(8) = 'Ratio'
  Ncol = 8
  un = ''
  if (type <= 6) then
    un(4) = 'b'
    un(5) = 'b'
  endif
! do isource = 1, 3
!   if (isource == 1) Nr = Ncomp
!   if (isource == 2) Nr = Nexp
!   if (isource == 3) Nr = Nlib
  call write_datablock(quantity,Ncol,Nres,col,un)
! call write_datablock(quantity,Ncol,Nr,col,un)
  do k = 1, Nres
    F = res_xs(k) / res_xs_sel
!   if (isource == 1 .and.  res_type(k) /= 'Compilation') cycle
!   if (isource == 2 .and.  res_type(k) /= 'EXFOR') cycle
!   if (isource == 3 .and.  res_type(k) /= 'NDL') cycle
    write(1, '(a30,a15,6x,i4,5x,2es15.6,3x,a12,es15.6)') res_author(k), res_type(k), res_year(k), res_xs(k), & 
 &    res_dxs(k), res_ref(k), F
  enddo
! enddo
  close(1)
  return
end subroutine writethermal
! Copyright A.J. Koning 2025
