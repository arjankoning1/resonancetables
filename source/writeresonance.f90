subroutine writeresonance(Z, A, Liso, type)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose: Write resonance data
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
  character(len=132) :: nucfile    ! nuclide file
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
  nucfile=trim(respath)//trim(react)//'/'//trim(targetnuclide)//'.res'
  write(*,*) Z, A, Liso, trim(nucfile), " ", Nres
  open (unit = 1, status = 'unknown', file = trim(nucfile))
  call write_header(topline,source,user,date,oformat)
  call write_target
  call write_reaction(react,0.D0,0.D0,0,0)
  write(1,'("# observables:")')
  call write_real(2,'selected value [eV]',res_xs_sel)
  call write_real(2,'selected value uncertainty [eV]',res_dxs_sel)
  call write_char(2,'selected value source',res_author_sel)
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
  call write_datablock(quantity,Ncol,Nres,col,un)
  do k = 1, Nres
    F = res_xs(k) / res_xs_sel
    write(1, '(a30,a15,6x,i4,5x,2es15.6,3x,a12,es15.6)') res_author(k), res_type(k), res_year(k), res_xs(k), res_dxs(k), &
 &    trim(adjustl(res_ref(k))), F
  enddo
  close(1)
  return
end subroutine writeresonance
! Copyright A.J. Koning 2025
