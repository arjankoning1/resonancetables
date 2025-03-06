subroutine writemacs(Z, A, Liso, Riso)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose: Write MACS data
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
  character(len=6)   :: dir
  character(len=20)  :: react      ! reaction
  character(len=15)  :: col(9)     ! header
  character(len=15)  :: un(9)      ! units
  character(len=80)  :: quantity   ! quantity
  character(len=132) :: topline    ! topline
  integer            :: Z          ! charge number
  integer            :: A          ! mass number
  integer            :: k          ! counter
  integer            :: Ncol       ! number of columns
  integer            :: Liso       ! target isomer
  integer            :: Riso       ! residual isomer
  real               :: F
!
! **************** Write databases for MACS *****
!
  if (.not.res_exist) return
  dir='ng/'
  Ztarget = Z
  Atarget = A
  quantity='MACS'
  react=reaction(4)
  topline=trim(targetnuclide)//trim(react)//' '//trim(quantity)
  nucfile=trim(macspath)//trim(dir)//trim(targetnuclide)//'.macs'
  write(*,*) Z, A, Liso, Riso, trim(nucfile), " ", Nres
  open (unit = 1, status = 'unknown', file = trim(nucfile))
  call write_header(topline,source,user,date,oformat)
  call write_target
  call write_reaction(react,0.D0,0.D0,0,0)
  write(1,'("# observables:")')
  call write_real(2,'selected value [b]',res_xs_sel)
  call write_real(2,'selected value uncertainty [b]',res_dxs_sel)
  call write_char(2,'selected value source',res_author_sel)
  col(1) = 'Author'
  col(2) = ''
  col(3) = 'Type'
  col(4) = 'Year'
  col(5) = 'Value'
  col(6) = 'dValue'
  col(7) = 'Reference'
  col(8) = 'Ratio'
  col(9) = 'Average'
  Ncol = 9
  un = ''
  un(5) = 'b'
  un(6) = 'b'
  if (Ncomp > 0) then
    quantity='Compilation'
    call write_quantity(quantity)
    call write_datablock(Ncol,Ncomp,col,un)
    do k = 1, Nres
      if (res_type(k) == 'Compilation' .and. res_av(k) == '') then
        F = res_xs(k) / res_xs_sel
        write(1, '(a30,a15,6x,i4,5x,2es15.6,3x,a12,es15.6,3x,a12)') res_author(k), res_type(k), res_year(k), res_xs(k), &
 &        res_dxs(k), res_ref(k), F, res_av(k)
      endif
    enddo
  endif
  if (Ncomp_av > 0) then
    quantity='Compilation spectrum-averaged'
    call write_quantity(quantity)
    call write_datablock(Ncol,Ncomp_av,col,un)
    do k = 1, Nres
      if (res_type(k) == 'Compilation' .and. res_av(k) /= '') then
        F = res_xs(k) / res_xs_sel
        write(1, '(a30,a15,6x,i4,5x,2es15.6,3x,a12,es15.6,3x,a12)') res_author(k), res_type(k), res_year(k), res_xs(k), &
 &        res_dxs(k), res_ref(k), F, res_av(k)
      endif
    enddo
  endif
  if (Nexp > 0) then
    quantity='EXFOR'
    call write_quantity(quantity)
    call write_datablock(Ncol,Nexp,col,un)
    do k = 1, Nres
      if (res_type(k) == 'EXFOR' .and. res_av(k) == '') then
        F = res_xs(k) / res_xs_sel
        write(1, '(a30,a15,6x,i4,5x,2es15.6,3x,a12,es15.6,3x,a12)') res_author(k), res_type(k), res_year(k), res_xs(k), &
 &        res_dxs(k), res_ref(k), F, res_av(k)
      endif
    enddo
  endif
  if (Nexp_av > 0) then
    quantity='EXFOR spectrum-averaged'
    call write_quantity(quantity)
    call write_datablock(Ncol,Nexp_av,col,un)
    do k = 1, Nres
      if (res_type(k) == 'EXFOR' .and. res_av(k) /= '') then
        F = res_xs(k) / res_xs_sel
        write(1, '(a30,a15,6x,i4,5x,2es15.6,3x,a12,es15.6,3x,a12)') res_author(k), res_type(k), res_year(k), res_xs(k), &
 &        res_dxs(k), res_ref(k), F, res_av(k)
      endif
    enddo
  endif
  if (Nlib > 0) then
    quantity='Nuclear data library'
    call write_quantity(quantity)
    call write_datablock(Ncol,Nlib,col,un)
    do k = 1, Nres
      if (res_type(k) == 'NDL') then
        F = res_xs(k) / res_xs_sel
        write(1, '(a30,a15,6x,i4,5x,2es15.6,3x,a12,es15.6,6x,a9)') res_author(k), res_type(k), res_year(k), res_xs(k), &
 &        res_dxs(k), res_ref(k), F, res_av(k)
      endif
    enddo
  endif
  close(1)
  return
end subroutine writemacs
! Copyright A.J. Koning 2025
