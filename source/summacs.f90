subroutine summacs
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose: Summarize MACS
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
  character(len=3)   :: Astring    ! mass string
  character(len=6)   :: dir
  character(len=6)   :: nuclide
  character(len=132) :: nucfile    ! nuclide file
  character(len=20)  :: react      ! reaction
  character(len=15)  :: col(9)     ! header
  character(len=15)  :: un(9)      ! units
  character(len=80)  :: quantity   ! quantity
  character(len=132) :: topline    ! topline
  logical            :: flagav
  integer            :: k          ! counter
  integer            :: N          ! counter
  integer            :: Ncol       ! number of columns
!
! **************** Write databases for thermal cross sections *****
!
  quantity='MACS'
  dir='ng/'
  react=reaction(4)
  topline=trim(react)//' '//trim(quantity)
  nucfile=trim(macspath)//trim(dir)//'ng.macs'
  write(*,*) " Writing to ", trim(nucfile) 
  open (unit = 1, status = 'unknown', file = trim(nucfile))
  call write_header(topline,source,user,date,oformat)
  call write_reaction(react,0.D0,0.D0,0,0)
  col(1) = 'Z'
  col(2) = 'A'
  col(3) = 'Liso'
  col(4) = 'Value'
  col(5) = 'dValue'
  col(6) = 'Reference'
  col(7) = '#Experiments'
  col(8) = 'Nuclide'
  col(9) = 'Average'
  Ncol = 9
  un = ''
  un(4) = 'b'
  un(5) = 'b'
  N = Nsave
  call write_quantity(quantity)
  call write_datablock(Ncol,N,col,un)
  do k = 1, N
    Astring='   '
    write(Astring(1:3),'(i3.3)') Asave(k)
    nuclide=trim(nuc(Zsave(k)))//Astring
    if (Lisosave(k) == 1) nuclide = trim(nuclide)//'m'
    if (Lisosave(k) == 2) nuclide = trim(nuclide)//'n'
    write(1, '(3(6x,i4,5x),2es15.6,2x,a15,4x,i4,11x,a6,3x,a9)') Zsave(k), Asave(k), Lisosave(k), xssave(k), dxssave(k), & 
 &     refsave(k), Nexpsave(k), nuclide, avsave(k)
  enddo
  close(1)
  return
end subroutine summacs
! Copyright A.J. Koning 2025
