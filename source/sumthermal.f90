subroutine sumthermal(Riso, type, flagav)
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose: Summarize thermal cross section data
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
  character(len=2)   :: iso        ! extension
  character(len=3)   :: Astring    ! mass string
  character(len=3)   :: exten
  character(len=6)   :: nuclide
  character(len=132) :: nucfile    ! nuclide file
  character(len=18)  :: rfile
  character(len=18)  :: react      ! reaction
  character(len=15)  :: col(8)     ! header
  character(len=15)  :: un(8)      ! units
  character(len=80)  :: quantity   ! quantity
  character(len=132) :: topline    ! topline
  integer            :: k          ! counter
  integer            :: N          ! counter
  integer            :: Ncol       ! number of columns
  integer            :: Riso       ! residual isomer
  integer            :: type       ! reaction type
!
! **************** Write databases for thermal cross sections *****
!
  if (flagav) then
    exten='_av'
  else
    exten=''
  endif
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
  topline=trim(react)//' '//trim(quantity)
  rfile=trim(reac(type))//trim(iso)//exten
  nucfile=trim(thermalpath)//trim(rfile)//'/thermal.'//trim(rfile)
  write(*,*) " Writing to ", trim(nucfile)
  open (unit = 1, status = 'unknown', file = trim(nucfile))
  call write_header(topline,source,user,date,oformat)
  call write_reaction(react,0.D0,0.D0,0,0)
  col(1) = 'Z'
  col(2) = 'A'
  col(3) = 'Liso'
  col(4) = 'xs'
  col(5) = 'dxs'
  col(6) = 'Reference'
  col(7) = '#Experiments'
  col(8) = 'Nuclide'
  Ncol = 8
  un = ''
  if (type <= 6) then
    un(4) = 'b'
    un(5) = 'b'
  endif
  N = Nsave
  call write_datablock(quantity,Ncol,N,col,un)
  do k = 1, N
    Astring='   '
    write(Astring(1:3),'(i3.3)') Asave(k)
    nuclide=trim(nuc(Zsave(k)))//Astring
    if (Lisosave(k) == 1) nuclide = trim(nuclide)//'m'
    if (Lisosave(k) == 2) nuclide = trim(nuclide)//'n'
    write(1, '(3(6x,i4,5x),2es15.6,2x,a15,4x,i4,11x,a6)') Zsave(k), Asave(k), Lisosave(k), xssave(k), dxssave(k), refsave(k), &
 &     Nexpsave(k), nuclide
  enddo
  close(1)
  return
end subroutine sumthermal
! Copyright A.J. Koning 2025
