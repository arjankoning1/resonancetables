subroutine writeresonance
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose: Write average resonance data
!
! Revision    Date      Author      Quality  Description
! ======================================================
!    1     01-01-2019   A.J. Koning    A     Original code
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
  character(len=3)   :: Astring       ! mass string
  character(len=4)   :: Zstring       ! Z string
! character(len=4)   :: year          ! year
  character(len=6)   :: nuclide       ! nuclide
  character(len=6), dimension(100)   :: filestatus    ! status of file
  character(len=40)  :: ref           ! reference
  character(len=40)  :: ref0          ! reference
! character(len=132) :: nucfile       ! nuclide file
  character(len=132) :: reacdir       ! reaction directory
  character(len=132) :: reacfile      ! reaction file
  character(len=132) :: reaclib       ! library file
  character(len=132) :: xsstring
  character(len=132) :: ratiostring
  integer            :: Z             ! charge number
  integer            :: A             ! mass number
  integer            :: i             ! counter
  integer            :: iUxs          ! unit for cross section
  integer            :: iUrat         ! unit for ratios
  integer            :: iUall         ! unit for all data
! integer            :: N             ! counter
  integer            :: lib           ! library
  integer            :: type          ! reaction type
  real(sgl)          :: ratio         ! ratio
  real(sgl)          :: xs            ! cross section
  real(sgl)          :: xs0           ! cross section
  real(sgl)          :: dxs           ! cross section uncertainty
!
! **************** Write databases for average resonance parameters *****
!
  xsstring='#  Z   A       Value           dValue     Ref                                  Nuc'
  ratiostring='#  Z   A       Ratio           R(this)      R0(final)     Ref                  Nuc'
  write(*, *) "Writing data to databases....."
  do type = 1, 5
    reacdir = trim(respath)//trim(restype(type))//'/'
    reacfile = trim(reacdir)//trim(restype(type))//'.'
    open (unit = 1, status = 'unknown', file = trim(reacfile)//'all')
    write(1, '(a)') trim(xsstring)
    do lib = 0, numlib
      filestatus = 'delete'
      if (lib == 2) cycle
      reaclib = trim(reacfile)//trim(ext(lib))
      open (unit = 11, status = 'unknown', file = trim(reaclib))
      write(11, '(a)') trim(xsstring)
      if (lib > 0) then
        open (unit = 21, status = 'unknown', file = trim(reaclib)//'.ratio')
        write(21, '(a)') trim(ratiostring)
      endif
      do Z = 1, numZ
        do A = 0, numA
          Astring='   '
          write(Astring(1:3),'(i3.3)') A
          nuclide=trim(nuc(Z))//Astring
          Zstring='    '
          write(Zstring(2:4),'(i3)') Z
          if (A == 0) Zstring(1:1) = '#'
!
! Compilations
!
          iUall = 1
          iUxs = 11
          iUrat = 21
          if (lib < 6) then
            if (res_R(lib, type, Z, A) > 0.) then
              xs = res_R(lib, type, Z, A)
              dxs = res_dR(lib, type, Z, A)
              ref = res_ref(lib, type, Z, A)
              ref0 = res_ref(lib, type, Z, A)
              write(iUxs, '(a4,i4,2es15.5,3x,a,t80,a6)') Zstring, A, xs ,dxs, trim(ref), nuclide
              filestatus(iUxs) = 'keep'
              if (lib > 0) then
                xs0 = res_R(0, type, Z, A)
                ratio = ratio_R(lib, type, Z, A)
                write(iUrat, '(a4,i4,3es15.5,3x,a,t80,a6)') Zstring, A, ratio, xs, xs0, trim(ref0), nuclide
                filestatus(iUrat) = 'keep'
              endif
            endif
          else
!
! EXFOR
!
!           N = Nres_R(type, Z, A)
!           if (N > 0) then
!             do i = 1, N
!               xs = Etherm_xs(type, Z, A, i)
!               dxs = Etherm_dxs(type, Z, A, i)
!               ref = Etherm_ref(type, Z, A, i)
!               write(iUxs, '(a4,i4,2es15.5,3x,a)') Zstring, A, xs ,dxs, trim(ref)
!               filestatus(iUxs) = 'keep'
!               if (lib > 0) then
!                 xs0 = res_R(0, type, Z, A)
!                 ratio = Eratio_R(type, Z, A, i)
!                 write(iUrat, '(a4,i4,3es15.5,3x,a)') Zstring, A, ratio, xs, xs0, trim(ref)
!                 filestatus(iUrat) = 'keep'
!               endif
!             enddo
!             nucfile=trim(respath)//'nuc/'//trim(nuclide)//'.'//reac(type)
!             open (unit = 7, status = 'unknown', file = trim(nucfile))
!             do i = 1, N
!               xs = Eres_R(type, Z, A,  i)
!               dxs = Eres_dR(type, Z, A,  i)
!               year = Eres_year(type, Z, A,  i)
!               ref = Eres_ref(type, Z, A,  i)
!               write(7, '(a4,2es15.5,3x,a)') year, xs ,dxs, trim(ref)
!             enddo
!             close (7)
!           endif
          endif
!
! All data
!
          if (res_R(0, type, Z, A) > 0.) then
            write(iUall, '(a4,i4,2es15.5,3x,a,t80,a6)') Zstring, A, xs ,dxs, trim(ref), nuclide
            do i = 1, numlib
              xs = res_R(i, type, Z, A)
              dxs = res_dR(i, type, Z, A)
              ref = res_ref(i, type, Z, A)
              if (xs > 0.) then
                write(iUall, '(4x,a12,2es15.5,3x,a)') ext(i), xs ,dxs, trim(ref)
                filestatus(iUall) = 'keep'
              endif
            enddo
!           N = Nres_R(type, Z, A)
!           do i = 1, N
!             xs = Eres_R(type, Z, A, , i)
!             dxs = Eres_dR(type, Z, A, , i)
!             ref = Eres_ref(type, Z, A, , i)
!             write(iUall, '(4x,a12,2es15.5,3x,a)') ext(numlib), xs ,dxs, trim(ref)
!             filestatus(iUall) = 'keep'
!           enddo
          endif
        enddo
      enddo
      close (11, status = filestatus(11))
      close (21, status = filestatus(21))
    enddo
    close (1)
    if (type == 4) then
      close (2)
      close (3)
      close (4)
    endif
  enddo
  return
end subroutine writeresonance
! Copyright A.J. Koning 2019
