subroutine writethermal
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose: Write thermal cross section data
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
  character(len=4)   :: year          ! year
  character(len=6)   :: nuclide       ! nuclide
  character(len=6), dimension(100)   :: filestatus    ! status of file
  character(len=40)  :: ref           ! reference
  character(len=40)  :: ref0          ! reference
  character(len=132) :: nucfile       ! nuclide file
  character(len=132) :: reacdir       ! reaction directory
  character(len=132) :: reacfile      ! reaction file
  character(len=132) :: reaclib       ! library file
  character(len=132) :: xsstring      ! 
  character(len=132) :: ratiostring   ! 
  integer            :: Z             ! charge number
  integer            :: A             ! mass number
  integer            :: i             ! counter
  integer            :: j             ! counter
  integer            :: iUxs          ! unit for cross section
  integer            :: iUrat         ! unit for ratios
  integer            :: iUall         ! unit for all data
  integer            :: N             ! counter
  integer            :: lib           ! library
  integer            :: isoT          ! target isomer
  integer            :: isoR          ! residual isomer
  integer            :: type          ! reaction type
  real(sgl)          :: ratio         ! ratio
  real(sgl)          :: xs            ! cross section
  real(sgl)          :: xs0           ! cross section
  real(sgl)          :: dxs           ! cross section uncertainty
!
! **************** Write databases for thermal cross sections *****
!
  xsstring='#  Z   A Tiso        xs             dxs       Ref                              Nuc'
  ratiostring='#  Z   A Tiso        Ratio       xs(this)      xs(final)    Ref                Nuc'
  write(*, *) "Writing data to databases....."
  do type = 1, numtype
    reacdir = trim(thermalpath)//trim(reac(type))//'/'
    reacfile = trim(reacdir)//'thermal_'//trim(reac(type))
    open (unit = 1, status = 'unknown', file = trim(reacfile)//'.all')
    write(1, '(a)') trim(xsstring)
    if (type == 4) then
      open (unit = 2, status = 'unknown', file = trim(reacfile)//'_g.all')
      open (unit = 3, status = 'unknown', file = trim(reacfile)//'_m.all')
      open (unit = 4, status = 'unknown', file = trim(reacfile)//'_n.all')
      write(2, '(a)') trim(xsstring)
      write(3, '(a)') trim(xsstring)
      write(4, '(a)') trim(xsstring)
    else
      if (type <= 6) open (unit = 8, status = 'unknown', file = trim(reacfile)//'.capratio')
    endif
    do lib = 0, numlib
      filestatus = 'delete'
      reaclib = trim(reacfile)//'.'//trim(ext(lib))
      open (unit = 11, status = 'unknown', file = trim(reaclib))
      write(11, '(a,"    #E")') trim(xsstring)
      if (type == 4) then
        open (unit = 12, status = 'unknown', file = trim(reacfile)//'_g.'//trim(ext(lib)))
        open (unit = 13, status = 'unknown', file = trim(reacfile)//'_m.'//trim(ext(lib)))
        open (unit = 14, status = 'unknown', file = trim(reacfile)//'_n.'//trim(ext(lib)))
        write(12, '(a,"    #E")') trim(xsstring)
        write(13, '(a,"    #E")') trim(xsstring)
        write(14, '(a,"    #E")') trim(xsstring)
      else
        if (lib == 0 .and. type <= 6) then
          write(8, '("#  Z   A Tiso        Ratio")')
        endif
      endif
      if (lib > 0 .and. type <= 6) then
        open (unit = 21, status = 'unknown', file = trim(reaclib)//'.ratio')
        write(21, '(a)') trim(ratiostring)
        if (type == 4) then
          open (unit = 22, status = 'unknown', file = trim(reacfile)//'_g.'//trim(ext(lib))//'.ratio')
          open (unit = 23, status = 'unknown', file = trim(reacfile)//'_m.'//trim(ext(lib))//'.ratio')
          open (unit = 24, status = 'unknown', file = trim(reacfile)//'_n.'//trim(ext(lib))//'.ratio')
          write(22, '(a)') trim(ratiostring)
          write(23, '(a)') trim(ratiostring)
          write(24, '(a)') trim(ratiostring)
        endif
      endif
      do Z = 1, numZ
        do A = 0, numA
          Astring='   '
          write(Astring(1:3),'(i3.3)') A
          nuclide=trim(nuc(Z))//Astring
          Zstring='    '
          write(Zstring(2:4),'(i3)') Z
          if (A == 0) Zstring(1:1) = '#'
          do isoT = 0, numisom
            if (isoT == 1) nuclide = trim(nuclide)//'m'
            if (isoT == 2) nuclide = trim(nuclide)//'n'
            do isoR = -1, numisom
!
! Compilations and ratios
!
              iUall = 2 + isoR
              iUxs = 12 + isoR
              iUrat = 22 + isoR
              if (lib < numlib) then
                if (therm_xs(lib, type, Z, A, isoT, isoR) > 0.) then
                  xs = therm_xs(lib, type, Z, A, isoT, isoR)
                  dxs = therm_dxs(lib, type, Z, A, isoT, isoR)
                  ref = therm_ref(lib, type, Z, A, isoT, isoR)
                  ref0 = therm_ref(0, type, Z, A, isoT, isoR)
                  N = Ntherm_xs(type, Z, A, isoT, isoR)
                  write(iUxs, '(a4,2i4,2es15.5,3x,a,t80,a6,i3)') Zstring, A, isoT, xs ,dxs, trim(ref), nuclide, N
                  filestatus(iUxs) = 'keep'
                  if (lib > 0 .and. type <= 6) then
                    xs0 = therm_xs(0, type, Z, A, isoT, isoR)
                    ratio = ratio_xs(lib, type, Z, A, isoT, isoR)
                    write(iUrat, '(a4,2i4,3es15.5,3x,a,t80,a6)') Zstring, A, isoT, ratio, xs, xs0, trim(ref0), nuclide
                    filestatus(iUrat) = 'keep'
                  endif
                endif
              else
!
! EXFOR
!
                N = Ntherm_xs(type, Z, A, isoT, isoR)
                if (N > 0) then
                  do i = 1, N
                    xs = Etherm_xs(type, Z, A, isoT, isoR, i)
                    dxs = Etherm_dxs(type, Z, A, isoT, isoR, i)
                    ref = Etherm_ref(type, Z, A, isoT, isoR, i)
                    write(iUxs, '(a4,2i4,2es15.5,3x,a,t80,a6,i4)') Zstring, A, isoT, xs ,dxs, trim(ref), nuclide, i
                    filestatus(iUxs) = 'keep'
                    if (lib > 0 .and. type <= 6) then
                      xs0 = therm_xs(0, type, Z, A, isoT, isoR)
                      ratio = Eratio_xs(type, Z, A, isoT, isoR, i)
                      write(iUrat, '(a4,2i4,3es15.5,3x,a,t95,a6)') Zstring, A, isoT, ratio, xs, xs0, trim(ref), nuclide
                      filestatus(iUrat) = 'keep'
                    endif
                  enddo
                  nucfile=trim(thermalpath)//'nuc/'//trim(nuclide)//'.'//reac(type)
                  open (unit = 7, status = 'unknown', file = trim(nucfile))
                  write(7, '("#year        xs             dxs       Ref")')
                  do i = 1, N
                    xs = Etherm_xs(type, Z, A, isoT, isoR, i)
                    dxs = Etherm_dxs(type, Z, A, isoT, isoR, i)
                    year = Etherm_year(type, Z, A, isoT, isoR, i)
                    ref = Etherm_ref(type, Z, A, isoT, isoR, i)
                    write(7, '(a4,2es15.5,3x,a)') year, xs ,dxs, trim(ref)
                  enddo
                  close (7)
                endif
              endif
!
! All data and ratio to capture
!
              if (lib == 0) then
                xs = therm_xs(0, type, Z, A, isoT, isoR)
                dxs = therm_dxs(0, type, Z, A, isoT, isoR)
                ref0 = therm_ref(0, type, Z, A, isoT, isoR)
                if (xs > 0.) then
                  write(iUall, '(a4,2i4,2es15.5,3x,a,t80,a6)') Zstring, A, isoT, xs ,dxs, trim(ref0), nuclide
                  filestatus(iUall) = 'keep'
                  do i = 1, numlib
                    xs = therm_xs(i, type, Z, A, isoT, isoR)
                    dxs = therm_dxs(i, type, Z, A, isoT, isoR)
                    ref = therm_ref(i, type, Z, A, isoT, isoR)
                    if (xs > 0.) then
                      write(iUall, '(4x,a12,2es15.5,3x,a)') ext(i), xs ,dxs, trim(ref)
                      filestatus(iUall) = 'keep'
                    endif
                  enddo
                  N = Ntherm_xs(type, Z, A, isoT, isoR)
                  do i = 1, N
                    xs = Etherm_xs(type, Z, A, isoT, isoR, i)
                    dxs = Etherm_dxs(type, Z, A, isoT, isoR, i)
                    ref = Etherm_ref(type, Z, A, isoT, isoR, i)
                    write(iUall, '(4x,a12,2es15.5,3x,a)') ext(numlib), xs ,dxs, trim(ref)
                    filestatus(iUall) = 'keep'
                  enddo
                  do j = 1, numndlib
                    xs = Ltherm_xs(j, type, Z, A, isoT, isoR)
                    ratio = Lratio_xs(j, type, Z, A, isoT, isoR)
                    write(iUall, '(4x,a10,2x,es15.5,11x,"CE=",es15.5)') ndlib(j), xs, ratio
                    filestatus(iUall) = 'keep'
                  enddo
                  if (type /= 4 .and. type <= 6) then
                    ratio = capratio_xs(type, Z, A, isoT, isoR)
                    if (ratio > 0.) then
                      write(8, '(a4,2i4,3es15.5,3x,a)') Zstring, A, isoT, ratio
                    endif
                  endif
                endif
              endif
            enddo
          enddo
        enddo
      enddo
      close (11, status = filestatus(11))
      close (21, status = filestatus(21))
      if (type == 4) then
        close (12, status = filestatus(12))
        close (13, status = filestatus(13))
        close (14, status = filestatus(14))
        close (22, status = filestatus(22))
        close (23, status = filestatus(23))
        close (24, status = filestatus(24))
      endif
    enddo
    close (1)
    if (type == 4) then
      close (2)
      close (3)
      close (4)
    else
      close (8)
    endif
  enddo
  return
end subroutine writethermal
! Copyright A.J. Koning 2019
