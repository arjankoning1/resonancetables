program resonancetables
!
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose: Unify databases for thermal cross sections, resonance integrals, average resonance parameters, and MACS
!
! Revision    Date      Author      Quality  Description
! ======================================================
!    1     2025-08-01   A.J. Koning    A     Original code
!-----------------------------------------------------------------------------------------------------------------------------------
!
!   |-------------------------------------------------------|
!   |                 Arjan Koning                          |
!   |                                                       |
!   | Email: A.Koning@@iaea.org                             |
!   |-------------------------------------------------------|
!
! MIT License
!
! Copyright (c) 2025 Arjan Koning
!
! Permission is hereby granted, free of charge, to any person obtaining a copy
! of this software and associated documentation files (the "Software"), to deal
! in the Software without restriction, including without limitation the rights
! to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
! copies of the Software, and to permit persons to whom the Software is
! furnished to do so, subject to the following conditions:
!
! The above copyright notice and this permission notice shall be included in all
! copies or substantial portions of the Software.
!
! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
! IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
! FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
! AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
! LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
! OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
! SOFTWARE.
!
! *** Use data from other modules
!
  use A0_resonancetables_mod
!
! *** Declaration of local data
!
  implicit none
!
! **************** Read data and process it into new database **********
!
  call machine
  call resonancetablesinitial
  call thermal
  call macs
  call resonance
end program resonancetables
! Copyright A.J. Koning 2025
