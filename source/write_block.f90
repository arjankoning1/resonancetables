  subroutine write_header(title,source,user,date,comment)
!
! +---------------------------------------------------------------------
! | Author: Arjan Koning
! | Date  : August 8, 2023
! | Task  : Write output block for title
! +---------------------------------------------------------------------
!
! ****************** Declarations **************************************
!
  implicit none
  character(len=*) title
  character(len=*) source
  character(len=*) user
  character(len=*) date
  character(len=*) comment
  integer          indent
!
! ************* Write output block *************************************
!
  write(1,'("# header:")')
  indent=2
  call write_char(indent,'title',trim(title))
  call write_char(indent,'source',trim(source))
  if (len_trim(user) > 0) call write_char(indent,'user',trim(user))
  call write_char(indent,'date',trim(date))
  if (len_trim(comment) > 0) call write_char(indent,'format',trim(comment))
  return
  end
!Copyright (C)  2023 A.J. Koning
  subroutine write_target
!
! +---------------------------------------------------------------------
! | Author: Arjan Koning
! | Date  : May 20, 2023
! | Task  : Write output block for target
! +---------------------------------------------------------------------
!
! *** Use data from other modules
!
  use A0_resonancetables_mod
!
! ****************** Declarations **************************************
!
  implicit none
! integer          Zixx
! integer          Nixx
! integer          ilev
! integer          P
! real             J
! real             Ee
! real             tau0
!
! ************* Write output block *************************************
!
  write(1,'("# target:")')
  call write_nuc(2,Ztarget,Atarget,targetnuclide)
! if (Liso > 0) then
!   Zixx=parZ(k0)
!   Nixx=parN(k0)
!   ilev=Lisomer(Zixx, Nixx, Liso)
!   Ee=edis(Zixx,Nixx,ilev)
!   J=jdis(Zixx,Nixx,ilev)
!   P=parlev(Zixx,Nixx,ilev)
!   tau0=tau(Zixx,Nixx,ilev)
!   call write_level(2,Liso,ilev,Ee,J,P,tau0)
! endif
  return
  end
!Copyright (C)  2023 A.J. Koning
  subroutine write_reaction(reaction,Q,E,MF,MT)
!
! +---------------------------------------------------------------------
! | Author: Arjan Koning
! | Date  : May 20, 2023
! | Task  : Write output block for reaction
! +---------------------------------------------------------------------
!
! ****************** Declarations **************************************
!
  implicit none
  character(len=*) reaction
  integer          indent
  integer          MF
  integer          MT
  real*8           Q
  real*8           E
!
! ************* Write output block *************************************
!
  write(1,'("# reaction:")')
  indent=2
  call write_char(indent,'type',trim(reaction))
  if (Q /= 0.D0 .or. E /= 0.D0) then
    call write_double(indent,'Q-value [MeV]',Q)
    call write_double(indent,'E-threshold [MeV]',E)
  endif
  return
  end
!Copyright (C)  2023 A.J. Koning
  subroutine write_nuc(indent,Z,A,nuc)
!
! +---------------------------------------------------------------------
! | Author: Arjan Koning
! | Date  : May 20, 2023
! | Task  : Write output block for nuclide
! +---------------------------------------------------------------------
!
! ****************** Declarations **************************************
!
  implicit none
  character(len=*) nuc
  integer          indent
  integer          Z
  integer          A
!
! ************* Write output block *************************************
!
  call write_integer(indent,'Z',Z)
  call write_integer(indent,'A',A)
  call write_char(indent,'nuclide',trim(nuc))
  return
  end
!Copyright (C)  2023 A.J. Koning
  subroutine write_quantity(quantity)
! 
! +---------------------------------------------------------------------
! | Author: Arjan Koning
! | Date  : February 23, 2025
! | Task  : Write output block for quantity
! +---------------------------------------------------------------------
! 
! ****************** Declarations **************************************
! 
  implicit none
  character(len=*) quantity
  integer          indent
!
! ************* Write output block *************************************
!
  write(1,'("# quantity:")')
  indent=2
  call write_char(indent,'type',trim(quantity))
  return
  end
!Copyright (C)  2025 A.J. Koning
  subroutine write_datablock(Nc,Ne,col,un)
!
! +---------------------------------------------------------------------
! | Author: Arjan Koning
! | Date  : August 8, 2023
! | Task  : Write output block for data block
! +---------------------------------------------------------------------
!
! ****************** Declarations **************************************
!
  implicit none
  integer            Ne
  integer            Nc
  integer            L
  integer            i
  integer            indent
  integer            k
  integer            ibeg
  integer            width
  character(len=*)   col(Nc)
  character(len=*)   un(Nc)
  character(len=15)  word
  character(len=3000) obs
  character(len=3000) units
!
! ************* Write output block *************************************
!
  write(1,'("# datablock:")')
  indent=2
  call write_integer(indent,'columns',Nc)
  call write_integer(indent,'entries',Ne)
  width=15
  obs="##"
  k=2
  do i=1,Nc
    word=''
    L=len_trim(col(i))
    ibeg=max(width/2-L/2,0)
    word(ibeg+1:ibeg+L)=col(i)(1:L)
    obs=obs(1:k)//word
    k=k+width
  enddo
  write(1, '(a)') trim(obs)
  units="##"
  k=2
  do i=1,Nc
    word=''
    L=len_trim(un(i))
    ibeg=max(width/2-L/2,0)
    word(ibeg+1:ibeg+L)=un(i)(1:L)
    word(ibeg:ibeg)='['
    word(ibeg+L+1:ibeg+L+1)=']'
    units=units(1:k)//word
    k=k+width
  enddo
  write(1, '(a)') trim(units)
  return
  end
!Copyright (C)  2023 A.J. Koning
  subroutine write_real(indent,key,x)
!
! +---------------------------------------------------------------------
! | Author: Arjan Koning
! | Date  : May 20, 2023
! | Task  : Write output block for real number
! +---------------------------------------------------------------------
!
! ****************** Declarations **************************************
!
  implicit none
  character(len=*)  key
  character(len=80) fmt
  integer           indent
  real              x
!
! ************* Write output block *************************************
!
  fmt = '("# ",    a,": ",es13.6)'
  if (indent > 0 ) then
    fmt = '("# ",  x,a,": ",es13.6)'
    write(fmt(7:8),'(i2)') indent
  endif
  write(1, fmt) trim(key),x
  return
  end
!Copyright (C)  2023 A.J. Koning
  subroutine write_double(indent,key,x)
!
! +---------------------------------------------------------------------
! | Author: Arjan Koning
! | Date  : May 20, 2023
! | Task  : Write output block for double precision number
! +---------------------------------------------------------------------
!
! ****************** Declarations **************************************
!
  implicit none
  character(len=*)  key
  character(len=80) fmt
  integer           indent
  real*8            x
!
! ************* Write output block *************************************
!
  fmt = '("# ",    a,": ",es13.6)'
  if (indent > 0 ) then
    fmt = '("# ",  x,a,": ",es13.6)'
    write(fmt(7:8),'(i2)') indent
  endif
  write(1, fmt) trim(key),x
  return
  end
!Copyright (C)  2023 A.J. Koning
  subroutine write_integer(indent,key,k)
!
! +---------------------------------------------------------------------
! | Author: Arjan Koning
! | Date  : May 20, 2023
! | Task  : Write output block for integer
! +---------------------------------------------------------------------
!
! ****************** Declarations **************************************
!
  implicit none
  character(len=*)  key
  character(len=80) fmt
  integer           indent
  integer           k
!
! ************* Write output block *************************************
!
  fmt = '("# ",    a,": ",i0)'
  if (indent > 0 ) then
    fmt = '("# ",  x,a,": ",i0)'
    write(fmt(7:8),'(i2)') indent
  endif
  write(1, fmt) trim(key),k
  return
  end
!Copyright (C)  2023 A.J. Koning
  subroutine write_char(indent,key,word)
!
! +---------------------------------------------------------------------
! | Author: Arjan Koning
! | Date  : May 20, 2023
! | Task  : Write output block for character
! +---------------------------------------------------------------------
!
! ****************** Declarations **************************************
!
  implicit none
  character(len=*)  key
  character(len=*)  word
  character(len=80) fmt
  integer           indent
!
! ************* Write output block *************************************
!
  fmt = '("# ",    a,": ",a)'
  if (indent > 0 ) then
    fmt = '("# ",  x,a,": ",a)'
    write(fmt(7:8),'(i2)') indent
  endif
  write(1, fmt) trim(key),trim(word)
  return
  end
!Copyright (C)  2023 A.J. Koning
