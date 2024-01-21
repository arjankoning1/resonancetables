module A0_resonancetables_mod
!                                                                                                                                   
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose: All global variables                                                                                                     
!                                                                                                                                   
! Revision    Date      Author           Description                                                                                
! ====================================================                                                                              
!    1     2023-12-29   A.J. Koning      Original code                                                                              
!-----------------------------------------------------------------------------------------------------------------------------------
!                                                                                                                                   
! *** Use data from other modules                                                                                                   
!                                                                                                                                   
  use A0_kinds_mod, only:     & ! Definition of single and double precision variables                                               
                 sgl            !  single precision kind                                                                            
!                                                                                                                                   
! *** Declaration of global data                                                                                                    
!                                                                                                                                   
!-----------------------------------------------------------------------------------------------------------------------------------
!                                                                                                                                   
! Array dimensions                                                                                                                  
!                                                                                                                                   
  integer, parameter :: numpar=7        ! maximum number of particles                                                               
  integer, parameter :: numZ=100        ! maximum number of elements                                                                
  integer, parameter :: numA=257        ! maximum number of masses                                                                  
  integer, parameter :: numisom=2       ! maximum number of isomers                                                                 
  integer, parameter :: numlib=6        ! maximum number of data libraries                                                          
  integer, parameter :: numndlib=5      ! maximum number of ND libraries                                                          
  integer, parameter :: numtype=9       ! number of reaction types
  integer, parameter :: numex=50        ! number of EXFOR subentries per reaction
!                                                                                                                                   
! machine                                                                                                                           
!                                                                                                                                   
  character(len=132) :: thermalpath ! directory containing data libraries
  character(len=132) :: macspath    ! directory containing data libraries
  character(len=132) :: respath     ! directory containing data libraries
  character(len=132) :: filespath   ! directory containing X4 and structure files to be read
  character(len=132) :: libspath    ! directory containing data from nuclear data libraries
!
! thermalbaseinitial
!
  character(len=2)   :: nuc(numZ)        ! symbol of nucleus                                                
  character(len=3)   :: reac(numtype)    ! reaction type
  character(len=6)   :: restype(numtype) ! resonance type
  character(len=20)  :: ext(0:numlib)    ! library
  character(len=10)  :: ndlib(numndlib)  ! ND library
!                                                                                                                                   
! readthermal                                                                                                                       
!                                                                                                                                   
  logical, dimension(0:numlib, numtype, -1:numisom)                                    :: therm_exist ! flag for existence
  logical, dimension(numtype, -1:numisom)                                              :: Etherm_exist ! flag for existence
  character(len=4), dimension(numtype, numZ, 0:numA, -1:numisom, -1:numisom, numex)    :: Etherm_year ! EXFOR thermal year
  character(len=40), dimension(0:numlib, numtype, numZ, 0:numA, -1:numisom, -1:numisom):: therm_ref ! reference
  character(len=40), dimension(numtype, numZ, 0:numA, -1:numisom, -1:numisom, numex)   :: Etherm_ref ! reference
  character(len=9), dimension(numtype, numZ, 0:numA, -1:numisom, -1:numisom, numex)    :: Etherm_subentry ! EXFOR thermal subentry
  character(len=24), dimension(numtype, numZ, 0:numA, -1:numisom, -1:numisom, numex)   :: Etherm_author ! EXFOR thermal author
  integer, dimension(numtype, numZ, 0:numA, -1:numisom, -1:numisom)                    :: Ntherm_xs ! number of EXFOR cases
  real, dimension(0:numlib, numtype, numZ, 0:numA, -1:numisom, -1:numisom)             :: therm_xs  ! thermal cross section
  real, dimension(0:numlib, numtype, numZ, 0:numA, -1:numisom, -1:numisom)             :: ratio_xs  ! ratio of thermal cross section
  real, dimension(numndlib, numtype, numZ, 0:numA, -1:numisom, -1:numisom)             :: Lratio_xs  ! NDL ratio of thermal cross section
  real, dimension(numtype, numZ, 0:numA, -1:numisom, -1:numisom)                       :: capratio_xs ! ratio to (n,g) thermal cross s.
  real, dimension(0:numlib, numtype, numZ, 0:numA, -1:numisom, -1:numisom)          :: therm_dxs ! thermal cross section uncertainty
  real, dimension(numtype, numZ, 0:numA, -1:numisom, -1:numisom, numex)                :: Etherm_xs ! EXFOR thermal cross section
  real, dimension(numtype, numZ, 0:numA, -1:numisom, -1:numisom, numex)      :: Etherm_dxs ! EXFOR thermal cross section uncertainty
  real, dimension(numtype, numZ, 0:numA, -1:numisom, -1:numisom, numex)                :: Eratio_xs ! EXFOR thermal ratio
  real, dimension(numndlib, numtype, numZ, 0:numA, -1:numisom, -1:numisom)             :: Ltherm_xs  ! NDL thermal cross section
!                                                                                                                                   
! readmacs                                                                                                                       
!                                                                                                                                   
  logical, dimension(0:numlib, numtype, -1:numisom)                                    :: macs_exist ! flag for existence
  logical, dimension(numtype, -1:numisom)                                              :: Emacs_exist ! flag for existence
  character(len=4), dimension(numtype, numZ, 0:numA, -1:numisom, -1:numisom, numex)    :: Emacs_year ! EXFOR MACS year
  character(len=40), dimension(0:numlib, numtype, numZ, 0:numA, -1:numisom, -1:numisom):: macs_ref ! reference
  character(len=40), dimension(numtype, numZ, 0:numA, -1:numisom, -1:numisom, numex)   :: Emacs_ref ! reference
  character(len=9), dimension(numtype, numZ, 0:numA, -1:numisom, -1:numisom, numex)    :: Emacs_subentry ! EXFOR MACS subentry
  character(len=24), dimension(numtype, numZ, 0:numA, -1:numisom, -1:numisom, numex)   :: Emacs_author ! EXFOR MACS author
  integer, dimension(numtype, numZ, 0:numA, -1:numisom, -1:numisom)                    :: Nmacs_xs ! number of EXFOR cases
  real, dimension(0:numlib, numtype, numZ, 0:numA, -1:numisom, -1:numisom)             :: macs_xs  ! MACS cross section
  real, dimension(0:numlib, numtype, numZ, 0:numA, -1:numisom, -1:numisom)             :: ratiomacs_xs  ! ratio of MACS cross section
  real, dimension(0:numlib, numtype, numZ, 0:numA, -1:numisom, -1:numisom)          :: macs_dxs ! MACS cross section uncertainty
  real, dimension(numtype, numZ, 0:numA, -1:numisom, -1:numisom, numex)                :: Emacs_xs ! EXFOR MACS cross section
  real, dimension(numtype, numZ, 0:numA, -1:numisom, -1:numisom, numex)      :: Emacs_dxs ! EXFOR MACS cross section uncertainty
  real, dimension(numtype, numZ, 0:numA, -1:numisom, -1:numisom, numex)                :: Eratiomacs_xs ! EXFOR MACS ratio
  real, dimension(numndlib, numtype, numZ, 0:numA, -1:numisom, -1:numisom)             :: Lmacs_xs  ! NDL MACS
  real, dimension(numndlib, numtype, numZ, 0:numA, -1:numisom, -1:numisom)             :: Lratiomacs_xs  ! NDL ratio of MACS
!                                                                                                                                   
! readresonance
!                                                                                                                                   
  logical, dimension(0:numlib, numtype)                         :: res_exist ! flag for existence
  logical, dimension(numtype)                                   :: Eres_exist ! flag for existence
  character(len=40), dimension(0:numlib, numtype, numZ, 0:numA) :: res_ref ! reference
  real, dimension(0:numlib, numtype, numZ, 0:numA)              :: res_R  ! 
  real, dimension(0:numlib, numtype, numZ, 0:numA)              :: res_dR ! 
  real, dimension(0:numlib, numtype, numZ, 0:numA)              :: ratio_R  ! 
end module A0_resonancetables_mod                                                                                          
!Copyright (C) 2019  A.J. Koning                                                                                                    
