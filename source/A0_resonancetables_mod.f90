module A0_resonancetables_mod
!                                                                                                                                   
!-----------------------------------------------------------------------------------------------------------------------------------
! Purpose: All global variables                                                                                                     
!                                                                                                                                   
! Revision    Date      Author           Description                                                                                
! ====================================================                                                                              
!    1     2025-02-27   A.J. Koning      Original code                                                                              
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
  integer, parameter :: numpar=7    ! maximum number of particles                                                               
  integer, parameter :: numZ=100    ! maximum number of elements                                                                
  integer, parameter :: numA=257    ! maximum number of masses                                                                  
  integer, parameter :: numisom=2   ! maximum number of isomers                                                                 
  integer, parameter :: numlib=6    ! maximum number of databases
  integer, parameter :: numndlib=5  ! maximum number of ND libraries                                                          
  integer, parameter :: numtype=9   ! number of reaction types
  integer, parameter :: numex=500   ! number of subentries per reaction
  integer, parameter :: numdat=1000 ! number of total Z, A data points
!                                                                                                                                   
! machine                                                                                                                           
!                                                                                                                                   
  character(len=132) :: filespath   ! directory containing structure files to be read
  character(len=132) :: libspath    ! directory containing data from nuclear data libraries
  character(len=132) :: exforpath   ! directory containing X4 data
  character(len=132) :: resbasepath ! directory containing TARES data
  character(len=132) :: thermalpath ! directory containing data libraries
  character(len=132) :: macspath    ! directory containing data libraries
  character(len=132) :: respath     ! directory containing data libraries
!
! thermalbaseinitial
!
  character(len=2)         :: nuc(numZ)         ! symbol of nucleus                                                
  character(len=3)         :: reac(numtype)     ! reaction type
  character(len=10)        :: restype(numtype)  ! resonance type
  character(len=20)        :: reaction(numtype) ! reaction
  character(len=10)        :: ndlib(numndlib)   ! ND library
  integer                  :: ndyear(numndlib)  ! year of ND library
  character(len=132)       :: source            ! source of data
  character(len=132)       :: oformat           ! format of data
  character(len=132)       :: user              ! user of data
  integer, dimension(numZ) :: heavy             ! heaviest isotope
  integer, dimension(numZ) :: light             ! lightest isotope
  character(len=10)        :: date              ! date
!                                                                                                                                   
! readresonance
!                                                                                                                                   
  logical                             :: res_exist ! flag for existence
  integer, dimension(numex)           :: res_year ! res year
  character(len=40), dimension(numex) :: res_ref ! reference
  character(len=40), dimension(numex) :: res_type ! type
  character(len=24), dimension(numex) :: res_author ! EXFOR res author
  integer                             :: Nres ! number of cases
  integer                             :: Nres_exp ! number of EXFOR cases
  real, dimension(numex)              :: res_xs  ! res parameter
  real, dimension(numex)              :: res_dxs ! res parameter uncertainty
!                                                                                                                                   
! procres                                                                                                                       
!                                                                                                                                   
  character(len=40)                   :: res_author_sel 
  real                                :: res_xs_sel  
  real                                :: res_dxs_sel 
!                                                                                                                                   
! readthermal                                                                                                                       
!                                                                                                                                   
  character(len=40), dimension(numdat) :: refsave
  integer, dimension(numdat)           :: Zsave
  integer, dimension(numdat)           :: Asave
  integer, dimension(numdat)           :: Lisosave
  integer, dimension(numdat)           :: Nexpsave
  real, dimension(numdat)              :: xssave
  real, dimension(numdat)              :: dxssave
  integer                              :: Nsave
!                                                                                                                                   
! writethermal                                                                                                                       
  integer                             :: Ztarget
  integer                             :: Atarget
  character(len=6)                    :: targetnuclide    
end module A0_resonancetables_mod                                                                                          
!Copyright (C) 2025  A.J. Koning                                                                                                    
