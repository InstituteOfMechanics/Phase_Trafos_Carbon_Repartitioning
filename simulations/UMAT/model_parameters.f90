! TODO: License etc

!=======================================================================
!
! This file contains a custom data type that holds all parameters for
! the material model, as well as a subroutine to set the values for
! these parameters according to the chosen material (100Cr6 or 
! 100CrMnSi6-4).
!
!=======================================================================

!=======================================================================
!
!        CUSTOM DATA TYPE TO HOLD CONSTANT MATERIAL PARAMETERS
!
!=======================================================================
module material_data_mod
    implicit none
    type, public :: material_data_type
        ! This data type is used to store all material model parameters
        ! Explanations for each value are given in the subroutine
        ! where the values are set
        double precision :: xC, xCB, density
        
        ! flag to disable evolution of martensite/bainite
        integer :: martensite_flag, bainite_flag
          
        ! KM model parameters
        double precision :: km_regularization_width
        double precision :: km_p0k, km_e1k, km_e2k, km_p0Ms, km_p1MS
          
        ! JMAK model parameters
        double precision :: jmak_nucleation_fraction
        double precision :: jmak_regularization_width
        double precision :: jmak_p0Bs, jmak_p1BS
        double precision :: jmak_p0logb, jmak_p1logb, jmak_p2logb
        double precision :: jmak_p0N, jmak_p1N, jmak_p2N
          
        ! mechanical parameters
        double precision :: mech_p0EM, mech_p1EM, mech_p0EB, &
                            mech_p1EB, mech_p0EA, mech_p1EA
		double precision :: mech_p0nuM, mech_p1nuM, mech_p0nuB, &
		                    mech_p1nuB, mech_p0nuA, mech_p1nuA
		double precision :: mech_p0alphaM, mech_p1alphaM, mech_p0alphaB, &
		                    mech_p1alphaB, mech_p0alphaA, mech_p1alphaA
		double precision :: mech_p0gammaM, mech_p1gammaM, &
		                    mech_p0gammaB, mech_p1gammaB
          
        ! latent heat parameters
        double precision :: lath_p0rhoDeltahAM, lath_p1rhoDeltahAM, &
                            lath_p0rhoDeltahAB, lath_p1rhoDeltahAB
          
        ! thermal parameters
        double precision :: thermal_p0cpM, thermal_p1cpM, &
                            thermal_p2cpM, thermal_p3cpM
		double precision :: thermal_p0cpB, thermal_p1cpB, &
		                    thermal_p2cpB, thermal_p3cpB
        double precision :: thermal_p0cpA, thermal_p1cpA, &
                            thermal_p2cpA, thermal_p3cpA
        double precision :: thermal_p0lambM, thermal_p1lambM, thermal_p2lambM
        double precision :: thermal_p0lambB, thermal_p1lambB, thermal_p2lambB
        double precision :: thermal_p0lambA, thermal_p1lambA, thermal_p2lambA
          
    end type material_data_type
end module material_data_mod

!=======================================================================
!
!    SUBROUTINE TO SET MATERIAL PARAMETERS BASED ON MATERIAL FLAG
!
!=======================================================================
subroutine get_material_data(material_data, material_flag, density, &
                             martensite_flag, bainite_flag)
	! Set a material_data_type container with all model parameters. 
	!
	! Input variables
	! ==================================================================
	! material_flag   : 0 for 100Cr6, 1 for 100CrMnSi6-4
	! density         : mass density (must be specified as a property)
	! martensite_flag : 1 if martensite evolution is suppressed
	! bainite_flag    : 1 if bainite evolution is suppressed
	!
	! Output variables
	! ==================================================================
	! material_data :   Hold all parameters for the material model
	
	use material_data_mod
	implicit none
	
	type(material_data_type), intent(out) :: material_data
	integer, intent(in) :: material_flag, martensite_flag, bainite_flag
	double precision, intent(in) :: density
	
	! flags
	material_data%martensite_flag = martensite_flag
	material_data%bainite_flag = bainite_flag
	
	! mass density (set based on definition in Abaqus)
	material_data%density = density
	
	! nominal carbon content
	material_data%xC = 1.0d-2
	
	! constant bainite phase carbon content
	if (material_flag == 0) then
		! 100Cr6
		material_data%xCB = material_data%xC
	else
		! 100CrMnSi6-4
		material_data%xCB = 3.0d-4
	end if 
	
	!===================================================================
	!               Parameters for martensite evolution
	!===================================================================
	
	! regularization width for activation function
	material_data%km_regularization_width = 4.0d0
    
    ! parametrization coefficients for Koistinen Marburger parameters
    if (material_flag == 0) then
		! 100Cr6
		material_data%km_p0Ms = 526.29d0
		material_data%km_p1MS = -29702.0d0
		
		material_data%km_p0k = 7.14d-3
		material_data%km_e1k = 0.0198d0
		material_data%km_e2k = -156.0d0
	else
		! 100CrMnSi6-4
		material_data%km_p0Ms = 492.71d0
		material_data%km_p1MS = -28290.6d0
		
		material_data%km_p0k = 6.88d-3
		material_data%km_e1k = 0.0198d0
		material_data%km_e2k = -156.0d0
	end if
    
    !===================================================================
	!               Parameters for bainite evolution
	!===================================================================
	
	! bainite fraction at which nucleation is considered complete
	material_data%jmak_nucleation_fraction = 1.0d-2
	
	! regularization width for activation function
	material_data%jmak_regularization_width = 4.0d0
    
    ! polynomial coefficients for JMAK parameters
    if (material_flag == 0) then
		! 100Cr6
		material_data%jmak_p0Bs = 681.51d0
		material_data%jmak_p1BS = -29618.0d0
	else
		! 100CrMnSi6-4
		material_data%jmak_p0Bs = 646.95d0
		material_data%jmak_p1BS = -28142.8d0
	end if
    
    
    material_data%jmak_p0logb = -1.8847d2
    material_data%jmak_p1logb = 8.1383d-1
    material_data%jmak_p2logb = -1.0303d-3
    
    material_data%jmak_p0N = 1.1313d1
    material_data%jmak_p1N = -4.5363d-2
    material_data%jmak_p2N = 5.8965d-5
          
    !===================================================================
	!               Parameters for stress computation
	!===================================================================
	
	! polynomial coefficients for Young's modulus (for each phase)
    material_data%mech_p0EM = 212000.0d0
    material_data%mech_p1EM = -52.7d0
    material_data%mech_p0EB = 212000.0d0
    material_data%mech_p1EB = -52.7d0
    material_data%mech_p0EA = 204000.0d0
    material_data%mech_p1EA = -91.0d0
    
	! polynomial coefficients for Poisson's ratio (for each phase)
    material_data%mech_p0nuM = 0.344d0
    material_data%mech_p1nuM = 0.00010d0
    material_data%mech_p0nuB = 0.344d0
    material_data%mech_p1nuB = 0.00010d0
    material_data%mech_p0nuA = 0.223d0
    material_data%mech_p1nuA = 0.00025d0
    
	! polynomial coefficients for linear heat expansion coefficient
	! (for each phase)
	material_data%mech_p0alphaM = 1.17d-5
	material_data%mech_p1alphaM = 0.0d0
	material_data%mech_p0alphaB = 1.17d-5
	material_data%mech_p1alphaB = 0.0d0
	material_data%mech_p0alphaA = 2.22d-5
	material_data%mech_p1alphaA = 0.0d0
	
	! polynomial coefficients for linear transfomation strain
	! (for each phase)
	material_data%mech_p0gammaM = 1.11d-2
	material_data%mech_p1gammaM = 0.0d0
	material_data%mech_p0gammaB = 5.0d-3
	material_data%mech_p1gammaB = 0.0d0
    
    
    !===================================================================
	!              Parameters for latent heat computation
	!===================================================================
	
	! polynomial coefficients for austenite -> martensite trafo
    material_data%lath_p0rhoDeltahAM = 6.4d2
    material_data%lath_p1rhoDeltahAM = 0.0d0
    
    ! polynomial coefficients for austenite -> bainite trafo
    material_data%lath_p0rhoDeltahAB = 1.56d3
    material_data%lath_p1rhoDeltahAB = -1.5d0
    
    !===================================================================
	!               Parameters for thermal model
	!===================================================================
	      
    ! polynomial coefficients for heat capacity (for all phases)
    material_data%thermal_p0cpM = 336.51d0
    material_data%thermal_p1cpM = 0.4875d0
    material_data%thermal_p2cpM = -6.16d-4
    material_data%thermal_p3cpM = 1.62d-6
    material_data%thermal_p0cpB = 336.51d0
    material_data%thermal_p1cpB = 0.4875d0
    material_data%thermal_p2cpB = -6.16d-4
    material_data%thermal_p3cpB = 1.62d-6
    material_data%thermal_p0cpA = 448.67d0
    material_data%thermal_p1cpA = 0.3007d0
    material_data%thermal_p2cpA = -2.49d-4
    material_data%thermal_p3cpA = 1.42d-7
    
    ! polynomial coefficients for thermal conductivity (for all phases)
    material_data%thermal_p0lambM = 3.93d-2
    material_data%thermal_p1lambM = 2.4d-5
    material_data%thermal_p2lambM = -6.43d-8
    material_data%thermal_p0lambB = 3.93d-2
    material_data%thermal_p1lambB = 2.4d-5
    material_data%thermal_p2lambB = -6.43d-8
    material_data%thermal_p0lambA = 1.68d-2
    material_data%thermal_p1lambA = 1.19d-5
    material_data%thermal_p2lambA = 0.0d0
	
end subroutine get_material_data 
