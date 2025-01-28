!=======================================================================
! This file contains the constitutive laws for the stresses and heat
! generation. All routines assume that the internal variables have
! already been updated.
!
! Voigt notation is used for the stresses/strains. Abaqus has a
! different order for the shear components, but since our material
! model is isotropic is does not change anything.
!=======================================================================
subroutine calculate_stresses(stress, DdstressDdstrain, DdstressDtemp, &
							  DdstressDkappa, dstrainTh, dstrainTr, &
							  ntens, ndi, nkappa, strain, strainTh, &
							  strainTr, dstrain, kappa, kappan, tempn, &
							  dtemp, mdata)
	! Calculate the Cauchy stress in Voigt notation, as well as the
	! partial derivatives w.r.t. strains, temperature, and SDVs.
	!
	! Variables passed in
	! ==================================================================
	! ntens             : Number of tensor components
	! ndi				: Number of direct tensor components
	! nkappa			: Number of components in the SDV vector kappa
	! strain			: Strain at beginning of increment (Voigt notation)
	! strainTh          : Linear thermal strain at beginning of increment
	! strainTr          : Linear transformation strain at beginning of increment
	! dstrain			: Strain increment
	! kappa				: SDVs at the end of the increment
	! kappan			: SDVs at the beginning of the increment
	! tempn				: Temperature at the beginning of the increment
	! dtemp				: Temperature increment
	! mdata             : All material model parameters
	!
	! Variables set by this subroutine
	! ==================================================================
	! stress            : Stress at the end of the increment
	! DdstressDdstrain  : Derivative of the stress increment w.r.t. strain increment
	! DdstressDtemp 	: Derivative of the stress increment w.r.t. temperature
	! DdstressDkappa	: Derivative of the stress increment w.r.t. SDVs
	! dstrainTh			: Increment of thermal strain (scalar)
	! dstrainTr			: Increment of transformation strain (scalar)
    use material_data_mod
    
    implicit none
    
    integer, intent(in) :: ntens, ndi, nkappa
    double precision, intent(in) :: strain(ntens), strainTh, strainTr, dstrain(ntens), kappa(nkappa), kappan(nkappa), tempn, dtemp
    double precision, intent(out) :: stress(ntens), DdstressDdstrain(ntens, ntens), DdstressDtemp(ntens), DdstressDkappa(ntens, nkappa), dstrainTh, dstrainTr
    type(material_data_type), intent(in) :: mdata
    
    double precision :: temp, strainElN(ntens), strainEl(ntens), dstrainEl(ntens), DdstrainElDtemp(ntens), DdstrainElDbetaM(ntens), DdstrainElDbetaB(ntens)
    double precision :: DdstrainThDtemp, DdstrainThDbetaM, DdstrainThDbetaB, DdstrainTrDtemp, DdstrainTrDbetaM, DdstrainTrDbetaB
    double precision :: CM(ntens,ntens), CB(ntens,ntens), CA(ntens,ntens), C(ntens,ntens)
    double precision :: dCMdE(ntens,ntens), dCMdnu(ntens,ntens), dCBdE(ntens,ntens), dCBdnu(ntens,ntens), dCAdE(ntens,ntens), dCAdnu(ntens,ntens)
    double precision :: dCMdtemp(ntens,ntens), dCBdtemp(ntens,ntens), dCAdtemp(ntens,ntens), dCdtemp(ntens,ntens)
    double precision :: betaM, betaB, betaA, betaMn, betaBn, betaAn
    double precision :: nuM, nuB, nuA, dnuMdtemp, dnuBdtemp, dnuAdtemp
    double precision :: EM, EB, EA, dEMdtemp, dEBdtemp, dEAdtemp
    double precision :: alphaM, alphaB, alphaA, alpha, dalphaMdtemp, dalphaBdtemp, dalphaAdtemp, dalphadtemp
    double precision :: gammaM, gammaB, dgammaMdtemp, dgammaBdtemp

    ! get the current state
    temp = tempn + dtemp
    
    betaM = kappa(1)
    betaB = kappa(2)
    betaA = 1.0d0 - betaM - betaB
    
    betaMn = kappan(1)
    betaBn = kappan(2)
    betaAn = 1.0d0 - betaMn - betaBn
    
    ! evaluate parameterization for temp-dependent material parameters
    ! and their derivatives
    EM = mdata%mech_p0EM + mdata%mech_p1EM*temp
    dEMdtemp = mdata%mech_p1EM
	EB = mdata%mech_p0EB + mdata%mech_p1EB*temp
    dEBdtemp = mdata%mech_p1EB
    EA = mdata%mech_p0EA + mdata%mech_p1EA*temp
    dEAdtemp = mdata%mech_p1EA

    nuM = mdata%mech_p0nuM + mdata%mech_p1nuM*temp
    dnuMdtemp = mdata%mech_p1nuM
    nuB = mdata%mech_p0nuB + mdata%mech_p1nuB*temp
    dnuBdtemp = mdata%mech_p1nuB
    nuA = mdata%mech_p0nuA + mdata%mech_p1nuA*temp
    dnuAdtemp = mdata%mech_p1nuA
	
	alphaM = mdata%mech_p0alphaM + mdata%mech_p1alphaM*temp
	dalphaMdtemp = mdata%mech_p1alphaM
	alphaB = mdata%mech_p0alphaB + mdata%mech_p1alphaB*temp
	dalphaBdtemp = mdata%mech_p1alphaB
	alphaA = mdata%mech_p0alphaA + mdata%mech_p1alphaA*temp
	dalphaAdtemp = mdata%mech_p1alphaA
	
	alpha = betaM*alphaM + betaB*alphaB + betaA*alphaA
	dalphadtemp = betaM*dalphaMdtemp + betaB*dalphaBdtemp + betaA*dalphaAdtemp
	
	gammaM = mdata%mech_p0gammaM + mdata%mech_p1gammaM*temp
	dgammaMdtemp = mdata%mech_p1gammaM
    gammaB = mdata%mech_p0gammaB + mdata%mech_p1gammaB*temp
    dgammaBdtemp = mdata%mech_p1gammaB
    
    ! construct stiffness tensor for each phase...
	call get_stiffness_tensor(CM, dCMdE, dCMdnu, ntens, ndi, EM, nuM)    
	call get_stiffness_tensor(CB, dCBdE, dCBdnu, ntens, ndi, EB, nuB)
	call get_stiffness_tensor(CA, dCAdE, dCAdnu, ntens, ndi, EA, nuA)
	
	! ...and their temperature drivatives
	DCMDtemp(:,:) = dCMdE(:,:)*dEMdtemp + dCMdnu(:,:)*dnuMdtemp
	DCBDtemp(:,:) = dCBdE(:,:)*dEBdtemp + dCBdnu(:,:)*dnuBdtemp
	DCADtemp(:,:) = dCAdE(:,:)*dEAdtemp + dCAdnu(:,:)*dnuAdtemp
	
	C(:,:) = betaM*CM(:,:) + betaB*CB(:,:) + betaA*CA(:,:)
	dCdtemp(:,:) = betaM*dCMdtemp(:,:) + betaB*dCBdtemp(:,:) + betaA*dCAdtemp(:,:)
	
    ! scalar thermal strain increment, and derivatives w.r.t. temperature and SDVs
    dstrainTh = alpha*dtemp
    
    DdstrainThDtemp = alpha + dalphadtemp*dtemp
    
    DdstrainThDbetaM = alphaM*dtemp
    DdstrainThDbetaB = alphaB*dtemp
    
    ! scalar transformation strain increment (austenite = 0), 
    ! and derivatives w.r.t. temperature and SDVs
    dstrainTr = gammaM*(betaM-betaMn) + gammaB*(betaB-betaBn)
    
    DdstrainTrDtemp = dgammaMdtemp*(betaM-betaMn) +  dgammaBdtemp*(betaB-betaBn)
    
    DdstrainTrdbetaM = gammaM
    DdstrainTrdbetaB = gammaB
    
    ! construct elastic strain at last increment
    strainElN(:) = strain(:)
    strainElN(1:ndi) = strainElN(1:ndi) - strainTh - strainTr
    
    ! construct elastic strain increment depsEl = deps - depsTh - depsTr,
    ! and derivatives w.r.t. temperature and SDVs
    dstrainEl(:) = dstrain(:)
	dstrainEl(1:ndi) = dstrainEl(1:ndi) - dstrainTh - dstrainTr
	
	DdstrainElDtemp(:) = 0.0d0
	DdstrainElDtemp(1:ndi) = - DdstrainThDtemp - DdstrainTrDtemp
	
	DdstrainElDbetaM(:) = 0.0d0
	DdstrainElDbetaM(1:ndi) = - DdstrainThdbetaM - DdstrainTrDbetaM
	DdstrainElDbetaB(:) = 0.0d0
	DdstrainElDbetaB(1:ndi) = - DdstrainThdbetaB - DdstrainTrDbetaB
	
	! calculate stress increment and update total stress
	strainEl(:) = strainElN(:) + dstrainEl(:)
	stress(:) = matmul(C, strainEl)
	
	! calculate partial derivatives for global tangent
	DdstressDkappa(:,:) = 0.0d0
	DdstressDkappa(:,1) = matmul(CM-CA, strainEl) + matmul(C, DdstrainelDbetaM)
	DdstressDkappa(:,2) = matmul(CB-CA, strainEl) + matmul(C, DdstrainelDbetaB)
	
    DdstressDdstrain(:,:) = C(:,:)
    Ddstressdtemp = matmul(dCdtemp, strainEl) + matmul(C, DdstrainEldtemp)
    
end subroutine calculate_stresses


subroutine calculate_heatgen(rpl, drpldstrain, drpldtemp, drpldkappa, &
                             ntens, nkappa, kappa, kappan, tempn, &
                             dtemp, dtime, mdata)
	! Calculate the rate of heat generation per unit mass, and its
	! partial derivatives w.r.t. strains, temperature, and SDVs.
	!
	! Variables passed in
	! ==================================================================
	! ntens         : Number of tensor components
	! nkappa		: Number of components in the SDV vector kappa
	! kappa			: SDVs at the end of the increment
	! kappan		: SDVs at the beginning of the increment
	! tempn			: Temperature at the beginning of the increment
	! dtemp			: Temperature increment
	! dtime         : Time increment
	! mdata         : All material model parameters
	!
	! Variables set by this subroutine
	! ==================================================================
	! rpl           : heat generation rate at the end of the increment
	! drpldstrain   : Derivative of rpl w.r.t. strain increment
	! Drpldtemp 	: Derivative of rpl w.r.t. temperature
	! drpldkappa	: Derivative of rpl w.r.t. SDVs
	use material_data_mod
	
	implicit none
    
    integer, intent(in) :: ntens, nkappa
    double precision, intent(in) :: kappa(nkappa), kappan(nkappa), tempn, dtemp, dtime
    double precision, intent(out) :: rpl, drpldstrain(ntens), drpldtemp, drpldkappa(nkappa)
    type(material_data_type), intent(in) :: mdata
    
    double precision :: rhoDeltahAM, rhoDeltahAB, drhoDeltahAMdtemp, drhoDeltahABdtemp, dbetaMdtime, dbetaBdtime, rM, rB, temp !, rhoDeltahAM0, rhoDeltahAB0, rhoDeltahAB1

    temp = tempn + dtemp
    
    ! approximate rate of volume fractions
    dbetaMdtime = (kappa(1) - kappan(1))/dtime
    dbetaBdtime = (kappa(2) - kappan(2))/dtime
    
    ! parametrized temperature-dependent latent heat
    rhoDeltahAM = mdata%lath_p0rhoDeltahAM + mdata%lath_p1rhoDeltahAM*temp
    rhoDeltahAB = mdata%lath_p0rhoDeltahAB + mdata%lath_p1rhoDeltahAB*temp
    
    drhoDeltahAMdtemp = mdata%lath_p1rhoDeltahAM
    drhoDeltahABdtemp = mdata%lath_p1rhoDeltahAB
    
    ! heat generation rate
    rM = rhoDeltahAM*dbetaMdtime
    rB = rhoDeltahAB*dbetaBdtime
    
    rpl = rM + rB
    
    ! Remark: Partial derivative of heat generation w.r.t. strain is 0 for this model
    drpldstrain = 0.0d0
    
    ! derivative of heat generation rate w.r.t. temperature
    ! this is easy enough to do it analytically
	drpldtemp = drhoDeltahAMdtemp*dbetaMdtime + drhoDeltahABdtemp*dbetaBdtime

	! derivative of heat generation rate w.r.t. sdvs
    ! this is easy enough to do it analytically
	drpldkappa(1) = rhoDeltahAM/dtime
	drpldkappa(2) = rhoDeltahAB/dtime
	drpldkappa(3) = 0.0d0

end subroutine calculate_heatgen


subroutine get_stiffness_tensor(C, dCdE, dCdnu, ntens, ndi, E, nu)
    ! Calculate the stiffness tensor in Voigt notation, and its derivatives
    ! w.r.t. Youngs modulus and Poisson's ratio
    !
    ! Variables passed in
	! ==================================================================
    ! ntens   : number of stress components
    ! ndi     : number of direct stress components
    ! E  	  : Young's modulus
    ! nu      : Poisson's ratio
    !
    ! Variables set by this subroutine
	! ==================================================================
    ! C       : Stiffness tensor
    ! dCdE    : derivative of C w.r.t. Young's modulus
    ! dCdnu   : derivative of C w.r.t. Poisson's Ratio
    implicit none
    integer, intent(in) :: ntens, ndi
    double precision, intent(in) :: E, nu
    double precision, intent(out) :: C(ntens,ntens), dCdE(ntens, ntens), dCdnu(ntens, ntens)
    
    integer i, j
    double precision :: la, mu, dladE, dladnu, dmudE, dmudnu
    
    la = E*nu/(1.0d0+nu)/(1.0d0-2.0d0*nu)
    dladE = nu/(1.0d0+nu)/(1.0d0-2.0d0*nu)
    dladnu = (2.0d0*nu**2+1.0d0)/(2.0d0*nu**2+nu-1.0d0)
    
    mu = E/2.0d0/(1.0d0+nu)
    dmudE = 0.5d0/(1.0d0+nu)
    dmudnu = -0.5d0/(1.0d0+nu)**2
    
    C(:,:) = 0.0d0
    forall (i=1:ndi, j=1:ndi) C(i,j) = la
    C(1:ndi, 1:ndi) = la
    forall (i=1:ndi) C(i,i) = C(i,i) + 2.0d0*mu
    forall (i=ndi+1:ntens) C(i,i) = mu

	dCdE(:,:) = 0.0d0
    forall (i=1:ndi, j=1:ndi) dCdE(i,j) = dladE
    forall (i=1:ndi) dCdE(i,i) = dCdE(i,i) + 2.0d0*dmudE
    forall (i=ndi+1:ntens) dCdE(i,i) = dmudE
    
    dCdnu(:,:) = 0.0d0
    forall (i=1:ndi, j=1:ndi) dCdnu(i,j) = dladnu
    forall (i=1:ndi) dCdnu(i,i) = dCdnu(i,i) + 2.0d0*dmudnu
    forall (i=ndi+1:ntens) dCdnu(i,i) = dmudnu

end subroutine get_stiffness_tensor
