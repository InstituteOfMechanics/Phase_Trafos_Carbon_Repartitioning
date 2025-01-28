!=======================================================================
! This file contains the main UMAT and UMATHT subroutines, and therefore
! defines the interface between Abaqus and the custom coding. 
!=======================================================================

! include separate files
include 'model_parameters.f90'
include 'local_evolution.f90'
include 'stresses_heatgen.f90'

!=======================================================================
!
!                        UMAT MAIN SUBROUTINE
!
!=======================================================================
SUBROUTINE UMAT(stress, statev, ddsdde, sse, spd, scd, rpl, ddsddt, &
                drplde, drpldt, stran, dstran, time, dtime, temp, dtemp, &
                predef, dpred, cmname, ndi, nshr, ntens, nstatv, props, &
                nprops, coords, drot, pnewdt, celent, dfgrd0, dfgrd1, &
                noel, npt, layer, kspt, jstep, kinc)
    ! The interface for the UMAT subroutine is documented in the Abaqus
    ! docs, so we do not repeat it here...
    !
    ! The following order is used for the state variables:
    ! ==================================================================
	! statev(1)  : martensite volume fraction
	! statev(2)  : bainite volume fraction
	! statev(3)  : austenite volume fraction
	! statev(4)  : carbon mass fraction in austenite phase
	! statev(5)  : accumulated linear thermal strain
    ! statev(6)  : accumulated linear transformation strain
	! statev(7)  : derivative of betaM w.r.t. temp (passed to UMATHT)
	! statev(8)  : derivative of betaB w.r.t. temp (passed to UMATHT)
	! statev(9)  : derivative of xCA w.r.t. temp (passed to UMATHT)
    ! statev(10) : partial bainite transformation factor f
	
    use material_data_mod
    
    implicit none
    
    CHARACTER*80 cmname ! Material name

    integer, intent(in) :: ntens, nprops, nstatv, ndi, nshr, jstep(4), &
                           kinc, kspt, noel, npt, layer

    double precision, intent(in) :: stran(ntens), dstran(ntens), &
                                    time(2), predef(1), dpred(1), &
                                    props(nprops), coords(3), drot(3,3), &
                                    dfgrd0(3,3), dfgrd1(3,3), &
                                    celent, dtime, temp, dtemp

    double precision, intent(out) :: ddsdde(ntens,ntens), ddsddt(ntens), &
                                     drplde(ntens), rpl, drpldt

    double precision, intent(inout) :: stress(ntens), sse, spd, scd, &
                                       pnewdt, statev(nstatv)
                                       
    ! nkappa : number of independent state variable in the metallurgical
    ! evolution equations
    integer, parameter :: nkappa = 3
    
    integer :: i, material_flag, martensite_flag, bainite_flag
    double precision :: dkappadtemp(nkappa), kappa(nkappa), kappan(nkappa)
	double precision :: dstressdkappa(ntens, nkappa), drpldkappa(nkappa)
    double precision :: strainTh, dstrainTh, strainTr, dstrainTr, density
    double precision :: fbainite
    
    type(material_data_type) :: mdata
    
    ! extract properties and combine them with parameters set in Fortran code
    
    ! mass density (specification in Abaqus needed, avoid double specification)
    density = props(1)
    
    ! material flag (0: 100Cr6,  1: 100CrMnSi6-4) from property array 
    material_flag = nint(props(2))
    
    ! switches for disabling martensite/bainite evolution
    martensite_flag = nint(props(3))
    bainite_flag = nint(props(4))
    
    ! get a struct with all model parameters to easily pass them to
    ! all subroutines
    call get_material_data(mdata, material_flag, density, &
                           martensite_flag, &
                           bainite_flag)
    
    ! Unique state variables (only those actually needed at the quadrature point level) :
	! kappa(1) : martensite volume fraction
	! kappa(2) : bainite volume fraction
	! kappa(3) : carbon mass fraction in austenite phase
	kappan(1) = statev(1)
	kappan(2) = statev(2)
	kappan(3) = statev(4)
	
    !===================================================================
    ! Step 1: Solve evolution equations for kappa and derivatives
    !         Derivative w.r.t. strains is 0 for this model.
    !===================================================================
    call local_evolution(kappa, dkappadtemp, pnewdt, fbainite, nkappa, &
                         kappan, dtime, temp, dtemp, mdata)

    ! if a smaller step is required, dont waste time
    if (pnewdt < 1.0d0) then
        return
    end if
    
    ! update the global state variable vector
    statev(1) = kappa(1)
    statev(2) = kappa(2)
    statev(3) = 1.0d0 - statev(1) - statev(2)
    statev(4) = kappa(3)
    
    statev(7:9) = dkappadtemp
    statev(10) = fbainite
    
    !===================================================================
    ! Step 2: Constitutive model and derivatives
    !===================================================================
    strainTh = statev(5)
    strainTr = statev(6)
    
    call calculate_stresses(stress, ddsdde, ddsddt, dstressdkappa, &
                            dstrainTh, dstrainTr, ntens, ndi, nkappa, &
                            stran, strainTh, strainTr, dstran, kappa, &
                            kappan, temp, dtemp, mdata)
                            
    statev(5) = statev(5) + dstrainTh
    statev(6) = statev(6) + dstrainTr
    
    call calculate_heatgen(rpl, drplde, drpldt, drpldkappa, ntens, &
                           nkappa, kappa, kappan, temp, dtemp, dtime, &
                           mdata)
    
    !===================================================================
    ! Step 4: Sensitivities
    !===================================================================
    ddsddt = ddsddt + matmul(dstressdkappa, dkappadtemp)
    drpldt = drpldt + dot_product(drpldkappa, dkappadtemp)
    
END SUBROUTINE UMAT


!=======================================================================
!
!                       UMATHT MAIN SUBROUTINE
!
!=======================================================================
SUBROUTINE UMATHT(u, dudt, dudg, flux, dfdt, dfdg, &
                  statev, tempn, dtemp, dtemdx, time, dtime, predef, &
                  dpred, cmname, ntgrd, nstatv, props, nprops, coords, &
                  pnewdt, noel, npt, layer, kspt, kstep, kinc)
    ! The interface for the UMATHT subroutine is documented in the Abaqus
    ! docs, so we do not repeat it here...
    !
    ! Variables that will be updated by this subroutine:
    ! ==================================================================
    ! U    : internal thermal energy per unit mass
    ! dudt : variation of U w.r.t. temperature
    ! dudg : variation of U w.r.t. temperature gradient (usually 0)
    ! flux : heat flux vector
    ! dfdt : variation of heat flux vector w.r.t. temperature
    ! dfdg : variation of heat flux vector w.r.t. temperature gradient
    use material_data_mod
    
    implicit none
    
    CHARACTER*80 cmname ! Material name

    integer, intent(in) :: ntgrd, nstatv, nprops, noel, npt, layer, &
                           kspt, kstep, kinc

    double precision, intent(inout) :: u, statev(nstatv), flux(ntgrd), pnewdt

    double precision, intent(out) :: dudt, dudg(ntgrd), &
                                     dfdt(ntgrd), dfdg(ntgrd,ntgrd)

    double precision, intent(in) :: tempn, dtemp, dtemdx(ntgrd), dtime, &
                                    time(2), predef(1), dpred(1), &
                                    props(nprops), coords(3)

    double precision :: temp, betaM, betaB, betaA
    double precision :: cpM, cpB, cpA, cp, lambM, lambB, lambA, lamb
    double precision :: dlambMdtemp, dlambBdtemp, dlambAdtemp, dlambdtemp
    double precision :: dbetaMdtemp, dbetaBdtemp, dbetaAdtemp
    
    type(material_data_type) :: mdata

    integer i
    
    ! get a struct with all model parameters
    call get_material_data(mdata, props(1), props(2))
    
    ! get the current state
    temp = tempn + dtemp
    
    betaM = statev(1)
    betaB = statev(2)
    betaA = 1.0d0 - betaM - betaB
    
    ! derivative of volume fractions w.r.t. temperature was computed in UMAT
    dbetaMdtemp = statev(7)
    dbetaBdtemp = statev(8)
    dbetaAdtemp = -dbetaBdtemp-dbetaMdtemp
    
    ! evaluate parametrization of temperature-dependent material parameters
    cpM = mdata%thermal_p0cpM + mdata%thermal_p1cpM*temp + mdata%thermal_p2cpM*temp**2
    cpB = mdata%thermal_p0cpB + mdata%thermal_p1cpB*temp + mdata%thermal_p2cpB*temp**2
    cpA = mdata%thermal_p0cpA + mdata%thermal_p1cpA*temp + mdata%thermal_p2cpA*temp**2
    
    lambM = mdata%thermal_p0lambM + mdata%thermal_p1lambM*temp + mdata%thermal_p2lambM*temp**2
    lambB = mdata%thermal_p0lambB + mdata%thermal_p1lambB*temp + mdata%thermal_p2lambB*temp**2
    lambA = mdata%thermal_p0lambA + mdata%thermal_p1lambA*temp + mdata%thermal_p2lambA*temp**2
    
    dlambMdtemp = mdata%thermal_p1lambM + 2.0d0*mdata%thermal_p2lambM
    dlambBdtemp = mdata%thermal_p1lambB + 2.0d0*mdata%thermal_p2lambB
    dlambAdtemp = mdata%thermal_p1lambA + 2.0d0*mdata%thermal_p2lambA
    
    cp = betaM*cpM + betaB*cpB + betaA*cpA
    lamb = betaM*lambM + betaB*lambB + betaA*lambA

    ! consider also the dependency of SDVS in temperature here
    dlambdtemp =   betaM*dlambMdtemp + betaB*dlambBdtemp + betaA*dlambAdtemp &
                 + dbetaMdtemp*lambM + dbetaBdtemp*lambB + dbetaAdtemp*lambA
    
    !===================================================================
    ! Update internal energy and heat flux
    ! The dependency on the internal variables was already included in
    ! the derivatives.
    !===================================================================
    u = u + cp*dtemp
    
    dudt = cp
    
    dudg(:) = 0.0d0
    

    flux(:) = -lamb*dtemdx(:)
    
    dfdt(:) = -dlambdtemp*dtemdx(:)
    
    dfdg = 0.0d0
    forall (i=1:ntgrd) dfdg(i,i) = -lamb
    
END SUBROUTINE UMATHT



