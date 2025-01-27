! TODO: LICENCE ETC

!=======================================================================
!
! This file contains all subroutines used for the solution of the
! local system of (nonlinear) evolution equations. The actual solving
! is implemented using Intel MKL, and especially the reverse 
! communication interface (RCI), where the solver routine is called in
! a loop and requests updates of the residual and Jacobian by a flag.
!
!=======================================================================

include 'mkl_rci.f90' ! needed for reverse communication in mkl solver
include 'mkl_lapack.f90' ! needed for linear equation solver

subroutine solve_linear_system(x, info, n, A, b)
	! Solve the linear system of equations A*x=b
	!
	! Input variables
	! ==================================================================
	! A: coefficient matrix
	! b: right hand side of the system
	! n: number of equations
	!
	! Output variables
	! ==================================================================
	! x: solution of the system
	
	implicit none
	
	double precision, intent(in) :: A(n,n), b(n)
	integer, intent(in) :: n
	double precision, intent(out) :: x(n)
	integer, intent(out) :: info

	integer ipiv(n)

	! copy b to x since ?getrs replaces the rhs array with the solution
	x = b

	call dgetrf(n,n,A,n,ipiv,info)

	if (info /= 0) then
		write (*,*) "Warning: Singular equation system encountered."
		x = 0.0d0
		return
	end if

	call dgetrs('N',n,1,A,n,ipiv,x,n,info)

end subroutine solve_linear_system

subroutine heaviside(x, x0, w, y, dydx, dydx0)
	! Regularized Heaviside function
	!
	! The Heaviside function is approximated as a tanh within a width of
	! 2*w around the jump. Actual values of 0 or 1 are returned outside
	! of this window.
	! 
	! input variables
	! ==================================================================
	! x     : Current coordinate
	! x0    : Coordinate where the Heaviside function jumps from 0 to 1
	! w     : Regularization width, h(x0-w/2) = 0.01, h(x0+w/2) = 0.99
	! 
	! output variables
	! ==================================================================
	! y     : Value of the regularized step function
	! dydx  : Derivative w.r.t. x
	! dydx0 : Derivative w.r.t. x0
	implicit none

	double precision, intent(in) :: x, x0, w
	double precision, intent(out) :: y, dydx, dydx0
	double precision :: a

	! this a guarantess that h(x0-w/2) = 0.01, h(x0+w/2) = 0.99
	a = 2.29756d0*2.0d0/w
	
	! if switch to prevent possible overflows of cosh()
	if (abs(x-x0) < w) then
		! use tanh near x0 where it is well defined
		y = 0.5d0*(tanh(a*(x - x0)) + 1.0d0)
		dydx = 0.5d0*a/cosh(a*(x - x0))**2
		dydx0 = -0.5d0*a/cosh(a*(x - x0))**2
	else
		! otherwise just use a hard switch
		if (x < x0) then
			y = 0.0d0
		else
			y = 1.0d0
		end if
		dydx = 0.0d0
		dydx0 = 0.0d0
	end if
	
end subroutine heaviside


subroutine max_bainite_fraction(betaBhat, dbetaBhatdT, dbetaBhatdbetaM, &
								dbetaBhatdbetaB, dbetaBhatdxCA, f, T, &
								betaM, betaB, xCA, mdata)
	! Calculate the maximum reachable bainite fraction
	!
	! Input variables
	! ==================================================================
	! T                 : Current temperature
	! betaM             : Martensite volume fraction
	! betaB             : Bainite volume fraction
	! xCA               : Austenite phase carbon fraction
	! mdata             : All material model parameters
	! 
	! Output variables
	! ==================================================================
	! betaBhat          : Maximum reachable bainite fraction
	! DbetaBhatDT       : Derivative w.r.t. temperature
	! DbetaBhatDbetaM   : Derivative w.r.t. martensite fraction
	! DbetaBhatDbetaB   : Derivative w.r.t. bainite fraction
	! DbetaBhatDxCA     : Derivative w.r.t. austenite carbon fraction
	! f                 : Factor for incomplete trafo 
	use material_data_mod
	
	implicit none
	
	type(material_data_type), intent(in) :: mdata
	
	double precision, intent(in) :: T, betaM, betaB, xCA
	double precision, intent(out) :: betaBhat, dbetaBhatdT, dbetaBhatdbetaM, &
	                                 dbetaBhatdbetaB, dbetaBhatdxCA, f
	                                 
	double precision, parameter :: p0Bs=681.51d0, &
								   p1Bs=-29618.0d0, &
								   xCB=3.0d-4
								   
	double precision :: betaA, xCAinf, dfdxCA, dfdxCAinf, dfdT, dxCAinfdT
	
	! calculate the carbon fraction for which T = Bs
	! T = Bs = p0 + p1*xCAinf
	betaA = 1.0d0 - betaM - betaB
	
	xCAinf = (T-mdata%jmak_p0Bs)/mdata%jmak_p1Bs
	
	f = (xCA-xCAinf)/(mdata%xCB-xCAinf)
	
	if (f < 0.0d0 ) then
		f = 0.0d0
	end if
	
	dfdxCA = -xCAinf/(mdata%xCB-xCAinf)
	dfdxCAinf = -1.0d0/(mdata%xCB-xCAinf) + (xCA-xCAinf)/(mdata%xCB-xCAinf)**2
	dxCAinfdT = -mdata%jmak_p0Bs/mdata%jmak_p1Bs
	dfdT = dfdxCAinf*dxCAinfdT
        
	betaBhat = f*betaA + betaB
	
	dbetaBhatdT = dfdT*betaA
	dbetaBhatdbetaM = -f
	dbetaBhatdbetaB = 1.0d0 - f
	dbetaBhatdxCA = betaA*dfdxCA

end subroutine max_bainite_fraction


subroutine martensite_evolution(RbetaM, DRbetaMDkappa, DRbetaMDtemp, &
                                request, nkappa, kappa, kappan, dtime, &
                                tempn, dtemp, mdata)
	! Calulate the residual for martensite evolution and its derivatives
	!
	! Input variables
	! ==================================================================
	! request         : 1: residual, 2: deriv. sdvs, 3: deriv. temp
	! nkappa          : Number of unique SDVs
	! kappa           : Unique SDVs
	! kappan          : Unique SDVS at last converged increment
	! dtime           : Time increment
	! tempn           : Temperature at last converged increment
	! dtemp           : Temperature increment
	! mdata           : All material model parameters
	!
	! Ouput variables
	! ==================================================================
	! RbetaM          : Residual for martensite evolution
	! DRbetaMDkappa   : Derivative w.r.t. SDVs
	! DRbetaMDtemp    : Derivative w.r.t. temperature
	 
	use material_data_mod

	implicit none
	
	integer, intent(in) :: nkappa, request
    double precision, intent(in) :: kappa(nkappa), kappan(nkappa)
    double precision, intent(in) :: dtime, tempn, dtemp 
	type(material_data_type), intent(in) :: mdata
	double precision, intent(out) :: RbetaM, DRbetaMDkappa(nkappa), DRbetaMDtemp

    double precision :: betaA, betaM, betaMn, betaB, xCA, temp
    double precision :: Ms, k, zetaM
    double precision :: dMsdxCA, dkdxCA, dzetaMdMs, dzetaMdxCA, dzetaMdT
    
    temp = tempn + dtemp

    betaM = kappa(1)
    betaB = kappa(2)
    xCA = kappa(3)

    betaMn = kappan(1)
    
    betaA = 1.0d0 - betaM - betaB
    
    ! Koistinen-Marburger parameters, and derivatives
    Ms = mdata%km_p0Ms + mdata%km_p1Ms*xCA
    dMsdxCA = mdata%km_p1Ms
        
    k = mdata%km_p0k + mdata%km_e1k*exp(mdata%km_e2k*xCA)
    dkdxCA = mdata%km_e1k*mdata%km_e2k*exp(mdata%km_e2k*xCA)
    
    ! KM activation function
    call heaviside(Ms, temp, mdata%km_regularization_width, zetaM, dzetaMdMs, dzetaMdT)
    dzetaMdxCA = dzetaMdMs*dMsdxCA
    
	select case (request)

        case (1)
			! do nothing if martensite evolution is disabled
			if (mdata%martensite_flag == 1) then
				RbetaM = betaM - betaMn
				return
			end if
        
			! calculate the residual
			if (dtemp <= 0.0d0) then
				RbetaM = betaM - betaMn + zetaM*betaA*k*dtemp
			else
				RbetaM = betaM - betaMn
			end if
			
			return
            
        case (2)
			! do nothing if martensite evolution is disabled
			if (mdata%martensite_flag == 1) then
				DRbetaMDkappa(1) = 1.0d0
				DRbetaMDkappa(2) = 0.0d0
				DRbetaMDkappa(3) = 0.0d0
				return
			end if
			
            ! calculate the derivative of the residual w.r.t. kappa
            if (dtemp <= 0.0d0) then
				DRbetaMDkappa(1) = 1.0d0 - zetaM*k*dtemp
				DRbetaMDkappa(2) = -zetaM*k*dtemp
				DRbetaMDkappa(3) =  betaA*dtemp*(zetaM*dkdxCA + dzetaMdxCA*k)
			else
				DRbetaMDkappa(1) = 1.0d0
				DRbetaMDkappa(2) = 0.0d0
				DRbetaMDkappa(3) = 0.0d0
			end if
			
            return
            
        case (3)
			! do nothing if martensite evolution is disabled
			if (mdata%martensite_flag == 1) then
				DRbetaMDtemp = 0.0d0
				return
			end if
			
			! calculate derivative of residual w.r.t. temperature
			if (dtemp < 0.0d0) then
				DRbetaMDtemp = betaA*k*(dzetaMdT*dtemp + zetaM)
			else
				DRbetaMDtemp = 0.0d0
			end if
			
        case default
            write (*,*) 'Wrong request in martensite evolution.'
            stop 1

    end select

end subroutine martensite_evolution



subroutine bainite_evolution(RbetaB, DRbetaBDkappa, DRbetaBDtemp, f, &
							 request, nkappa, kappa, kappan, dtime, &
							 tempn, dtemp, mdata)
	! Calulate the residual for bainite evolution and its derivatives
	!
	! Input variables
	! ==================================================================
	! request         : 1: residual, 2: deriv. sdvs, 3: deriv. temp
	! nkappa          : Number of unique SDVs
	! kappa           : Unique SDVs
	! kappan          : Unique SDVS at last converged increment
	! dtime           : Time increment
	! tempn           : Temperature at last converged increment
	! dtemp           : Temperature increment
	! mdata           : All material model parameters
	!
	! Ouput variables
	! ==================================================================
	! RbetaB          : Residual for bainite evolution
	! DRbetaBDkappa   : Derivative w.r.t. SDVs
	! DRbetaBDtemp    : Derivative w.r.t. temperature
	! f               : Factor for incomplete bainite transformation
	use material_data_mod
							 
	implicit none

    integer, intent(in) :: nkappa, request
    double precision, intent(in) :: kappa(nkappa), kappan(nkappa)
    double precision, intent(in) :: dtime, tempn, dtemp
    type(material_data_type), intent(in) :: mdata
    double precision, intent(out) :: RbetaB, DRbetaBDkappa(nkappa), &
                                     DRbetaBDtemp, f
	
	
	double precision :: betaA, betaM, betaB, betaBn, xCA, temp
    double precision :: Bs, N, b, zetaB, dNdT, dbdT
    double precision :: dBsdxCA, dzetaBdBs, dzetaBdxCA, dzetaBdT
    double precision :: betaBhat, dbetaBhatdT, dbetaBhatdbetaM, &
                        dbetaBhatdbetaB, dbetaBhatdxCA
    double precision :: tau, dtaudb, dtaudN, dtaudT
    double precision :: AA
    
    double precision, parameter :: tol_transformation_stop = 1.0d-4
	
    temp = tempn + dtemp

    betaM = kappa(1)
    betaB = kappa(2)
    xCA = kappa(3)

    betaBn = kappan(2)

    betaA = 1.0d0 - betaM - betaB

    ! JMAK parameters (no derivatives w.r.t. sdvs)
    Bs = mdata%jmak_p0Bs + mdata%jmak_p1Bs*xCA
    dBsdxCA = mdata%jmak_p1Bs

    b = exp(mdata%jmak_p0logb + mdata%jmak_p1logb*temp + mdata%jmak_p2logb*temp**2)
    dbdT = (mdata%jmak_p1logb + 2.0d0*mdata%jmak_p2logb*temp)*b
    
    N = mdata%jmak_p0N + mdata%jmak_p1N*temp + mdata%jmak_p2N*temp**2
    dNdT = mdata%jmak_p1N + 2.0d0*mdata%jmak_p2N*temp
    
    call max_bainite_fraction(betaBhat, dbetaBhatdT, dbetaBhatdbetaM, &
							  dbetaBhatdbetaB, dbetaBhatdxCA, f, &
							  temp, betaM, betaB, xCA, mdata)

    ! JMAK activation function
    call heaviside(Bs, temp, mdata%jmak_regularization_width, zetaB, &
                   dzetaBdBs, dzetaBdT)
	dzetaBdxCA = dzetaBdBs*dBsdxCA
	
	! JMAK nucleation time
	tau = (1/b*log(1.0d0/0.99d0))**(1.0d0/N)
	
    select case (request)

        case (1)
			! calculate the residual
			
			! do nothing if bainite evolution is disabled
			if (mdata%bainite_flag == 1) then
				RbetaB = betaB - betaBn
				return
			end if

			if (betaBn < mdata%jmak_nucleation_fraction) then
				! nucleation phase
				RbetaB = betaB - betaBn - mdata%jmak_nucleation_fraction*zetaB*dtime/tau
			else
				! JMAK law after nucleation is complete
				if (betaBhat - betaB > tol_transformation_stop) then
					! transformation still ongoing
					AA = log(betaBhat) - log(betaBhat - betaB)
					RbetaB = betaB - betaBn - dtime*zetab*N*b**(1.0d0/N)*(betaBhat - betaB)*AA**(1.0d0-1.0d0/N)
				else
					! close to 0 transformable austenite left, no more evolution
					RbetaB = betaB - betaBn
				end if
				 
			end if
		
			return
            
        case (2)
            ! calculate the derivative of the residual w.r.t. kappa
			
			! do nothing if bainite evolution is disabled
			if (mdata%bainite_flag == 1) then
				DRbetaBDkappa(1) = 0.0d0
				DRbetaBDkappa(2) = 1.0d0
				DRbetaBDkappa(3) = 0.0d0
				return
			end if
			
			if (betaBn < mdata%jmak_nucleation_fraction) then
				! nucleation phase
				DRbetaBDkappa(1) = 0.0d0
				DRbetaBDkappa(2) = 1.0d0
				DRbetaBDkappa(3) = -mdata%jmak_nucleation_fraction*dtime/tau*dzetaBdxCA
			else
				! JMAK law after nucleation is complete
				if (betaBhat - betaB > tol_transformation_stop) then
					! transformation still ongoing
					AA = log(betaBhat) - log(betaBhat - betaB)
					DRbetaBDkappa(1) = -dtime*zetaB*(b/AA)**(1.0d0/N)*(N*AA - (N-1)*betaB/betaBhat)*dbetaBhatdbetaM
					DRbetaBDkappa(2) = 1.0d0 - dtime*zetaB*(b/AA)**(1.0d0/N)* &
							   (N*AA*(dbetaBhatdbetaB - 1.0d0) + (N - 1.0d0)*(1.0d0 - betaB/betaBhat*dbetaBhatdbetaB))
					DRbetaBDkappa(3) = -dtime*zetaB*(b/AA)**(1.0d0/N)*(N*AA - (N - 1.0d0)*betaB/betaBhat)*dbetaBhatdxCA & 
					           -dtime*N*b**(1.0d0/N)*(betaBhat-betaB)*AA**(1.0d0-1.0d0/N)*dzetaBdxCA
				else
					! close to 0 transformable austenite left, no more evolution
					DRbetaBDkappa(1) = 0.0d0
					DRbetaBDkappa(2) = 1.0d0
					DRbetaBDkappa(3) = 0.0d0
				end if
			end if
			
            return
            
        case (3)
			! calculate derivative of residual w.r.t. temperature
			
			! do nothing if bainite evolution is disabled
			if (mdata%bainite_flag == 1) then
				DRbetaBDTemp = 0.0d0
				return
			end if
			
			if (betaBn < mdata%jmak_nucleation_fraction) then
				! nucleation phase
				dtaudb = -tau/b/N
				dtaudN = -tau/N**2*log(1.0d0/b*log(1.0d0/0.99d0))
				dtaudT = dtaudb*dbdT + dtaudN*dNdT
				DRbetaBDtemp = mdata%jmak_nucleation_fraction*dtime*(zetaB/tau**2*dtaudT + dzetaBdT/tau)
			else
				! JMAK law after nucleation is complete
				if (betaBhat - betaB > tol_transformation_stop) then
					! transformation still ongoing
					AA = log(betaBhat) - log(betaBhat - betaB)
					DRbetaBDtemp = -dtime*AA**(-1.0d0/N)*(  &
						  N*AA*b**(1.0d0/N)*(betaBhat-betaB)*dzetaBdT &
					    + zetaB*AA*b**(1.0d0/N)*(betaBhat-betaB)*(1.0d0-log(b)/N+log(AA)/N)*dNdT &
					    + zetaB*AA*(betaBhat-betaB)*b**(1.0d0/N-1.0d0)*dbdT &
					    + zetaB*b**(1.0d0/N)*(AA*N-(N-1.0d0)*betaB/betaBhat)*dbetaBhatdT )
				else
					! close to 0 transformable austenite left, no more evolution
					DRbetaBDtemp = 0.0d0
				end if
			end if
		
        case default
            ! request should only ever be 1 or 2
            write (*,*) 'Wrong request in bainite evolution.'
            stop 1

    end select

end subroutine bainite_evolution


subroutine carbon_evolution(RxCA, DRxCADkappa, DRxCADtemp, &
							request, nkappa, kappa, kappan, dtime, &
							tempn, dtemp, mdata)
	! Calulate the residual for evolution of austenite phase carbon
	! fraction and its derivatives
	!
	! Input variables
	! ==================================================================
	! request         : 1: residual, 2: deriv. sdvs, 3: deriv. temp
	! nkappa          : Number of unique SDVs
	! kappa           : Unique SDVs
	! kappan          : Unique SDVS at last converged increment
	! dtime           : Time increment
	! tempn           : Temperature at last converged increment
	! dtemp           : Temperature increment
	! mdata           : All material model parameters
	!
	! Ouput variables
	! ==================================================================
	! RxCA            : Residual for carbon fraction evolution
	! DRxCAMDkappa    : Derivative w.r.t. SDVs
	! DRxCADtemp      : Derivative w.r.t. temperature

	use material_data_mod
	
	implicit none
	
	integer, intent(in) :: nkappa, request
    double precision, intent(in) :: kappa(nkappa), kappan(nkappa)
    double precision, intent(in) :: dtime, tempn, dtemp
    type(material_data_type), intent(in) :: mdata
    double precision, intent(out) :: RxCA, DRxCADkappa(nkappa), DRxCADtemp
	
	
	double precision :: betaA, betaM, betaB, betaBn, xCA, xCAn
	
	double precision, parameter :: tol_transformation_stop = 1.0d-4
	
    betaM = kappa(1)
    betaB = kappa(2)
    xCA = kappa(3)

    betaBn = kappan(2)
    xCAn = kappan(3)

    betaA = 1.0d0 - betaM - betaB
	
	select case (request)

        case (1)
			! calculate the residual
			
			if (betaA > tol_transformation_stop) then
			    ! still some austenite left
				RxCA = betaA*(xCA - xCAn) - (betaB - betaBn)*(xCA - mdata%xCB)
			else
				! if no austenite is left, calculating an increase in 
				! carbon content is not well defined
				RxCA = xCA - xCAn
			end if
        
			return
            
        case (2)
            ! calculate derivative of residual w.r.t. kappa
			if (betaA > tol_transformation_stop) then
				DRxCADkappa(1) = -(xCA - xCAn)
				DRxCADkappa(2) = -(xCA - xCAn) - (xCA - mdata%xCB)
				DRxCADkappa(3) = betaA - (betaB - betaBn)
			else
				DRxCADkappa(1) = 0.0d0
				DRxCADkappa(2) = -(xCA - mdata%xCB)
				DRxCADkappa(3) = -(betaB - betaBn)
			end if
			
            return
            
        case (3)
			! calculate derivative of residual w.r.t. temperature			
			DRxCADtemp = 0.0d0
		
        case default
            ! request should only ever be 1 or 2
            write (*,*) 'Wrong request in austenite phase carbon evolution.'
            stop 1

    end select

end subroutine carbon_evolution


subroutine local_residual(res, jac, DresDtemp, fbainite, request, nkappa, & 
						  kappa, kappan, dtime, tempn, dtemp, mdata)
	! Update local residual and its derivatives.
	!
    ! Depending on the request flag, do one of the following:
	! request == 1:
	!     calculate the residual, i.e.
	!     R = [RbetaM, RbetaB, RxCA]
	!
	! request == 2:
	! 	  calculate the Jacobian, i.e.
	!         | dR1/dx1   dR1/dx2   dR1/dx3 |
	!     J = | dR2/dx1   dR2/dx2   dR2/dx3 |
	!         | dR3/dx1   dR3/dx2   dR3/dx3 |
	!
	! request == 3:
	!     calculate the temperature derivative, i.e.
	!     dR/dtemp = [dRbetaMdtemp, dRbetaBdtemp, dRxCAdtemp]
	!
	! Input variables:
	! ==================================================================
	! request         : 1: residual, 2: deriv. sdvs, 3: deriv. temp
	! nkappa          : Number of unique SDVs
	! kappa           : Unique SDVs
	! kappan          : Unique SDVS at last converged increment
	! dtime           : Time increment
	! tempn           : Temperature at last converged increment
	! dtemp           : Temperature increment
	! mdata           : All material model parameters
	! 
	! Output variables:
	! ==================================================================
	! res             : Residual of evolution equations
	! jac             : Derivative w.r.t. kappa
	! DresDtemp       : Derivative w.r.t. temperature
	! fbainite        : Factor for incomplete bainite trafo
	
	use material_data_mod
	
    implicit none

    integer, intent(in) :: nkappa, request
    double precision, intent(in) :: kappa(nkappa), kappan(nkappa)
    double precision, intent(in) :: dtime, tempn, dtemp
    type(material_data_type), intent(in) :: mdata
    double precision, intent(out) :: res(nkappa), jac(nkappa, nkappa), &
                                     dresdtemp(nkappa), fbainite
	
	call martensite_evolution(res(1), jac(1,:), DresDtemp(1), &
                              request, nkappa, kappa, kappan, dtime, &
                              tempn, dtemp, mdata)
	
	call bainite_evolution(res(2), jac(2,:), DresDtemp(2), fbainite, &
						   request, nkappa, kappa, kappan, dtime, &
						   tempn, dtemp, mdata)
	
    call carbon_evolution(res(3), jac(3,:), DresDtemp(3), &
						  request, nkappa, kappa, kappan, dtime, &
						  tempn, dtemp, mdata)

end subroutine local_residual


subroutine local_evolution(kappa, DkappaDtemp, pnewdt, fbainite, &
						   nkappa, kappan, dtime, tempn, dtemp, mdata)
	! Solve the local system of evolution equations, and calculate the
	! derivative of the SDVs w.r.t. temperature
	!
	! Input variables
	! ==================================================================
	! nkappa       : Number of unique SDVs
	! kappan       : Unique SDVS at last converged increment
	! dtime        : Time increment
	! tempn        : Temperature at last converged increment
	! dtemp        : Temperature increment
	! mdata        : All material model parameters
	!
	! Output variables
	! ==================================================================
	! kappa        : Solution for unique SDVs for the current increment
	! DkappaDtemp  : Derivative w.r.t. temperature
	! pnewdt       : Factor for automatic time incrementation
	! fbainite     : Factor for incomplete bainite trafo
	
    ! MKL Reverse Communication Interface
    use MKL_RCI
    use MKL_RCI_TYPE
    use material_data_mod

    implicit none

    ! Supplementary routine to track/free memory
    external MKL_FREE_BUFFERS

    ! objective function
    ! the implementation of DJACOBIX is complicated with type checking,
    ! therefore we declare local_residual as external
    external local_residual

    integer, intent(in) :: nkappa
    double precision, intent(in) :: dtime, tempn, dtemp, kappan(nkappa)
    type(material_data_type), intent(in) :: mdata
    
    double precision, intent(out) :: kappa(nkappa), dkappadtemp(nkappa), &
                                     fbainite
                                     
    double precision, intent(inout) :: pnewdt
    
    ! function value vector, derivative w.r.t. SDVs
    double precision res(nkappa), jac(nkappa, nkappa)

    ! increments of volume fractions
    double precision dbetaM, dbetaB

    ! threshold values for volume increments to accept step
    double precision, parameter :: dbetaMmax=5.0d-1, dbetaBmax= 5.0d-1

    ! stopping criteria:
    ! 1: trust region area < eps(1)
    ! 2: function norm < eps(2)
    ! 3: jacobian norm < eps(3)
    ! 4: trial step norm < eps(4)
    ! 5: function change norm < eps(5)
    ! 6: trial step precision < eps(6)
    double precision, parameter, dimension(6) :: eps = &
					[ 1.0d-8, 1.0d-6, 1.0d-8, 1.0d-8, 1.0d-6, 1.0d-8]

	! iter1: maximum number of iterations
    ! iter2: maximum number of trial step iterations
	integer, parameter :: iter1=1000, iter2=100

    ! initial step bound / radius
    double precision, parameter :: rs = 1.0d0
    
    ! identifier for RCI
    integer rci_request

    ! iter: iteration counter
    ! st_cr: number of stopping criterion
    ! successful: whether a solution was found
    integer iter, st_cr, successful
	
    ! initial / final residual
    double precision r1, r2

    ! probably an mkl type for the handle...
    TYPE(handle_tr) :: handle

    ! counters
    integer i, j

    ! results of input parameter checks
    integer INFO(6)
    
	! for central difference approximation of partial derivatives
    double precision :: minusdresdtemp(nkappa)
	
	! initialize solution vector for updated sdvs and fixed parameters
    kappa(:) = kappan(:)

    ! initialize fvec and jjac
    res(:) = 0.0d0
    jac(:,:) = 0.0d0

    !** INITIALIZE SOLVER (ALLOCATE MEMORY, SET INITIAL VALUES)
    !**     HANDLE    IN/OUT: TR SOLVER HANDLE
    !**     N         IN:     NUMBER OF FUNCTION VARIABLES
    !**     M         IN:     DIMENSION OF FUNCTION VALUE
    !**     X         IN:     SOLUTION VECTOR. CONTAINS VALUES X FOR F(X)
    !**     EPS       IN:     PRECISIONS FOR STOP-CRITERIA
    !**     ITER1     IN:     MAXIMUM NUMBER OF ITERATIONS
    !**     ITER2     IN:     MAXIMUM NUMBER OF ITERATIONS OF CALCULATION OF TRIAL-STEP
    !**     RS        IN:     INITIAL STEP BOUND
    if (DTRNLSP_INIT(handle, nkappa, nkappa, kappa, eps, iter1, iter2, rs) /= TR_SUCCESS) THEN
        call MKL_FREE_BUFFERS

        write (*,*) '| ERROR IN DTRNLSP_INIT'
        stop 1
    end if


    !** CHECKS THE CORRECTNESS OF HANDLE AND ARRAYS CONTAINING JACOBIAN MATRIX,
    !** OBJECTIVE FUNCTION, LOWER AND UPPER BOUNDS, AND STOPPING CRITERIA.
    if (DTRNLSP_CHECK(handle, nkappa, nkappa, jac, res, eps, info) /= TR_SUCCESS) THEN
        call MKL_FREE_BUFFERS

        write (*,*) '| ERROR IN DTRNLSPBC_INIT'
        stop 1
    else
        ! info = validity for [handle, fjac, fvec, eps, ?, ?], 0 means valid
        if ( info(1) /= 0 .or. info(2) /= 0 .or. info(3) /= 0 .or. info(4) /= 0 ) THEN
            call MKL_FREE_BUFFERS

            write (*,*) '| INPUT PARAMETERS ARE NOT VALID'
            stop 1
        end if
    end if


    ! RCI cycle variables
    rci_request = 0
    successful = 0

    ! RCI cycle
    do while (successful == 0)
        !** CALL TR SOLVER
        !**   HANDLE        IN/OUT: TR SOLVER HANDLE
        !**   FVEC          IN:     VECTOR
        !**   FJAC          IN:     JACOBI MATRIX
        !**   RCI_REQUEST   IN/OUT: RETURN NUMBER WHICH DENOTE NEXT STEP FOR PERFORMING
        if (DTRNLSP_SOLVE (handle, res, jac, rci_request) /= TR_SUCCESS) THEN
            call MKL_FREE_BUFFERS

            write (*,*) '| ERROR IN DTRNLSP_SOLVE'
            stop 1
        end if

        ! next step depends on the rci_request:
        ! rci_request = -6, ... -1 : Stopping criterion reached
        ! rci_request = 1 : fvec must be recalculated
        ! rci_request = 2 : fjac must be updated
        ! rci_request = 0 : successful iteration on radius, but x may or may not have been changed
        select case (rci_request)

            case (-1, -2, -3, -4, -5, -6)
                successful = 1
                
            case (1, 2)
                ! update residual vector or jacobian
                call local_residual(res, jac, minusdresdtemp, fbainite, &
                                    rci_request, nkappa, kappa, kappan, &
                                    dtime, tempn, dtemp, mdata)
                
            case (0)
				! do nothing, just go to the next iteration
                
            case default
				write (*,*) 'Wrong request during solution of evolution equations.'
            stop 1
				
        end select
    end do
    
    ! handle: tr solver handle
    ! iter: number of iterations
    ! st_cr: stopping criterion triggered
    !        1:  number of iterations reached
    !        2:  trust region area < eps(1)
    !        3:  function norm < eps(2)
    !        4:  Jacobian is singular (some norm < eps(3))
    !        5:  trial step norm < eps(4)
    !        6:  change in function norm < eps(5)
    ! r1: initial residual
    ! r2: final residual
    if (DTRNLSP_GET(handle, iter, st_cr, r1, r2) /= TR_SUCCESS) THEN
        call MKL_FREE_BUFFERS
        write (*,*) '| ERROR IN DTRNLSP_GET'
        stop 1
    end if

    !** FREE HANDLE MEMORY
    if (DTRNLSP_DELETE(handle) /= TR_SUCCESS) then
        call MKL_FREE_BUFFERS

        write (*,*) '| ERROR IN DTRNLSP_DELETE'
        stop 1
    end if

    ! at this point we have extracted all information and can free any memory still used by mkl
    call MKL_FREE_BUFFERS

    if (st_cr == 1) then
        write (*,*) "No solution found after ", iter, " iterations, reducing time step."
        pnewdt = 0.5
        return
    end if
    
    ! check the residual explicitly
    call local_residual(res, jac, minusdresdtemp, fbainite, 1, nkappa, kappa, &
						kappan, dtime, tempn, dtemp, mdata)
    
    if (sqrt(dot_product(res, res)) > eps(2)) then
		pnewdt = 0.5d0
		return
	end if
    
    ! update the jacobian at the solution point
    call local_residual(res, jac, minusdresdtemp, fbainite, 2, nkappa, kappa, &
						kappan, dtime, tempn, dtemp, mdata)

    ! check if the change of phase is small enough
    dbetaM = kappa(1) - kappan(1)
    dbetaB = kappa(2) - kappan(2)

    if (abs(dbetaM) > dbetaMmax) then
        write (*,*) "Change in betaM was ", dbetaM, ", reducing time step."
        pnewdt = 0.2
        return
    end if

    if (abs(dbetaB) > dbetaBmax) then
        write (*,*) "Change in betaB was ", dbetaB, ", reducing time step."
        pnewdt = 0.2
        return
    end if

	! implicit function theorem: solve the system of equations
	!    dRdkappa * dkappadT = - dRdT
	! to obtain the derivative of the SDVs w.r.t. temperature
    call local_residual(res, jac, minusdresdtemp, fbainite, 3, nkappa, kappa, &
                        kappan, dtime, tempn, dtemp, mdata)
    minusdresdtemp(:) = -minusdresdtemp(:)
    
    call solve_linear_system(dkappadTemp, info, nkappa, jac, minusdresdtemp)
    
    
end subroutine

