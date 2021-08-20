{
    // Solve the energy equation
    volTensorField gradU = fvc::grad(U);
    tmp <volScalarField> nu = laminarTransport.nu();
    // Equation with viscous dissipation 
    tmp <volTensorField> tau = -nu * rho * (gradU + gradU.T()) + (2.0/3.0 * nu * rho * fvc::div(U)) * I;
    tmp <volScalarField> tauGradU = tau && gradU;
	
    fvScalarMatrix TEqn
    (
        fvm::div(phi, T)
      - fvm::laplacian((lambda/(cp*rho)), T)
      //+ (- laminarTransport.nu()*rho*twoSymm(fvc::grad(U)) && (fvc::grad(U)))/(rho*cp)
      + tauGradU/(cp*rho)
    );
	
    // Energy solution
// +MAZ
    TEqn.relax();
    TEqn.solve();
// -MAZ
}


// #########################################################################
// KTM
// #########################################################################
	// Equation with viscous dissipation complete formulation inkl. div(U)
/*
    tmp<fvScalarMatrix> TEqn
    (
        fvm::div(phi, T)
      - fvm::laplacian((lambda/(cp*rho)), T)
      + ((- laminarTransport.nu()*rho*twoSymm(fvc::grad(U))+
    		(2.0/3.0*laminarTransport.nu()*rho*fvc::div(U))*I)
			 && (fvc::grad(U)))/(rho*cp)
	);
 */

	// Equation with viscous dissipation for incompressible fluids: div(U)=0
/*     tmp<fvScalarMatrix> TEqn
    (
        fvm::div(phi, T)
      - fvm::laplacian((lambda/(cp*rho)), T)
      - (laminarTransport.nu()*rho*twoSymm(fvc::grad(U)) && (fvc::grad(U)))/(rho*cp)
	); */
