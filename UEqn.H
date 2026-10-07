// Solve the Momentum equation

tmp<fvVectorMatrix> tUEqn
(
    fvm::ddt(U) 
	+ fvm::div(phi,U)
	+ turbulence->divDevReff(U)
	- pc::e*(nP-nN)*E/rho0
	==
    fvOptions(U)
);
fvVectorMatrix& UEqn = tUEqn.ref();

UEqn.relax();

fvOptions.constrain(UEqn);

if (pimple.momentumPredictor())
{
    solve(UEqn == -fvc::grad(p));

    fvOptions.correct(U);
}
