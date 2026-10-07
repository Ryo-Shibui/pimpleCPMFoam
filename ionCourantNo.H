scalar ionCoNum = 0.0;
scalar meanIonCoNum = 0.0;

surfaceVectorField Ef = fvc::interpolate(E);
surfaceScalarField nPf = fvc::interpolate(nP);
surfaceScalarField nNf = fvc::interpolate(nN);

surfaceScalarField snGradP = fvc::snGrad(nP);
surfaceScalarField snGradN = fvc::snGrad(nN);

dimensionedScalar epsP("epsP", nPf.dimensions(), SMALL);
dimensionedScalar epsN("epsN", nNf.dimensions(), SMALL);
    
surfaceScalarField driftMagnitudeP = mag(muP * nPf * (Ef & mesh.Sf()));
surfaceScalarField phiDriftP = muP * (Ef & mesh.Sf()) * pos(driftMagnitudeP);
surfaceScalarField driftMagnitudeN = mag(muN * nNf * (Ef & mesh.Sf()));
surfaceScalarField phiDriftN = -muN * (Ef & mesh.Sf()) * pos(driftMagnitudeN);
    
//surfaceScalarField phiDiffusionP = -DP * snGradP * mesh.magSf() * pos(nPf) / max(nPf, epsP);
//surfaceScalarField phiDiffusionN = -DN * snGradN * mesh.magSf()* pos(nNf) / max(nNf, epsN);

surfaceScalarField phiTotalP = phi + phiDriftP;
surfaceScalarField phiTotalN = phi + phiDriftN;

scalarField sumPhiP = fvc::surfaceSum(mag(phiTotalP))().primitiveField();
scalarField sumPhiN = fvc::surfaceSum(mag(phiTotalN))().primitiveField();

scalarField maxSumPhi = Foam::max(sumPhiP, sumPhiN);
    
ionCoNum = 0.5 * gMax(maxSumPhi / mesh.V().field()) * runTime.deltaTValue();
meanIonCoNum = 0.5 * (gSum(maxSumPhi) / gSum(mesh.V().field())) * runTime.deltaTValue();


Info << "Ion Courant Number mean: " << meanIonCoNum << " max: " << ionCoNum << endl;
