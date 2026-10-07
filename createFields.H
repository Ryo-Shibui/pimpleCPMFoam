#include "createRDeltaT.H"
#include "physicalConstants.H"

Info<< "Reading physicalProperties\n" << endl;

IOdictionary physicalProperties
(
    IOobject
    (
        "physicalProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

dimensionedScalar T
(
    "T",
    dimensionSet(0, 0, 0, 1, 0, 0, 0),
    physicalProperties
);

dimensionedScalar rho0
(
    "rho0",
    dimensionSet(1, -3, 0, 0, 0, 0, 0),
    physicalProperties
);

dimensionedScalar nu
(
    "nu",
    dimensionSet(0, 2, -1, 0, 0, 0, 0),
    physicalProperties
);

// Capacitance of charge plate surface
dimensionedScalar Ccap 
(
    "Ccap",
    dimensionSet(-1, -2, 4, 0, 0, 2, 0),
    physicalProperties
);

// ion mobility  [m2 V‑1 s‑1]  →  kg^-1 s^2 A   (dim: -1 0 2 0 0 1 0)
dimensionedScalar muP
(
    "muP",
    dimensionSet(-1, 0, 2, 0, 0, 1, 0),
    physicalProperties
);

dimensionedScalar muN
(
    "muN",
    dimensionSet(-1, 0, 2, 0, 0, 1, 0),
    physicalProperties
);

// ion diffusion coeff. [m2 s‑1]  →  (0 2 -1 0 0 0 0)
dimensionedScalar DP
(
    "DP",
    dimensionSet(0, 2, -1, 0, 0, 0, 0),
    physicalProperties
);

dimensionedScalar DN
(
    "DN",
    dimensionSet(0, 2, -1, 0, 0, 0, 0),
    physicalProperties
);

// ion‑ion recombination coeff. β [m3 s‑1] → (0 3 -1 0 0 0 0)
dimensionedScalar beta
(
	"beta",
	dimensionSet(0, 3, -1, 0, 0, 0, 0),
	physicalProperties
);


Info<< "Reading field p\n" << endl;
volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

#include "createPhi.H"


Info<< "Reading field electrostatic potential phiE\n" << endl;
volScalarField phiE
(
    IOobject
    (
        "phiE",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading/calculating field E (electric field)\n" << endl;
volVectorField E
(
    IOobject
    (
        "E",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    -fvc::grad(phiE)
);

Info<< "Reading field positive ion density nP\n" << endl;
volScalarField nP
(
    IOobject
    (
        "nP",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field negative ion density nN\n" << endl;
volScalarField nN
(
    IOobject
    (
        "nN",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);


Info<< "Reading/calculating face flux field of charge field rho\n" << endl;
volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    pc::e*(nN-nP)
);

Info<< "Reading/calculating face flux field of electric field phiEF\n" << endl; 
surfaceScalarField phiEF
(
    IOobject
    (
       "phiEF",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
 	fvc::flux(E)
);

label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell(p, pimple.dict(), pRefCell, pRefValue);
mesh.setFluxRequired(p.name());


singlePhaseTransportModel laminarTransport(U, phi);

autoPtr<incompressible::turbulenceModel> turbulence
(
    incompressible::turbulenceModel::New(U, phi, laminarTransport)
);

#include "createMRF.H"
#include "createFvOptions.H"
