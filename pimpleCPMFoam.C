/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2019 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    pimpleFoam.C

Group
    grpIncompressibleSolvers

Description
    Transient solver for incompressible, turbulent flow of Newtonian fluids
    on a moving mesh.

    \heading Solver details
    The solver uses the PIMPLE (merged PISO-SIMPLE) algorithm to solve the
    continuity equation:

        \f[
            \div \vec{U} = 0
        \f]

    and momentum equation:

        \f[
            \ddt{\vec{U}} + \div \left( \vec{U} \vec{U} \right) - \div \gvec{R}
          = - \grad p + \vec{S}_U
        \f]

    Where:
    \vartable
        \vec{U} | Velocity
        p       | Pressure
        \vec{R} | Stress tensor
        \vec{S}_U | Momentum source
    \endvartable

    Sub-models include:
    - turbulence modelling, i.e. laminar, RAS or LES
    - run-time selectable MRF and finite volume options, e.g. explicit porosity

    \heading Required fields
    \plaintable
        U       | Velocity [m/s]
        p       | Kinematic pressure, p/rho [m2/s2]
        \<turbulence fields\> | As required by user selection
    \endplaintable

Note
   The motion frequency of this solver can be influenced by the presence
   of "updateControl" and "updateInterval" in the dynamicMeshDict.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "physicalConstants.H"
#include "dynamicFvMesh.H"	// it potential unnesessary
#include "singlePhaseTransportModel.H"	
#include "turbulentTransportModel.H"
#include "pimpleControl.H"
#include "CorrectPhi.H"
#include "fvOptions.H"
#include "localEulerDdtScheme.H"
#include "fvcSmooth.H"
//#include "createBoundaryValue.H"
//#include "createIonPatchData.H"
//#include "computeIonFlux.H"
//#include "updateIonConcentration.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Transient solver for incompressible, turbulent flow"
        " of Newtonian fluids on a moving mesh."
    );

    #include "postProcess.H"

    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createDynamicFvMesh.H"
    #include "initContinuityErrs.H"
    #include "createDyMControls.H"
    #include "createFields.H"
    #include "createUfIfPresent.H"                           
    #include "CourantNo.H"
    #include "ionCourantNo.H"             // CoNum based on ion transport velocity
    //#include "setInitialDeltaT.H"       // If you want to set deltaT based on ionCourantNo.H, this is unnecessary
    #include "ionSetInitialDeltaT.H"      // If you want to set deltaT based on CourantNo.H, this is unnecessary

    const label CPM = mesh.boundaryMesh().findPatchID("CPM");
    const bool hasCPM = (CPM >= 0);

    if (!hasCPM)
    {
        WarningInFunction
            << "Patch 'CPM' not found. Skipping ion trapping/current accounting and phiE update on CPM."
            << nl << "Available patches: " << mesh.boundaryMesh().names() << nl << endl;
    }

    const label walls = mesh.boundaryMesh().findPatchID("walls");
    const bool hasWalls = (walls >= 0);

    if (!hasWalls)
    {
        WarningInFunction
            << "Patch 'walls' not found. Skipping ion trapping/current accounting on walls."
            << nl << "Available patches: " << mesh.boundaryMesh().names() << nl << endl;
    }

    scalar phiE_CPM_accum = 0.0;

    if (hasCPM)
    {
        const fvPatchField<scalar>& phiPatch = phiE.boundaryField()[CPM];
        if (phiPatch.type() == "fixedValue")
        {
            const fixedValueFvPatchScalarField& fixedPatch =
                refCast<const fixedValueFvPatchScalarField>(phiPatch);

            if (fixedPatch.size() > 0)
            {
                phiE_CPM_accum = average(fixedPatch);
            }

            reduce(phiE_CPM_accum, maxOp<scalar>(), UPstream::worldComm);
            Info << "Initial phiE_CPM_accum = " << phiE_CPM_accum << " V" << endl;
        }
        else
        {
            WarningInFunction
                << "Patch 'CPM' exists but phiE is not fixedValue on it (type="
                << phiPatch.type() << "). Skipping initialization of phiE_CPM_accum." << nl << endl;
        }
    }


    turbulence->validate();

    word adjustTimeType = runTime.controlDict().lookupOrDefault<word>("adjustTimeType", "phi");

    if (!LTS)
    {
        #include "CourantNo.H"
        #include "ionCourantNo.H"
        
        if (adjustTimeType == "ionBased")
        {
            #include "ionSetInitialDeltaT.H" 
        }
        else if (adjustTimeType == "phi")
        {
            #include "setInitialDeltaT.H"
        }
    }
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "readDyMControls.H"

        if (LTS)
        {
            //#include "setRDeltaT.H"          // If you want to set deltaT based on ionCourantNo.H, this is unnecessary
            #include "ionSetRDeltaT.H"         // If you want to set deltaT based on CourantNo.H, this is unnecessary
        }
        else if (adjustTimeStep)
        {
            #include "CourantNo.H"
            #include "ionCourantNo.H"

            if (adjustTimeType == "ionBased")
            {
                #include "ionSetDeltaT.H"
            }
            else if (adjustTimeType == "phi")
            {
                #include "setDeltaT.H"
            }
        }
        else
        {
            #include "CourantNo.H"
            #include "ionCourantNo.H"
            #include "setDeltaT.H"
        }


        ++runTime;

        Info<< "Time = " << runTime.timeName() << nl << endl;

		// --- Pressure-velocity PIMPLE corrector loop
		while (pimple.loop())
		{

			// --- Poisson's equation
		    solve
		    (
		        fvm::laplacian(phiE) 
                + pc::e*(nP-nN)/pc::epsilon0
		    );
			E = (-1.0)*fvc::grad(phiE);

			// --- Momentum equation
			#include "UEqn.H" // solving momentum equation for U

			// --- Pressure corrector loop
			while (pimple.correct())
			{
				#include "pEqn.H"	// continuity equation for p
			}

			// --- Charge conservation
            phiEF = fvc::flux(E);
            fvScalarMatrix nPEqn
            (
                fvm::ddt(nP)
                + fvm::div(phi, nP)
                + fvm::div(muP*phiEF, nP)
                - fvm::laplacian(DP, nP)
                + beta * nP * nN
            );
            nPEqn.relax();
            nPEqn.solve();
            nP = Foam::max(nP, dimensionedScalar("zero", nP.dimensions(), 0.0));

            fvScalarMatrix nNEqn
            (
                fvm::ddt(nN)
                + fvm::div(phi, nN)
                - fvm::div(muN*phiEF, nN)
                - fvm::laplacian(DN, nN)
                + beta * nP * nN
            );
            nNEqn.relax();
            nNEqn.solve();
            nN = Foam::max(nN, dimensionedScalar("zero", nN.dimensions(), 0.0));

            rho = pc::e*(nP-nN);
            rho.correctBoundaryConditions();

			if (pimple.turbCorr())
            {
                laminarTransport.correct();
                turbulence->correct();
            }
		}

        if (hasCPM || hasWalls)
        {
            const scalar dt = mesh.time().deltaTValue();
            const scalar eVal   = pc::e.value();
            const scalar muPVal = muP.value();
            const scalar muNVal = muN.value();
            const scalar DPVal  = DP.value();
            const scalar DNVal  = DN.value();

            tmp<surfaceScalarField> tSnGradNP = fvc::snGrad(nP);
            tmp<surfaceScalarField> tSnGradNN = fvc::snGrad(nN);
            const surfaceScalarField& snGradNP = tSnGradNP();
            const surfaceScalarField& snGradNN = tSnGradNN();

            // ---- CPM 処理（存在するときのみ）----
            if (hasCPM)
            {
                const fvsPatchScalarField& phiEFbf = phiEF.boundaryField()[CPM];
                const fvsPatchScalarField& phibf   = phi.boundaryField()[CPM];
                const scalarField& magA            = mesh.magSf().boundaryField()[CPM];
                const scalarField& snGradPbf       = snGradNP.boundaryField()[CPM];
                const scalarField& snGradNbf       = snGradNN.boundaryField()[CPM];
                const labelList& faceCells         = mesh.boundaryMesh()[CPM].faceCells();

                scalarField nPcells(faceCells.size());
                scalarField nNcells(faceCells.size());
                forAll(faceCells, i)
                {
                    nPcells[i] = nP[faceCells[i]];
                    nNcells[i] = nN[faceCells[i]];
                }

                // convective [/s]
                scalarField jConvP = nPcells * phibf;
                scalarField jConvN = nNcells * phibf;

                // drift [/s]
                scalarField jDriftP = muPVal * phiEFbf * nPcells;
                scalarField jDriftN = (-1) * muNVal * phiEFbf * nNcells;

                // diffusive [/s]
                scalarField jDiffP = (-1) * DPVal * snGradPbf * magA;
                scalarField jDiffN = (-1) * DNVal * snGradNbf * magA;

                // total current [C/s]
                scalarField Jtotal = eVal * (jDriftP + jConvP + jDiffP - jDriftN - jConvN - jDiffN);

                scalar Inet = sum(Jtotal);
                reduce(Inet, sumOp<scalar>(), UPstream::worldComm);

                // Δφ = +Δt * I / C  （符号はあなたのモデル定義に従う）
                scalar dPhi = dt * Inet / Ccap.value();
                if (Pstream::master())
                {
                    phiE_CPM_accum += dPhi;
                }
                Pstream::broadcast(phiE_CPM_accum);

                // CPM面への“トラップ反映”（セルから減算）と集計
                scalar totalNP = 0.0;
                scalar totalNN = 0.0;

                forAll(faceCells, i)
                {
                    const label cellI = faceCells[i];

                    scalar dNP = (jDriftP[i] + jConvP[i] + jDiffP[i]) * dt;
                    scalar dNN = (jDriftN[i] + jConvN[i] + jDiffN[i]) * dt;

                    scalar dnP = dNP / mesh.V()[cellI];
                    scalar dnN = dNN / mesh.V()[cellI];

                    if (dnP > 0)
                    {
                        nP[cellI] = max(nP[cellI] - dnP, 0.0);
                        totalNP += dNP;
                    }
                    if (dnN > 0)
                    {
                        nN[cellI] = max(nN[cellI] - dnN, 0.0);
                        totalNN += dNN;
                    }
                }

                reduce(totalNP, sumOp<scalar>(), UPstream::worldComm);
                reduce(totalNN, sumOp<scalar>(), UPstream::worldComm);

                Info<< "  [CPM] phiE increment due to ion deposition (dPhi) = " << dPhi << " V" << endl;
                Info<< "  [CPM] Number of positive ions deposited = " << max(totalNP, 0.0) << endl;
                Info<< "  [CPM] Number of negative ions deposited = " << max(totalNN, 0.0) << endl;

                // CPM の phiE 更新（fixedValue の場合のみ）
                fvPatchField<scalar>& phiPatchRef = phiE.boundaryFieldRef()[CPM];
                if (phiPatchRef.type() == "fixedValue")
                {
                    fixedValueFvPatchScalarField& fixedPatch =
                        refCast<fixedValueFvPatchScalarField>(phiPatchRef);

                    fixedPatch == scalarField(fixedPatch.size(), phiE_CPM_accum);

                    Info<< "  [CPM] Updated surface potential (phiE) = " << phiE_CPM_accum << " V" << endl;
                }
                else
                {
                    WarningInFunction
                        << "Patch 'CPM' exists but phiE is not fixedValue (type="
                        << phiPatchRef.type() << "). Skipping phiE update." << nl << endl;
                }
            }

            // ---- walls 処理（存在するときのみ）----
            if (hasWalls)
            {
                const fvsPatchScalarField& wallphiEFbf = phiEF.boundaryField()[walls];
                const fvsPatchScalarField& wallphibf   = phi.boundaryField()[walls];
                const scalarField& wallmagA            = mesh.magSf().boundaryField()[walls];
                const scalarField& wallsnGradPbf       = snGradNP.boundaryField()[walls];
                const scalarField& wallsnGradNbf       = snGradNN.boundaryField()[walls];
                const labelList& wallFaceCells         = mesh.boundaryMesh()[walls].faceCells();

                scalarField wallnPcells(wallFaceCells.size());
                scalarField wallnNcells(wallFaceCells.size());
                forAll(wallFaceCells, i)
                {
                    wallnPcells[i] = nP[wallFaceCells[i]];
                    wallnNcells[i] = nN[wallFaceCells[i]];
                }

                scalarField walljConvP = wallnPcells * wallphibf;
                scalarField walljConvN = wallnNcells * wallphibf;

                scalarField walljDriftP = muPVal * wallphiEFbf * wallnPcells;
                scalarField walljDriftN = (-1) * muNVal * wallphiEFbf * wallnNcells;

                scalarField walljDiffP = (-1) * DPVal * wallsnGradPbf * wallmagA;
                scalarField walljDiffN = (-1) * DNVal * wallsnGradNbf * wallmagA;

                scalar walltotalNP = 0.0;
                scalar walltotalNN = 0.0;

                forAll(wallFaceCells, i)
                {
                    const label cellI = wallFaceCells[i];

                    scalar dNP = (walljDriftP[i] + walljConvP[i] + walljDiffP[i]) * dt;
                    scalar dNN = (walljDriftN[i] + walljConvN[i] + walljDiffN[i]) * dt;

                    scalar dnP = dNP / mesh.V()[cellI];
                    scalar dnN = dNN / mesh.V()[cellI];

                    if (dnP > 0)
                    {
                        nP[cellI] = max(nP[cellI] - dnP, 0.0);
                        walltotalNP += dNP;
                    }
                    if (dnN > 0)
                    {
                        nN[cellI] = max(nN[cellI] - dnN, 0.0);
                        walltotalNN += dNN;
                    }
                }

                reduce(walltotalNP, sumOp<scalar>(), UPstream::worldComm);
                reduce(walltotalNN, sumOp<scalar>(), UPstream::worldComm);

                Info<< "  [walls] Number of positive ions deposited = " << max(walltotalNP, 0.0) << endl;
                Info<< "  [walls] Number of negative ions deposited = " << max(walltotalNN, 0.0) << endl;
            }
        }
		runTime.write();
		runTime.printExecutionTime(Info);

    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
