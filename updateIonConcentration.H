#ifndef UPDATE_ION_CONCENTRATION_H
#define UPDATE_ION_CONCENTRATION_H

#include "createIonPatchData.H"

inline void updateIonConcentration
(
    const IonPatchData& patch,
    const fvMesh& mesh,
    const scalar dt,
    const scalarField& jDriftP,
    const scalarField& jDriftN,
    const scalarField& jConvP,
    const scalarField& jConvN,
    const scalarField& jDiffP,
    const scalarField& jDiffN,
    volScalarField& nP,
    volScalarField& nN,
    scalar& totalNP,
    scalar& totalNN
)
{
    totalNP = 0.0;
    totalNN = 0.0;

    forAll(patch.faceCells, i)
    {
        label cellI = patch.faceCells[i];

        scalar dNP = (jDriftP[i] + jConvP[i] + jDiffP[i]) * dt;
        scalar dNN = (jDriftN[i] + jConvN[i] + jDiffN[i]) * dt;

        scalar dnP = dNP / mesh.V()[cellI];
        scalar dnN = dNN / mesh.V()[cellI];

        totalNP += dNP;
        totalNN += dNN;

        nP[cellI] = max(nP[cellI] - dnP, 0.0);
        nN[cellI] = max(nN[cellI] - dnN, 0.0);
    }
}

#endif
