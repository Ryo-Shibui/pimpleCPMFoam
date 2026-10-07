#ifndef COMPUTE_ION_CURRENT_H
#define COMPUTE_ION_CURRENT_H

#include "createIonPatchData.H"

inline void computeIonCurrent
(
    const IonPatchData& patch,
    const scalar muPVal,
    const scalar muNVal,
    const scalar DPVal,
    const scalar DNVal,
    scalarField& jDriftP,
    scalarField& jDriftN,
    scalarField& jConvP,
    scalarField& jConvN,
    scalarField& jDiffP,
    scalarField& jDiffN
)
{
    jConvP = patch.nPcells * patch.phibf;
    jConvN = patch.nNcells * patch.phibf;

    jDriftP = muPVal * patch.phiEFbf * patch.nPcells;
    jDriftN = (-1) * muNVal * patch.phiEFbf * patch.nNcells;

    jDiffP = (-1) * DPVal * patch.snGradPbf * patch.magA;
    jDiffN = (-1) * DNVal * patch.snGradNbf * patch.magA;
}

#endif
