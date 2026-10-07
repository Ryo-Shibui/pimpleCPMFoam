#ifndef CREATE_ION_PATCH_DATA_H
#define CREATE_ION_PATCH_DATA_H

#include "createBoundaryValue.H"

struct IonPatchData
{
    word name;
    tmp<surfaceScalarField> tGradNP;
    tmp<surfaceScalarField> tGradNN;
    labelList faceCells;
    scalarField phiEFbf;
    scalarField phibf;
    scalarField magA;
    scalarField snGradPbf;
    scalarField snGradNbf;
    scalarField nPcells;
    scalarField nNcells;

    IonPatchData(const word& patchName, const volScalarField& nP, const volScalarField& nN)
    :
        name(patchName),
        tGradNP(fvc::snGrad(nP)),
        tGradNN(fvc::snGrad(nN)),
        faceCells(boundaryValue::faceCells(patchName)),
        phiEFbf(boundaryValue::phiEF(patchName).primitiveField()),
        phibf(boundaryValue::phi(patchName).primitiveField()),
        magA(boundaryValue::magA(patchName)),
        snGradPbf(boundaryValue::snGrad(tGradNP(), patchName)),
        snGradNbf(boundaryValue::snGrad(tGradNN(), patchName)),
        nPcells(boundaryValue::cells(nP, faceCells)),
        nNcells(boundaryValue::cells(nN, faceCells))
    {}
};

#endif
