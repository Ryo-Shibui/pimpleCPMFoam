#ifndef CREATE_BOUNDARY_VALUE_H
#define CREATE_BOUNDARY_VALUE_H

#include "fvMesh.H"
#include "volFields.H"
#include "surfaceFields.H"
#include <map>

namespace boundaryValue
{
    static std::map<word, label> patchIDMap;

    static const fvMesh* meshPtr = nullptr;
    static const surfaceScalarField* phiEFPtr = nullptr;
    static const surfaceScalarField* phiPtr = nullptr;

    inline void initPatchIDs(const fvMesh& mesh)
    {
        const polyBoundaryMesh& patches = mesh.boundaryMesh();
        patchIDMap.clear();
        forAll(patches, i)
        {
            patchIDMap[patches[i].name()] = i;
        }
        meshPtr = &mesh;
    }

    inline void setPhiEF(const surfaceScalarField& phiEF)
    {
        phiEFPtr = &phiEF;
    }

    inline void setPhi(const surfaceScalarField& phi)
    {
        phiPtr = &phi;
    }

    inline label patchID(const word& patchName)
    {
        auto it = patchIDMap.find(patchName);
        if (it == patchIDMap.end())
        {
            FatalErrorIn("boundaryValue::patchID")
                << "Patch name \"" << patchName << "\" not found." << exit(FatalError);
        }
        return it->second;
    }

    inline const fvsPatchScalarField& phiEF(const word& patchName)
    {
        return phiEFPtr->boundaryField()[patchID(patchName)];
    }

    inline const fvsPatchScalarField& phi(const word& patchName)
    {
        return phiPtr->boundaryField()[patchID(patchName)];
    }

    inline const scalarField& magA(const word& patchName)
    {
        return meshPtr->magSf().boundaryField()[patchID(patchName)];
    }

    inline const labelList& faceCells(const word& patchName)
    {
        return meshPtr->boundaryMesh()[patchID(patchName)].faceCells();
    }

    inline const scalarField& snGrad(const surfaceScalarField& snGradField, const word& patchName)
    {
        return snGradField.boundaryField()[patchID(patchName)];
    }

    inline scalarField cells(const volScalarField& field, const labelList& faceCells)
    {
        scalarField vals(faceCells.size());
        forAll(faceCells, i)
        {
            vals[i] = field[faceCells[i]];
        }
        return vals;
    }
}

#endif
