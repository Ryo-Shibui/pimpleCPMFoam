{
    volScalarField& rDeltaT = trDeltaT.ref();

    const dictionary& pimpleDict = pimple.dict();

    scalar maxCo
    (
        pimpleDict.getOrDefault<scalar>("maxCo", 0.8)
    );

    scalar rDeltaTSmoothingCoeff
    (
        pimpleDict.getOrDefault<scalar>("rDeltaTSmoothingCoeff", 0.02)
    );

    scalar rDeltaTDampingCoeff
    (
        pimpleDict.getOrDefault<scalar>("rDeltaTDampingCoeff", 1.0)
    );

    scalar maxDeltaT
    (
        pimpleDict.getOrDefault<scalar>("maxDeltaT", GREAT)
    );

    volScalarField rDeltaT0("rDeltaT0", rDeltaT);

    surfaceScalarField phiIon = Foam::max(mag(phiTotalP), mag(phiTotalN));

    // rDeltaTの計算（= 1/deltaT）: 2*maxCo は内部面補正
    rDeltaT.ref() = max
    (
        1 / dimensionedScalar("maxDeltaT", dimTime, maxDeltaT),
        fvc::surfaceSum(phiIon)()() / ((2 * maxCo) * mesh.V())
    );

    rDeltaT.correctBoundaryConditions();

    Info<< "Ion-based time scale min/max = "
        << gMin(1 / rDeltaT.primitiveField())
        << ", " << gMax(1 / rDeltaT.primitiveField()) << endl;

    // スムージング
    if (rDeltaTSmoothingCoeff < 1.0)
    {
        fvc::smooth(rDeltaT, rDeltaTSmoothingCoeff);
    }

    Info<< "Smoothed ion time scale min/max = "
        << gMin(1 / rDeltaT.primitiveField())
        << ", " << gMax(1 / rDeltaT.primitiveField()) << endl;

    // ダンピング（急激な時間拡大を抑制）
    if
    (
        rDeltaTDampingCoeff < 1.0
     && runTime.timeIndex() > runTime.startTimeIndex() + 1
    )
    {
        rDeltaT =
            rDeltaT0
           *max(rDeltaT / rDeltaT0, scalar(1) - rDeltaTDampingCoeff);

        Info<< "Damped ion time scale min/max = "
            << gMin(1 / rDeltaT.primitiveField())
            << ", " << gMax(1 / rDeltaT.primitiveField()) << endl;
    }
}
