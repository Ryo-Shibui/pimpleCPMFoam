if (adjustTimeStep)
{
    scalar maxDeltaTFact = maxCo / (ionCoNum + SMALL);
    scalar deltaTFact = min(min(maxDeltaTFact, 1.0 + 0.1*maxDeltaTFact), 1.2);

    runTime.setDeltaT
    (
        min
        (
            deltaTFact * runTime.deltaTValue(),
            maxDeltaT
        )
    );

    Info<< "deltaT (ion-based) = " << runTime.deltaTValue() << endl;
}
