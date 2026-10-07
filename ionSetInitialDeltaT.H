if (adjustTimeStep)
{
    if ((runTime.timeIndex() == 0) && (ionCoNum > SMALL))
    {
        runTime.setDeltaT
        (
            min
            (
                maxCo * runTime.deltaTValue() / ionCoNum,
                min(runTime.deltaTValue(), maxDeltaT)
            )
        );
    }
}
