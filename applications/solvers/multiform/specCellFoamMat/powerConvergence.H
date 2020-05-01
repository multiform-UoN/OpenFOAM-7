// - Check convergence
{
    currPowerIter++;

    powerConverged =
        (
            (
                lambdaRes < eigenValTol
                //||
                //Required for convergence with Neumann BCs
            //    lambda.value() < scalar(1e-6)
            )
            &&
            (
              psiRes < eigenFunTol
            &&
              psiAdjRes < eigenFunTol
            )
        )
        ||
        (
            currPowerIter > maxPowerIter
        );

}
