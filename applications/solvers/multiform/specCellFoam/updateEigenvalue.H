scalar lambdaOld = lambda.value();
scalar lambdaAdjOld = lambdaAdj.value();

{
    dimensionedScalar one("one",dimless/dimTime,1.);
    dimensionedScalar smallV
    (
        "smallV",
        psi.dimensions()*psi.dimensions()*dimVol,
        small
    );


    if(implicitRayleigh)
    {
        volScalarField LpsiAdj
        (
            psiAdjEqn.A()*psiAdj - psiAdjEqn.H()
        );

        volScalarField Lpsi(psiEqn.A()*psi - psiEqn.H());

        lambda =
        (
            fvc::domainIntegrate(Lpsi*psi)
            /
            fvc::domainIntegrate(psi*psi)
        );


        lambdaAdj =
        (
            fvc::domainIntegrate(LpsiAdj*psiAdj)
            /
            fvc::domainIntegrate(psiAdj*psiAdj)
        );

    }
    else
    {
        lambda =
        (
            lambda*fvc::domainIntegrate(psi.prevIter()*psi)
            /
            fvc::domainIntegrate(psi*psi)
        );


        lambdaAdj =
        (
            lambdaAdj*fvc::domainIntegrate(psiAdj.prevIter()*psiAdj)
            /
            fvc::domainIntegrate(psiAdj*psiAdj)
        );

    }

    // //- Evaluate lambda as the average eigenvalue
    // lambda = 0.9*lambda + 0.1* (lambda + lambdaAdj)/scalar(2.);

    //- Calculate new residual
    lambdaRes = mag(lambdaOld - lambda.value())/mag(lambdaOld);

}
    //- Evaluate error on eigenvalues
    dimensionedScalar eigenErr = mag(lambdaAdj -lambda)/(mag(lambda));
