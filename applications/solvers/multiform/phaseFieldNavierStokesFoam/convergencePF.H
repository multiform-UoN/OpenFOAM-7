//- Check convergence
{
    //- Residuals on mu
    scalar resMu
    (
        fvc::domainIntegrate(mag(mu-mu.prevIter())).value()
        /fvc::domainIntegrate(mag(mu)).value()
    );

    //- Residuals on alpha
    scalar resAlpha
    (
        fvc::domainIntegrate(mag(alpha-alpha.prevIter())).value()
        /fvc::domainIntegrate(mag(alpha)).value()
    );

    Info<< "Residulas on mu: " << resMu
        << "    Residulas on alpha: " << resAlpha << endl;

    if(resMu < muTol && resAlpha < alphaTol )
    {
        converged = true;
        Info<<"Phase field corrector converged in "
            << corr << " iterations" << endl;
    }
    else if(corr == pfMaxCorr)
    {
        converged = true;
        Info<<"\nPhase field corrector NOT converged in "
            << pfMaxCorr << " iterations" << endl;
    }

    Info<<endl;

}
