Info<< "Avg Energy : " <<
    (
        fvc::domainIntegrate
        (
            sqr(sqr(alpha)-scalar(1.))/scalar(4.)
            + epsSq*magSqr(fvc::grad(alpha))/scalar(2.)
        )/vol
    ).value() << endl;

Info<< "Avg mass : " << (fvc::domainIntegrate(alpha)/vol).value() << endl;

Info << "Avg (mu-pot) : " << (fvc::domainIntegrate(mu-pot)/vol).value() << endl;
