Info<< "Energy : " <<
    (
        fvc::domainIntegrate
        (
            sqr(sqr(alpha)-scalar(1.))/scalar(4.)
            + sqr(epsilon)*magSqr(fvc::grad(alpha))/scalar(2.)
        )
    ).value() << endl;

Info<< "Total mass : " << (fvc::domainIntegrate(alpha)).value() << endl;
