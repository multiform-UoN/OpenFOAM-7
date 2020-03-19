Info<< "Energy = " <<
      (
        fvc::domainIntegrate
        (
        pow(pow(alpha,2)-scalar(1),2)/scalar(4)
        + pow(epsilon,2)*magSqr(fvc::grad(alpha))/scalar(2)
        )
      ).value() << endl;

Info<< "Total mass = " <<
      (
        fvc::domainIntegrate
        (
        alpha
        )
      ).value() << endl;
