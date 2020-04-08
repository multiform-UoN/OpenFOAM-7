const scalar sumS((fvc::domainIntegrate(- pot + mu)).value());

mu -= sumS;
mu.correctBoundaryConditions();
