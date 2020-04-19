mu = pot - epsSq * fvc::laplacian(alpha);
mu -= fvc::domainIntegrate(mu-pot)/vol;
mu.correctBoundaryConditions();
