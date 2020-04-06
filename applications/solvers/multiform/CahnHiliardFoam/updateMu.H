mu = pot - sqr(epsilon) * fvc::laplacian(alpha);
mu -= fvc::domainIntegrate(mu)/vol;
mu.correctBoundaryConditions();
