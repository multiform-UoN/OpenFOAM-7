mu = pot - sqr(epsilon) * fvc::laplacian(alpha);
mu.correctBoundaryConditions();
