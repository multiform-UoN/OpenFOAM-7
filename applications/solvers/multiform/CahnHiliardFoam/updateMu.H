mu = pot - epsilon*epsilon * fvc::laplacian(alpha);
mu.correctBoundaryConditions();
