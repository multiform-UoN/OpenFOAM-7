// Concentration Evolution Equation
while (pimple.correctNonOrthogonal())
{
      fvScalarMatrix cEqn
    (

         fvm::ddt(c)
         + fvm::div(phic,c)
         - fvm::laplacian(Dc,c)
    );

    cEqn.relax();
    cEqn.solve();

    c.relax();
    c.correctBoundaryConditions();
}
