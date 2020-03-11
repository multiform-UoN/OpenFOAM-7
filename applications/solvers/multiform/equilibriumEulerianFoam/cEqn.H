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
    // Guarantees boundedness (I'd rather use a bounded scheme)
    c = 0.5*(c+mag(c));
    // Guarantee that it is lower than cm-0.1 (so it dont reach a singularity):
    //c = min(c,cm-0.1); (do we really need this?)
    c.correctBoundaryConditions();
}
