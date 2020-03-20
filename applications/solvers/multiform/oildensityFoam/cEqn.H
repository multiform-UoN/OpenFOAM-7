//Concentration Evolution Equation

  while (pimple.correctNonOrthogonal())
  {
   fvScalarMatrix cEqn
   (
        fvm::ddt(w,c)
      + fvm::div(phic,c)
      - fvm::laplacian(Dc*w,c)
      ==
      - fvm::Sp(w*ScImp,c)
      - w*ScExp
    );
   cEqn.relax();
   cEqn.solve();

  }
