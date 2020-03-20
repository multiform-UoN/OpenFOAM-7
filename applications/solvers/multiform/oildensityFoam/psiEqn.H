  //Psi Evolution Equation
  while (pimple.correctNonOrthogonal())
  {
      fvScalarMatrix psiEqn
    (

           fvm::ddt(psi)
         + fvm::div(phipsi,psi)
         - fvm::laplacian(Dpsi,psi)
         ==
				 - fvm::Sp(SpsiImp,psi)
				 - SpsiExp
    );

    psiEqn.relax();
    psiEqn.solve();
    //psi = 0.5*(psi+mag(psi));
    //psi.correctBoundaryConditions();
  }
