while (pimple.correctNonOrthogonal())
{
    const volVectorField Uslip0
          ( U - (St
            *
            (I-St*Gr.T())   // this correction might not be compatible with this extended slip model
            &
            (
              fvc::ddt(U)
              +
              (U & fvc::grad(U))
            ))
          );

      // Momentum predictor
      tmp<fvVectorMatrix> tUslipEqn
      (
         fvm::ddt(Uslip)
         +
         fvm::div(phi, Uslip)
         -
         fvm::laplacian(psi*nu, Uslip)
         -
         fvc::div(psi*nu*dev2(fvc::grad(Uslip)().T()))
         ==
         -
         fvm::Sp(invSt,Uslip)
      );

      fvVectorMatrix& UslipEqn = tUslipEqn.ref();

      UslipEqn.relax();

      solve(UslipEqn == (Uslip0*invSt) );
}

Uslip = Uslip - U;
phislip = fvc::flux(Uslip);
