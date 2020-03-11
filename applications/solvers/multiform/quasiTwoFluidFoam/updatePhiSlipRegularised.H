while (pimple.correctNonOrthogonal())
{
    const volVectorField Uslip0
          ( U - (St
            *
            (
              fvc::ddt(U)
              +
              (U & fvc::grad(U))
            ))
          );

    // Momentum predictor
    tmp<fvVectorMatrix> tUslipEqn
    (
       // Stokes approximation
       // fvm::ddt(Uslip)
       // +
       // fvm::div(phi, Uslip)
       -
       fvm::laplacian(c*nu, Uslip)
       -
       fvc::div(c*nu*dev2(fvc::grad(Uslip)().T()))
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
