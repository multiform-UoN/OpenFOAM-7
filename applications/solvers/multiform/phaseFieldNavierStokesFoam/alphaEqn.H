{
    //- Create matrix for mu (positive definite)
    fvScalarMatrix muEqn
    (
        - fvm::laplacian(M,mu)
    );

    //- Get diagonal (positive)
    volScalarField Amu(muEqn.A());

    //- Initialize loop controls
    label corr(0);
    bool  converged(false);
    surfaceScalarField phi(fvc::flux(U));
    volScalarField divU(fvc::div(phi));

    do
    {
        corr++;

        //- Get off-diagonal terms
        volScalarField Hmu(muEqn.H());
        alpha.storePrevIter();

        while(pimple.correctNonOrthogonal())
        {
            //- Assemble matrix for alpha (positive definite)
            fvScalarMatrix alphaEqn
            (
                fvm::ddt(alpha)
              + fvm::div(phi,alpha,"div(phi,alpha)")
              - fvc::Sp(divU,alpha)
              - Amu*
                (
                    fvm::laplacian(epsSq,alpha)
                )
            );

            alphaEqn.relax();

            //- Solve for alpha
            solve(alphaEqn == Hmu - Amu*fPrime);

        }
        //- Relax alpha
        alpha.relax();

        //- Update fPrime
        fPrime = fea*pow(alpha,3.) - feb*(alpha);

        mu.storePrevIter();
        //- Update potential
        mu = -fvc::laplacian(epsSq,alpha) + fPrime;

        mu.relax();
        mu.correctBoundaryConditions();

        #include "convergencePF.H"


    }
    while
    (
        !converged
    );

    #include "writeEnergyAndMass.H"

    //- Update real indicator field
    alpha1 = (alpha +1.0)/0.5;
    alpha1.correctBoundaryConditions();
}
