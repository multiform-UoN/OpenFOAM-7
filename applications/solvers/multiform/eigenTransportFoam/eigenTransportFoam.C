/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    eigenTransportFoam

Description
    Power method to solve an eigenvalue problem L(psi) = lambda psi, where 
    L is a transport (advection-diffusion) operator with velocity U and 
    diffusion coefficient Gamma.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "fvOptions.H"
#include "simpleControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCaseLists.H"

    #include "createTime.H"
    #include "createMesh.H"

    simpleControl pwrctrl(mesh,"PowerControl");

    #include "createFields.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nCalculating spectral radius and corresponding"
        << " eigenfunction \n" << endl;

    while (pwrctrl.loop(runTime))
    {
        Info<< "Power iteration = " << runTime.timeName() << nl << endl;

        //- Store previous iteration eigenvalue
        psi.storePrevIter();
        psiAdj.storePrevIter();
        dimensionedScalar lambdaOld(lambda);
        dimensionedScalar lambdaAdjOld(lambdaAdj);

        while (pwrctrl.correctNonOrthogonal())
        {
            //- Direct eigenvalue problem
            {
                //- Solve eigenproblem with previous eigenfunction to 
                //  help convergence

                fvScalarMatrix psiEqn
                (
                    fvm::div(phi,psi)
                 ==
                    fvm::laplacian(Gamma,psi)
                );

                psiEqn.relax();

                solve(psiEqn == lambda*psi.prevIter());

                //- Update the eigenvalue only at the last iteration
                if(pwrctrl.finalNonOrthogonalIter())
                {

                    //- Normalize psi
                    psi *= vol/fvc::domainIntegrate(psi);
                    
                    //- Use Rayleight quotient
                    volScalarField Lpsi(psiEqn.A()*psi - psiEqn.H());

                    lambda = 
                    (
                        fvc::domainIntegrate(Lpsi*psi)
                        /
                        fvc::domainIntegrate(psi*psi)
                    );

                    Info<<"eigenvalue = " << lambda.value() << endl; 

                    //- Note that this corresponds to an additional "step"
                    //  in the power method. 

                }

            }

            //- Adjoint eigenvalue problem
            {
                //- Solve eigenproblem with previous eigenfunction to 
                //  help convergence

                fvScalarMatrix psiAdjEqn
                (
                    fvm::div(phiAdj,psiAdj)
                 ==
                    fvm::laplacian(Gamma,psiAdj)
                );

                psiAdjEqn.relax();

                solve(psiAdjEqn == lambdaAdj*psiAdj.prevIter());

                //- Update the eigenvalue only at the last iteration
                if(pwrctrl.finalNonOrthogonalIter())
                {

                    //- Normalize psi
                    psiAdj *= vol/fvc::domainIntegrate(psiAdj);
                    
                    //- Use Rayleight quotient
                    volScalarField LpsiAdj
                    (
                        psiAdjEqn.A()*psiAdj - psiAdjEqn.H()
                    );

                    lambdaAdj = 
                    (
                        fvc::domainIntegrate(LpsiAdj*psiAdj)
                        /
                        fvc::domainIntegrate(psiAdj*psiAdj)
                    );

                    Info<<"adjoint eigenvalue = " << lambdaAdj.value() << endl; 

                    //- Note that this corresponds to an additional "step"
                    //  in the power method. 

                }

            }            
        }   

        //- Check convergence
        scalar lambdaRes
        (
            sameEigen
            ?
                mag
                (
                    (lambda.value() - lambdaAdj.value())
                    /
                    lambda.value()
                )         
            :
                max
                (
                    mag
                    (
                        (lambdaOld.value() - lambda.value())
                        /
                        lambda.value()
                    ),

                    mag
                    (
                        (lambdaAdjOld.value() - lambdaAdj.value())
                        /
                        lambdaAdj.value()
                    )
                )
        );


        Info<<"Residual on eigenvalue = " << lambdaRes << endl;

        if(lambdaRes < eigenTol)
        {
            Info<<"Power iterations converged!" <<endl;
            runTime.writeAndEnd();
        }

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
