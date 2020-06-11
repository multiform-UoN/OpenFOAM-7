/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2020 OpenFOAM Foundation
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
    CahnHilliardFoam

Description
    Solves the Cahn-Hilliard equation with a fourth order free energy.
    This is a segregated implementation for the mixed formulation, where
    just one equation (for alpha) is solved and where the fourth order operator
    is splitted in a semi-implicit way.

Authors
    Federico Municchi, University of Nottingham (2020)
    Matteo Icardi, University of Nottingham (2020)

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "pimpleControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createFields.H"

    pimpleControl pimple(mesh);

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    while (runTime.run())
    {
        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        //- Create matrix for mu (positive definite)
        fvScalarMatrix muEqn
        (
            - fvm::laplacian(M,mu)
        );

        //- Get diagonal (positive)
        volScalarField Amu(muEqn.A());

        while(pimple.loop())
        {

            //- Get off-diagonal terms
            volScalarField Hmu(muEqn.H());

            //- Assemble matrix for alpha (positive definite)
            fvScalarMatrix alphaEqn
            (
                fvm::ddt(alpha)
              - Amu*
                (
                    fvm::laplacian(epsSq,alpha)
                )
            );

            alphaEqn.relax();

            //- Solve for alpha
            solve(alphaEqn == Hmu - Amu*fPrime);

            //- Relax alpha
            alpha.relax();
            alpha.correctBoundaryConditions();

            //- Update fPrime
            fPrime = fea*pow(alpha,3.) - feb*alpha;

            //- Update potential
            mu = -fvc::laplacian(epsSq,alpha) + fPrime;

            mu.relax();
            mu.correctBoundaryConditions();


        }

        #include "writeEnergyAndMass.H"

        runTime.write();

        Info<< nl << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
