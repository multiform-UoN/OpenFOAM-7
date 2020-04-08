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
CahnHiliardFoam

Description
Solves Cahn-Hilliard equation

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "fvOptions.H"
#include "pimpleControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCaseLists.H"

    #include "createTime.H"
    #include "createMesh.H"

    pimpleControl pimple(mesh);

    #include "createFields.H"


    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nCalculating Cahn-Hiliard solution\n" << endl;

    label muRefCell = 0;
    scalar muRefValue = 0.0;

    Info << "Initial Energy calculation" << endl;
    #include "computeEnergy.H"

    Info << "Start time loop" << endl;

    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        Info << "Start internal loop" << endl;

        #include "updatePot.H"
        #include "updateMu.H"

        // Newton's iteration
        volScalarField alpha0(alpha);
        while (pimple.loop())
        {

            volScalarField alphan(alpha);
            while (pimple.correctNonOrthogonal())
            {
                fvScalarMatrix alphaEqn
                (
                    - fvm::laplacian(epsilon*epsilon, alpha)
                    + fvm::Sp(pot_imp,alpha)
                    ==
                      fvc::laplacian(epsilon*epsilon*alphan)
                    - pot + mu
                );

                alphaEqn.relax();
                alphaEqn.solve();
            }

            alpha += alphan;
            alpha.relax();
            alpha.correctBoundaryConditions();
            #include "updatePot.H"

            while (pimple.correctNonOrthogonal())
            {
                fvScalarMatrix muEqn
                (
                    fvc::ddt(alpha) - fvm::laplacian(M, mu)
                    //==
                    //fvOptions(mu)
                    //M/scalar(2)*fvc::laplacian(mu.prevIter())
                );

                //fvOptions.constrain(muEqn);
                muEqn.setReference(muRefCell, muRefValue);
                muEqn.relax();
                muEqn.solve();
                //fvOptions.correct(mu);
            }

            mu.relax();
            //#include "computeEnergy.H"

        }

        #include "computeEnergy.H"

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
        << "  ClockTime = " << runTime.elapsedClockTime() << " s"
        << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
