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
    constrainedLaplacianFoam

Description
    Solves a simple Laplace equation,
    with a Lagrange multiplier to set the average

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

    simpleControl simple(mesh);

    #include "createFields.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nCalculating temperature distribution\n" << endl;

    while (simple.loop(runTime))
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        const dimensionedScalar
              lambda(
                      (fvc::domainIntegrate(T)/vol)
                    );
        Info<< "Average = " << lambda << endl;

        while (simple.correctNonOrthogonal())
        {

            fvScalarMatrix TEqn
            (
                - fvm::laplacian(DT, T)
                // + fvm::SuSp(-coeff,T)
             ==
                source
                - fvm::Sp(rate*lambda,T)
                + rate*avgT*T
                // - fvm::Sp(rate*posPart(lambda),T)
                // + fvm::Sp(rate*negPart(avgT),T)
                // + rate*posPart(avgT)*T
                // - rate*negPart(lambda)*T
            );

            // -- THIS CAN BE ADDED TO fvMatrix
            //  TEqn.constrain(weight);
            //  which computes lambda as
            //  fvc::domainIntegrate(weight*T)/vol
            //  and adds these to the matrix
            //   - fvm::Sp(rate*weight*lambda,T) + rate*weight*avgT*T
            //  rate can be estimated based on the diagonal of the matrix


            TEqn.relax();
            TEqn.solve();
            T.relax();

            Info<< "Residuals = " <<  gMax((fvc::laplacian(DT,T)+rate*(avgT-lambda)*T+source)->primitiveField()) << endl;
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
