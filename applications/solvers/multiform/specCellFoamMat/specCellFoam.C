/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015-2019
     \\/     M anipulation  | Matteo Icardi, Federico Municchi
-------------------------------------------------------------------------------
License
    This file is derivative work of OpenFOAM.

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
    specCellFoam

Description
    Compute leading order and first correction of the homogenised advection
    diffusion equation with reactive internal boundaries by solving the
    cell problem with spectral decomposition.

Developers
    Federico Municchi, Nottingham (2019)
    Matteo Icardi,  Nottingham (2019)
\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "wallDist.H"
#include "simpleControl.H"


int main(int argc, char *argv[]) {

    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createCellProblemControls.H"

    #include "createFields.H"

    // Since solver contains no time loop it would never execute
    // function objects so do it ourselves
    runTime.functionObjects().start();

    //- Power iterations (i.e., spectral problem)
    Info << "Power iterations for the spectral problem" << endl;
    do
    {

        // Manually advance the time index
        runTime++;

        #include "powerIterSettings.H"

        Info<<"\nPower iteration " << currPowerIter << endl;

        // DIRECT
        while(pwrctrl.correctNonOrthogonal())
        {
            #include    "psiEqn.H"
        }

        #include    "updateEigenvalue.H"
        #include    "normalisePsi.H"

        // ADJOINT
        while(pwrctrl.correctNonOrthogonal())
        {

            #include    "psiAdjEqn.H"

        }

        #include    "updateEigenvalueAdj.H"
        #include    "normalisePsiAdj.H"

        #include    "aitkenRelaxEigenfunctions.H"
        #include    "powerConvergence.H"

        runTime.write();

    } while (
        !powerConverged
    );

    //- Calculate modified velocity and beta
    #include    "calculateBeta.H"
    #include    "calculateUplus.H"


    //- Solve corrector problem
    Info << "Solving for scalar transport corrector field" << endl;
    do {

        // Manually advance the time index
        runTime++;

        #include    "cellIterSettings.H"

        while(pwrctrl.correctNonOrthogonal())
        {
            #include    "XEqn.H"
        }

        #include "rescaleX.H"
        #include "aitkenRelaxation.H"
        #include "cellConvergence.H"

        runTime.write();

    } while (
        !cellConverged
    );

    #include    "rescaleX.H"
    #include    "calculateEffectiveParam.H"

    runTime.functionObjects().end();
    runTime.writeAndEnd();

    #include    "writeEffectiveParam.H"

    Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
        << "  ClockTime = " << runTime.elapsedClockTime() << " s"
        << nl << endl;

    Info<< "End" << endl;

    return 0;
}
