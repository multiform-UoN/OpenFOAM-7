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

  Info << "Initial Energy calculation" << endl;

  #include "computeEnergy.H"
  #include "updateMu.H"


  while (runTime.loop())
  {
    Info<< endl << "Time = " << runTime.timeName() << endl;

    // update potential
    #include "updatePot.H"

    // -- DEBUG
    if (debugCH)
    {
      mu.write();
      alpha.write();
      runTime++;
    }
    // -- DEBUG

    mu.storePrevIter();
    alpha.storePrevIter();

    while (pimple.loop())
    {

      // Shall we update the non-linearity here?
      // #include "updatePot.H"

      #include "alphaEqn.H"

      mu.storePrevIter();
      #include "muEqn.H"

      #include "computeEnergy.H"

      // -- DEBUG
      if (debugCH)
      {
        mu.write();
        alpha.write();
        runTime++;
      }
      // -- DEBUG

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
