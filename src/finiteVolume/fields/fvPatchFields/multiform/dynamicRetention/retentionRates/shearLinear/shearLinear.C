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


Authors:
    Federico Municchi, Nottingham (2019)
\*---------------------------------------------------------------------------*/
#include "shearLinear.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{
namespace retentionRates
{
    defineTypeNameAndDebug(shearLinear, 0);

    addToRunTimeSelectionTable
    (
        retentionRate,
        shearLinear,
        dictionary
    );
}
}

using namespace Foam;

/*------------------------  Constructors  ------------------------------------*/
retentionRates::shearLinear::shearLinear
(
    const dictionary& dict
)
:
    retentionRate(dict),
    K_(readScalar(dict.lookup("K")))
{}
/*------------------------  Destructors  -------------------------------------*/
retentionRates::shearLinear::~shearLinear()
{}
/*------------------------  Member functions ---------------------------------*/
tmp<scalarField>
retentionRates::shearLinear::value(const scalarField& S)
{
    tmp<scalarField> tf (

        new scalarField(S*scalar(K_))

    );

    return tf;
}

tmp<scalarField>
retentionRates::shearLinear::ddS(const scalarField& S)
{
    //- Zero, because here S is the shear rate
    tmp<scalarField> tf (

        new scalarField(S.size(),scalar(0))

    );

    return tf;
}
