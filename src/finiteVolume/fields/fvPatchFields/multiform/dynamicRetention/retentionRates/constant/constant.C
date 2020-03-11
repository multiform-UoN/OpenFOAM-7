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
#include "constant.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{
namespace retentionRates
{
    defineTypeNameAndDebug(constant, 0);

    addToRunTimeSelectionTable
    (
        retentionRate,
        constant,
        dictionary
    );
}
}

using namespace Foam;

/*------------------------  Constructors  ------------------------------------*/
retentionRates::constant::constant
(
    const dictionary& dict
)
:
    retentionRate(dict)
{}
/*------------------------  Destructors  -------------------------------------*/
retentionRates::constant::~constant()
{}
/*------------------------  Member functions ---------------------------------*/
tmp<scalarField>
retentionRates::constant::value(const scalarField& S)
{
    tmp<scalarField> tf (

        new scalarField(S.size(), scalar(1))

    );

    return tf;
}

tmp<scalarField>
retentionRates::constant::ddS(const scalarField& S)
{
    tmp<scalarField> tf (

        new scalarField(S.size(),scalar(0))

    );

    return tf;
}
