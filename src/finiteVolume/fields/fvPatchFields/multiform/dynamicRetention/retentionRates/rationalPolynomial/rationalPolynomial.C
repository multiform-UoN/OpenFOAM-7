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
#include "rationalPolynomial.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{
namespace retentionRates
{
    defineTypeNameAndDebug(rationalPolynomial, 0);

    addToRunTimeSelectionTable
    (
        retentionRate,
        rationalPolynomial,
        dictionary
    );
}
}

using namespace Foam;

/*------------------------  Constructors  ------------------------------------*/
retentionRates::rationalPolynomial::rationalPolynomial
(
    const dictionary& dict
)
:
    retentionRate(dict),
    n_(readScalar(dict.lookup("n"))),
    S0_(readScalar(dict.lookup("S0"))),
    Kref_(readScalar(dict.lookup("Kref")))
{}
/*------------------------  Destructors  -------------------------------------*/
retentionRates::rationalPolynomial::~rationalPolynomial()
{}
/*------------------------  Member functions ---------------------------------*/
tmp<scalarField>
retentionRates::rationalPolynomial::value(const scalarField& S)
{
    tmp<scalarField> tf (

        new scalarField(S.size())

    );

    scalarField& f = tf.ref();

    f = ( 1 + Foam::pow(S/S0_,n_) ) /
          ( 1 + Kref_* Foam::pow(S/S0_,n_));

    return tf;
}

tmp<scalarField>
retentionRates::rationalPolynomial::ddS(const scalarField& S)
{
    tmp<scalarField> tf (

        new scalarField(S.size())

    );

    scalarField& f = tf.ref();

    f = (n_/S0_) * Foam::pow(S/S0_,n_-scalar(1))
        *
        ( scalar(1) - Kref_*(this->value(S)) )
        /
        ( scalar(1) + Kref_*Foam::pow(S/S0_,n_) );
    
    return tf;
}
