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


Class
    Foam::retentionRates::shearLinear

Description
    Linear shear rate

Authors:
    Federico Municchi, Nottingham (2019)
\*---------------------------------------------------------------------------*/

#ifndef RETENTION_SHEARLINEAR_H
#define RETENTION_SHEARLINEAR_H

#include "retentionRate.H"

namespace Foam
{
namespace retentionRates
{
class shearLinear
:
    public retentionRate
{
private:

    //- Scaling coefficient
    scalar K_;

public:

    //- Runtime type information
    TypeName("shearLinear");

    //- Constructor
        shearLinear( const dictionary& dict);

    //- Destructor
        virtual ~shearLinear();

    //- Member functions

        //- Calculate retention rate
        virtual tmp<scalarField> value(const scalarField& S);

        //- Return derivative w.r.t. S
        virtual tmp<scalarField> ddS(const scalarField& S);



};
//--------------- End namespace Foam ------------------------------------------/
}
//----------------End namespace constants --------------------------------/
}

#endif
