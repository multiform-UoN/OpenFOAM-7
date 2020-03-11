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
    Foam::retentionRates::retentionRate

Description
    Base abstract class for retention rates

Authors:
    Federico Municchi, Nottingham (2019)
\*---------------------------------------------------------------------------*/

#ifndef RETENTION_RATE_H
#define RETENTION_RATE_H

#include "dictionary.H"
#include "runTimeSelectionTables.H"
#include "fvCFD.H"

namespace Foam
{
namespace retentionRates
{
class retentionRate
{
private:

    //- Dictionary with settings
    dictionary  dict_;

public:

    //- Runtime type information
    TypeName("retentionRate");

    // Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        retentionRate,
        dictionary,
        (
            const dictionary& dict
        ),
        (dict)
    );

    //- Constructor
        retentionRate( const dictionary& dict );

    //- Destructor
        virtual ~retentionRate();

    //- Selector
        static autoPtr<retentionRate> New
        (
            const dictionary& dict
        );

    //- Member functions

        //- Calculate retention rate
        virtual tmp<scalarField> value(const scalarField& S) = 0;

        //- Return derivative w.r.t. S
        virtual tmp<scalarField> ddS(const scalarField& S) = 0;

        //- Return dictionary
        dictionary& dict() {return dict_;};
	
	//- Write 
	virtual void write(Ostream& o) const;

};
//--------------- End namespace Foam ------------------------------------------/
}
//----------------End namespace retentionRates --------------------------------/
}

#endif
