/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "mappedDirichletFvPatchField.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::mappedDirichletFvPatchField<Type>::mappedDirichletFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    fixedValueFvPatchField<Type>(p, iF),
    newMappedPatchFieldBase<Type>(this->mapper(p, iF), *this)
{}


template<class Type>
Foam::mappedDirichletFvPatchField<Type>::mappedDirichletFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<Type>(p, iF, dict),
    newMappedPatchFieldBase<Type>(this->mapper(p, iF), *this, dict)
{}


template<class Type>
Foam::mappedDirichletFvPatchField<Type>::mappedDirichletFvPatchField
(
    const mappedDirichletFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<Type>(ptf, p, iF, mapper),
    newMappedPatchFieldBase<Type>(this->mapper(p, iF), *this, ptf)
{}


template<class Type>
Foam::mappedDirichletFvPatchField<Type>::mappedDirichletFvPatchField
(
    const mappedDirichletFvPatchField<Type>& ptf
)
:
    fixedValueFvPatchField<Type>(ptf),
    newMappedPatchFieldBase<Type>(ptf)
{}


template<class Type>
Foam::mappedDirichletFvPatchField<Type>::mappedDirichletFvPatchField
(
    const mappedDirichletFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    fixedValueFvPatchField<Type>(ptf, iF),
    newMappedPatchFieldBase<Type>(this->mapper(this->patch(), iF), *this, ptf)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
const Foam::mappedPatchBase& Foam::mappedDirichletFvPatchField<Type>::mapper
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
{
    if (!isA<mappedPatchBase>(p.patch()))
    {
        FatalErrorInFunction
            << "' not type '" << mappedPatchBase::typeName << "'"
            << "\n    for patch " << p.patch().name()
            << " of field " << iF.name()
            << " in file " << iF.objectPath()
            << exit(FatalError);
    }
    return refCast<const mappedPatchBase>(p.patch());
}


template<class Type>
void Foam::mappedDirichletFvPatchField<Type>::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    this->operator==(this->newMappedField());

    //if (debug)
    {
        Info<< "mapped on field:"
            << this->internalField().name()
            << " patch:" << this->patch().name()
            << " grad:" << gAverage(this->newMappedGrad())
            << "  avg:" << gAverage(*this)
            << "  min:" << gMin(*this)
            << "  max:" << gMax(*this)
            << endl;
    }

    fixedValueFvPatchField<Type>::updateCoeffs();
}


template<class Type>
void Foam::mappedDirichletFvPatchField<Type>::write(Ostream& os) const
{
    fvPatchField<Type>::write(os);
    newMappedPatchFieldBase<Type>::write(os);
    writeEntry(os, "value", *this);
}


// ************************************************************************* //
