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


\*---------------------------------------------------------------------------*/

#include "mappedContinuitySolidFvPatchScalarField.H"
#include "fvPatchFieldMapper.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::mappedContinuitySolidFvPatchScalarField::mappedContinuitySolidFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    RobinFvPatchScalarField(p, iF),
    newMappedPatchFieldBase<scalar>(this->mapper(p, iF), *this),
    phiName_("phi"),
    RobinKeff_(p.size()),
    RobinFeff_(p.size()),
    Dfluid_(p.size())
{}


Foam::mappedContinuitySolidFvPatchScalarField::mappedContinuitySolidFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    RobinFvPatchScalarField(p, iF, dict),
    newMappedPatchFieldBase<scalar>(this->mapper(p, iF), *this, dict),
    phiName_(dict.lookupOrDefault<word>("phi", "phi")),
    RobinKeff_(p.size()),
    RobinFeff_(p.size()),
    Dfluid_("Dfluid", dict, p.size())
{
    //updateCoeffs();
}


Foam::mappedContinuitySolidFvPatchScalarField::mappedContinuitySolidFvPatchScalarField
(
    const mappedContinuitySolidFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    RobinFvPatchScalarField(p, iF),
    newMappedPatchFieldBase<scalar>(this->mapper(p, iF), *this, ptf),
    phiName_(ptf.phiName_),
    RobinKeff_(mapper(ptf.RobinKeff_)),
    RobinFeff_(mapper(ptf.RobinFeff_)),
    Dfluid_(mapper(ptf.Dfluid_))
{}


Foam::mappedContinuitySolidFvPatchScalarField::mappedContinuitySolidFvPatchScalarField
(
    const mappedContinuitySolidFvPatchScalarField& ptf
)
:
    RobinFvPatchScalarField(ptf),
    newMappedPatchFieldBase<scalar>(ptf),
    phiName_(ptf.phiName_),
    RobinKeff_(ptf.RobinKeff_),
    RobinFeff_(ptf.RobinFeff_),
    Dfluid_(ptf.Dfluid_)
{}


Foam::mappedContinuitySolidFvPatchScalarField::mappedContinuitySolidFvPatchScalarField
(
    const mappedContinuitySolidFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    RobinFvPatchScalarField(ptf, iF),
    newMappedPatchFieldBase<scalar>(this->mapper(this->patch(), iF), *this, ptf),
    phiName_(ptf.phiName_),
    RobinKeff_(ptf.RobinKeff_),
    RobinFeff_(ptf.RobinFeff_),
    Dfluid_(ptf.Dfluid_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::mappedPatchBase& Foam::mappedContinuitySolidFvPatchScalarField::mapper
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
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

// void Foam::mappedContinuitySolidFvPatchScalarField::autoMap
// (
//     const fvPatchFieldMapper&   m
// )
// {
//     RobinFvPatchScalarField::autoMap(m);
// //    m(phiName_,phiName_);
// //    m(RobinKeff_,RobinKeff_);
// }
//
// void Foam::mappedContinuitySolidFvPatchScalarField::rmap
// (
//     const fvPatchField<scalar>& ptf,
//     const labelList& addr
// )
// {
//     RobinFvPatchScalarField::rmap(ptf,addr);
//
//     const mappedContinuitySolidFvPatchScalarField& mptf =
//         refCast<const mappedContinuitySolidFvPatchScalarField>(ptf);
//
// //    phiName_.rmap(mptf.phiName_,addr);
// //    RobinKeff_.rmap(mptf.RobinKeff_,addr);
// }

void Foam::mappedContinuitySolidFvPatchScalarField::write(Ostream& os) const
{
    RobinFvPatchScalarField::write(os);
    newMappedPatchFieldBase::write(os);
    writeEntry(os, "phi", phiName_);
//    writeEntry(os, "RobinKeff", RobinKeff_);
//    writeEntry(os, "RobinFeff", RobinFeff_);
    writeEntry(os, "Dfluid", Dfluid_);
}

// void Foam::mappedContinuitySolidFvPatchScalarField::evaluate
// (
//     const Pstream::commsTypes commsType
// )
// {
//
//     const scalarField& RobinK = RobinFvPatchScalarField::RobinK();
//     const scalarField& RobinF = RobinFvPatchScalarField::RobinF();
//
//     //- Calculate effective Robin coefficient
//     RobinKeff_ =
//                   - this->newMappedDelta()*Dfluid_
//                   + RobinK;
//     RobinFeff_ = this->newMappedDelta()*this->newMappedInternalField()*Dfluid_
//                   + RobinF;
//
//     //- Evaluate Robin boundary condition
//     RobinFvPatchScalarField::evaluate();
//
// }

void Foam::mappedContinuitySolidFvPatchScalarField::updateCoeffs()
{
//    Info << "update " << this->updated() <<  endl;

    if (this->updated())
    {
        return;
    }

    const scalarField& RobinK = RobinFvPatchScalarField::RobinK();
    const scalarField& RobinF = RobinFvPatchScalarField::RobinF();

    //- Calculate effective Robin coefficient
    RobinKeff_ =  mappedSurField<scalar>(phiName_)/patch().magSf()
                  - this->newMappedDelta()*Dfluid_
                  + RobinK;
    RobinFeff_ = this->newMappedDelta()*this->newMappedInternalField()*Dfluid_
                  + RobinF;

    //- Evaluate Robin boundary condition
    RobinFvPatchScalarField::updateCoeffs();

}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        mappedContinuitySolidFvPatchScalarField
    );
}

// ************************************************************************* //
