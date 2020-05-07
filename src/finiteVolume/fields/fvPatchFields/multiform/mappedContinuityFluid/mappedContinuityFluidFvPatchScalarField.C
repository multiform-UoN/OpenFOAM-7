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

#include "mappedContinuityFluidFvPatchScalarField.H"
#include "fvPatchFieldMapper.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::mappedContinuityFluidFvPatchScalarField::mappedContinuityFluidFvPatchScalarField
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
    Dsolid_(p.size())
{}


Foam::mappedContinuityFluidFvPatchScalarField::mappedContinuityFluidFvPatchScalarField
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
    Dsolid_(p.size())
{
  if (dict.found("Dsolid"))
  {
    Dsolid_ = scalarField("Dsolid",dict,p.size());
  }
  else
  {
    Dsolid_ = scalarField(p.size(),1.0);
  }

  if (dict.found("Dfluid"))
  {
    RobinD(scalarField("Dfluid",dict,p.size()));
  }

    //updateCoeffs();
}


Foam::mappedContinuityFluidFvPatchScalarField::mappedContinuityFluidFvPatchScalarField
(
    const mappedContinuityFluidFvPatchScalarField& ptf,
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
    Dsolid_(mapper(ptf.Dsolid_))
{}


Foam::mappedContinuityFluidFvPatchScalarField::mappedContinuityFluidFvPatchScalarField
(
    const mappedContinuityFluidFvPatchScalarField& ptf
)
:
    RobinFvPatchScalarField(ptf),
    newMappedPatchFieldBase<scalar>(ptf),
    phiName_(ptf.phiName_),
    RobinKeff_(ptf.RobinKeff_),
    RobinFeff_(ptf.RobinFeff_),
    Dsolid_(ptf.Dsolid_)
{}


Foam::mappedContinuityFluidFvPatchScalarField::mappedContinuityFluidFvPatchScalarField
(
    const mappedContinuityFluidFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    RobinFvPatchScalarField(ptf, iF),
    newMappedPatchFieldBase<scalar>(this->mapper(this->patch(), iF), *this, ptf),
    phiName_(ptf.phiName_),
    RobinKeff_(ptf.RobinKeff_),
    RobinFeff_(ptf.RobinFeff_),
    Dsolid_(ptf.Dsolid_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::mappedPatchBase& Foam::mappedContinuityFluidFvPatchScalarField::mapper
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

// void Foam::mappedContinuityFluidFvPatchScalarField::autoMap
// (
//     const fvPatchFieldMapper&   m
// )
// {
//     RobinFvPatchScalarField::autoMap(m);
// //    m(phiName_,phiName_);
// //    m(RobinKeff_,RobinKeff_);
// }
//
// void Foam::mappedContinuityFluidFvPatchScalarField::rmap
// (
//     const fvPatchField<scalar>& ptf,
//     const labelList& addr
// )
// {
//     RobinFvPatchScalarField::rmap(ptf,addr);
//
//     const mappedContinuityFluidFvPatchScalarField& mptf =
//         refCast<const mappedContinuityFluidFvPatchScalarField>(ptf);
//
// //    phiName_.rmap(mptf.phiName_,addr);
// //    RobinKeff_.rmap(mptf.RobinKeff_,addr);
// }

void Foam::mappedContinuityFluidFvPatchScalarField::write(Ostream& os) const
{
    RobinFvPatchScalarField::write(os);
    newMappedPatchFieldBase::write(os);
    writeEntry(os, "phiName", phiName_);
//    writeEntry(os, "RobinKeff", RobinKeff_);
//    writeEntry(os, "RobinFeff", RobinFeff_);
    writeEntry(os, "Dsolid", Dsolid_);
}

// void Foam::mappedContinuityFluidFvPatchScalarField::evaluate
// (
//     const Pstream::commsTypes commsType
// )
// {
//
//     const scalarField& RobinK = RobinFvPatchScalarField::RobinK();
//     const scalarField& RobinF = RobinFvPatchScalarField::RobinF();
//
//     const fvsPatchField<scalar>& phip =
//         patch().lookupPatchField<surfaceScalarField, scalar>(phiName_);
//
//     //- Calculate effective Robin coefficient
//     RobinKeff_ = phip/patch().magSf()
//                   - this->newMappedDelta()*Dsolid_
//                   + RobinK;
//     RobinFeff_ = this->newMappedDelta()*this->newMappedInternalField()*Dsolid_
//                   + RobinF;
//
//     //- Evaluate Robin boundary condition
//     RobinFvPatchScalarField::evaluate();
//
// }


void Foam::mappedContinuityFluidFvPatchScalarField::updateCoeffs()
{
//    Info << "update " << this->updated() <<  endl;

    if (this->updated())
    {
        return;
    }

    const scalarField& RobinK = RobinFvPatchScalarField::RobinK();
    const scalarField& RobinF = RobinFvPatchScalarField::RobinF();

    const fvsPatchField<scalar>& phip =
        patch().lookupPatchField<surfaceScalarField, scalar>(phiName_);

    //- Calculate effective Robin coefficient
    RobinKeff_ = phip/patch().magSf()
                  - this->newMappedDelta()*Dsolid_
                  + RobinK;
    RobinFeff_ = this->newMappedDelta()*this->newMappedInternalField()*Dsolid_
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
        mappedContinuityFluidFvPatchScalarField
    );
}

// ************************************************************************* //
