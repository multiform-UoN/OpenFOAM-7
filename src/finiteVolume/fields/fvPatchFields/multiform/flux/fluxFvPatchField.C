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
#include "fluxFvPatchField.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::fluxFvPatchField<Type>::fluxFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    fvPatchField<Type>(p, iF),
    fluxICoeffs_(p.size()),
    fluxBCoeffs_(p.size()),
    snGrad_(p.size()),
    fluxCorrected_(true)
{
}


template<class Type>
Foam::fluxFvPatchField<Type>::fluxFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    fvPatchField<Type>(p, iF, dict, false),
    fluxICoeffs_("fluxICoeffs", dict, p.size()),
    fluxBCoeffs_("fluxBCoeffs", dict, p.size()),
    snGrad_("snGrad",dict,p.size()),
    fluxCorrected_(false)
{
    evaluate();
}


template<class Type>
Foam::fluxFvPatchField<Type>::fluxFvPatchField
(
    const fluxFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper,
    const bool mappingRequired
)
:
    fvPatchField<Type>(ptf, p, iF, mapper, mappingRequired),
    fluxICoeffs_(mapper(ptf.fluxICoeffs_)),
    fluxBCoeffs_(mapper(ptf.fluxBCoeffs_)),
    snGrad_(mapper(ptf.snGrad_)),
    fluxCorrected_(ptf.fluxCorrected_)
{
    if (mappingRequired && notNull(iF) && mapper.hasUnmapped())
    {
        WarningInFunction
            << "On field " << iF.name() << " patch " << p.name()
            << " patchField " << this->type()
            << " : mapper does not map all values." << nl
            << "    To avoid this warning fully specify the mapping in derived"
            << " patch fields." << endl;
    }
}


template<class Type>
Foam::fluxFvPatchField<Type>::fluxFvPatchField
(
    const fluxFvPatchField<Type>& ptf
)
:
    fvPatchField<Type>(ptf),
    fluxICoeffs_(ptf.fluxICoeffs_),
    fluxBCoeffs_(ptf.fluxBCoeffs_),
    snGrad_(ptf.snGrad_),
    fluxCorrected_(ptf.fluxCorrected_)
{}


template<class Type>
Foam::fluxFvPatchField<Type>::fluxFvPatchField
(
    const fluxFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    fvPatchField<Type>(ptf, iF),
    fluxICoeffs_(ptf.fluxICoeffs_),
    fluxBCoeffs_(ptf.fluxBCoeffs_),
    snGrad_(ptf.snGrad_),
    fluxCorrected_(ptf.fluxCorrected_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::fluxFvPatchField<Type>::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fvPatchField<Type>::autoMap(m);
    m(fluxICoeffs_,fluxICoeffs_);
    m(fluxBCoeffs_,fluxBCoeffs_);
    m(snGrad_,snGrad_);
}

template<class Type>
void Foam::fluxFvPatchField<Type>::rmap
(
    const fvPatchField<Type>& ptf,
    const labelList& addr
)
{
    fvPatchField<Type>::rmap(ptf, addr);

    const fluxFvPatchField<Type>& mptf =
        refCast<const fluxFvPatchField<Type>>(ptf);

    fluxICoeffs_.rmap(mptf.fluxICoeffs_, addr);
    fluxBCoeffs_.rmap(mptf.fluxBCoeffs_, addr);
    snGrad_.rmap(mptf.snGrad_, addr);
}

template<class Type>
void Foam::fluxFvPatchField<Type>::evaluate(const Pstream::commsTypes)
{
    if (!this->updated())
    {
        this->updateCoeffs();
    }

    fvPatchField<Type>::evaluate();
}

template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::fluxFvPatchField<Type>::snGrad() const
{
    return snGrad_;
}

template<class Type>
void Foam::fluxFvPatchField<Type>::write(Ostream& os) const
{
    fvPatchField<Type>::write(os);
    writeEntry(os,"fluxICoeffs",fluxICoeffs_);
    writeEntry(os,"fluxBCoeffs",fluxBCoeffs_);
    writeEntry(os,"snGrad",fluxBCoeffs_);
    writeEntry(os,"value",*this);
}



template<class Type>
void Foam::fluxFvPatchField<Type>::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    //- This still does not work. Users need to be careful with this BC.
    // if(
    //     (
    //         this->db().time().value()
    //         >
    //         (
    //             this->db().time().startTime().value()
    //             + this->db().time().deltaT().value()
    //         )
    //     )
    //     &&
    //     !fluxCorrected_
    // )
    // {
    //     FatalErrorInFunction
    //     << "constrainFluxes was not called for this patch and therefore\n"
    //     << "this BC has been used improperly."
    //     << abort(FatalError);
    // }

    fvPatchField<Type>::updateCoeffs();
}

// ************************************************************************* //
