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

#include "dynamicRetentionFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"

#define MAX_NEWTON_ITER 1000
#define NEWTON_TOLERANCE 1e-13
#define STEADY_TIMESTEP 1e10

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dynamicRetentionFvPatchScalarField::dynamicRetentionFvPatchScalarField
(
  const fvPatch& p,
  const DimensionedField<scalar, volMesh>& iF
)
:
RobinPhiFvPatchScalarField(p, iF),
S_(p.size()),
S0_(p.size()),
RobinKeff_(p.size()),
Kd_(p.size()),
// shearDetachment_(false),
RobinFeff_(p.size()),
timeIndex_(-1),
rR_(NULL),
dR_(NULL)
{

}


Foam::dynamicRetentionFvPatchScalarField::dynamicRetentionFvPatchScalarField
(
  const fvPatch& p,
  const DimensionedField<scalar, volMesh>& iF,
  const dictionary& dict
)
:
RobinPhiFvPatchScalarField(p, iF, dict),
S_("S",dict,p.size()),
S0_(S_),
RobinKeff_(p.size(),scalar(0)),
Kd_("Kd",dict,p.size()),
// shearDetachment_(dict.lookupOrDefault<bool>("shearDetachment",false)),
RobinFeff_(p.size(),scalar(0)),
timeIndex_(iF.mesh().time().timeIndex()),
rR_
(
  retentionRates::retentionRate::New
  (
    dict.subDict("forwardReaction")
  )
),
dR_
(
  retentionRates::retentionRate::New
  (
    dict.subDict("backwardReaction")
  )
)
{
  // if(shearDetachment_)
  // {
  //   sdR_.set
  //   (
  //     retentionRates::retentionRate::New
  //     (
  //       dict.subDict("shearDetachmentRate")
  //     ).ptr()
  //   );
  //
  // }
}


Foam::dynamicRetentionFvPatchScalarField::dynamicRetentionFvPatchScalarField
(
  const dynamicRetentionFvPatchScalarField& ptf,
  const fvPatch& p,
  const DimensionedField<scalar, volMesh>& iF,
  const fvPatchFieldMapper& mapper
)
:
RobinPhiFvPatchScalarField(p, iF),
S_(mapper(ptf.S_)),//,mapper),
S0_(mapper(ptf.S0_)),//,mapper),
RobinKeff_(mapper(ptf.RobinKeff_)),//,mapper),
Kd_(mapper(ptf.Kd_)),//,mapper),
RobinFeff_(mapper(ptf.RobinFeff_)),//,mapper),
timeIndex_(ptf.timeIndex_)
{

}


Foam::dynamicRetentionFvPatchScalarField::dynamicRetentionFvPatchScalarField
(
  const dynamicRetentionFvPatchScalarField& ptf
)
:
RobinPhiFvPatchScalarField(ptf),
S_(ptf.S_),
S0_(ptf.S0_),
RobinKeff_(ptf.RobinKeff_),
Kd_(ptf.Kd_),
// shearDetachment_(ptf.shearDetachment_),
RobinFeff_(ptf.RobinFeff_),
timeIndex_(ptf.timeIndex_)
{

}


Foam::dynamicRetentionFvPatchScalarField::dynamicRetentionFvPatchScalarField
(
  const dynamicRetentionFvPatchScalarField& ptf,
  const DimensionedField<scalar, volMesh>& iF
)
:
RobinPhiFvPatchScalarField(ptf, iF),
S_(ptf.S_),
S0_(ptf.S0_),
RobinKeff_(ptf.RobinKeff_),
Kd_(ptf.Kd_),
// shearDetachment_(ptf.shearDetachment_),
RobinFeff_(ptf.RobinFeff_),
timeIndex_(ptf.timeIndex_)
{

}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
void Foam::dynamicRetentionFvPatchScalarField::autoMap
(
  const fvPatchFieldMapper&   m
)
{
  RobinPhiFvPatchScalarField::autoMap(m);
  m(S_,S_);//S_.autoMap(mapper);
  m(S0_,S0_);//S0_.autoMap(mapper);
  m(RobinKeff_,RobinKeff_);//.autoMap(mapper);
  m(Kd_,Kd_);//.autoMap(mapper);
  m(RobinFeff_,RobinFeff_);//.autoMap(mapper);
}

void Foam::dynamicRetentionFvPatchScalarField::rmap
(
  const fvPatchField<scalar>& ptf,
  const labelList& addr
)
{
  RobinPhiFvPatchScalarField::rmap(ptf,addr);
  const dynamicRetentionFvPatchScalarField& mptf =
  refCast<const dynamicRetentionFvPatchScalarField>(ptf);
  S_.rmap(mptf.S_,addr);
  S0_.rmap(mptf.S0_,addr);
  RobinKeff_.rmap(mptf.RobinKeff_,addr);
  Kd_.rmap(mptf.Kd_,addr);
  RobinFeff_.rmap(mptf.RobinFeff_,addr);
  rR_ = mptf.rR_;
  dR_ = mptf.dR_;
}

void Foam::dynamicRetentionFvPatchScalarField::write(Ostream& os) const
{
  RobinPhiFvPatchScalarField::write(os);
  writeEntry(os, "Kd", Kd_);
  // writeEntry(os, "shearDetachment", shearDetachment_);
  writeEntry(os, "S", S_);
  writeEntry(os, "RobinKeff2", RobinKeff_);
  writeEntry(os, "RobinFeff", RobinFeff_);

  if(!rR_.valid() || !dR_.valid())
  {
    return;
  }
  // Write retention and detachment dicts
  os<<"\n        forwardReaction";
  rR_().write(os);
  os<<"\n        backwardReaction	";
  dR_().write(os);
  // if(shearDetachment_)
  // {
  //   os<<"\n        shearDetachmentRate	";
  //   sdR_().write(os);
  // }

}

void Foam::dynamicRetentionFvPatchScalarField::updateCoeffs()
// void Foam::dynamicRetentionFvPatchScalarField::evaluate
// (
//   const Pstream::commsTypes commsType
// )
{

  RobinPhiFvPatchScalarField::updateCoeffs();

  const fvMesh& mesh = this->internalField().mesh();

  word ddtScheme
  (
      mesh.ddtScheme("dynamicRetention")
  );

  scalar deltaT = this->db().time().deltaTValue();

  if (ddtScheme=="steadyState")
  {
    deltaT = STEADY_TIMESTEP;
  }

  scalarField& C(*this);
  const scalarField& RobinKorig = RobinFvPatchScalarField::RobinK();
  const scalarField& RobinK0 = RobinPhiFvPatchScalarField::RobinK();
  const scalarField& RobinF0 = RobinPhiFvPatchScalarField::RobinF();
  const scalar area(gSum(this->patch().magSf()));

  // //- Evaluate shear rate (always computed)
  // const volVectorField& U(this->db().lookupObject<volVectorField>("U"));
  // const scalarField srb
  // (
  //   mag
  //   (
  //     ( tensor::I - (this->patch().nf()*this->patch().nf()) )
  //     & (fvc::grad(U))->boundaryField()[this->patch().index()]
  //     & this->patch().nf()
  //   )
  // );

  //- Get old time concentration
  if(timeIndex_!=this->db().time().timeIndex())
  {
    S0_ = S_;
    timeIndex_ = this->db().time().timeIndex();

    // - Evolution equation for S (linearised source term)
    // - With Newton below, this is now just the first guess
    S_ =
    (
      - RobinKorig*C*
      (
        rR_().value(S_) -
        rR_().ddS(S0_)*S0_
      )
      + S0_/deltaT
    )
    /
    (
      scalar(1)/deltaT
      + RobinKorig*C*rR_().ddS(S0_)
      - (
        // shearDetachment_
        // ?
        // Kd_*(dR_().value((*this))) - sdR_().value(srb)
        // :
        Kd_*(dR_().value((*this)))
      )
    );

    //- Firt guess for Robin coefficients
    RobinKeff_ = RobinK0-RobinKorig
               + RobinKorig*rR_().value(S_)
               - Kd_*S_*dR_().ddS(C);

    RobinFeff_ =
    RobinF0 +
    - S_*
    (
      (
        // shearDetachment_
        // ?
        // Kd_*(dR_().value(C)) - sdR_().value(srb)
        // :
        Kd_*(dR_().value(C))
      )
      - Kd_*dR_().ddS(C)*C
    );

  }

  //- Evolution equation for S (Newton's method)
  label n = 0;
  scalarField Sn(S_);
  do
  {

    if (n>MAX_NEWTON_ITER)
    {
      FatalErrorInFunction
      << "Max number of Newton iterations reached for "
      << "dynamicRetentionFvPatchScalarField: "
      << this->patch().name()
      << exit(FatalError);
    }

    RobinFvPatchScalarField::evaluate();

    Sn = S_;
    S_ = Sn -
    (
      RobinKorig*C*rR_().value(Sn)
      - (
        // shearDetachment_
        // ?
        // Kd_*(dR_().value(C)) - sdR_().value(srb)
        // :
        Kd_*(dR_().value(C))
        )*Sn
      + (Sn-S0_)/deltaT
    )
    /
    (
      scalar(1)/deltaT
      + RobinKorig*C*rR_().ddS(Sn)
      - (
        // shearDetachment_
        // ?
        // Kd_*(dR_().value(C)) - sdR_().value(srb)
        // :
        Kd_*(dR_().value(C))
        )
    );
    n++;

    //- Calculate effective Robin coefficient
    RobinKeff_ = RobinK0-RobinKorig
               + RobinKorig*rR_().value(S_)
               - Kd_*S_*dR_().ddS(C);

    RobinFeff_ = RobinF0 +
                - S_ *
                (
                  (
                    // shearDetachment_
                    // ?
                    // Kd_*(dR_().value(C)) - sdR_().value(srb)
                    // :
                    Kd_*(dR_().value(C))
                  )
                  - Kd_*dR_().ddS(C)*C
                );


  }
  while ( (Foam::sqrt(gSum(pow(S_-Sn,2)))) > NEWTON_TOLERANCE );

  Info<<"  dynamicRetention " << this->patch().name()
  <<": converged in " << n << " iterations, S="
  << gSum(this->patch().magSf()*S_)/area
  << " Cp=" << gSum(this->patch().magSf()*(this->patchInternalField()))/area
  << " Cf=" << gSum(this->patch().magSf()*C)/gSum(this->patch().magSf())
  << " R=" << gSum(this->patch().magSf()*(RobinKeff_*C + RobinFeff_))/area
  << " dSdt=" << gSum(this->patch().magSf()*(S_-S0_)/deltaT)/area
  << " dCdn=" << gSum(this->patch().magSf()*this->snGrad())/area << endl;
  // << " Shear=" << gSum(srb*this->patch().magSf())/area << endl;


  //    RobinFvPatchScalarField::evaluate();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
  makePatchTypeField
  (
    fvPatchScalarField,
    dynamicRetentionFvPatchScalarField
  );
}

// ************************************************************************* //
