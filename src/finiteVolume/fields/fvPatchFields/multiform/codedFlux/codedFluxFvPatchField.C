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

#include "codedFluxFvPatchField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "dynamicCode.H"
#include "dynamicCodeContext.H"
#include "stringOps.H"

// * * * * * * * * * * * * Private Static Data Members * * * * * * * * * * * //

template<class Type>
const Foam::wordList Foam::codedFluxFvPatchField<Type>::codeKeys_ =
    {"code", "codeInclude", "localCode"};


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<class Type>
const Foam::word Foam::codedFluxFvPatchField<Type>::codeTemplateC =
    "fluxFvPatchFieldTemplate.C";

template<class Type>
const Foam::word Foam::codedFluxFvPatchField<Type>::codeTemplateH =
    "fluxFvPatchFieldTemplate.H";


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

template<class Type>
void Foam::codedFluxFvPatchField<Type>::setFieldTemplates
(
    dynamicCode& dynCode
)
{
    word fieldType(pTraits<Type>::typeName);

    // template type for fvPatchField
    dynCode.setFilterVariable("TemplateType", fieldType);

    // Name for fvPatchField - eg, ScalarField, VectorField, ...
    fieldType[0] = toupper(fieldType[0]);
    dynCode.setFilterVariable("FieldType", fieldType + "Field");
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
const Foam::IOdictionary& Foam::codedFluxFvPatchField<Type>::dict() const
{
    const objectRegistry& obr = this->db();

    if (obr.foundObject<IOdictionary>("codeDict"))
    {
        return obr.lookupObject<IOdictionary>("codeDict");
    }
    else
    {
        return obr.store
        (
            new IOdictionary
            (
                IOobject
                (
                    "codeDict",
                    this->db().time().system(),
                    this->db(),
                    IOobject::MUST_READ_IF_MODIFIED,
                    IOobject::NO_WRITE
                )
            )
        );
    }
}


template<class Type>
Foam::dlLibraryTable& Foam::codedFluxFvPatchField<Type>::libs() const
{
    return const_cast<dlLibraryTable&>(this->db().time().libs());
}


template<class Type>
void Foam::codedFluxFvPatchField<Type>::prepare
(
    dynamicCode& dynCode,
    const dynamicCodeContext& context
) const
{
    // take no chances - typeName must be identical to name_
    dynCode.setFilterVariable("typeName", name_);

    // set TemplateType and FieldType filter variables
    // (for fvPatchField)
    setFieldTemplates(dynCode);

    // compile filtered C template
    dynCode.addCompileFile(codeTemplateC);

    // copy filtered H template
    dynCode.addCopyFile(codeTemplateH);


    // debugging: make BC verbose
    //  dynCode.setFilterVariable("verbose", "true");
    //  Info<<"compile " << name_ << " sha1: "
    //      << context.sha1() << endl;

    // define Make/options
    dynCode.setMakeOptions
        (
            "EXE_INC = -g \\\n"
            "-I$(LIB_SRC)/finiteVolume/lnInclude \\\n"
            "-I$(ELECTROFOAM)/src/lnInclude \\\n"
            + context.options()
            + "\n\nLIB_LIBS = \\\n"
            + "    -lOpenFOAM \\\n"
            + "    -lfiniteVolume \\\n"
            + "     -L$(FOAM_USER_LIBBIN) \\\n"
            + "     -lelectrokineticFoam \\\n"
            + context.libs()
        );
}


template<class Type>
const Foam::dictionary& Foam::codedFluxFvPatchField<Type>::codeDict()
const
{
    // use system/codeDict or in-line
    return
    (
        dict_.found("code")
      ? dict_
      : this->dict().subDict(name_)
    );
}


template<class Type>
const Foam::wordList& Foam::codedFluxFvPatchField<Type>::codeKeys() const
{
    return codeKeys_;
}


template<class Type>
Foam::string Foam::codedFluxFvPatchField<Type>::description() const
{
    return
        "patch "
      + this->patch().name()
      + " on field "
      + this->internalField().name();
}


template<class Type>
void Foam::codedFluxFvPatchField<Type>::clearRedirect() const
{
    // remove instantiation of fvPatchField provided by library
    redirectPatchFieldPtr_.clear();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::codedFluxFvPatchField<Type>::codedFluxFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    fluxFvPatchField<Type>(p, iF),
    codedBase(),
    redirectPatchFieldPtr_()
{}


template<class Type>
Foam::codedFluxFvPatchField<Type>::codedFluxFvPatchField
(
    const codedFluxFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fluxFvPatchField<Type>(ptf, p, iF, mapper),
    codedBase(),
    dict_(ptf.dict_),
    name_(ptf.name_),
    redirectPatchFieldPtr_()
{}


template<class Type>
Foam::codedFluxFvPatchField<Type>::codedFluxFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    fluxFvPatchField<Type>(p, iF, dict),
    codedBase(),
    dict_(dict),
    name_
    (
        dict.found("redirectType")
      ? dict.lookup("redirectType")
      : dict.lookup("name")
    ),
    redirectPatchFieldPtr_()
{
    updateLibrary(name_);
}


template<class Type>
Foam::codedFluxFvPatchField<Type>::codedFluxFvPatchField
(
    const codedFluxFvPatchField<Type>& ptf
)
:
    fluxFvPatchField<Type>(ptf),
    codedBase(),
    dict_(ptf.dict_),
    name_(ptf.name_),
    redirectPatchFieldPtr_()
{}


template<class Type>
Foam::codedFluxFvPatchField<Type>::codedFluxFvPatchField
(
    const codedFluxFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    fluxFvPatchField<Type>(ptf, iF),
    codedBase(),
    dict_(ptf.dict_),
    name_(ptf.name_),
    redirectPatchFieldPtr_()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
const Foam::fluxFvPatchField<Type>&
Foam::codedFluxFvPatchField<Type>::redirectPatchField() const
{
    if (!redirectPatchFieldPtr_.valid())
    {
        // Construct a patch
        // Make sure to construct the patchfield with up-to-date value

        // Write the data from the flux b.c.
        OStringStream os;
        fluxFvPatchField<Type>::write(os);
        IStringStream is(os.str());
        // Construct dictionary from it.
        dictionary dict(is);

        // Override the type to enforce the fvPatchField::New constructor
        // to choose our type
        dict.set("type", name_);

        redirectPatchFieldPtr_.set
        (
            dynamic_cast<fluxFvPatchField<Type>*>
            (
                fvPatchField<Type>::New
                (
                    this->patch(),
                    this->internalField(),
                    dict
                ).ptr()
            )
        );
    }
    return redirectPatchFieldPtr_();
}


template<class Type>
void Foam::codedFluxFvPatchField<Type>::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    // Make sure library containing user-defined fvPatchField is up-to-date
    updateLibrary(name_);

    const fluxFvPatchField<Type>& fvp = redirectPatchField();

    const_cast<fluxFvPatchField<Type>&>(fvp).updateCoeffs();

    // Copy through coefficients
    this->fluxInternalCoeffs() = fvp.fluxInternalCoeffs();
    this->fluxBoundaryCoeffs() = fvp.fluxBoundaryCoeffs();
    this->snGradRef() = fvp.snGradRef();

    // Copy through value
    this->operator==(fvp);

    fluxFvPatchField<Type>::updateCoeffs();
}


template<class Type>
void Foam::codedFluxFvPatchField<Type>::evaluate
(
    const Pstream::commsTypes commsType
)
{
    // Make sure library containing user-defined fvPatchField is up-to-date
    updateLibrary(name_);

    const fluxFvPatchField<Type>& fvp = redirectPatchField();

    // - updates the value of fvp (though not used)
    // - resets the updated() flag
    const_cast<fluxFvPatchField<Type>&>(fvp).evaluate(commsType);

    // Update the value (using the coefficients) locally
    fluxFvPatchField<Type>::evaluate(commsType);
}


template<class Type>
void Foam::codedFluxFvPatchField<Type>::write(Ostream& os) const
{
    fluxFvPatchField<Type>::write(os);
    os.writeKeyword("name") << name_
        << token::END_STATEMENT << nl;

    if (dict_.found("codeInclude"))
    {
        os.writeKeyword("codeInclude")
            << token::HASH << token::BEGIN_BLOCK;

        os.writeQuoted(string(dict_["codeInclude"]), false)
            << token::HASH << token::END_BLOCK
            << token::END_STATEMENT << nl;
    }

    if (dict_.found("localCode"))
    {
        os.writeKeyword("localCode")
            << token::HASH << token::BEGIN_BLOCK;

        os.writeQuoted(string(dict_["localCode"]), false)
            << token::HASH << token::END_BLOCK
            << token::END_STATEMENT << nl;
    }

    if (dict_.found("code"))
    {
        os.writeKeyword("code")
            << token::HASH << token::BEGIN_BLOCK;

        os.writeQuoted(string(dict_["code"]), false)
            << token::HASH << token::END_BLOCK
            << token::END_STATEMENT << nl;
    }

    if (dict_.found("codeOptions"))
    {
        os.writeKeyword("codeOptions")
            << token::HASH << token::BEGIN_BLOCK;

        os.writeQuoted(string(dict_["codeOptions"]), false)
            << token::HASH << token::END_BLOCK
            << token::END_STATEMENT << nl;
    }

    if (dict_.found("codeLibs"))
    {
        os.writeKeyword("codeLibs")
            << token::HASH << token::BEGIN_BLOCK;

        os.writeQuoted(string(dict_["codeLibs"]), false)
            << token::HASH << token::END_BLOCK
            << token::END_STATEMENT << nl;
    }
}


// ************************************************************************* //
