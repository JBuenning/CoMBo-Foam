#include "meltFlowVelocityFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "fieldMapper.H"
#include "vectorField.H"
#include "volFields.H"
#include "pointFields.H"
#include "PrimitivePatchInterpolation.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::meltFlowVelocityFvPatchVectorField::
meltFlowVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF),
    rhoRatio_(dict.lookup<scalar>("rhoRatio")),
    URef_(dict.lookupOrDefault<vector>("URef", Zero))
{
    // is problematic for paraview
    // fixedValueFvPatchVectorField::evaluate();

    fvPatchVectorField::operator=
    (
        vectorField("value", iF.dimensions(), dict, p.size())
    );
}


Foam::meltFlowVelocityFvPatchVectorField::
meltFlowVelocityFvPatchVectorField
(
    const meltFlowVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper),
    rhoRatio_(ptf.rhoRatio_),
    URef_(ptf.URef_)
{}


Foam::meltFlowVelocityFvPatchVectorField::
meltFlowVelocityFvPatchVectorField
(
    const meltFlowVelocityFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(ptf, iF),
    rhoRatio_(ptf.rhoRatio_),
    URef_(ptf.URef_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::meltFlowVelocityFvPatchVectorField::map
(
    const fvPatchVectorField& ptf,
    const fieldMapper& mapper
)
{
    fixedValueFvPatchVectorField::map(ptf, mapper);
}


void Foam::meltFlowVelocityFvPatchVectorField::reset
(
    const fvPatchVectorField& ptf
)
{
    fixedValueFvPatchVectorField::reset(ptf);
}


void Foam::meltFlowVelocityFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const polyPatch& pp = patch().patch();
    const PrimitivePatchInterpolation<polyPatch> linearInterp(pp);

    const pointVectorField& pointMotionU = db().lookupObject<pointVectorField>("pointMotionU");
    const pointPatchField<vector>& pointPatchVelociy = pointMotionU.boundaryField()[pp.index()];

    const vectorField faceVelocity
    (
        linearInterp.pointToFaceInterpolate(pointPatchVelociy.patchInternalField())
    );

    fixedValueFvPatchVectorField::operator==
    (
        (1 - rhoRatio_) * ((faceVelocity - URef_) & pp.faceNormals()) * pp.faceNormals()
        +
        URef_
    );
    fixedValueFvPatchVectorField::updateCoeffs();
}


void Foam::meltFlowVelocityFvPatchVectorField::write
(
    Ostream& os
) const
{
    fvPatchVectorField::write(os);
    writeEntry(os, "rhoRatio", rhoRatio_);
    writeEntry(os, "URef", URef_);
    writeEntry(os, "value", *this);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        meltFlowVelocityFvPatchVectorField
    );
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
