#include "kinematicSurfaceVelocityPointPatchVectorField.H"
#include "surfaceVelocityTools.H"
#include "pointPatchFields.H"
#include "addToRunTimeSelectionTable.H"
#include "polyMesh.H"
#include "vectorField.H"
#include "volFields.H"


namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

kinematicSurfaceVelocityPointPatchVectorField::
kinematicSurfaceVelocityPointPatchVectorField
(
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF,
    const dictionary& dict
)
:
    fixedValuePointPatchField<vector>(p, iF, dict),
    surfaceVelocityTools_(p, dict),
    UName_(dict.lookupOrDefault<word>("UName", "U"))
{}


kinematicSurfaceVelocityPointPatchVectorField::
kinematicSurfaceVelocityPointPatchVectorField
(
    const kinematicSurfaceVelocityPointPatchVectorField& ptf,
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF,
    const pointPatchFieldMapper& mapper
)
:
    fixedValuePointPatchField<vector>(ptf, p, iF, mapper),
    surfaceVelocityTools_(ptf.surfaceVelocityTools_, p),
    UName_(ptf.UName_)
{}


kinematicSurfaceVelocityPointPatchVectorField::
kinematicSurfaceVelocityPointPatchVectorField
(
    const kinematicSurfaceVelocityPointPatchVectorField& ptf,
    const DimensionedField<vector, pointMesh>& iF
)
:
    fixedValuePointPatchField<vector>(ptf, iF),
    surfaceVelocityTools_(ptf.surfaceVelocityTools_),
    UName_(ptf.UName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void kinematicSurfaceVelocityPointPatchVectorField::map
(
    const pointPatchField<vector>& ptf,
    const pointPatchFieldMapper& mapper
)
{
    const kinematicSurfaceVelocityPointPatchVectorField& oVptf =
        refCast<const kinematicSurfaceVelocityPointPatchVectorField>(ptf);

    fixedValuePointPatchField<vector>::map(oVptf, mapper);
    surfaceVelocityTools_.map(oVptf.surfaceVelocityTools_, mapper);
}


void kinematicSurfaceVelocityPointPatchVectorField::reset
(
    const pointPatchField<vector>& ptf
)
{
    const kinematicSurfaceVelocityPointPatchVectorField& oVptf =
        refCast<const kinematicSurfaceVelocityPointPatchVectorField>(ptf);

    fixedValuePointPatchField<vector>::reset(oVptf);
    surfaceVelocityTools_.reset(oVptf.surfaceVelocityTools_);
}


void kinematicSurfaceVelocityPointPatchVectorField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    const volVectorField& UField = this->db().lookupObject<volVectorField>(UName_);
    const fvPatchVectorField& UPatchField = UField.boundaryField()[this->patch().index()];
    const polyPatch& pp = UPatchField.patch().patch();

    const vectorField pointVelocity =
        surfaceVelocityTools_.pointVelocityFromFaceVelocity(UPatchField & pp.faceNormals(), pp);

    Field<vector>::operator=
    (
        std::move(pointVelocity)
    );
    fixedValuePointPatchField<vector>::updateCoeffs();
}


void kinematicSurfaceVelocityPointPatchVectorField::write(Ostream& os) const
{
    pointPatchField<vector>::write(os);
    surfaceVelocityTools_.write(os);
    writeEntry(os, "UName", UName_);
    writeEntry(os, "value", *this);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePointPatchTypeField
(
    pointPatchVectorField,
    kinematicSurfaceVelocityPointPatchVectorField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam
