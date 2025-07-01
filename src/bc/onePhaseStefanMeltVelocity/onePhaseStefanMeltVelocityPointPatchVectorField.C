#include "onePhaseStefanMeltVelocityPointPatchVectorField.H"
#include "fixedValuePointPatchField.H"
#include "scalarField.H"
#include "surfaceVelocityTools.H"
#include "pointPatchFields.H"
#include "addToRunTimeSelectionTable.H"
#include "vectorField.H"
#include "volFields.H"


namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

onePhaseStefanMeltVelocityPointPatchVectorField::
onePhaseStefanMeltVelocityPointPatchVectorField
(
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF,
    const dictionary& dict
)
:
    fixedValuePointPatchField<vector>(p, iF, dict),
    surfaceVelocityTools_(p, dict),
    TName_(dict.lookupOrDefault<word>("TName", "T")),
    kappaOverRhoH_(dict.lookup<scalar>("kappaOverRhoH")),
    URef_(dict.lookupOrDefault<vector>("URef", Zero))
{}


onePhaseStefanMeltVelocityPointPatchVectorField::
onePhaseStefanMeltVelocityPointPatchVectorField
(
    const onePhaseStefanMeltVelocityPointPatchVectorField& ptf,
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF,
    const pointPatchFieldMapper& mapper
)
:
    fixedValuePointPatchField<vector>(ptf, p, iF, mapper),
    surfaceVelocityTools_(ptf.surfaceVelocityTools_, p),
    TName_(ptf.TName_),
    kappaOverRhoH_(ptf.kappaOverRhoH_),
    URef_(ptf.URef_)
{}


onePhaseStefanMeltVelocityPointPatchVectorField::
onePhaseStefanMeltVelocityPointPatchVectorField
(
    const onePhaseStefanMeltVelocityPointPatchVectorField& ptf,
    const DimensionedField<vector, pointMesh>& iF
)
:
    fixedValuePointPatchField<vector>(ptf, iF),
    surfaceVelocityTools_(ptf.surfaceVelocityTools_),
    TName_(ptf.TName_),
    kappaOverRhoH_(ptf.kappaOverRhoH_),
    URef_(ptf.URef_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void onePhaseStefanMeltVelocityPointPatchVectorField::map
(
    const pointPatchField<vector>& ptf,
    const pointPatchFieldMapper& mapper
)
{
    const onePhaseStefanMeltVelocityPointPatchVectorField& oVptf =
        refCast<const onePhaseStefanMeltVelocityPointPatchVectorField>(ptf);

    fixedValuePointPatchField<vector>::map(oVptf, mapper);
    surfaceVelocityTools_.map(oVptf.surfaceVelocityTools_, mapper);
}


void onePhaseStefanMeltVelocityPointPatchVectorField::reset
(
    const pointPatchField<vector>& ptf
)
{
    const onePhaseStefanMeltVelocityPointPatchVectorField& oVptf =
        refCast<const onePhaseStefanMeltVelocityPointPatchVectorField>(ptf);

    fixedValuePointPatchField<vector>::reset(oVptf);
    surfaceVelocityTools_.reset(oVptf.surfaceVelocityTools_);
}


void onePhaseStefanMeltVelocityPointPatchVectorField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    const volScalarField& TField = this->db().lookupObject<volScalarField>(TName_);
    const fvPatchScalarField& TPatchField = TField.boundaryField()[this->patch().index()];
    const polyPatch& pp = TPatchField.patch().patch();

    const scalarField faceNormalVelocity
    (
        -kappaOverRhoH_ * TPatchField.snGrad() + (pp.faceNormals() & URef_)
    );

    const vectorField pointVelocity
        = surfaceVelocityTools_.pointVelocityFromFaceVelocity(faceNormalVelocity, pp);

    Field<vector>::operator=
    (
        std::move(pointVelocity)
    );
    fixedValuePointPatchField<vector>::updateCoeffs();
}


void onePhaseStefanMeltVelocityPointPatchVectorField::write(Ostream& os) const
{
    pointPatchField<vector>::write(os);
    surfaceVelocityTools_.write(os);
    writeEntry(os, "TName", TName_);
    writeEntry(os, "kappaOverRhoH", kappaOverRhoH_);
    writeEntry(os, "URef", URef_);
    writeEntry(os, "value", *this);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePointPatchTypeField
(
    pointPatchVectorField,
    onePhaseStefanMeltVelocityPointPatchVectorField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam
