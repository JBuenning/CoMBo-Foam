#include "scalarField.H"
#include "twoPhaseStefanMeltVelocityPointPatchVectorField.H"
#include "surfaceVelocityTools.H"
#include "pointPatchFields.H"
#include "pointPatch.H"
#include "addToRunTimeSelectionTable.H"
#include "polyMesh.H"
#include "vectorField.H"
#include "volFields.H"
#include "mappedPatchBase.H"
#include "coupledPointPatchFieldMap.H"


namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

twoPhaseStefanMeltVelocityPointPatchVectorField::
twoPhaseStefanMeltVelocityPointPatchVectorField
(
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF,
    const dictionary& dict
)
:
    fixedValuePointPatchField<vector>(p, iF, dict),
    surfaceVelocityTools_(p, dict),
    TName_(dict.lookupOrDefault<word>("TName", "T")),
    kappa_(dict.lookup<scalar>("kappa")),
    rhoTimesH_(dict.lookup<scalar>("rhoTimesH")),
    URef_(dict.lookupOrDefault<vector>("URef", Zero)),
    externally_updated_(false)
{}


twoPhaseStefanMeltVelocityPointPatchVectorField::
twoPhaseStefanMeltVelocityPointPatchVectorField
(
    const twoPhaseStefanMeltVelocityPointPatchVectorField& ptf,
    const pointPatch& p,
    const DimensionedField<vector, pointMesh>& iF,
    const fieldMapper& mapper
)
:
    fixedValuePointPatchField<vector>(ptf, p, iF, mapper),
    surfaceVelocityTools_(ptf.surfaceVelocityTools_, p),
    TName_(ptf.TName_),
    kappa_(ptf.kappa_),
    rhoTimesH_(ptf.rhoTimesH_),
    URef_(ptf.URef_),
    externally_updated_(false)
{}


twoPhaseStefanMeltVelocityPointPatchVectorField::
twoPhaseStefanMeltVelocityPointPatchVectorField
(
    const twoPhaseStefanMeltVelocityPointPatchVectorField& ptf,
    const DimensionedField<vector, pointMesh>& iF
)
:
    fixedValuePointPatchField<vector>(ptf, iF),
    surfaceVelocityTools_(ptf.surfaceVelocityTools_),
    TName_(ptf.TName_),
    kappa_(ptf.kappa_),
    rhoTimesH_(ptf.rhoTimesH_),
    URef_(ptf.URef_),
    externally_updated_(false)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void twoPhaseStefanMeltVelocityPointPatchVectorField::map
(
    const pointPatchField<vector>& ptf,
    const fieldMapper& mapper
)
{
    const twoPhaseStefanMeltVelocityPointPatchVectorField& oVptf =
        refCast<const twoPhaseStefanMeltVelocityPointPatchVectorField>(ptf);

    fixedValuePointPatchField<vector>::map(oVptf, mapper);
    surfaceVelocityTools_.map(oVptf.surfaceVelocityTools_, mapper);
}


void twoPhaseStefanMeltVelocityPointPatchVectorField::reset
(
    const pointPatchField<vector>& ptf
)
{
    const twoPhaseStefanMeltVelocityPointPatchVectorField& oVptf =
        refCast<const twoPhaseStefanMeltVelocityPointPatchVectorField>(ptf);

    fixedValuePointPatchField<vector>::reset(oVptf);
    surfaceVelocityTools_.reset(oVptf.surfaceVelocityTools_);
}

void twoPhaseStefanMeltVelocityPointPatchVectorField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }
    if (externally_updated_)
    {
        fixedValuePointPatchField<vector>::updateCoeffs();
        return;
    }

    // (following two lines copied from coupledTemperatureFvPatchScalarField
    // I hope it's enough to account for everything)
    // Since we're inside initEvaluate/evaluate there might be processor
    // comms underway. Change the tag we use.
    int oldTag = UPstream::msgType();
    UPstream::msgType() = oldTag + 1;

    const volScalarField& TField = this->db().lookupObject<volScalarField>(TName_);
    const fvPatchScalarField& TPatchField = TField.boundaryField()[this->patch().index()];
    const polyPatch& pp = TPatchField.patch().poly();

    const mappedPatchBase& mpp = mappedPatchBase::getMap(pp);
    const label patchINbr = mpp.nbrPolyPatch().index();

    pointVectorField& nbrField = mpp.nbrMesh().lookupObjectRef<pointVectorField>(this->internalField().name());
    pointPatchField<vector>& nbrPatchField = nbrField.boundaryFieldRef()[patchINbr];
    if (!isA<twoPhaseStefanMeltVelocityPointPatchVectorField>(nbrPatchField))
    {
        FatalErrorInFunction
            << "Patch field for " << internalField().name() << " on "
            << this->patch().name() << " is of type "
            << twoPhaseStefanMeltVelocityPointPatchVectorField::typeName
            << endl << "The neighbouring patch field on"
            << nbrPatchField.patch().name() << " is required to be the same, but is "
            << "currently of type " << nbrPatchField.type() << exit(FatalError);
    }

    twoPhaseStefanMeltVelocityPointPatchVectorField& twoPhaseStefanNbr =
        refCast<twoPhaseStefanMeltVelocityPointPatchVectorField>(nbrPatchField);

    const fvPatch& fvPatchNbr =
        refCast<const fvMesh>(mpp.nbrMesh()).boundary()[patchINbr];

    const fvPatchScalarField& TPatchFieldNbr =
        fvPatchNbr.lookupPatchField<volScalarField, scalar>(twoPhaseStefanNbr.TName_);

    const scalarField faceNormalVelocity
    (
        (-kappa_ * TPatchField.snGrad() -
        twoPhaseStefanNbr.kappa_ * mpp.fromNeighbour(TPatchFieldNbr.snGrad())
        ) / rhoTimesH_
        +
        (URef_ & pp.faceNormals())
    );

    vectorField pointVelocity =
        surfaceVelocityTools_.pointVelocityFromFaceVelocity(faceNormalVelocity, pp);

    // update neighbour BC
    const coupledPointPatchFieldMap pointMap(pp);
    twoPhaseStefanNbr.Field<vector>::operator=
    (
        pointMap.toNeighbour(pointVelocity)
    );
    twoPhaseStefanNbr.externally_updated_ = true;

    Field<vector>::operator=
    (
        std::move(pointVelocity)
    );
    fixedValuePointPatchField<vector>::updateCoeffs();
}


void twoPhaseStefanMeltVelocityPointPatchVectorField::write(Ostream& os) const
{
    pointPatchField<vector>::write(os);
    surfaceVelocityTools_.write(os);
    writeEntry(os, "TName", TName_);
    writeEntry(os, "kappa", kappa_);
    writeEntry(os, "rhoTimesH", rhoTimesH_);
    writeEntry(os, "URef", URef_);
    writeEntry(os, "value", *this);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePointPatchTypeField
(
    pointPatchVectorField,
    twoPhaseStefanMeltVelocityPointPatchVectorField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam
