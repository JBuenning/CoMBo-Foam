#include "surfaceVelocityTools.H"
#include "UList.H"
#include "fieldMapper.H"
#include "scalarField.H"
#include "polyMesh.H"
#include "Time.H"
#include "pointPatch.H"
#include "polyPatch.H"
#include "vector.H"
#include "vectorField.H"
#include "polyPatch.H"
#include "addToRunTimeSelectionTable.H"
#include "vectorField.H"
#include "PrimitivePatchInterpolation.H"
#include "surfaceVelocityInterpolation.H"


namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

surfaceVelocityTools::
surfaceVelocityTools
(
    const pointPatch& p,
    const dictionary& dict
)
:
    linearUpwindBlendingFactor_(dict.lookupOrDefault<scalar>("linearUpwindBlendingFactor", 0)),
    laplaceSmoothing_(dict.lookupOrDefault<bool>("laplaceSmoothing", false)),
    timeIndex_(-1),
    laplacePointCorrection_(p.size(), vector(0, 0, 0)),
    upwindInterp_(nullptr)
{}


surfaceVelocityTools::
surfaceVelocityTools
(
    const surfaceVelocityTools& other,
    const pointPatch& p
)
:
    linearUpwindBlendingFactor_(other.linearUpwindBlendingFactor_),
    laplaceSmoothing_(other.laplaceSmoothing_),
    timeIndex_(-1),
    laplacePointCorrection_(p.size(), vector(0, 0, 0)),
    upwindInterp_(nullptr)
{}


surfaceVelocityTools::
surfaceVelocityTools
(
    const surfaceVelocityTools& other
)
:
    linearUpwindBlendingFactor_(other.linearUpwindBlendingFactor_),
    laplaceSmoothing_(other.laplaceSmoothing_),
    timeIndex_(-1),
    laplacePointCorrection_(other.laplacePointCorrection_.size(), vector(0, 0, 0)),
    upwindInterp_(nullptr)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void surfaceVelocityTools::updateLaplaceCorrection(const polyPatch& pp)
{
    List<labelHashSet> pointPatchAddressing(pp.nPoints(), labelHashSet());
    forAll(pp.meshPoints(), localPointIdx)
    {
        const label& globalPointIdx = pp.meshPoints()[localPointIdx];

        forAll(pp.boundaryMesh(), patchI)
        {
            const polyPatch& other_pp = pp.boundaryMesh()[patchI];
            const label otherLocalIdx = other_pp.whichPoint(globalPointIdx);
            if (otherLocalIdx != -1) // point is both on this and on other patch
            {
                pointPatchAddressing[localPointIdx].insert(patchI);
            }
        }
    }

    const vectorField& points = pp.localPoints();
    List<label> nNeighbourPoints(pp.nPoints(), 0);
    laplacePointCorrection_ = vectorField(pp.nPoints(), vector(0, 0, 0));
    forAll(pp.edges(), edgeI)
    {
        const edge& e = pp.edges()[edgeI];
        const label& p1 = e[0];
        const label& p2 = e[1];

        const labelHashSet edgePatchAddressing(pointPatchAddressing[p1] & pointPatchAddressing[p2]);

        if (pointPatchAddressing[p1] == edgePatchAddressing)
        {
            laplacePointCorrection_ [p1] += points[p2] - points[p1];
            nNeighbourPoints[p1]++;
        }

        if (pointPatchAddressing[p2] == edgePatchAddressing)
        {
            laplacePointCorrection_ [p2] += points[p1] - points[p2];
            nNeighbourPoints[p2]++;
        }
    }

    forAll(pp.meshPoints(), localPointIdx)
    {
        vector& correction = laplacePointCorrection_ [localPointIdx];
        const label& globalPointIdx = pp.meshPoints()[localPointIdx];
        List<vector> restrictions;

        forAllConstIter(labelHashSet, pointPatchAddressing[localPointIdx], patchIdxIter)
        {
            const polyPatch& other_pp = pp.boundaryMesh()[patchIdxIter.key()];
            const label otherLocalIdx = other_pp.whichPoint(globalPointIdx);
            vector newRestriction = other_pp.pointNormals()[otherLocalIdx];
            forAll(restrictions, restrictionI) //Gram-Schmidt proc.
            {
                newRestriction -= (restrictions[restrictionI] & newRestriction) * restrictions[restrictionI];
            }
            newRestriction = normalised(newRestriction);
            correction -= (correction & newRestriction) * newRestriction;
            restrictions.append(std::move(newRestriction));
        }
        // factor of 2 to make laplace iteration stable, for both implicit and explicit mesh update
        correction /= 2. * max(1., nNeighbourPoints[localPointIdx]);
    }
}


vectorField surfaceVelocityTools::boundaryConformingPointNormals(const polyPatch &pp) const
{
    vectorField n = pp.pointNormals();

    forAll(pp.boundaryPoints(), boundaryPointIdx)
    {
        const label& localPointIdx = pp.boundaryPoints()[boundaryPointIdx];
        const label& globalPointIdx = pp.meshPoints()[localPointIdx];
        List<vector> restrictions;

        forAll(pp.boundaryMesh(), patchI)
        {
            const polyPatch& other_pp = pp.boundaryMesh()[patchI];

            if (other_pp.index() == pp.index()) continue; //skip this patch

            const label otherLocalIdx = other_pp.whichPoint(globalPointIdx);
            if (otherLocalIdx != -1) // point is both on this and on other patch
            {
                vector newRestriction = other_pp.pointNormals()[otherLocalIdx];
                forAll(restrictions, restrictionI) //Gram-Schmidt proc.
                {
                    newRestriction -= (restrictions[restrictionI] & newRestriction) * restrictions[restrictionI];
                }
                newRestriction = normalised(newRestriction);
                n[localPointIdx] -= (n[localPointIdx] & newRestriction) * newRestriction;
                restrictions.append(std::move(newRestriction));
            }
        }
        n[localPointIdx] = normalised(n[localPointIdx]);
    }
    return n;
}


void surfaceVelocityTools::map
(
    const surfaceVelocityTools& other,
    const fieldMapper& mapper
)
{
    mapper(laplacePointCorrection_, other.laplacePointCorrection_);
}


void surfaceVelocityTools::reset
(
    const surfaceVelocityTools& other
)
{
    laplacePointCorrection_.reset(other.laplacePointCorrection_);
}


vectorField surfaceVelocityTools::pointVelocityFromFaceVelocity
(
    const scalarField& faceNormalVelocity,
    const polyPatch& pp
)
{
    const Time& time = pp.boundaryMesh().db().time();
    if (time.timeIndex() != timeIndex_)
    {
        timeIndex_ = time.timeIndex();

        if(laplaceSmoothing_)
        {
            updateLaplaceCorrection(pp);
        }

        const PrimitivePatchInterpolation<polyPatch> naiveInterp(pp);
        const vectorField upwindFlux
        (
            naiveInterp.faceToPointInterpolate(faceNormalVelocity * pp.faceNormals())
        );
        upwindInterp_ = new UpwindSurfaceVelocityInterpolation(upwindFlux, pp);
    }

    const vectorField pointNormals = boundaryConformingPointNormals(pp);

    vectorField pointVelocity
    (
        linearUpwindBlendingFactor_ *
        linearVelocityFaceToPointInterpolate(faceNormalVelocity, pointNormals, pp) * pointNormals
        +
        (1-linearUpwindBlendingFactor_) *
        upwindInterp_->velocityFaceToPointInterpolate(faceNormalVelocity, pp.faceNormals(), pointNormals) * pointNormals
    );

    if (laplaceSmoothing_)
    {
        pointVelocity += laplacePointCorrection_ / time.deltaTValue();
    }
    return pointVelocity;
}

void surfaceVelocityTools::overwriteNeighborVeolocities
(
    const polyPatch& pp,
    pointVectorField::Boundary& velocityField
) const
{
    const pointPatchField<vector>& thisVelocityField = velocityField[pp.index()];

    forAll(pp.boundaryMesh(), patchI)
    {
        const polyPatch& other_pp = pp.boundaryMesh()[patchI];

        if (other_pp.index() == pp.index()) continue; //skip own patch

        pointPatchField<vector>& otherGenericVelocityField = velocityField[patchI];
        if (!isA<valuePointPatchField<vector>>(otherGenericVelocityField))
        {
            continue; // The patch field does not expose writable patch-local values
        }
        valuePointPatchField<vector>& otherVelocityField =
            refCast<valuePointPatchField<vector>>(otherGenericVelocityField);

        // maybe dangerous
        otherVelocityField.updateCoeffs();

        forAll(pp.boundaryPoints(), boundaryPointIdx)
        {
            const label& localPointIdx = pp.boundaryPoints()[boundaryPointIdx];
            const label& globalPointIdx = pp.meshPoints()[localPointIdx];

            const label otherLocalIdx = other_pp.whichPoint(globalPointIdx);
            if (otherLocalIdx != -1) // point is both on this and on other patch
            {
                otherVelocityField[otherLocalIdx] = thisVelocityField.internalField()[localPointIdx];
            }
        }
    }
}


void surfaceVelocityTools::write(Ostream& os) const
{
    writeEntry(os, "linearUpwindBlendingFactor", linearUpwindBlendingFactor_);
    writeEntry(os, "laplaceSmoothing", laplaceSmoothing_);
}

} // End namespace Foam
