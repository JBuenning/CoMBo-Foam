#ifndef surfaceVelocityInterpolation_H
#define surfaceVelocityInterpolation_H

#include "labelList.H"
#include "polyPatch.H"
#include "scalarField.H"
#include "vectorField.H"
#include "zero.H"
#include <limits>

namespace Foam
{

// linear interpolation of surface velocity from faces to points with
// inverse distance weighting
inline scalarField linearVelocityFaceToPointInterpolate
(
    const scalarField& faceNormalVelocity,
    const vectorField& pointNormals,
    const polyPatch& pp
)
{
    scalarField pointNormalVelocity(pp.nPoints(), Zero);

    const labelListList& pointFaceAddressing = pp.pointFaces();
    const vectorField facePositions(pp.faceCellCentres());

    forAll(pointFaceAddressing, pointIdx)
    {
        const vector& pointPosition = pp.localPoints()[pointIdx];
        const labelList& nbrFaces = pointFaceAddressing[pointIdx];
        scalar totWeight = 0;
        forAll(nbrFaces, nbrFaceI)
        {
            const label& faceIdx = nbrFaces[nbrFaceI];
            const scalar w = 1. / mag(pointPosition - facePositions[faceIdx]);
            totWeight += w;
            pointNormalVelocity[pointIdx]
                += w * faceNormalVelocity[faceIdx] / (pointNormals[pointIdx] & pp.faceNormals()[faceIdx]);
        }
        pointNormalVelocity[pointIdx] /= totWeight;
    }

    return pointNormalVelocity;
}


class UpwindSurfaceVelocityInterpolation
{
    labelList pointFaceMap_;

public:
    UpwindSurfaceVelocityInterpolation(const vectorField& flux, const polyPatch& pp)
    : pointFaceMap_(pp.nPoints())
    {
        const labelListList& pointFaceAddressing = pp.pointFaces();
        const vectorField facePositions(pp.faceCellCentres());
        forAll(pointFaceAddressing, pointIdx)
        {
            const vector& pointPosition = pp.localPoints()[pointIdx];
            const labelList& nbrFaces = pointFaceAddressing[pointIdx];
            label bestMatchFaceIdx = -1;
            scalar bestMatchValue = std::numeric_limits<scalar>::lowest();
            forAll(nbrFaces, nbrFaceI)
            {
                const label& faceIdx = nbrFaces[nbrFaceI];
                const vector cellMidToPointDirection =
                    normalised(pointPosition - facePositions[faceIdx]);
                const scalar candidate = cellMidToPointDirection & flux[pointIdx];
                if (candidate > bestMatchValue)
                {
                    bestMatchValue = candidate;
                    bestMatchFaceIdx = faceIdx;
                }
            }
            pointFaceMap_[pointIdx] = bestMatchFaceIdx;
        }
    }

    // template<typename T>
    // Field<T> faceToPointInterpolate(const Field<T>& faceField) const
    // {
    //     Field<T> result(pointFaceMap_.size());
    //     forAll(result, pointIdx)
    //     {
    //         result[pointIdx] = faceField[pointFaceMap_[pointIdx]];
    //     }
    //     return result;
    // }

    scalarField velocityFaceToPointInterpolate
    (
        const scalarField& faceNormalVelocity,
        const vectorField& faceNormals,
        const vectorField& pointNormals
    )
    {
        scalarField pointVelocity(pointFaceMap_.size());
        forAll(pointVelocity, pointIdx)
        {
            pointVelocity[pointIdx]
                = faceNormalVelocity[pointFaceMap_[pointIdx]]
                / (pointNormals[pointIdx] & faceNormals[pointFaceMap_[pointIdx]]);
        }
        return pointVelocity;
    }
};

} // namespace Foam

#endif
