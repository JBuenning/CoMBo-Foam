#ifndef coupledPointPatchFieldMap_H
#define coupledPointPatchFieldMap_H

#include "polyMesh.H"
#include "mappedPatchBase.H"

namespace Foam
{

class coupledPointPatchFieldMap
{
private:
    labelList pointToPointAddressing_;

public:
    coupledPointPatchFieldMap(const polyPatch& pp)
    {
        pointToPointAddressing_ = labelList(pp.nPoints(), -1);
        const mappedPatchBase& mpp = mappedPatchBase::getMap(pp);
        const faceList& nbrLocalFaces = mpp.nbrPolyPatch().localFaces();
        const faceList& thisLocalFaces = pp.localFaces();
        forAll(thisLocalFaces, patchFacei)
        {
            const face& fA = thisLocalFaces[patchFacei];
            const face& fB = nbrLocalFaces[patchFacei];

            forAll(fA, indexA)
            {
                label patchPointA = fA[indexA];

                if (pointToPointAddressing_[patchPointA] == -1)
                {
                    label indexB = (fB.size() - indexA) % fB.size();
                    pointToPointAddressing_[patchPointA] = fB[indexB];
                }
            }
        }
    }

    template<typename T>
    tmp<Field<T>> fromNeighbour(const Field<T>& nbrField) const
    {
        tmp<Field<T>> tmappedField(new Field<T>(nbrField.size()));
        Field<T>& mappedField = tmappedField.ref();
        forAll(mappedField, i)
        {
            mappedField[i] = nbrField[pointToPointAddressing_[i]];
        }
        return tmappedField;
    }

    template<typename T>
    tmp<Field<T>> fromNeighbour(const tmp<Field<T>>& tNbrField) const
    {
        const Field<T>& nbrField = tNbrField();
        return fromNeighbour(nbrField);
    }

    template<typename T>
    tmp<Field<T>> toNeighbour(const Field<T>& localField) const
    {
        tmp<Field<T>> tmappedField(new Field<T>(localField.size()));
        Field<T>& mappedField = tmappedField.ref();
        forAll(localField, i)
        {
            mappedField[pointToPointAddressing_[i]] = localField[i];
        }
        return tmappedField;
    }

    template<typename T>
    tmp<Field<T>> toNeighbour(const tmp<Field<T>>& tLocalField) const
    {
        const Field<T>& localField = tLocalField();
        return toNeighbour(localField);
    }

    ~coupledPointPatchFieldMap() = default;
};

} // namespace Foam

#endif
