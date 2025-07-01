#include "explicitImplicitVelocityLaplacianFvMotionSolver.H"
#include "motionDiffusivity.H"
#include "fvmLaplacian.H"
#include "addToRunTimeSelectionTable.H"
#include "volPointInterpolation.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(explicitImplicitVelocityLaplacianFvMotionSolver, 0);

    addToRunTimeSelectionTable
    (
        motionSolver,
        explicitImplicitVelocityLaplacianFvMotionSolver,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::explicitImplicitVelocityLaplacianFvMotionSolver::explicitImplicitVelocityLaplacianFvMotionSolver
(
    const word& name,
    const polyMesh& mesh,
    const dictionary& dict
)
:
    velocityMotionSolver(name, mesh, dict, typeName),
    fvMotionSolver(mesh),
    cellMotionU_
    (
        IOobject
        (
            "cellMotionU",
            mesh.time().name(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        fvMesh_,
        dimensionedVector
        (
            "cellMotionU",
            pointMotionU_.dimensions(),
            Zero
        ),
        cellMotionBoundaryTypes<vector>(pointMotionU_.boundaryField())
    ),
    diffusivityPtr_
    (
        motionDiffusivity::New(fvMesh_, coeffDict().lookup("diffusivity"))
    ),
    timeIndex_(fvMesh_.time().timeIndex()),
    theta_(coeffDict().lookup<scalar>("theta"))
{
    // make sure old time is saved and available in first iteration
    pointMotionU_.oldTime();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::explicitImplicitVelocityLaplacianFvMotionSolver::~explicitImplicitVelocityLaplacianFvMotionSolver()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::pointField>
Foam::explicitImplicitVelocityLaplacianFvMotionSolver::curPoints() const
{
    volPointInterpolation::New(fvMesh_).interpolate
    (
        cellMotionU_,
        pointMotionU_
    );

    tmp<pointField> tcurPoints;

    if (timeIndex_ != fvMesh_.time().timeIndex())
    {
        timeIndex_ = fvMesh_.time().timeIndex();
        tcurPoints = tmp<pointField>
        (
            fvMesh_.points()
            + fvMesh_.time().deltaTValue() *
            (theta_ * pointMotionU_.primitiveField() + (1. - theta_) * pointMotionU_.oldTime().primitiveField())
        );
    }
    else
    {
        tcurPoints = tmp<pointField>
        (
            fvMesh_.oldPoints()
            + fvMesh_.time().deltaTValue() *
            (theta_ * pointMotionU_.primitiveField() + (1. - theta_) * pointMotionU_.oldTime().primitiveField())
        );
    }

    twoDCorrectPoints(tcurPoints.ref());

    return tcurPoints;
}


void Foam::explicitImplicitVelocityLaplacianFvMotionSolver::solve()
{

    // The points have moved so before interpolation update
    // the fvMotionSolver accordingly
    movePoints(fvMesh_.points());

    diffusivityPtr_->correct();

    pointMotionU_.boundaryFieldRef().updateCoeffs();

    Foam::solve
    (
        fvm::laplacian
        (
            diffusivityPtr_->operator()(),
            cellMotionU_,
            "laplacian(diffusivity,cellMotionU)"
        )
    );
}


//void Foam::explicitImplicitVelocityLaplacianFvMotionSolver::movePoints(const pointField& p)
//{
//    // Movement of pointMesh and volPointInterpolation already
//    // done by polyMesh,fvMesh
//}


void Foam::explicitImplicitVelocityLaplacianFvMotionSolver::topoChange
(
    const polyTopoChangeMap& map
)
{
    velocityMotionSolver::topoChange(map);

    // Update diffusivity. Note two stage to make sure old one is de-registered
    // before creating/registering new one.
    diffusivityPtr_.reset(nullptr);
    diffusivityPtr_ = motionDiffusivity::New
    (
        fvMesh_,
        coeffDict().lookup("diffusivity")
    );
    theta_ = coeffDict().lookup<scalar>("theta"); //maybe unnecessary
}


void Foam::explicitImplicitVelocityLaplacianFvMotionSolver::mapMesh
(
    const polyMeshMap& map
)
{
    velocityMotionSolver::mapMesh(map);

    // Update diffusivity. Note two stage to make sure old one is de-registered
    // before creating/registering new one.
    diffusivityPtr_.reset(nullptr);
    diffusivityPtr_ = motionDiffusivity::New
    (
        fvMesh_,
        coeffDict().lookup("diffusivity")
    );
    // theta_ = coeffDict().lookup<scalar>("theta"); //maybe unnecessary
}
