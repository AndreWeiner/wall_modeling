/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2019-2020 OpenCFD Ltd.
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

#include "nutUIntegralSpaldingWallFunctionFvPatchScalarField.H"
#include "turbulenceModel.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

Foam::tmp<Foam::scalarField>
Foam::nutUIntegralSpaldingWallFunctionFvPatchScalarField::calcNut() const
{
    const label patchi = patch().index();

    const turbulenceModel& turbModel = db().lookupObject<turbulenceModel>
    (
        IOobject::groupName
        (
            turbulenceModel::propertiesName,
            internalField().group()
        )
    );
    const fvPatchVectorField& Uw = U(turbModel).boundaryField()[patchi];
    const scalarField magGradU(mag(Uw.snGrad()));
    const tmp<scalarField> tnuw = turbModel.nu(patchi);
    const scalarField& nuw = tnuw();


    // Calculate new viscosity
    tmp<scalarField> tnutw
    (
        max
        (
            scalar(0),
            sqr(calcIntegralUTau(magGradU))/(magGradU + ROOTVSMALL) - nuw
        )
    );
    /*
    if (tolerance_ != 0.01)
    {
        // User-specified tolerance. Find out if current nut already satisfies
        // eqns.

        // Run ut for one iteration
        scalarField err;
        tmp<scalarField> UTau(calcUTau(magGradU, 1, err));

        // Preserve nutw if the initial error (err) already lower than the
        // tolerance.

        scalarField& nutw = tnutw.ref();
        forAll(err, facei)
        {
            if (err[facei] < tolerance_)
            {
                nutw[facei] = this->operator[](facei);
            }
        }
    }*/
    return tnutw;
}


Foam::tmp<Foam::scalarField>
Foam::nutUIntegralSpaldingWallFunctionFvPatchScalarField::calcUTau
(
    const scalarField& magGradU
) const
{
    scalarField err;
    return calcUTau(magGradU, maxIter_, err);
}


Foam::tmp<Foam::scalarField>
Foam::nutUIntegralSpaldingWallFunctionFvPatchScalarField::calcUTau
(
    const scalarField& magGradU,
    const label maxIter,
    scalarField& err
) const
{
    const label patchi = patch().index();

    const turbulenceModel& turbModel = db().lookupObject<turbulenceModel>
    (
        IOobject::groupName
        (
            turbulenceModel::propertiesName,
            internalField().group()
        )
    );
    const scalarField& y = turbModel.y()[patchi];

    const fvPatchVectorField& Uw = U(turbModel).boundaryField()[patchi];
    const scalarField magUp(mag(Uw.patchInternalField() - Uw));

    const tmp<scalarField> tnuw = turbModel.nu(patchi);
    const scalarField& nuw = tnuw();

    const scalarField& nutw = *this;

    tmp<scalarField> tuTau(new scalarField(patch().size(), Zero));
    scalarField& uTau = tuTau.ref();

    err.setSize(uTau.size());
    err = 0.0;

    forAll(uTau, facei)
    {
        scalar ut = sqrt((nutw[facei] + nuw[facei])*magGradU[facei]);
        // Note: for exact restart seed with laminar viscosity only:
        //scalar ut = sqrt(nuw[facei]*magGradU[facei]);

        if (ROOTVSMALL < ut)
        {
            int iter = 0;

            do
            {
                scalar kUu = min(kappa_*magUp[facei]/ut, 50);
                scalar fkUu = exp(kUu) - 1 - kUu*(1 + 0.5*kUu);

                scalar f =
                    - ut*y[facei]/nuw[facei]
                    + magUp[facei]/ut
                    + 1/E_*(fkUu - 1.0/6.0*kUu*sqr(kUu));

                scalar df =
                    y[facei]/nuw[facei]
                  + magUp[facei]/sqr(ut)
                  + 1/E_*kUu*fkUu/ut;

                scalar uTauNew = ut + f/df;
                err[facei] = mag((ut - uTauNew)/ut);
                ut = uTauNew;

                //iterations_++;

            } while
            (
                ut > ROOTVSMALL
             && err[facei] > tolerance_
             && ++iter < maxIter
            );

            uTau[facei] = max(0.0, ut);

            //invocations_++;
            //if (iter > 1)
            //{
            //    nontrivial_++;
            //}
            //if (iter >= maxIter_)
            //{
            //    nonconvergence_++;
            //}
        }
    }

    return tuTau;
}

Foam::tmp<Foam::scalarField>
Foam::nutUIntegralSpaldingWallFunctionFvPatchScalarField::calcIntegralUTau
(
    const scalarField& magGradU
) const
{
    scalarField err;
    return calcIntegralUTau(magGradU, maxIter_, err);
}

Foam::tmp<Foam::scalarField>
Foam::nutUIntegralSpaldingWallFunctionFvPatchScalarField::calcIntegralUTau
(
    const scalarField& magGradU,
    const label maxIter,
    scalarField& err
) const
{
    const label patchi = patch().index();

    const turbulenceModel& turbModel = db().lookupObject<turbulenceModel>
    (
        IOobject::groupName
        (
            turbulenceModel::propertiesName,
            internalField().group()
        )
    );
    const scalarField& y = turbModel.y()[patchi];

    const fvPatchVectorField& Uw = U(turbModel).boundaryField()[patchi];
    const scalarField magUp(mag(Uw.patchInternalField() - Uw));

    const tmp<scalarField> tnuw = turbModel.nu(patchi);
    const scalarField& nuw = tnuw();

    const scalarField& nutw = *this;

    tmp<scalarField> tuTau(new scalarField(patch().size(), Zero));
    scalarField& uTau = tuTau.ref();

    /*********************************************/
    // Test to access mesh methods    
    const fvMesh& mesh = patch().boundaryMesh().mesh();
    // Access polyPatch
    const polyPatch& bouPatch = patch().patch();

    // Get global IDs of faces in boundary patch
    labelList patchFaceIDs(bouPatch.size());
    forAll (bouPatch,faceI)
    {
        patchFaceIDs[faceI] = bouPatch.start()+faceI;
    }

    // Get IDs of primary adjacent cells
    const labelList& adjacentCellIDs = bouPatch.faceCells();

    // Get IDs of faces of the cell that are opposite to the patch face
    labelList oppFaceIDs(patchFaceIDs.size());

    forAll (patchFaceIDs, faceI)
    {
        oppFaceIDs[faceI] =
            mesh.cells()[adjacentCellIDs[faceI]].opposingFaceLabel
            (
                patchFaceIDs[faceI],mesh.faces()
            );
    }
    
    // Get IDs of secondary adjacent cells
    labelList secAdjacentCellIDs(adjacentCellIDs.size());

    label globFace = -1;
    forAll (oppFaceIDs, faceI)
    {
        globFace = oppFaceIDs[faceI];

        if (mesh.owner()[globFace] == adjacentCellIDs[faceI])
        {
            secAdjacentCellIDs[faceI] = mesh.neighbour()[globFace];
        }
        else
        {
            secAdjacentCellIDs[faceI] = mesh.owner()[globFace];
        }
    }

    // Get IDs of faces of the cell that are opposite to the opponent face
    labelList secOppFaceIDs(oppFaceIDs.size());

    forAll (oppFaceIDs, faceI)
    {
        secOppFaceIDs[faceI] =
            mesh.cells()[secAdjacentCellIDs[faceI]].opposingFaceLabel
            (
                oppFaceIDs[faceI],mesh.faces()
            );
    }

    // compute distance between boundary face and first cell face opposite
    scalarList distF(bouPatch.size(), 0.0);
    forAll(bouPatch, faceI)
    {
        distF[faceI] = mag(
                       patch().Cf()[faceI]
                     - patch().Cf()[oppFaceIDs[faceI]]
                    );
    }

    // area of boundary faces
    const scalarField magSf(mag(patch().Sf()));
        
    // print velocity for faces of the cell that are opposite to the patch face
    for (int i = 0; i < oppFaceIDs.size(); i++)
    {
        Info << "Index = " << oppFaceIDs[i] << nl
        << "Velocity = " << U(turbModel).internalField()[oppFaceIDs[i]] << endl;
    }

    // print velocity for secondary adjacent cells
    for (int i = 0; i < secOppFaceIDs.size(); i++)
    {
        Info << "Index = " << secOppFaceIDs[i] << nl
        << "Velocity = " << U(turbModel).internalField()[secOppFaceIDs[i]] << endl;
    }
    /*********************************************/

    err.setSize(uTau.size());
    err = 0.0;

    forAll(uTau, facei)
    {
        scalar ut = sqrt((nutw[facei] + nuw[facei])*magGradU[facei]);
        // Note: for exact restart seed with laminar viscosity only:
        //scalar ut = sqrt(nuw[facei]*magGradU[facei]);

        // Number of integral points
        scalar num_points = 10;
        // Initial uTau values
        scalar ut_old1 = ut;
        scalar ut_old2 = 20.0;

        // diff = u_avg[facei] - magUp[facei]
        // Find the solution of the eqution "u_avg[facei] = magUp[facei]"
        scalar diff = 10.0;
        
        if (y[facei]*2 < 5e-4) // Newton's method
        {
            scalar diff_old1 = average_velocity(y[facei]*2, ut_old1, nuw[facei], "NEWTON", E_, kappa_, num_points) - magUp[facei];
            scalar diff_old2 = average_velocity(y[facei]*2, ut_old2, nuw[facei], "NEWTON", E_, kappa_, num_points) - magUp[facei];
            if (ROOTVSMALL < ut)
            {
                int iter = 0;

                do
                {
                    // By Secant method, we are able to find u_tau without using Spalding's function.
                    // This method needs two previous status.
                    ut = ut_old1 - diff_old1*(ut_old1 - ut_old2)/(diff_old1 - diff_old2 + ROOTVSMALL);
                    diff = average_velocity(y[facei]*2, ut, nuw[facei], "NEWTON", E_, kappa_, num_points) - magUp[facei];

                    err[facei] = abs(diff)/magUp[facei];

                    diff_old2 = diff_old1;
                    ut_old2 = ut_old1;
                    diff_old1 = diff;
                    ut_old1 = ut;

                } while
                (
                    ut > ROOTVSMALL
                && err[facei] > tolerance_
                && ++iter < maxIter
                );

                uTau[facei] = max(0.0, ut);
            }        
        }

        else // Bisection method
        {
            scalar diff_old1 = average_velocity(y[facei]*2, ut_old1, nuw[facei], "BISECTION", E_, kappa_, num_points) - magUp[facei];
            scalar diff_old2 = average_velocity(y[facei]*2, ut_old2, nuw[facei], "BISECTION", E_, kappa_, num_points) - magUp[facei];
            if (ROOTVSMALL < ut)
            {
                int iter = 0;

                do
                {
                    // By Secant method, we are able to find u_tau without using Spalding's function.
                    // This method needs two previous status.
                    ut = ut_old1 - diff_old1*(ut_old1 - ut_old2)/(diff_old1 - diff_old2 + ROOTVSMALL);
                    diff = average_velocity(y[facei]*2, ut, nuw[facei], "BISECTION", E_, kappa_, num_points) - magUp[facei];

                    err[facei] = abs(diff)/magUp[facei];

                    diff_old2 = diff_old1;
                    ut_old2 = ut_old1;
                    diff_old1 = diff;
                    ut_old1 = ut;

                } while
                (
                    ut > ROOTVSMALL
                && err[facei] > tolerance_
                && ++iter < maxIter
                );

                uTau[facei] = max(0.0, ut);
            }        
        }
    }

    return tuTau;
}

Foam::scalar
Foam::nutUIntegralSpaldingWallFunctionFvPatchScalarField::spaldings_law
(
    const scalar u_plus,
    const scalar E_,
    const scalar kappa_
) const
{
    return u_plus + 1.0/E_ * (exp(kappa_*u_plus) - 1.0 - kappa_*u_plus - (kappa_*u_plus)*(kappa_*u_plus)/2.0 - (kappa_*u_plus)*(kappa_*u_plus)*(kappa_*u_plus)/6.0);
}

Foam::scalar
Foam::nutUIntegralSpaldingWallFunctionFvPatchScalarField::spalding_velocity
(
    const scalar y_facei,
    scalar ut,
    const scalar nuw_facei,
    string method,
    const scalar E_,
    const scalar kappa_
) const
{
    scalar magUp_facei = ut;
    scalar n_iter = 0;
    scalar residual = 1.0;

    if (method == "NEWTON") // Newton's method
    {
        do
        {            
            scalar f = spaldings_law(magUp_facei/ut, E_, kappa_) - y_facei*ut/nuw_facei;
            //scalar f = -(y_facei*ut/nuw_facei) + magUp_facei/ut + 1/E_*(exp(kappa_*magUp_facei/ut) - 1 - kappa_*magUp_facei/ut - (1/2)*(kappa_*magUp_facei/ut)*(kappa_*magUp_facei/ut) - (1/6)*(kappa_*magUp_facei/ut)*(kappa_*magUp_facei/ut)*(kappa_*magUp_facei/ut));
            // Derivative of f with respect to magUp_facei
            scalar df = 1/ut + 1/E_*(kappa_/ut*exp(kappa_*magUp_facei/ut) - kappa_/ut - (kappa_*kappa_/(ut*ut))*magUp_facei - 1/2*(kappa_*kappa_*kappa_/(ut*ut*ut))*(magUp_facei*magUp_facei));

            scalar u_new = magUp_facei - f/df;
            residual = abs(u_new - magUp_facei);

            magUp_facei = u_new;
            n_iter += 1;

        } while (residual > ROOTVSMALL && n_iter < maxIter_);        
    }

    else if (method == "BISECTION") // Bisection method
    {   
        scalar u_upper = ut * 100.0;
        scalar u_lower = 0.0;
        do
        {
            magUp_facei = 0.5 * (u_upper + u_lower);
            scalar f_upper = spaldings_law(u_upper/ut, E_, kappa_) - y_facei*ut/nuw_facei;
            scalar f_center = spaldings_law(magUp_facei/ut, E_, kappa_) - y_facei*ut/nuw_facei;
            if (f_center * f_upper > 0.0)
            {
                u_upper = magUp_facei;
            }
            else
            {
                u_lower = magUp_facei;
            }
            residual = abs(f_center);
            n_iter += 1;
        } while (residual > ROOTVSMALL && n_iter < maxIter_);        
    }

    return magUp_facei;
}

Foam::scalar
Foam::nutUIntegralSpaldingWallFunctionFvPatchScalarField::average_velocity
(
    const scalar yf_facei,
    scalar ut,
    const scalar nuw_facei,
    string method,
    const scalar E_,
    const scalar kappa_,
    scalar num_points
) const
{
    scalar trap_sum = 0;
    scalar u_avg_facei = 0;

    if (method == "NEWTON")
    {
        for (int i = 0; i < num_points; i++)
        {
            // Use Trapezoidal Method for numerical integration
            trap_sum += spalding_velocity((yf_facei/num_points)*i, ut, nuw_facei, "NEWTON", E_, kappa_) + spalding_velocity((yf_facei/num_points)*(i+1), ut, nuw_facei, "NEWTON", E_, kappa_);
        }
        // LHS and RHS are already divided by the cell height.
        u_avg_facei = trap_sum*0.5*(1/num_points);
    }

    else if (method == "BISECTION")
    {
        for (int i = 0; i < num_points; i++)
        {
            // Use Trapezoidal Method for numerical integration
            trap_sum += spalding_velocity((yf_facei/num_points)*i, ut, nuw_facei, "BISECTION", E_, kappa_) + spalding_velocity((yf_facei/num_points)*(i+1), ut, nuw_facei, "BISECTION", E_, kappa_);
        }
        // LHS and RHS are already divided by the cell height.
        u_avg_facei = trap_sum*0.5*(1/num_points);
    }

    return u_avg_facei;
}

void Foam::nutUIntegralSpaldingWallFunctionFvPatchScalarField::writeLocalEntries
(
    Ostream& os
) const
{
    nutWallFunctionFvPatchScalarField::writeLocalEntries(os);

    os.writeEntryIfDifferent<label>("maxIter", 10, maxIter_);
    os.writeEntryIfDifferent<scalar>("tolerance", 0.01, tolerance_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::nutUIntegralSpaldingWallFunctionFvPatchScalarField::
nutUIntegralSpaldingWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    nutWallFunctionFvPatchScalarField(p, iF),
    maxIter_(10),
    tolerance_(0.01)
    //invocations_(0),
    //nontrivial_(0),
    //nonconvergence_(0),
    //iterations_(0)
{}


Foam::nutUIntegralSpaldingWallFunctionFvPatchScalarField::
nutUIntegralSpaldingWallFunctionFvPatchScalarField
(
    const nutUIntegralSpaldingWallFunctionFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    nutWallFunctionFvPatchScalarField(ptf, p, iF, mapper),
    maxIter_(ptf.maxIter_),
    tolerance_(ptf.tolerance_)
    //invocations_(0),
    //nontrivial_(0),
    //nonconvergence_(0),
    //iterations_(0)
{}


Foam::nutUIntegralSpaldingWallFunctionFvPatchScalarField::
nutUIntegralSpaldingWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    nutWallFunctionFvPatchScalarField(p, iF, dict),
    maxIter_(dict.getOrDefault<label>("maxIter", 10)),
    tolerance_(dict.getOrDefault<scalar>("tolerance", 0.01))
    //invocations_(0),
    //nontrivial_(0),
    //nonconvergence_(0),
    //iterations_(0)
{}


Foam::nutUIntegralSpaldingWallFunctionFvPatchScalarField::
nutUIntegralSpaldingWallFunctionFvPatchScalarField
(
    const nutUIntegralSpaldingWallFunctionFvPatchScalarField& wfpsf
)
:
    nutWallFunctionFvPatchScalarField(wfpsf),
    maxIter_(wfpsf.maxIter_),
    tolerance_(wfpsf.tolerance_)
    //invocations_(wfpsf.invocations_),
    //nontrivial_(wfpsf.nontrivial_),
    //nonconvergence_(wfpsf.nonconvergence_),
    //iterations_(wfpsf.iterations_)
{}


Foam::nutUIntegralSpaldingWallFunctionFvPatchScalarField::
nutUIntegralSpaldingWallFunctionFvPatchScalarField
(
    const nutUIntegralSpaldingWallFunctionFvPatchScalarField& wfpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    nutWallFunctionFvPatchScalarField(wfpsf, iF),
    maxIter_(wfpsf.maxIter_),
    tolerance_(wfpsf.tolerance_)
    //invocations_(0),
    //nontrivial_(0),
    //nonconvergence_(0),
    //iterations_(0)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::nutUIntegralSpaldingWallFunctionFvPatchScalarField::
~nutUIntegralSpaldingWallFunctionFvPatchScalarField()
{
    //if (debug)
    //{
    //    Info<< "nutUIntegralSpaldingWallFunctionFvPatchScalarField :"
    //        << " total invocations:"
    //        << returnReduce(invocations_, sumOp<label>())
    //        << " total iterations:"
    //        << returnReduce(iterations_, sumOp<label>())
    //        << " total non-convergence:"
    //        << returnReduce(nonconvergence_, sumOp<label>())
    //        << " total non-trivial:"
    //        << returnReduce(nontrivial_, sumOp<label>())
    //        << endl;
    //}
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::scalarField>
Foam::nutUIntegralSpaldingWallFunctionFvPatchScalarField::yPlus() const
{
    const label patchi = patch().index();

    const turbulenceModel& turbModel = db().lookupObject<turbulenceModel>
    (
        IOobject::groupName
        (
            turbulenceModel::propertiesName,
            internalField().group()
        )
    );
    const scalarField& y = turbModel.y()[patchi];
    const fvPatchVectorField& Uw = U(turbModel).boundaryField()[patchi];
    const tmp<scalarField> tnuw = turbModel.nu(patchi);
    const scalarField& nuw = tnuw();

    return y*calcUTau(mag(Uw.snGrad()))/nuw;
}


void Foam::nutUIntegralSpaldingWallFunctionFvPatchScalarField::write
(
    Ostream& os
) const
{
    fvPatchField<scalar>::write(os);
    writeLocalEntries(os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        nutUIntegralSpaldingWallFunctionFvPatchScalarField
    );
}


// ************************************************************************* //
