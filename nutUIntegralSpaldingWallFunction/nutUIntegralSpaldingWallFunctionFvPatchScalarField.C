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

        // error = u_p_int[facei] - magUp[facei]
        // Find the solution of the eqution "u_p_int[facei] = magUp[facei]" by "errorIntegralFunc"
        scalar error = 10.0;
        scalar error_old1 = errorIntegralFunc(E_, kappa_, magUp[facei], ut_old1, y[facei]*2, nuw[facei], num_points, 1.0);
        scalar error_old2 = errorIntegralFunc(E_, kappa_, magUp[facei], ut_old2, y[facei]*2, nuw[facei], num_points, 1.0);

        if (ROOTVSMALL < ut)
        {
            int iter = 0;

            do
            {
                // By Secant method, we are able to find u_tau without using Spalding's function.
                // This method needs two previous status.
                ut = ut_old1 - error_old1*(ut_old1 - ut_old2)/(error_old1 - error_old2 + ROOTVSMALL);
                error = errorIntegralFunc(E_, kappa_, magUp[facei], ut, y[facei]*2, nuw[facei], num_points, error_old1);

                err[facei] = abs(error)/magUp[facei];

                error_old2 = error_old1;
                ut_old2 = ut_old1;
                error_old1 = error;
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

    return tuTau;
}

Foam::scalar
Foam::nutUIntegralSpaldingWallFunctionFvPatchScalarField::velNewton
(
    const scalar E_,
    const scalar kappa_,
    scalar magUp_facei,
    scalar ut,
    const scalar y_facei,
    const scalar nuw_facei
) const
{
    scalar error = 1.0;

    do
    {
        scalar f = -(y_facei*ut/nuw_facei) + magUp_facei/ut + 1/E_*(exp(kappa_*magUp_facei/ut) - 1 - kappa_*magUp_facei/ut - (1/2)*(kappa_*magUp_facei/ut)*(kappa_*magUp_facei/ut) - (1/6)*(kappa_*magUp_facei/ut)*(kappa_*magUp_facei/ut)*(kappa_*magUp_facei/ut));
        // Derivative of f with respect to magUp_facei
        scalar df = 1/ut + 1/E_*(kappa_/ut*exp(kappa_*magUp_facei/ut) - kappa_/ut - (kappa_*kappa_/(ut*ut))*magUp_facei - 1/2*(kappa_*kappa_*kappa_/(ut*ut*ut))*(magUp_facei*magUp_facei));

        scalar u_new = magUp_facei - f/df;
        error = abs(u_new - magUp_facei);

        magUp_facei = u_new;

    } while (error > ROOTVSMALL);
    
    return magUp_facei;
}

Foam::scalar
Foam::nutUIntegralSpaldingWallFunctionFvPatchScalarField::velBisection
(
    const scalar E_,
    const scalar kappa_,    
    scalar ut,
    const scalar y_facei,
    const scalar nuw_facei
) const
{
    // Initial interval for Bisection Method
    scalar a = -1e3;
    scalar b = 1e6;
    scalar c = 0.0;

    do
    {
        // Left limit of interval
        scalar f_a = -(y_facei*ut/nuw_facei) + a/ut + 1/E_*(exp(kappa_*a/ut) - 1 - kappa_*a/ut - (1/2)*(kappa_*a/ut)*(kappa_*a/ut) - (1/6)*(kappa_*a/ut)*(kappa_*a/ut)*(kappa_*a/ut));
        // Right limit of interval
        //scalar f_b = -(y_facei*ut/nuw_facei) + b/ut + 1/E_*(exp(kappa_*b/ut) - 1 - kappa_*b/ut - (1/2)*(kappa_*b/ut)*(kappa_*b/ut) - (1/6)*(kappa_*b/ut)*(kappa_*b/ut)*(kappa_*b/ut));

        c = (a + b)/2;
        scalar f_c = -(y_facei*ut/nuw_facei) + c/ut + 1/E_*(exp(kappa_*c/ut) - 1 - kappa_*c/ut - (1/2)*(kappa_*c/ut)*(kappa_*c/ut) - (1/6)*(kappa_*c/ut)*(kappa_*c/ut)*(kappa_*c/ut));
        // If the solution is found before the interval is smaller than "ROOTVSMALL", "break" is applied.
        if (f_c == 0)
            break;

        if ((f_c >= 0 && f_a >= 0) || (f_c < 0 && f_a < 0))
            a = c;
        else
            b = c;
    } while ((b - a)/2 > ROOTVSMALL);

    scalar magUp_facei = c;

    return magUp_facei;   
}

Foam::scalar
Foam::nutUIntegralSpaldingWallFunctionFvPatchScalarField::integralNewton
(
    const scalar E_,
    const scalar kappa_,
    const scalar magUp_facei,
    scalar ut,
    const scalar yf_facei,
    const scalar nuw_facei,
    scalar num_points
) const
{
    scalar trap_sum = 0;

    for(int i = 0; i < num_points; i++)
    {
        // Use Trapezoidal Method for numerical integration
        trap_sum += velNewton(E_, kappa_, magUp_facei, ut, (yf_facei/num_points)*i, nuw_facei) + velNewton(E_, kappa_, magUp_facei, ut, (yf_facei/num_points)*(i+1), nuw_facei);
    }
    // LHS and RHS are already divided by the cell height.
    scalar u_p_int_facei = trap_sum*0.5*(1/num_points);

    return u_p_int_facei;
}

Foam::scalar
Foam::nutUIntegralSpaldingWallFunctionFvPatchScalarField::integralBisection
(
    const scalar E_,
    const scalar kappa_,    
    scalar ut,
    const scalar yf_facei,
    const scalar nuw_facei,
    scalar num_points
) const
{
    scalar trap_sum = 0;

    for(int i = 0; i < num_points; i++)
    {
        // Use Trapezoidal Method for numerical integration
        trap_sum += velBisection(E_, kappa_, ut, (yf_facei/num_points)*i, nuw_facei) + velBisection(E_, kappa_, ut, (yf_facei/num_points)*(i+1), nuw_facei);
    }
    // LHS and RHS are already divided by the cell height.
    scalar u_p_int_facei = trap_sum*0.5*(1/num_points);

    return u_p_int_facei;
}

Foam::scalar
Foam::nutUIntegralSpaldingWallFunctionFvPatchScalarField::integralCombined
(
    const scalar E_,
    const scalar kappa_,
    const scalar magUp_facei,
    scalar ut,
    const scalar yf_facei,
    const scalar nuw_facei,
    scalar num_points,
    scalar u_p_int
) const
{
    // Combine Newton and Bisection methods depending on u_p_int value
    if (u_p_int > 1e6)
        u_p_int = integralBisection(E_, kappa_, ut, yf_facei, nuw_facei, num_points);
    else
        u_p_int = integralNewton(E_, kappa_, magUp_facei, ut, yf_facei, nuw_facei, num_points);

    return u_p_int;
}

Foam::scalar
Foam::nutUIntegralSpaldingWallFunctionFvPatchScalarField::errorIntegralFunc
(
    const scalar E_,
    const scalar kappa_,
    const scalar magUp_facei,
    scalar ut,
    const scalar yf_facei,
    const scalar nuw_facei,
    scalar num_points,
    scalar error
) const
{
    scalar u_p_int = error + magUp_facei;
    // Use this function to find the solution of "u_p_int = magUp_facei"
    scalar func = integralCombined(E_, kappa_, magUp_facei, ut, yf_facei, nuw_facei, num_points, u_p_int) - magUp_facei;

    return func;
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
