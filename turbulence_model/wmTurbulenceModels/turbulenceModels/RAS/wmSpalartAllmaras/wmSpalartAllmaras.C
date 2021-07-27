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

#include "wmSpalartAllmaras.H"
#include "fvOptions.H"
#include "bound.H"
#include "wallDist.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class BasicTurbulenceModel>
tmp<volScalarField> wmSpalartAllmaras<BasicTurbulenceModel>::chi() const
{
    return nuTilda_/this->nu();
}


template<class BasicTurbulenceModel>
tmp<volScalarField> wmSpalartAllmaras<BasicTurbulenceModel>::fv1
(
    const volScalarField& chi
) const
{
    const volScalarField chi3(pow3(chi));
    return chi3/(chi3 + pow3(Cv1_));
}


template<class BasicTurbulenceModel>
tmp<volScalarField> wmSpalartAllmaras<BasicTurbulenceModel>::fv2
(
    const volScalarField& chi,
    const volScalarField& fv1
) const
{
    return 1.0 - chi/(1.0 + chi*fv1);
}


template<class BasicTurbulenceModel>
tmp<volScalarField> wmSpalartAllmaras<BasicTurbulenceModel>::Stilda
(
    const volScalarField& chi,
    const volScalarField& fv1
) const
{
    volScalarField Omega(::sqrt(2.0)*mag(skew(fvc::grad(this->U_))));

    return
    (
        max
        (
            Omega
          + fv2(chi, fv1)*nuTilda_/sqr(kappa_*y_),
            Cs_*Omega
        )
    );
}


template<class BasicTurbulenceModel>
tmp<volScalarField> wmSpalartAllmaras<BasicTurbulenceModel>::fw
(
    const volScalarField& Stilda
) const
{
    volScalarField r
    (
        min
        (
            nuTilda_
           /(
               max
               (
                   Stilda,
                   dimensionedScalar("SMALL", Stilda.dimensions(), SMALL)
               )
              *sqr(kappa_*y_)
            ),
            scalar(10)
        )
    );
    r.boundaryFieldRef() == 0.0;

    const volScalarField g(r + Cw2_*(pow6(r) - r));

    return g*pow((1.0 + pow6(Cw3_))/(pow6(g) + pow6(Cw3_)), 1.0/6.0);
}


template<class BasicTurbulenceModel>
void wmSpalartAllmaras<BasicTurbulenceModel>::correctNut
(
    const volScalarField& fv1
)
{
    this->nut_ = nuTilda_*fv1;
    /*---------------------- Newly added code ------------------------*/
    // Required geometrical information
    label surfaceID(-1);
    forAll (this->mesh_.boundary(), patchI)
    {
        if (this->mesh_.boundary()[patchI].name() == "bottomWall")
        {
            surfaceID = patchI;
        }
    }    

    const polyPatch& bouPatch = this->mesh_.Cf().boundaryField()[surfaceID].patch().patch();

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
            this->mesh_.cells()[adjacentCellIDs[faceI]].opposingFaceLabel
            (
                patchFaceIDs[faceI],this->mesh_.faces()
            );
    }

    // Get IDs of secondary adjacent cells
    labelList secAdjacentCellIDs(adjacentCellIDs.size());

    label globFace = -1;
    forAll (oppFaceIDs, faceI)
    {
        globFace = oppFaceIDs[faceI];

        if (this->mesh_.owner()[globFace] == adjacentCellIDs[faceI])
        {
            secAdjacentCellIDs[faceI] = this->mesh_.neighbour()[globFace];
        }
        else
        {
            secAdjacentCellIDs[faceI] = this->mesh_.owner()[globFace];
        }
    }

    // Get IDs of faces of the cell that are opposite to the opponent face
    labelList secOppFaceIDs(oppFaceIDs.size());

    forAll (oppFaceIDs, faceI)
    {
        secOppFaceIDs[faceI] =
            this->mesh_.cells()[secAdjacentCellIDs[faceI]].opposingFaceLabel
            (
                oppFaceIDs[faceI],this->mesh_.faces()
            );
    }

    // Area of boundary faces
    const scalarField magSf(mag(this->mesh_.Sf().boundaryField()[surfaceID]));

    const surfaceScalarField y_face = fvc::interpolate(this->y_);
    // Interpolate face values for cell center velocities
    const surfaceVectorField U_face = fvc::interpolate(this->U_);
    const surfaceVectorField U_sngrad = fvc::snGrad(this->U_);
    
        
    // For temporary sanity check at the middle of the plate
    Info << "yfirst = " << this->y_[adjacentCellIDs[200]] << "\tUfirst = "
    << this->U_[adjacentCellIDs[200]] << "\tnutfirst = " << this->nut_[adjacentCellIDs[200]]
    << nl << "yface200 = " << y_face[oppFaceIDs[200]] << "\tUface200 = "
    << U_face[oppFaceIDs[200]] << "\tnutface200 = " << this->nut_face[oppFaceIDs[200]] << nl
    << "ysecond = " << this->y_[secAdjacentCellIDs[200]] << "\tUsecond = "
    << this->U_[secAdjacentCellIDs[200]] << "\tnutsecond = " << this->nut_[secAdjacentCellIDs[200]]
    << endl;


    // The nut calculation result is saved to this->nut_ directly
    // This calculation is executed instead of wall functions
    // Additionally, nut at the first cell face opposite to the wall is also calculated
    calcNut
    (
        surfaceID, adjacentCellIDs, oppFaceIDs, secAdjacentCellIDs, secOppFaceIDs,
        U_face, U_sngrad
    );    
    /*----------------------------------------------------------------*/

    this->nut_.correctBoundaryConditions();
    fv::options::New(this->mesh_).correct(this->nut_);
    
    BasicTurbulenceModel::correctNut();
}

/*--------------------------- Newly added code -----------------------------*/
template<class BasicTurbulenceModel>
void wmSpalartAllmaras<BasicTurbulenceModel>::calcNut
(
    label patchi, 
    const labelList& adjacentCellIDs, 
    labelList oppFaceIDs, 
    labelList secAdjacentCellIDs,
    labelList secOppFaceIDs,
    const surfaceVectorField U_face,
    const surfaceVectorField U_sngrad
)
{
    // Gradient of the velocity at the wall patch
    const scalarField magGradU(mag(U_sngrad.boundaryField()[patchi]));
    
    // nu at the wall
    tmp<scalarField> tnuw = this->nu(patchi);
    const scalarField& nuw = tnuw();

    // Calculate new viscosity at the wall
    tmp<scalarField> tnutw
    (
        max
        (
            scalar(0),
            sqr(calcIntegralUTau
            (
                magGradU, patchi, adjacentCellIDs, oppFaceIDs, 
                secAdjacentCellIDs, U_face, "Wall"
            ))/(magGradU + ROOTVSMALL) - nuw
        )
    );

    // A reference variable nut_wall for nut_ at the bottomWall patch
    scalarField& nut_wall = this->nut_.boundaryFieldRef()[patchi];

    // Save tnutw value to nut_wall (automatically applied to nut_ at the wall)    
    nut_wall = tnutw;

    
    // Gradient of the velocity at the first cell face opposite to the wall
    scalarField tmp_magGradUf(mag(U_sngrad.boundaryField()[patchi]));    
    forAll (oppFaceIDs, faceI)
    {
        // Set an arbitrary value in order to avoid zero denominator at the variable tnuf and 
        // overflow of the exponential term in the function spaldings_law for the first time step
        if (mag(U_sngrad.internalField()[oppFaceIDs[faceI]]) == 0)
        {
            tmp_magGradUf[faceI] = 1e5;
        }        
        else
        {
            tmp_magGradUf[faceI] = mag(U_sngrad.internalField()[oppFaceIDs[faceI]]);
        }        
    }
    
    const scalarField magGradUf(tmp_magGradUf);


    // Sanity check for magGradU at the middle of the plate
    Info << "magGradU[200] : " << magGradU[200] << nl
    << "magGradUf[200] : " << magGradUf[200] << endl;

    
    // nu at the first cell face opposite to the wall    
    scalarField nuf = this->nu(patchi);
    const volScalarField& tnu = this->nu();
    forAll (oppFaceIDs, faceI)
    {
        nuf[faceI] = tnu[oppFaceIDs[faceI]];
    }    

    // Calculate new viscosity at the first cell face opposite to the wall
    tmp<scalarField> tnutf
    (
        max
        (
            scalar(0),
            sqr(calcIntegralUTau
            (
                magGradUf, patchi, adjacentCellIDs, oppFaceIDs, 
                secAdjacentCellIDs, U_face, "Face"
            ))/(magGradUf + ROOTVSMALL) - nuf
        )
    );

    // After updating the wall value for every time step, nut_face is interpolated from nut_
    this->nut_face = fvc::interpolate(this->nut_);

    // Save tnutf value to the first cell face opposite to the wall
    scalarField tnut = tnutf;
    forAll (adjacentCellIDs, faceI)
    {
        // Average nut_ between the wall and the first cell face
        this->nut_[adjacentCellIDs[faceI]]
        = (nut_wall[faceI] + tnut[faceI])/2;
        // Average nut_ between the first cell face and the second cell face
        this->nut_[secAdjacentCellIDs[faceI]]
        = (tnut[faceI] + this->nut_face[secOppFaceIDs[faceI]])/2;
        // Save nut_ at the first cell face for later use
        this->nut_face[oppFaceIDs[faceI]] = tnut[faceI];
    }    
}

template<class BasicTurbulenceModel>
tmp<scalarField> wmSpalartAllmaras<BasicTurbulenceModel>::calcIntegralUTau
(
    const scalarField& magGradU,
    label patchi,
    const labelList& adjacentCellIDs,
    labelList oppFaceIDs, 
    labelList secAdjacentCellIDs,    
    const surfaceVectorField U_face,
    string field
)
{
    scalarField err;
    return calcIntegralUTau
    (
        magGradU, this->maxIter_, err, patchi, adjacentCellIDs, oppFaceIDs, 
        secAdjacentCellIDs, U_face, field
    );
}

template<class BasicTurbulenceModel>
tmp<scalarField> wmSpalartAllmaras<BasicTurbulenceModel>::calcIntegralUTau
(
    const scalarField& magGradU,
    const label maxIter,
    scalarField& err,
    label patchi,
    const labelList& adjacentCellIDs,
    labelList oppFaceIDs, 
    labelList secAdjacentCellIDs,
    const surfaceVectorField U_face,
    string field
)
{
    // Employ a variable kappa in order to use scalar instead of dimensionedScalar
    const scalar kappa = this->kappa_.value();

    // Initialize a variable nut_wf (wall and face)
    scalarField nut_wf(this->nut_.boundaryField()[patchi].size());    

    // Distance between the first cell center normal to the wall and the wall
    scalarField dist_wall = mag(this->mesh_.Cf().boundaryField()[patchi]);
    forAll (adjacentCellIDs, faceI)
    {
        dist_wall[faceI] = mag(this->mesh_.C()[adjacentCellIDs[faceI]] 
                        - this->mesh_.Cf().boundaryField()[patchi][faceI]);
    }    

    // Distance between the first cell face opposite to the wall and the wall
    scalarField dist_fface = mag(this->mesh_.Cf().boundaryField()[patchi]);
    forAll (oppFaceIDs, faceI)
    {
        dist_fface[faceI] = mag(this->mesh_.Cf()[oppFaceIDs[faceI]] 
                        - this->mesh_.Cf().boundaryField()[patchi][faceI]);
    }    
    
    // Initialization of variables for the wall patch
    // Save the first cell center height as a default value (for the wall patch)
    scalarField dist = dist_wall;    
    
    // Use the velocity information of the first cell center normal to the wall
    vectorField Uc = this->U_.boundaryField()[patchi];
    forAll (adjacentCellIDs, faceI)
    {
        Uc[faceI] = this->U_.internalField()[adjacentCellIDs[faceI]];
    }
    scalarField magUp(mag(this->U_.boundaryField()[patchi] - Uc));
    scalarField nu_wf = this->nu(patchi);    

    if (field == "Wall")
    {        
        // magUp, dist and nu_wf are already defined above
    
        // nut at the wall
        nut_wf = this->nut_.boundaryField()[patchi];
    }

    else if (field == "Face")
    {        
        dist = dist_fface;
        
        // Use the velocity information of the first cell face opposite to the wall
        vectorField Uf = this->U_.boundaryField()[patchi];
        forAll (oppFaceIDs, faceI)
        {            
            Uf[faceI] = U_face[oppFaceIDs[faceI]] - this->U_.boundaryField()[patchi][faceI];
        }        

        // magUp is newly determined for the first cell face opposite to the wall
        magUp = mag(Uf);
    
        // nu at the first cell face opposite to the wall
        const volScalarField& tnu = this->nu();
        forAll (oppFaceIDs, faceI)
        {
            nu_wf[faceI] = tnu[oppFaceIDs[faceI]];
        }        

        // nut at the first cell face opposite to the wall
        forAll (oppFaceIDs, faceI)
        {
            nut_wf[faceI] = this->nut_face[oppFaceIDs[faceI]];
        }
    }

    tmp<scalarField> tuTau(new scalarField(this->mesh_.Cf().boundaryField()[patchi].patch().size(), Zero));
    scalarField& uTau = tuTau.ref();
    err.setSize(uTau.size());
    err = 0.0;

    forAll(uTau, facei)
    {
        scalar ut = sqrt((nut_wf[facei] + nu_wf[facei])*magGradU[facei]);        
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
        
        if (dist[facei]*2 < 5e-4) // Newton's method
        {
            scalar diff_old1 = average_velocity(dist[facei]*2, ut_old1, nu_wf[facei], "NEWTON", kappa, num_points) - magUp[facei];
            scalar diff_old2 = average_velocity(dist[facei]*2, ut_old2, nu_wf[facei], "NEWTON", kappa, num_points) - magUp[facei];
            if (ROOTVSMALL < ut)
            {
                int iter = 0;

                do
                {
                    // By Secant method, we are able to find u_tau without using Spalding's function.
                    // This method needs two previous status.
                    ut = ut_old1 - diff_old1*(ut_old1 - ut_old2)/(diff_old1 - diff_old2 + ROOTVSMALL);
                    diff = average_velocity(dist[facei]*2, ut, nu_wf[facei], "NEWTON", kappa, num_points) - magUp[facei];

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
            scalar diff_old1 = average_velocity(dist[facei]*2, ut_old1, nu_wf[facei], "BISECTION", kappa, num_points) - magUp[facei];
            scalar diff_old2 = average_velocity(dist[facei]*2, ut_old2, nu_wf[facei], "BISECTION", kappa, num_points) - magUp[facei];
            if (ROOTVSMALL < ut)
            {
                int iter = 0;

                do
                {
                    // By Secant method, we are able to find u_tau without using Spalding's function.
                    // This method needs two previous status.
                    ut = ut_old1 - diff_old1*(ut_old1 - ut_old2)/(diff_old1 - diff_old2 + ROOTVSMALL);
                    diff = average_velocity(dist[facei]*2, ut, nu_wf[facei], "BISECTION", kappa, num_points) - magUp[facei];

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

template<class BasicTurbulenceModel>
scalar wmSpalartAllmaras<BasicTurbulenceModel>::spaldings_law
(
    const scalar u_plus,
    const scalar kappa
)
{
    return u_plus + 1.0/E_ * (exp(kappa*u_plus) - 1.0 - kappa*u_plus - (kappa*u_plus)*(kappa*u_plus)/2.0 - (kappa*u_plus)*(kappa*u_plus)*(kappa*u_plus)/6.0);
}

template<class BasicTurbulenceModel>
scalar wmSpalartAllmaras<BasicTurbulenceModel>::spalding_velocity
(
    const scalar dist_facei,
    scalar ut,
    const scalar nuwf_facei,
    string method,
    const scalar kappa
)
{
    scalar magUp_facei = ut;
    scalar n_iter = 0;
    scalar residual = 1.0;

    if (method == "NEWTON") // Newton's method
    {
        do
        {            
            scalar f = spaldings_law(magUp_facei/ut, kappa) - dist_facei*ut/nuwf_facei;            
            // Derivative of f with respect to magUp_facei
            scalar df = 1/ut + 1/E_*(kappa/ut*exp(kappa*magUp_facei/ut) - kappa/ut - (kappa*kappa/(ut*ut))*magUp_facei - 1/2*(kappa*kappa*kappa/(ut*ut*ut))*(magUp_facei*magUp_facei));

            scalar u_new = magUp_facei - f/df;            
            residual = abs(u_new - magUp_facei);

            magUp_facei = u_new;
            n_iter += 1;

        } while (residual > ROOTVSMALL && n_iter < this->maxIter_);
    }

    else if (method == "BISECTION") // Bisection method
    {   
        scalar u_upper = ut * 100.0;
        scalar u_lower = 0.0;
        do
        {
            magUp_facei = 0.5 * (u_upper + u_lower);
            
            scalar f_upper = spaldings_law(u_upper/ut, kappa) - dist_facei*ut/nuwf_facei;
            scalar f_center = spaldings_law(magUp_facei/ut, kappa) - dist_facei*ut/nuwf_facei;
            
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
        } while (residual > ROOTVSMALL && n_iter < this->maxIter_);        
    }

    return magUp_facei;
}

template<class BasicTurbulenceModel>
scalar wmSpalartAllmaras<BasicTurbulenceModel>::average_velocity
(
    const scalar distf_facei,
    scalar ut,
    const scalar nuwf_facei,
    string method,
    const scalar kappa,
    scalar num_points
)
{
    scalar trap_sum = 0;
    scalar u_avg_facei = 0;

    if (method == "NEWTON")
    {
        for (int i = 0; i < num_points; i++)
        {
            // Use Trapezoidal Method for numerical integration
            trap_sum += spalding_velocity((distf_facei/num_points)*i, ut, nuwf_facei, "NEWTON", kappa) + spalding_velocity((distf_facei/num_points)*(i+1), ut, nuwf_facei, "NEWTON", kappa);
        }
        // LHS and RHS are already divided by the cell height.
        u_avg_facei = trap_sum*0.5*(1/num_points);
    }

    else if (method == "BISECTION")
    {
        for (int i = 0; i < num_points; i++)
        {
            // Use Trapezoidal Method for numerical integration
            trap_sum += spalding_velocity((distf_facei/num_points)*i, ut, nuwf_facei, "BISECTION", kappa) + spalding_velocity((distf_facei/num_points)*(i+1), ut, nuwf_facei, "BISECTION", kappa);
        }
        // LHS and RHS are already divided by the cell height.
        u_avg_facei = trap_sum*0.5*(1/num_points);
    }

    return u_avg_facei;
}
/*--------------------------------------------------------------------------*/

template<class BasicTurbulenceModel>
void wmSpalartAllmaras<BasicTurbulenceModel>::correctNut()
{
    correctNut(fv1(this->chi()));
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
wmSpalartAllmaras<BasicTurbulenceModel>::wmSpalartAllmaras
(
    const alphaField& alpha,
    const rhoField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const transportModel& transport,
    const word& propertiesName,
    const word& type
)
:
    eddyViscosity<RASModel<BasicTurbulenceModel>>
    (
        type,
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        transport,
        propertiesName
    ),

    sigmaNut_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "sigmaNut",
            this->coeffDict_,
            0.66666
        )
    ),
    kappa_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "kappa",
            this->coeffDict_,
            0.41
        )
    ),

    Cb1_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "Cb1",
            this->coeffDict_,
            0.1355
        )
    ),
    Cb2_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "Cb2",
            this->coeffDict_,
            0.622
        )
    ),
    Cw1_(Cb1_/sqr(kappa_) + (1.0 + Cb2_)/sigmaNut_),
    Cw2_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "Cw2",
            this->coeffDict_,
            0.3
        )
    ),
    Cw3_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "Cw3",
            this->coeffDict_,
            2.0
        )
    ),
    Cv1_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "Cv1",
            this->coeffDict_,
            7.1
        )
    ),
    Cs_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "Cs",
            this->coeffDict_,
            0.3
        )
    ),

    nuTilda_
    (
        IOobject
        (
            "nuTilda",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),

    y_(wallDist::New(this->mesh_).y())
{
    if (type == typeName)
    {
        this->printCoeffs(type);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
bool wmSpalartAllmaras<BasicTurbulenceModel>::read()
{
    if (eddyViscosity<RASModel<BasicTurbulenceModel>>::read())
    {
        sigmaNut_.readIfPresent(this->coeffDict());
        kappa_.readIfPresent(this->coeffDict());

        Cb1_.readIfPresent(this->coeffDict());
        Cb2_.readIfPresent(this->coeffDict());
        Cw1_ = Cb1_/sqr(kappa_) + (1.0 + Cb2_)/sigmaNut_;
        Cw2_.readIfPresent(this->coeffDict());
        Cw3_.readIfPresent(this->coeffDict());
        Cv1_.readIfPresent(this->coeffDict());
        Cs_.readIfPresent(this->coeffDict());

        return true;
    }

    return false;
}


template<class BasicTurbulenceModel>
tmp<volScalarField> wmSpalartAllmaras<BasicTurbulenceModel>::DnuTildaEff() const
{
    return tmp<volScalarField>
    (
        new volScalarField("DnuTildaEff", (nuTilda_ + this->nu())/sigmaNut_)
    );
}


template<class BasicTurbulenceModel>
tmp<volScalarField> wmSpalartAllmaras<BasicTurbulenceModel>::k() const
{
    WarningInFunction
        << "Turbulence kinetic energy not defined for "
        << "Spalart-Allmaras model. Returning zero field"
        << endl;

    return tmp<volScalarField>::New
    (
        IOobject
        (
            "k",
            this->runTime_.timeName(),
            this->mesh_
        ),
        this->mesh_,
        dimensionedScalar(sqr(dimLength)/sqr(dimTime), Zero),
        zeroGradientFvPatchField<scalar>::typeName
    );
}


template<class BasicTurbulenceModel>
tmp<volScalarField> wmSpalartAllmaras<BasicTurbulenceModel>::epsilon() const
{
    WarningInFunction
        << "Turbulence kinetic energy dissipation rate not defined for "
        << "Spalart-Allmaras model. Returning zero field"
        << endl;

    return tmp<volScalarField>::New
    (
        IOobject
        (
            "epsilon",
            this->runTime_.timeName(),
            this->mesh_
        ),
        this->mesh_,
        dimensionedScalar(sqr(dimLength)/pow3(dimTime), Zero),
        zeroGradientFvPatchField<scalar>::typeName
    );
}


template<class BasicTurbulenceModel>
tmp<volScalarField> wmSpalartAllmaras<BasicTurbulenceModel>::omega() const
{
    WarningInFunction
        << "Specific dissipation rate not defined for "
        << "Spalart-Allmaras model. Returning zero field"
        << endl;

    return tmp<volScalarField>::New
    (
        IOobject
        (
            IOobject::groupName("omega", this->alphaRhoPhi_.group()),
            this->runTime_.timeName(),
            this->mesh_
        ),
        this->mesh_,
        dimensionedScalar(dimless/dimTime, Zero)
    );
}


template<class BasicTurbulenceModel>
void wmSpalartAllmaras<BasicTurbulenceModel>::correct()
{
    if (!this->turbulence_)
    {
        return;
    }

    {
        // Local references
        const alphaField& alpha = this->alpha_;
        const rhoField& rho = this->rho_;
        const surfaceScalarField& alphaRhoPhi = this->alphaRhoPhi_;
        fv::options& fvOptions(fv::options::New(this->mesh_));

        eddyViscosity<RASModel<BasicTurbulenceModel>>::correct();

        const volScalarField chi(this->chi());
        const volScalarField fv1(this->fv1(chi));

        const volScalarField Stilda(this->Stilda(chi, fv1));

        tmp<fvScalarMatrix> nuTildaEqn
        (
            fvm::ddt(alpha, rho, nuTilda_)
          + fvm::div(alphaRhoPhi, nuTilda_)
          - fvm::laplacian(alpha*rho*DnuTildaEff(), nuTilda_)
          - Cb2_/sigmaNut_*alpha*rho*magSqr(fvc::grad(nuTilda_))
         ==
            Cb1_*alpha*rho*Stilda*nuTilda_
          - fvm::Sp(Cw1_*alpha*rho*fw(Stilda)*nuTilda_/sqr(y_), nuTilda_)
          + fvOptions(alpha, rho, nuTilda_)
        );

        nuTildaEqn.ref().relax();
        fvOptions.constrain(nuTildaEqn.ref());
        solve(nuTildaEqn);
        fvOptions.correct(nuTilda_);
        bound(nuTilda_, dimensionedScalar(nuTilda_.dimensions(), Zero));
        nuTilda_.correctBoundaryConditions();
    }

    // Update nut with latest available k,epsilon
    correctNut();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// ************************************************************************* //
