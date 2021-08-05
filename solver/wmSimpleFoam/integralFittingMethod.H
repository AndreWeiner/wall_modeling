// Employ a variable nu_ in order to use scalar instead of dimensionedScalar
const scalar nu_ = nu.value();
// Constant nu value is inserted to nu_patch (scalarField)
scalarField nu_patch(nuEff.boundaryField()[surfaceID].size());
forAll (adjacentCellIDs, faceI)
{
    nu_patch[faceI] = nu_;
}

// Interpolate face values for cell center velocities
const surfaceVectorField U_face = fvc::interpolate(U);
const surfaceVectorField U_sngrad = fvc::snGrad(U);

//***********************************************************************/
// Wall Correction
//***********************************************************************/

// Gradient of the velocity at the wall patch
const scalarField magGradU(mag(U_sngrad.boundaryField()[surfaceID]));

// Initialize a variable nut_w (for the wall)
scalarField nut_w(nuEff.boundaryField()[surfaceID].size());

// Distance between the first cell center normal to the wall and the wall
scalarField dist_wall = mag(mesh.Cf().boundaryField()[surfaceID]);
forAll (adjacentCellIDs, faceI)
{
    dist_wall[faceI] = mag(mesh.C()[adjacentCellIDs[faceI]] 
                    - mesh.Cf().boundaryField()[surfaceID][faceI]);
}
    
// Initialization of variables for the wall patch
    
// Use the velocity information of the first cell center normal to the wall
vectorField Uc = U.boundaryField()[surfaceID];
forAll (adjacentCellIDs, faceI)
{
    Uc[faceI] = U.internalField()[adjacentCellIDs[faceI]];
}
scalarField magUp(mag(U.boundaryField()[surfaceID] - Uc));

// nut at the wall
nut_w = nuEff.boundaryField()[surfaceID] - nu_patch;

tmp<scalarField> tuTau(new scalarField(mesh.Cf().boundaryField()[surfaceID].patch().size(), Zero));
scalarField& uTau = tuTau.ref();
scalarField err;
err.setSize(uTau.size());
err = 0.0;
// Direct calculation of uTau for the wall
forAll(uTau, facei)
{
    scalar ut = Foam::sqrt((nut_w[facei] + nu_)*magGradU[facei]);
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
        
    if (dist_wall[facei]*2 < 5e-4) // Newton's method
    {
        scalar diff_old1 = average_velocity(dist_wall[facei]*2, ut_old1, nu_, "NEWTON", kappa, num_points) - magUp[facei];
        scalar diff_old2 = average_velocity(dist_wall[facei]*2, ut_old2, nu_, "NEWTON", kappa, num_points) - magUp[facei];
        if (ROOTVSMALL < ut)
        {
            int iter = 0;

            do
            {
                // By Secant method, we are able to find u_tau without using Spalding's function.
                // This method needs two previous status.
                ut = ut_old1 - diff_old1*(ut_old1 - ut_old2)/(diff_old1 - diff_old2 + ROOTVSMALL);
                diff = average_velocity(dist_wall[facei]*2, ut, nu_, "NEWTON", kappa, num_points) - magUp[facei];

                err[facei] = abs(diff)/magUp[facei];

                diff_old2 = diff_old1;
                ut_old2 = ut_old1;
                diff_old1 = diff;
                ut_old1 = ut;

            } while
            (
                ut > ROOTVSMALL
            && err[facei] > tolerance_
            && ++iter < maxIter_
            );

            uTau[facei] = max(0.0, ut);
        }        
    }

    else // Bisection method
    {
        scalar diff_old1 = average_velocity(dist_wall[facei]*2, ut_old1, nu_, "BISECTION", kappa, num_points) - magUp[facei];
        scalar diff_old2 = average_velocity(dist_wall[facei]*2, ut_old2, nu_, "BISECTION", kappa, num_points) - magUp[facei];
        if (ROOTVSMALL < ut)
        {
            int iter = 0;

            do
            {
                // By Secant method, we are able to find u_tau without using Spalding's function.
                // This method needs two previous status.
                ut = ut_old1 - diff_old1*(ut_old1 - ut_old2)/(diff_old1 - diff_old2 + ROOTVSMALL);
                diff = average_velocity(dist_wall[facei]*2, ut, nu_, "BISECTION", kappa, num_points) - magUp[facei];

                err[facei] = abs(diff)/magUp[facei];

                diff_old2 = diff_old1;
                ut_old2 = ut_old1;
                diff_old1 = diff;
                ut_old1 = ut;

            } while
            (
                ut > ROOTVSMALL
            && err[facei] > tolerance_
            && ++iter < maxIter_
            );

            uTau[facei] = max(0.0, ut);
        }        
    }
}

// Calculate new viscosity at the wall
tmp<scalarField> tnutw
(
    max
    (
        scalar(0),
        sqr(tuTau)/(magGradU + ROOTVSMALL) - nu_patch
    )
);

// A reference variable nut_wall for nut_ at the bottomWall patch
scalarField& nuEff_wall = nuEff.boundaryFieldRef()[surfaceID];

// Save tnutw value to nut_wall (automatically applied to nuEff at the wall)
nuEff_wall = tnutw + nu_patch;


//***********************************************************************/
// First Cell Face Correction
//***********************************************************************/

// Initialize a variable nut_f (for the first cell face)
scalarField nut_f(nuEff.boundaryField()[surfaceID].size());

// Gradient of the velocity at the first cell face opposite to the wall
scalarField tmp_magGradUf(mag(U_sngrad.boundaryField()[surfaceID]));    
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

// Distance between the first cell face opposite to the wall and the wall
scalarField dist_fface = mag(mesh.Cf().boundaryField()[surfaceID]);
forAll (oppFaceIDs, faceI)
{
    dist_fface[faceI] = mag(mesh.Cf()[oppFaceIDs[faceI]] 
                    - mesh.Cf().boundaryField()[surfaceID][faceI]);
}  
        
// Use the velocity information of the first cell face opposite to the wall
vectorField Uf = U.boundaryField()[surfaceID];
forAll (oppFaceIDs, faceI)
{            
    Uf[faceI] = U_face[oppFaceIDs[faceI]] - U.boundaryField()[surfaceID][faceI];
}        

// magUpf is newly determined for the first cell face opposite to the wall
scalarField magUpf(mag(U.boundaryField()[surfaceID]));
magUpf = mag(Uf);

// nut at the first cell face opposite to the wall
forAll (oppFaceIDs, faceI)
{
    nut_f[faceI] = nuEff[oppFaceIDs[faceI]] - nu_patch[faceI];
}

tmp<scalarField> tuTauf(new scalarField(mesh.Cf().boundaryField()[surfaceID].patch().size(), Zero));
scalarField& uTauf = tuTauf.ref();
scalarField errf;
errf.setSize(uTauf.size());
errf = 0.0;
// Direct calculation of uTau for the first cell face opposite to the wall
forAll(uTauf, facei)
{
    scalar ut = Foam::sqrt((nut_f[facei] + nu_)*magGradUf[facei]);
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
        
    if (dist_fface[facei]*2 < 5e-4) // Newton's method
    {
        scalar diff_old1 = average_velocity(dist_fface[facei]*2, ut_old1, nu_, "NEWTON", kappa, num_points) - magUpf[facei];
        scalar diff_old2 = average_velocity(dist_fface[facei]*2, ut_old2, nu_, "NEWTON", kappa, num_points) - magUpf[facei];
        if (ROOTVSMALL < ut)
        {
            int iter = 0;

            do
            {
                // By Secant method, we are able to find u_tau without using Spalding's function.
                // This method needs two previous status.
                ut = ut_old1 - diff_old1*(ut_old1 - ut_old2)/(diff_old1 - diff_old2 + ROOTVSMALL);
                diff = average_velocity(dist_fface[facei]*2, ut, nu_, "NEWTON", kappa, num_points) - magUpf[facei];

                errf[facei] = abs(diff)/magUpf[facei];

                diff_old2 = diff_old1;
                ut_old2 = ut_old1;
                diff_old1 = diff;
                ut_old1 = ut;

            } while
            (
                ut > ROOTVSMALL
            && errf[facei] > tolerance_
            && ++iter < maxIter_
            );

            uTauf[facei] = max(0.0, ut);
        }        
    }

    else // Bisection method
    {
        scalar diff_old1 = average_velocity(dist_fface[facei]*2, ut_old1, nu_, "BISECTION", kappa, num_points) - magUpf[facei];
        scalar diff_old2 = average_velocity(dist_fface[facei]*2, ut_old2, nu_, "BISECTION", kappa, num_points) - magUpf[facei];
        if (ROOTVSMALL < ut)
        {
            int iter = 0;

            do
            {
                // By Secant method, we are able to find u_tau without using Spalding's function.
                // This method needs two previous status.
                ut = ut_old1 - diff_old1*(ut_old1 - ut_old2)/(diff_old1 - diff_old2 + ROOTVSMALL);
                diff = average_velocity(dist_fface[facei]*2, ut, nu_, "BISECTION", kappa, num_points) - magUpf[facei];

                errf[facei] = abs(diff)/magUpf[facei];

                diff_old2 = diff_old1;
                ut_old2 = ut_old1;
                diff_old1 = diff;
                ut_old1 = ut;

            } while
            (
                ut > ROOTVSMALL
            && errf[facei] > tolerance_
            && ++iter < maxIter_
            );

            uTauf[facei] = max(0.0, ut);
        }        
    }
}

// Calculate new viscosity at the first cell face opposite to the wall
tmp<scalarField> tnutf
(
    max
    (
        scalar(0),
        sqr(tuTauf)/(magGradUf + ROOTVSMALL) - nu_patch
    )
);

// Save tnutf value to the first cell face opposite to the wall
scalarField nutf = tnutf;
forAll (oppFaceIDs, faceI)
{
    // Save nut_ at the first cell face
    nuEff[oppFaceIDs[faceI]] = nutf[faceI] + nu_patch[faceI];
}


// Temp use
scalarField dist_second = mag(mesh.Cf().boundaryField()[surfaceID]);
forAll (secAdjacentCellIDs, faceI)
{
    dist_second[faceI] = mag(mesh.C()[secAdjacentCellIDs[faceI]] 
                    - mesh.Cf().boundaryField()[surfaceID][faceI]);
}
// For temporary sanity check at the middle of the plate
Info << "nuEff wall = " << nuEff.boundaryField()[surfaceID][200] << nl
<< "yfirst = " << dist_wall[200] << "\tUfirst = "
<< U[adjacentCellIDs[200]]
<< nl << "yfface = " << dist_fface[200] << "\tUface200 = "
<< U_face[oppFaceIDs[200]] << "\tnuEffface200 = " << nuEff[oppFaceIDs[200]] << nl
<< "ysecond = " << dist_second[200] << "\tUsecond = "
<< U[secAdjacentCellIDs[200]]
<< endl;
