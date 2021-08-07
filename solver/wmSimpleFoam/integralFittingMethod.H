// Employ a variable nu_ in order to use scalar instead of dimensionedScalar
const scalar nu_ = nu.value();

// Interpolate face values for cell center velocities
surfaceVectorField U_face = fvc::interpolate(U);
surfaceVectorField U_sngrad = fvc::snGrad(U);

//***********************************************************************/
// Wall Correction
//***********************************************************************/

// Gradient of the velocity at the wall patch
const scalarField magGradU(mag(U_sngrad.boundaryField()[surfaceID]));

// Distance between the first cell center normal to the wall and the wall
scalarField dist_wall = mag(mesh.Cf().boundaryField()[surfaceID]);
forAll (adjacentCellIDs, faceI)
{
    dist_wall[faceI] = mag(mesh.C()[adjacentCellIDs[faceI]] 
                    - mesh.Cf().boundaryField()[surfaceID][faceI]);
}

// Distance between the first cell face opposite to the wall and the wall
scalarField dist_fface = mag(mesh.Cf().boundaryField()[surfaceID]);
forAll (oppFaceIDs, faceI)
{
    dist_fface[faceI] = mag(mesh.Cf()[oppFaceIDs[faceI]] 
                    - mesh.Cf().boundaryField()[surfaceID][faceI]);
}
    
// Use the velocity information of the first cell center normal to the wall
vectorField Uc = U.boundaryField()[surfaceID];
forAll (adjacentCellIDs, faceI)
{
    Uc[faceI] = U.internalField()[adjacentCellIDs[faceI]];
}
scalarField magUp(mag(U.boundaryField()[surfaceID] - Uc));

tmp<scalarField> tuTau(new scalarField(mesh.Cf().boundaryField()[surfaceID].patch().size(), Zero));
scalarField& uTau = tuTau.ref();
scalarField err;
err.setSize(uTau.size());
err = 0.0;
// Direct calculation of uTau for the wall
forAll(uTau, faceI)
{
    scalar ut = Foam::sqrt(nuEff.boundaryField()[surfaceID][faceI]*magGradU[faceI]);
    // Note: for exact restart seed with laminar viscosity only:
    //scalar ut = sqrt(nuw[faceI]*magGradU[faceI]);

    // Number of integral points
    scalar num_points = 10;
    // Initial uTau values
    scalar ut_old1 = ut;
    scalar ut_old2 = 20.0;

    // diff = u_avg[faceI] - magUp[faceI]
    // Find the solution of the eqution "u_avg[faceI] = magUp[faceI]"
    scalar diff = 10.0;
        
    if (dist_wall[faceI]*2 < 5e-4) // Newton's method
    {
        scalar diff_old1 = average_velocity(dist_wall[faceI]*2, ut_old1, nu_, "NEWTON", kappa, num_points) - magUp[faceI];
        scalar diff_old2 = average_velocity(dist_wall[faceI]*2, ut_old2, nu_, "NEWTON", kappa, num_points) - magUp[faceI];
        if (ROOTVSMALL < ut)
        {
            int iter = 0;

            do
            {
                // By Secant method, we are able to find u_tau without using Spalding's function.
                // This method needs two previous status.
                ut = ut_old1 - diff_old1*(ut_old1 - ut_old2)/(diff_old1 - diff_old2 + ROOTVSMALL);
                diff = average_velocity(dist_wall[faceI]*2, ut, nu_, "NEWTON", kappa, num_points) - magUp[faceI];

                err[faceI] = abs(diff)/magUp[faceI];

                diff_old2 = diff_old1;
                ut_old2 = ut_old1;
                diff_old1 = diff;
                ut_old1 = ut;

            } while
            (
                ut > ROOTVSMALL
            && err[faceI] > tolerance_
            && ++iter < maxIter_
            );

            uTau[faceI] = max(0.0, ut);
        }        
    }

    else // Bisection method
    {
        scalar diff_old1 = average_velocity(dist_wall[faceI]*2, ut_old1, nu_, "BISECTION", kappa, num_points) - magUp[faceI];
        scalar diff_old2 = average_velocity(dist_wall[faceI]*2, ut_old2, nu_, "BISECTION", kappa, num_points) - magUp[faceI];
        if (ROOTVSMALL < ut)
        {
            int iter = 0;

            do
            {
                // By Secant method, we are able to find u_tau without using Spalding's function.
                // This method needs two previous status.
                ut = ut_old1 - diff_old1*(ut_old1 - ut_old2)/(diff_old1 - diff_old2 + ROOTVSMALL);
                diff = average_velocity(dist_wall[faceI]*2, ut, nu_, "BISECTION", kappa, num_points) - magUp[faceI];

                err[faceI] = abs(diff)/magUp[faceI];

                diff_old2 = diff_old1;
                ut_old2 = ut_old1;
                diff_old1 = diff;
                ut_old1 = ut;

            } while
            (
                ut > ROOTVSMALL
            && err[faceI] > tolerance_
            && ++iter < maxIter_
            );

            uTau[faceI] = max(0.0, ut);
        }        
    }
}
/*
// Calculate new viscosity at the wall
scalarField tnuEffw = sqr(uTau)/(magGradU + ROOTVSMALL);
scalarField& nuEff_wall = nuEff.boundaryFieldRef()[surfaceID];
nuEff_wall = tnuEffw;
*/

// A small distance eps
scalar eps = dist_fface[0]/1000.0;

// Calculate dU/dy at the wall
scalarField dUdy_w_spalding(nuEff.boundaryField()[surfaceID].size());
forAll (patchFaceIDs, faceI)
{
    dUdy_w_spalding[faceI] = (spalding_velocity(0, uTau[faceI], nu_, "NEWTON", kappa)
    - spalding_velocity(eps, uTau[faceI], nu_, "NEWTON", kappa))/eps;
}

// Correction of nuEff at the wall
scalarField& nuEff_wall = nuEff.boundaryFieldRef()[surfaceID];
forAll (patchFaceIDs, faceI)
{
    nuEff_wall[faceI] = max(nu_, nu_*abs(dUdy_w_spalding[faceI])/(magGradU[faceI] + ROOTVSMALL));
}


//***********************************************************************/
// First Cell Face Correction
//***********************************************************************/

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

// Use the velocity information of the first cell face opposite to the wall
vectorField Uf = U.boundaryField()[surfaceID];
forAll (oppFaceIDs, faceI)
{            
    Uf[faceI] = U_face[oppFaceIDs[faceI]] - U.boundaryField()[surfaceID][faceI];
}        

// magUpf is newly determined for the first cell face opposite to the wall
scalarField magUpf(mag(U.boundaryField()[surfaceID]));
magUpf = mag(Uf);

tmp<scalarField> tuTauf(new scalarField(mesh.Cf().boundaryField()[surfaceID].patch().size(), Zero));
scalarField& uTauf = tuTauf.ref();
scalarField errf;
errf.setSize(uTauf.size());
errf = 0.0;
// Direct calculation of uTau for the first cell face opposite to the wall
forAll(uTauf, faceI)
{
    scalar ut = Foam::sqrt(nuEff[oppFaceIDs[faceI]]*magGradUf[faceI]);
    // Note: for exact restart seed with laminar viscosity only:
    //scalar ut = sqrt(nuw[faceI]*magGradU[faceI]);

    // Number of integral points
    scalar num_points = 10;
    // Initial uTau values
    scalar ut_old1 = ut;
    scalar ut_old2 = 20.0;

    // diff = u_avg[faceI] - magUp[faceI]
    // Find the solution of the eqution "u_avg[faceI] = magUp[faceI]"
    scalar diff = 10.0;
        
    if (dist_fface[faceI]*2 < 5e-4) // Newton's method
    {
        scalar diff_old1 = average_velocity(dist_fface[faceI]*2, ut_old1, nu_, "NEWTON", kappa, num_points) - magUpf[faceI];
        scalar diff_old2 = average_velocity(dist_fface[faceI]*2, ut_old2, nu_, "NEWTON", kappa, num_points) - magUpf[faceI];
        if (ROOTVSMALL < ut)
        {
            int iter = 0;

            do
            {
                // By Secant method, we are able to find u_tau without using Spalding's function.
                // This method needs two previous status.
                ut = ut_old1 - diff_old1*(ut_old1 - ut_old2)/(diff_old1 - diff_old2 + ROOTVSMALL);
                diff = average_velocity(dist_fface[faceI]*2, ut, nu_, "NEWTON", kappa, num_points) - magUpf[faceI];

                errf[faceI] = abs(diff)/magUpf[faceI];

                diff_old2 = diff_old1;
                ut_old2 = ut_old1;
                diff_old1 = diff;
                ut_old1 = ut;

            } while
            (
                ut > ROOTVSMALL
            && errf[faceI] > tolerance_
            && ++iter < maxIter_
            );

            uTauf[faceI] = max(0.0, ut);
        }        
    }

    else // Bisection method
    {
        scalar diff_old1 = average_velocity(dist_fface[faceI]*2, ut_old1, nu_, "BISECTION", kappa, num_points) - magUpf[faceI];
        scalar diff_old2 = average_velocity(dist_fface[faceI]*2, ut_old2, nu_, "BISECTION", kappa, num_points) - magUpf[faceI];
        if (ROOTVSMALL < ut)
        {
            int iter = 0;

            do
            {
                // By Secant method, we are able to find u_tau without using Spalding's function.
                // This method needs two previous status.
                ut = ut_old1 - diff_old1*(ut_old1 - ut_old2)/(diff_old1 - diff_old2 + ROOTVSMALL);
                diff = average_velocity(dist_fface[faceI]*2, ut, nu_, "BISECTION", kappa, num_points) - magUpf[faceI];

                errf[faceI] = abs(diff)/magUpf[faceI];

                diff_old2 = diff_old1;
                ut_old2 = ut_old1;
                diff_old1 = diff;
                ut_old1 = ut;

            } while
            (
                ut > ROOTVSMALL
            && errf[faceI] > tolerance_
            && ++iter < maxIter_
            );

            uTauf[faceI] = max(0.0, ut);
        }        
    }
}

// Calculate dU/dy at the first cell face opposite to the wall
scalarField dUdy_f_spalding(nuEff.boundaryField()[surfaceID].size());
forAll (patchFaceIDs, faceI)
{
    dUdy_f_spalding[faceI] = (spalding_velocity(dist_fface[faceI] - eps, uTauf[faceI], nu_, "NEWTON", kappa)
    - spalding_velocity(dist_fface[faceI] + eps, uTauf[faceI], nu_, "NEWTON", kappa))/(2*eps);
}

// Correction of nuEff at the first cell face opposite to the wall
forAll (patchFaceIDs, faceI)
{
    nuEff[oppFaceIDs[faceI]] = max(nu_, nu_*abs(dUdy_f_spalding[faceI])/(magGradUf[faceI] + ROOTVSMALL));
}



// Sanity check for magGradU at the middle of the plate
Info << "magGradU[200] : " << magGradU[200] << nl
<< "magGradUf[200] : " << magGradUf[200] << endl;
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

const vectorField Cf(mesh.Cf().boundaryField()[surfaceID]);
// write surface data to file
//OFstream outputFile(runTime.path()/runTime.timeName()/"surfaceFields.csv");
OFstream outputFile("surfaceFields.csv");
outputFile.precision(12);
outputFile << "ccx,tauw,nuEff_wall,magGradU,uTau_w,nuEff_face,magGradUf,uTau_f";

forAll(patchFaceIDs, faceI)
{
    outputFile << "\n";        
    outputFile << Cf[faceI].x() << ", "
                << (uTau[faceI]*uTau[faceI]) << ", "
                << nuEff.boundaryField()[surfaceID][faceI] << ", "
                << magGradU[faceI] << ", "
                << uTau[faceI] << ", "
                << nuEff[oppFaceIDs[faceI]] << ", "
                << magGradUf[faceI] << ", "
                << uTauf[faceI];
}
