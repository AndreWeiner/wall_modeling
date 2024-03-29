// Employ a variable nu_ in order to use scalar instead of dimensionedScalar
const scalar nu_ = nu.value();
// Velocity at wall
scalarField U_wall(mag(U.boundaryField()[surfaceID]));

// Interpolate face values for cell center velocities
surfaceVectorField U_face = fvc::interpolate(U);
surfaceVectorField U_sngrad = fvc::snGrad(U);
surfaceScalarField p_face = fvc::interpolate(p);

//***********************************************************************/
// Wall Correction
//***********************************************************************/

// Gradient of the velocity at the wall patch
const scalarField magGradU(mag(U_sngrad.boundaryField()[surfaceID]));

// Distance between the first cell face opposite to the wall and the wall
scalarField dist_fface = mag(mesh.Cf().boundaryField()[surfaceID]);
forAll (oppFaceIDs, faceI)
{
    dist_fface[faceI] = mag(mesh.Cf()[oppFaceIDs[faceI]] 
                    - mesh.Cf().boundaryField()[surfaceID][faceI]);
}
    
// Use the velocity information of the first cell center normal to the wall
scalarField magUp(mag(U.boundaryField()[surfaceID]));
forAll (adjacentCellIDs, faceI)
{
    magUp[faceI] = mag(U.internalField()[adjacentCellIDs[faceI]] - U.boundaryField()[surfaceID][faceI]);
}

// To check convergence of calculating uTau
boolField conv_uTau(mesh.Cf().boundaryField()[surfaceID].patch().size());

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

    // Newton's method
    scalar diff_old1 = average_velocity(dist_fface[faceI], ut_old1, nu_, "NEWTON", kappa, num_points) - magUp[faceI];
    scalar diff_old2 = average_velocity(dist_fface[faceI], ut_old2, nu_, "NEWTON", kappa, num_points) - magUp[faceI];
    if (ROOTVSMALL < ut)
    {
        int iter = 0;

        do
        {
            // By Secant method, we are able to find u_tau without using Spalding's function.
            // This method needs two previous status.
            ut = ut_old1 - diff_old1*(ut_old1 - ut_old2)/(diff_old1 - diff_old2 + ROOTVSMALL);
            diff = average_velocity(dist_fface[faceI], ut, nu_, "NEWTON", kappa, num_points) - magUp[faceI];

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

        // To check convergence of calculating uTau
        if (err[faceI] <= tolerance_)
        {
            conv_uTau[faceI] = 1;
        }
        else
        {
            conv_uTau[faceI] = 0;
        }

        uTau[faceI] = max(0.0, ut);
    }        
}

// yPlus at face
scalarField yPlusFace(nuEff.boundaryField()[surfaceID].size());
forAll (oppFaceIDs, faceI)
{    
    yPlusFace[faceI] = dist_fface[faceI]*uTau[faceI]/nu_;
}

// A small distance eps
scalar eps = dist_fface[0]/100.0;

// Calculate dU/dy at the wall
scalarField dUdy_w_spalding(nuEff.boundaryField()[surfaceID].size());
forAll (patchFaceIDs, faceI)
{
    dUdy_w_spalding[faceI] = (spalding_velocity(0, uTau[faceI], nu_, "NEWTON", kappa)
    - spalding_velocity(eps, uTau[faceI], nu_, "NEWTON", kappa))/eps;
}

// Normalization of velocity for blending
scalar alpha = 1.0;
scalarField U_tilde(mesh.Cf().boundaryField()[surfaceID].patch().size());
scalarField wall_b(mesh.Cf().boundaryField()[surfaceID].patch().size());
forAll (patchFaceIDs, faceI)
{
    // Between 0 and 1
    U_tilde[faceI] = (mag(U.internalField()[adjacentCellIDs[faceI]]) - U_wall[faceI])/(U_ref - U_wall[faceI]);
    // Blending function for wall : If y+ is larger, more portion of the slopes of Spalding are applied to the numerator.
    wall_b[faceI] = Foam::pow(U_tilde[faceI], alpha)*abs(dUdy_w_spalding[faceI]) + Foam::pow((1 - U_tilde[faceI]), alpha)*magGradU[faceI];
}

// Correction of nuEff at the wall
scalarField& nuEff_wall = nuEff.boundaryFieldRef()[surfaceID];

// Version 1 : blended wall correction
forAll (patchFaceIDs, faceI)
{
    nuEff_wall[faceI] = nu_*wall_b[faceI]/(magGradU[faceI] + ROOTVSMALL);
}

// Version 2 : fixed wall correction
forAll (patchFaceIDs, faceI)
{
    //nuEff_wall[faceI] = nu_*abs(dUdy_w_spalding[faceI])/(magGradU[faceI] + ROOTVSMALL);
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
        tmp_magGradUf[faceI] = 1e2;
    }        
    else
    {
        tmp_magGradUf[faceI] = mag(U_sngrad.internalField()[oppFaceIDs[faceI]]);
    }        
}
    
const scalarField magGradUf(tmp_magGradUf);

// Calculate dU/dy at the first cell face opposite to the wall
scalarField dUdy_f_spalding(nuEff.boundaryField()[surfaceID].size());
forAll (patchFaceIDs, faceI)
{
    dUdy_f_spalding[faceI] = (spalding_velocity(dist_fface[faceI] - eps, uTau[faceI], nu_, "NEWTON", kappa)
    - spalding_velocity(dist_fface[faceI] + eps, uTau[faceI], nu_, "NEWTON", kappa))/(2*eps);
    // To compare the velocity between -eps and +eps
    //Info << "spalding vel - eps diff " << faceI << " = " << spalding_velocity(dist_fface[faceI], uTau[faceI], nu_, "NEWTON", kappa) - spalding_velocity(dist_fface[faceI] - eps, uTau[faceI], nu_, "NEWTON", kappa) << nl
    //<< "spalding vel + eps diff " << faceI << " = " << spalding_velocity(dist_fface[faceI] + eps, uTau[faceI], nu_, "NEWTON", kappa) - spalding_velocity(dist_fface[faceI], uTau[faceI], nu_, "NEWTON", kappa) << nl << endl;
}

scalarField face_b(mesh.Cf().boundaryField()[surfaceID].patch().size());
forAll (patchFaceIDs, faceI)
{
    // Reversed blending function for 1st face : If y+ is smaller, more portion of the slopes of Spalding are applied to the denominator
    // in order to avoid exceedingly small magGradUf values. (ex. U_tilde = 0 -> no face correction is applied.)
    face_b[faceI] = Foam::pow(U_tilde[faceI], alpha)*magGradUf[faceI] + Foam::pow((1 - U_tilde[faceI]), alpha)*abs(dUdy_f_spalding[faceI]);
}

// Correction of nuEff at the first cell face opposite to the wall
forAll (oppFaceIDs, faceI)
{
    nuEff[oppFaceIDs[faceI]] *= abs(dUdy_f_spalding[faceI])/(face_b[faceI] + ROOTVSMALL);
}


// write surface data to file
const vectorField Cf(mesh.Cf().boundaryField()[surfaceID]);

string filename1 = getenv("yp");
//string filename2 = "yp_surfaceFields";
string filename2 = "yp_surfaceFields_iter_";
string filename3 = std::to_string(count);
string filename4 = ".csv";
filename1 += filename2 + filename3 + filename4;
//filename1 += filename2 + filename4;
string path1 = "yplus_";
string path2 = getenv("yp");
path1 += path2;

OFstream outputFile(path1/filename1);
outputFile.precision(12);
outputFile << "ccx,nuEff_wall_norm,nuEff_face_norm,nuEff_wall,nuEff_1st_face,dUdy_sp_wall,dUdy_sp_face,yPlusFace,nut_face,uTau,Conv_uTau,magGradU_wall,magGradU_face,U_tilde,wall_blending,face_blending,Ux_1st_cell,Uy_1st_cell,Ux_1st_face,Uy_1st_face,p_1st_cell,p_1st_face";

forAll(patchFaceIDs, faceI)
{
    outputFile << "\n";        
    outputFile << Cf[faceI].x() << ", "
                << nuEff.boundaryField()[surfaceID][faceI]/nu_ << ", "
                << nuEff[oppFaceIDs[faceI]]/nu_ << ", "
                << nuEff.boundaryField()[surfaceID][faceI] << ", "
                << nuEff[oppFaceIDs[faceI]] << ", "
                << dUdy_w_spalding[faceI] << ", "
                << dUdy_f_spalding[faceI] << ", "
                << yPlusFace[faceI] << ", "
                << nuEff[oppFaceIDs[faceI]] - nu_ << ", "
                << uTau[faceI] << ", "
                << conv_uTau[faceI] << ", "
                << magGradU[faceI] << ", "
                << magGradUf[faceI] << ", "
                << U_tilde[faceI] << ", "
                << wall_b[faceI] << ", "
                << face_b[faceI] << ", "
                << U[adjacentCellIDs[faceI]][0] << ", "
                << U[adjacentCellIDs[faceI]][1] << ", "
                << U_face[oppFaceIDs[faceI]][0] << ", "
                << U_face[oppFaceIDs[faceI]][1] << ", "
                << p[adjacentCellIDs[faceI]] << ", "
                << p_face[oppFaceIDs[faceI]];
}

count++;
