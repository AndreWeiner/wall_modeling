// Hardcoded coefficients
scalar E_ = 9.8;
label maxIter_ = 100;
scalar tolerance_ = 1e-8;
scalar kappa = 0.41;

// Functions for the new solver
scalar spaldings_law
(
    const scalar u_plus,
    const scalar kappa
)
{
    return u_plus + 1.0/E_ * (Foam::exp(kappa*u_plus) - 1.0 - kappa*u_plus - (kappa*u_plus)*(kappa*u_plus)/2.0 - (kappa*u_plus)*(kappa*u_plus)*(kappa*u_plus)/6.0);
}

scalar spalding_velocity
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
            scalar df = 1/ut + 1/E_*(kappa/ut*Foam::exp(kappa*magUp_facei/ut) - kappa/ut - (kappa*kappa/(ut*ut))*magUp_facei - 1/2*(kappa*kappa*kappa/(ut*ut*ut))*(magUp_facei*magUp_facei));

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
        } while (residual > ROOTVSMALL && n_iter < maxIter_);        
    }

    return magUp_facei;
}

scalar average_velocity
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
