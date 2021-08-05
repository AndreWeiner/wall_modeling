    // Momentum predictor

    MRF.correctBoundaryVelocity(U);
    surfaceScalarField nuEff = fvc::interpolate(turbulence->nuEff());
    #include "integralFittingMethod.H"
    tmp<fvVectorMatrix> tUEqn
    (
        fvm::div(phi, U)
      + MRF.DDt(U)
      // - fvc::div(turbulence->nuEff() * dev2(T(fvc::grad(U))))
      - fvm::laplacian(nuEff, U)
     ==
        fvOptions(U)
    );
    fvVectorMatrix& UEqn = tUEqn.ref();

    UEqn.relax();

    fvOptions.constrain(UEqn);

    if (simple.momentumPredictor())
    {
        solve(UEqn == -fvc::grad(p));

        fvOptions.correct(U);
    }
