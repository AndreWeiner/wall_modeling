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

Class
    Foam::RASModels::wmSpalartAllmaras

Group
    grpRASTurbulence

Description
    Spalart-Allmaras one-eqn mixing-length model for incompressible and
    compressible external flows.

    Reference:
    \verbatim
        Spalart, P.R., & Allmaras, S.R. (1994).
        A one-equation turbulence model for aerodynamic flows.
        La Recherche Aerospatiale, 1, 5-21.
    \endverbatim

    The model is implemented without the trip-term and hence the ft2 term is
    not needed.

    It is necessary to limit the Stilda generation term as the model generates
    unphysical results if this term becomes negative which occurs for complex
    flow.  Several approaches have been proposed to limit Stilda but it is not
    clear which is the most appropriate.  Here the limiter proposed by Spalart
    is implemented in which Stilda is clipped at Cs*Omega with the default value
    of Cs = 0.3.

    The default model coefficients are
    \verbatim
        wmSpalartAllmarasCoeffs
        {
            Cb1         0.1355;
            Cb2         0.622;
            Cw2         0.3;
            Cw3         2.0;
            Cv1         7.1;
            Cs          0.3;
            sigmaNut    0.66666;
            kappa       0.41;
        }
    \endverbatim

SourceFiles
    wmSpalartAllmaras.C

\*---------------------------------------------------------------------------*/

#ifndef wmSpalartAllmaras_H
#define wmSpalartAllmaras_H

#include "RASModel.H"
#include "eddyViscosity.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
                      Class wmSpalartAllmaras Declaration
\*---------------------------------------------------------------------------*/

template<class BasicTurbulenceModel>
class wmSpalartAllmaras
:
    public eddyViscosity<RASModel<BasicTurbulenceModel>>
{
    // Private Member Functions

        //- No copy construct
        wmSpalartAllmaras(const wmSpalartAllmaras&) = delete;

        //- No copy assignment
        void operator=(const wmSpalartAllmaras&) = delete;


protected:

    // Protected data

        // Model coefficients

            dimensionedScalar sigmaNut_;
            dimensionedScalar kappa_;

            dimensionedScalar Cb1_;
            dimensionedScalar Cb2_;
            dimensionedScalar Cw1_;
            dimensionedScalar Cw2_;
            dimensionedScalar Cw3_;
            dimensionedScalar Cv1_;
            dimensionedScalar Cs_;
        /*---------------------- Newly added code ------------------------*/
        // Additional coefficients
            scalar E_ = 9.8;
            label maxIter_ = 100;
            scalar tolerance_ = 1e-8;
        /*----------------------------------------------------------------*/
        
        // Fields

            volScalarField nuTilda_;

            //- Wall distance
            //  Note: different to wall distance in parent RASModel
            //  which is for near-wall cells only
            const volScalarField& y_;


    // Protected Member Functions

        tmp<volScalarField> chi() const;

        tmp<volScalarField> fv1(const volScalarField& chi) const;

        tmp<volScalarField> fv2
        (
            const volScalarField& chi,
            const volScalarField& fv1
        ) const;

        tmp<volScalarField> Stilda
        (
            const volScalarField& chi,
            const volScalarField& fv1
        ) const;

        tmp<volScalarField> fw(const volScalarField& Stilda) const;

        void correctNut(const volScalarField& fv1);
        /*---------------------- Newly added code ------------------------*/
        // Calculate the turbulent viscosity
        scalarField calcNut(label patchi, const labelList& adjacentCellIDs);

        // Calculate uTau by the integral method
        tmp<scalarField> calcIntegralUTau(const scalarField& magGradU, label patchi, const labelList& adjacentCellIDs);

        tmp<scalarField> calcIntegralUTau
        (
            const scalarField& magGradU,
            const label maxIter,
            scalarField& err,
            label patchi,
            const labelList& adjacentCellIDs
        );

        // Spalding's law of the wall
        scalar spaldings_law
        (
            const scalar u_plus,
            const scalar kappa
        );

        // Invert Spalding's law to compute physical velocity
        scalar spalding_velocity
        (
            const scalar y_facei,
            scalar ut,
            const scalar nuw_facei,
            string method,
            const scalar kappa
        );

        // Execute a numerical integration method to calculate a new average velocity
        scalar average_velocity
        (
            const scalar yf_facei,
            scalar ut,
            const scalar nuw_facei,
            string method,
            const scalar kappa,
            scalar num_points
        );
        /*----------------------------------------------------------------*/
        virtual void correctNut();


public:

    typedef typename BasicTurbulenceModel::alphaField alphaField;
    typedef typename BasicTurbulenceModel::rhoField rhoField;
    typedef typename BasicTurbulenceModel::transportModel transportModel;


    //- Runtime type information
    TypeName("wmSpalartAllmaras");


    // Constructors

        //- Construct from components
        wmSpalartAllmaras
        (
            const alphaField& alpha,
            const rhoField& rho,
            const volVectorField& U,
            const surfaceScalarField& alphaRhoPhi,
            const surfaceScalarField& phi,
            const transportModel& transport,
            const word& propertiesName = turbulenceModel::propertiesName,
            const word& type = typeName
        );


    //- Destructor
    virtual ~wmSpalartAllmaras() = default;


    // Member Functions

        //- Re-read model coefficients if they have changed
        virtual bool read();

        //- Return the effective diffusivity for nuTilda
        tmp<volScalarField> DnuTildaEff() const;

        //- Return the turbulence kinetic energy
        virtual tmp<volScalarField> k() const;

        //- Return the turbulence kinetic energy dissipation rate
        virtual tmp<volScalarField> epsilon() const;

        //- Return the (estimated) specific dissipation rate
        virtual tmp<volScalarField> omega() const;

        //- Solve the turbulence equations and correct the turbulence viscosity
        virtual void correct();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "wmSpalartAllmaras.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
