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
    Foam::nutUIntegralSpaldingWallFunctionFvPatchScalarField

Group
    grpWallFunctions

Description
    This boundary condition provides a wall constraint on the turbulent
    viscosity, i.e. \c nut, based on velocity, i.e. \c U. Using Spalding's
    law gives a continuous \c nut profile to the wall.

        \f[
            y^+ = u^+ + \frac{1}{E} \left[exp(\kappa u^+) - 1 - \kappa u^+\,
                - 0.5 (\kappa u^+)^2 - \frac{1}{6} (\kappa u^+)^3\right]
        \f]

    where
    \vartable
        y^+     | Non-dimensional position
        u^+     | Non-dimensional velocity
        \kappa  | von Kármán constant
    \endvartable

Usage
    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        // Mandatory entries (unmodifiable)
        type            nutUIntegralSpaldingWallFunction;

        // Optional entries (unmodifiable)
        maxIter         10;
        tolerance       0.0001;

        // Optional (inherited) entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property  | Description                         | Type   | Req'd  | Dflt
      type      | Type name: nutUBlendedWallFunction  | word   | yes    | -
      maxIter   | Number of Newton-Raphson iterations | label  | no     | 10
      tolerance | Convergence tolerance               | scalar | no     | 0.0001
    \endtable

    The inherited entries are elaborated in:
      - \link nutWallFunctionFvPatchScalarField.H \endlink

Note
    - Suffers from non-exact restart since \c correctNut() (called through
    \c turbulence->validate) returns a slightly different value every time
    it is called. This is since the seed for the Newton-Raphson iteration
    uses the current value of \c *this (\c =nut ).
    - This can be avoided by overriding the tolerance. This also switches on
    a pre-detection whether the current nut already satisfies the turbulence
    conditions and if so does not change it at all. This means that the nut
    only changes if it 'has really changed'. This probably should be used with
    a tight tolerance, to make sure to kick every iteration, e.g.
        maxIter     100;
        tolerance   1e-7;

SourceFiles
    nutUIntegralSpaldingWallFunctionFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef nutUIntegralSpaldingWallFunctionFvPatchScalarField_H
#define nutUIntegralSpaldingWallFunctionFvPatchScalarField_H

#include "nutWallFunctionFvPatchScalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
          Class nutUIntegralSpaldingWallFunctionFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class nutUIntegralSpaldingWallFunctionFvPatchScalarField
:
    public nutWallFunctionFvPatchScalarField
{
protected:

    // Protected Data

        //- Max iterations in calcNut
        const label maxIter_;

        //- Convergence tolerance
        const scalar tolerance_;

        //- Uncomment in case of intrumentation
        //mutable uint64_t invocations_;
        //mutable uint64_t nontrivial_;
        //mutable uint64_t nonconvergence_;
        //mutable uint64_t iterations_;


    // Protected Member Functions

        //- Calculate the turbulent viscosity
        virtual tmp<scalarField> calcNut() const;

        //- Calculate the friction velocity
        virtual tmp<scalarField> calcUTau(const scalarField& magGradU) const;

        //- Calculate the friction velocity and number of iterations for
        //- convergence
        virtual tmp<scalarField> calcUTau
        (
            const scalarField& magGradU,
            const label maxIter,
            scalarField& err
        ) const;

        // Calculate uTau by the integral method
        virtual tmp<scalarField> calcIntegralUTau(const scalarField& magGradU) const;

        virtual tmp<scalarField> calcIntegralUTau
        (
            const scalarField& magGradU,
            const label maxIter,
            scalarField& err
        ) const;

        // Spalding's law of the wall
        virtual scalar spaldings_law
        (
            const scalar u_plus,
            const scalar E_,
            const scalar kappa_
        ) const;

        // Invert Spalding's law to compute physical velocity
        virtual scalar spalding_velocity
        (
            const scalar y_facei,
            scalar ut,
            const scalar nuw_facei,
            string method,
            const scalar E_,
            const scalar kappa_
        ) const;

        // Execute a numerical integration method to calculate a new average velocity
        virtual scalar average_velocity
        (
            const scalar yf_facei,
            scalar ut,
            const scalar nuw_facei,
            string method,
            const scalar E_,
            const scalar kappa_,
            scalar num_points
        ) const;
        
        //- Write local wall function variables
        virtual void writeLocalEntries(Ostream&) const;


public:

    //- Runtime type information
    TypeName("nutUIntegralSpaldingWallFunction");


    // Constructors

        //- Construct from patch and internal field
        nutUIntegralSpaldingWallFunctionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        nutUIntegralSpaldingWallFunctionFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //- nutUIntegralSpaldingWallFunctionFvPatchScalarField
        //- onto a new patch
        nutUIntegralSpaldingWallFunctionFvPatchScalarField
        (
            const nutUIntegralSpaldingWallFunctionFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        nutUIntegralSpaldingWallFunctionFvPatchScalarField
        (
            const nutUIntegralSpaldingWallFunctionFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new nutUIntegralSpaldingWallFunctionFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        nutUIntegralSpaldingWallFunctionFvPatchScalarField
        (
            const nutUIntegralSpaldingWallFunctionFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new nutUIntegralSpaldingWallFunctionFvPatchScalarField(*this, iF)
            );
        }


    //- Destructor
    virtual ~nutUIntegralSpaldingWallFunctionFvPatchScalarField();


    // Member Functions

        // Evaluation functions

            //- Calculate and return the yPlus at the boundary
            virtual tmp<scalarField> yPlus() const;


        // I-O

            //- Write
            virtual void write(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
