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
    Foam::RASModels::negSpalartAllmaras

Group
    grpRASTurbulence

Description
    Negative Spalart-Allmaras one-eqn mixing-length model for incompressible and
    compressible external flows.

    Reference:
    \verbatim
        Spalart, P.R., & Allmaras, S.R. (1994).
        A one-equation turbulence model for aerodynamic flows.
        La Recherche Aerospatiale, 1, 5-21.

        Allmaras, S.R., Johnson, F.T., & Spalart, P.R. (2012)
        Modifications and Clarifications for the Implementation of the Spalart-Allmaras
        Turbulence Model
        7th International Conference on Computational Fluid Dynamics

        Anderson, W.K., Woos, S.L., & Allmaras, S.R. (2019)
        An initial exploration of Improved Numerics within  the Guidelines of the 
        Negative Spalart-Allmaras Turbulence Model
        NASA-TM-2019-220429 
    \endverbatim

    The model is implemented without the trip-term ft1 and hence the ft2 term is
    not needed.

    The coefficient alphaNegSA in the destruction in the PDE when nuTilda is negative 
    is set to 10 as a default value, but in the original paper of Allmaras et al. its 
    value is one. According to experiments carried out by Anderson et al., including
    the alphaNegSA coefficient enhances the convergence of the model. The value of 
    alphaNegSA is between 10 and 50. The coefficient Ct3 may be set to zero in some
    occasions to improve the stability of the model, especially at transonic regimes.

    The default model coefficients are
    \verbatim
        negSpalartAllmarasCoeffs
        {
            Cb1         0.1355;
            Cb2         0.622;
            Cw2         0.3;
            Cw3         2.0;
            Cv1         7.1;
            Cv2         0.7;
            Cv3         0.9
            Cn1         16;
            Ct3         1.2;
            alphaNegSA  10;
            sigmaNut    0.66666;
            kappa       0.41;
        }
    \endverbatim

SourceFiles
    negSpalartAllmaras.C

\*---------------------------------------------------------------------------*/

#ifndef negSpalartAllmaras_H
#define negSpalartAllmaras_H

#include "RASModel.H"
#include "eddyViscosity.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
                      Class negSpalartAllmaras Declaration
\*---------------------------------------------------------------------------*/

template<class BasicTurbulenceModel>
class negSpalartAllmaras
:
    public eddyViscosity<RASModel<BasicTurbulenceModel>>
{
    // Private Member Functions

        //- No copy construct
        negSpalartAllmaras(const negSpalartAllmaras&) = delete;

        //- No copy assignment
        void operator=(const negSpalartAllmaras&) = delete;


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
            dimensionedScalar Cv2_;
            dimensionedScalar Cv3_;
            dimensionedScalar Cn1_;
            dimensionedScalar Ct3_;
            dimensionedScalar alphaNegSA_;
            //dimensionedScalar Cs_;


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
        
        //  New function for diffussion term 

        tmp<volScalarField> fn(const volScalarField& chi) const;

        void correctNut(const volScalarField& fv1);
        virtual void correctNut();


public:

    typedef typename BasicTurbulenceModel::alphaField alphaField;
    typedef typename BasicTurbulenceModel::rhoField rhoField;
    typedef typename BasicTurbulenceModel::transportModel transportModel;


    //- Runtime type information
    TypeName("negSpalartAllmaras");


    // Constructors

        //- Construct from components
        negSpalartAllmaras
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
    virtual ~negSpalartAllmaras() = default;


    // Member Functions

        //- Re-read model coefficients if they have changed
        virtual bool read();

        //- Return the effective diffusivity for nuTilda
        tmp<volScalarField> DnuTildaEff(const volScalarField& chi) const;

        /*
        //- Return the effective diffusivity for nuTilda when negative
        tmp<volScalarField> DnuTildaEffNeg(const volScalarField& chi) const;
        */
       
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
    #include "negSpalartAllmaras.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
