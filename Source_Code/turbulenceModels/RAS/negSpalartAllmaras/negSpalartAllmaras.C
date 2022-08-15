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

#include "negSpalartAllmaras.H"
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
tmp<volScalarField> negSpalartAllmaras<BasicTurbulenceModel>::chi() const
{
    return nuTilda_/this->nu();
}


template<class BasicTurbulenceModel>
tmp<volScalarField> negSpalartAllmaras<BasicTurbulenceModel>::fv1
(
    const volScalarField& chi
) const
{
    const volScalarField chi3(pow3(chi));
    return chi3/(chi3 + pow3(Cv1_));
}


template<class BasicTurbulenceModel>
tmp<volScalarField> negSpalartAllmaras<BasicTurbulenceModel>::fv2
(
    const volScalarField& chi,
    const volScalarField& fv1
) const
{
    return 1.0 - chi/(1.0 + chi*fv1);
}

// New definition of Stilda

template<class BasicTurbulenceModel>
tmp<volScalarField> negSpalartAllmaras<BasicTurbulenceModel>::Stilda
(
    const volScalarField& chi,
    const volScalarField& fv1
) const
{
    volScalarField Omega(::sqrt(2.0)*mag(skew(fvc::grad(this->U_))));
    volScalarField Sbar(nuTilda_*fv2(chi, fv1)/sqr(kappa_*y_));

    if(  Sbar >= -Cv2_*Omega  )
    {
        return (Omega + Sbar);
    }
    else
    {
        return
        (
            Omega + (Omega * (Cv2_*Cv2_*Omega + Cv3_*Sbar))/((Cv3_-2*Cv2_)*Omega - Sbar)
        );
    }

}


template<class BasicTurbulenceModel>
tmp<volScalarField> negSpalartAllmaras<BasicTurbulenceModel>::fw
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

    const volScalarField g
    (
       min
       (
          r + Cw2_*(pow6(r) - r),
          dimensionedScalar("GREAT", r.dimensions(), GREAT)
       ) 
    );

    return g*pow((1.0 + pow6(Cw3_))/(pow6(g) + pow6(Cw3_)), 1.0/6.0);
}

//New function for diffussion term

template<class BasicTurbulenceModel>
tmp<volScalarField> negSpalartAllmaras<BasicTurbulenceModel>::fn
(
    const volScalarField& chi
) const
{
    const volScalarField chi3(pow3(chi));
    return ( (Cn1_+ chi3)/(Cn1_ - chi3) );
}

template<class BasicTurbulenceModel>
void negSpalartAllmaras<BasicTurbulenceModel>::correctNut
(
    const volScalarField& fv1
)
{
    const dimensionedScalar zero ("zero", this->nut_.dimensions(), 0.0);

    if(pos0(nuTilda_))
    {
        this->nut_ = nuTilda_*fv1;
    }
    else
    {
        this->nut_ = zero;
    }

    this->nut_.correctBoundaryConditions();
    fv::options::New(this->mesh_).correct(this->nut_);

    BasicTurbulenceModel::correctNut();
}


template<class BasicTurbulenceModel>
void negSpalartAllmaras<BasicTurbulenceModel>::correctNut()
{
    correctNut(fv1(this->chi()));
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
negSpalartAllmaras<BasicTurbulenceModel>::negSpalartAllmaras
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
    //1st Coefficient for modified Stilda
    Cv2_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "Cv2",
            this->coeffDict_,
            0.7
        )
    ),
    //2nd Coefficient for modified Stilda
    Cv3_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "Cv3",
            this->coeffDict_,
            0.9
        )
    ),
    //Coefficient for new diffusion coefficient function (fn)
    Cn1_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "Cn1",
            this->coeffDict_,
            16
        )
    ),
    //Coefficient for production term in negative equation
    Ct3_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "Ct3",
            this->coeffDict_,
            1.2
        )
    ),
    //Coefficient for destruction term
    alphaNegSA_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "alphaNegSA",
            this->coeffDict_,
            10
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
bool negSpalartAllmaras<BasicTurbulenceModel>::read()
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
        Cv2_.readIfPresent(this->coeffDict());
        Cv3_.readIfPresent(this->coeffDict());
        Cn1_.readIfPresent(this->coeffDict());
        Ct3_.readIfPresent(this->coeffDict());
        alphaNegSA_.readIfPresent(this->coeffDict());

        return true;
    }

    return false;
}


template<class BasicTurbulenceModel>
tmp<volScalarField> negSpalartAllmaras<BasicTurbulenceModel>::DnuTildaEff
(
    const volScalarField& chi
) const
{
    if(pos0(nuTilda_))
    {
        return tmp<volScalarField>
        (
            new volScalarField("DnuTildaEff", (nuTilda_ + this->nu())/sigmaNut_)
        );
    }
    else
    {
        return tmp<volScalarField>
        (
            new volScalarField("DnuTildaEff", (nuTilda_ + (this->nu())*fn(chi))/sigmaNut_)
        );
    }
    
}

template<class BasicTurbulenceModel>
tmp<volScalarField> negSpalartAllmaras<BasicTurbulenceModel>::k() const
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
tmp<volScalarField> negSpalartAllmaras<BasicTurbulenceModel>::epsilon() const
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
tmp<volScalarField> negSpalartAllmaras<BasicTurbulenceModel>::omega() const
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
void negSpalartAllmaras<BasicTurbulenceModel>::correct()
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
        const volScalarField Omega(::sqrt(2.0)*mag(skew(fvc::grad(this->U_))));

        tmp<fvScalarMatrix> nuTildaEqn
        (
            fvm::ddt(alpha, rho, nuTilda_)
        + fvm::div(alphaRhoPhi, nuTilda_)
        - fvm::laplacian(alpha*rho*DnuTildaEff(chi), nuTilda_)
        - Cb2_/sigmaNut_*alpha*rho*magSqr(fvc::grad(nuTilda_))
        ==
          pos0(nuTilda_) * ( Cb1_*alpha*rho*Stilda*nuTilda_
                             - fvm::Sp(Cw1_*alpha*rho*fw(Stilda)*nuTilda_/sqr(y_), nuTilda_)
                             + fvOptions(alpha, rho, nuTilda_) )
        + neg(nuTilda_) * (  fvm::Sp(Cb1_*rho*(1.0 - Ct3_)*Omega, nuTilda_)
                             + alphaNegSA_*Cw1_*alpha*rho*nuTilda_*nuTilda_/sqr(y_)
                             + fvOptions(alpha, rho, nuTilda_) )
        );

        
        nuTildaEqn.ref().relax();
        fvOptions.constrain(nuTildaEqn.ref());
        solve(nuTildaEqn);
        fvOptions.correct(nuTilda_);
        nuTilda_.correctBoundaryConditions();   
        
    }

    // Update nut with latest available k,epsilon
    correctNut();

}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// ************************************************************************* //
