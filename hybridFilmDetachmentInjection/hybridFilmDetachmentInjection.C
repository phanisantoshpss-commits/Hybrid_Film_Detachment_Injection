/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2017 OpenFOAM Foundation
    Copyright (C) 2020 OpenCFD Ltd.
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

#include "hybridFilmDetachmentInjection.H"
#include "addToRunTimeSelectionTable.H"
#include "kinematicSingleLayer.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace surfaceFilmModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(hybridFilmDetachmentInjection, 0);
addToRunTimeSelectionTable(injectionModel, hybridFilmDetachmentInjection, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

hybridFilmDetachmentInjection::hybridFilmDetachmentInjection
(
    surfaceFilmRegionModel& film,
    const dictionary& dict
)
:
    injectionModel(type(), film, dict),
    bondNumCrit_(coeffDict_.getOrDefault<scalar>("bondNumCrit", 1)),
    weberNumCrit_(coeffDict_.getOrDefault<scalar>("weberNumCrit", 5)),
    cm_(coeffDict_.getOrDefault<scalar>("cm", 0.25)),
    dCoeff_(coeffDict_.getOrDefault<scalar>("dCoeff", 3.3))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

hybridFilmDetachmentInjection::~hybridFilmDetachmentInjection()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void hybridFilmDetachmentInjection::correct
(
    scalarField& availableMass,
    scalarField& massToInject,
    scalarField& diameterToInject
)
{
    const kinematicSingleLayer& film =
        refCast<const kinematicSingleLayer>(this->film());

    // Calculate available dripping mass
    tmp<volScalarField> tgNorm(film.gNorm());
    const scalarField& gNorm = tgNorm();
    const scalar magg = mag(film.g().value());
    const volVectorField& nHat = film.nHat();

    // Properties of liquid film
    const scalarField& delta = film.delta();
    const scalarField& rho = film.rho();
    const scalarField& sigma = film.sigma();
    const volVectorField& UFilm = film.U();
    // Properties of carrier or gas field
    const volVectorField& Ugas = film.UPrimary();
    const scalarField& rhoGas = film.rhoPrimary();

    const fvMesh& mesh = film.regionMesh();
    const scalar deltaT = mesh.time().deltaTValue();


    forAll(delta, cellI)
    {
        bool dripOrStrip = false;
        const vector UrelVec = Ugas[cellI] - UFilm[cellI];
        const vector UrelTan = UrelVec - (UrelVec & nHat[cellI])*nHat[cellI];
        const scalar Urel = mag(UrelTan);
        const scalar magUg = mag(Ugas[cellI]);

        scalar bondNum = 0;

        if(gNorm[cellI] > SMALL)
        {
            //Bond number
            bondNum = ((rho[cellI] - rhoGas[cellI]) * gNorm[cellI] * pow(delta[cellI],2)) / sigma[cellI];       
        }

        //Weber number
        scalar weberNum = (rhoGas[cellI] * pow(Urel, 2) * delta[cellI])/sigma[cellI];

        const scalar FBo = max(bondNum/bondNumCrit_ - 1.0, 0.0);
        const scalar FWe = max(weberNum/weberNumCrit_ - 1.0, 0.0);

        //Determining Exceedence
        const scalar E = FBo + FWe;

        //Local capillary timescale
        scalar tau = max(pow(((rho[cellI] * pow(delta[cellI], 3))/sigma[cellI]), 0.5), SMALL);

        //Lambda determination
        scalar lambda = (cm_ * E) / tau;

        //Fraction of the locally available film mass
        scalar fDet = min(max((1 - (exp(-lambda * deltaT))), 0.0), 1.0);

        //Available mass for detachment is 
        scalar massLeave = fDet * availableMass[cellI];

        if(massLeave > 0)
        {    
            massToInject[cellI] += massLeave;
            availableMass[cellI] -= massLeave;
            addToInjectedMass(massLeave);
            scalar lc = sqrt(sigma[cellI]/(rho[cellI] * magg));
            scalar dBo = dCoeff_*lc;
            scalar dWe = 0.67 * ((pow(sigma[cellI], 0.75)) / (pow(magUg, 1.57)) );

            if(FBo > 0 && FWe == 0)
            { 
                scalar diam = dBo;
                diameterToInject[cellI] = diam;
            }

            else if(FBo == 0 && FWe > 0)
            {
                scalar diam = dWe;
                diameterToInject[cellI] = diam;
            }

            else if(FBo > 0 && FWe > 0)
            {
                //Weights of gravity dripping Wbo and shear stripping Wwe
                const scalar Wbo = (FBo/(FBo + FWe));
                const scalar Wwe = (FWe/(FBo + FWe));

                scalar diam = 1/((Wbo / dBo) + (Wwe / dWe));
                diameterToInject[cellI] = diam;
            }

            dripOrStrip = true;
        }

        if (!dripOrStrip)
        {
            diameterToInject[cellI] = 0;
            massToInject[cellI] = 0;
        }
    }

    injectionModel::correct();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceFilmModels
} // End namespace regionModels
} 