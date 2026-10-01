/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2017 OpenFOAM Foundation
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
    Foam::regionModels::surfaceFilmModels::hybridFilmDetachmentInjection

Description
 Hybrid surface-film detachment and cloud-injection model combining
 gravity-driven dripping and aerodynamic gas-shear stripping.

 A local Bond number characterises gravity detachment
 Bo = (rho_l - rho_g)*g_n*delta^2/sigma
 and shear stripping by a film Weber number based on the tangential
 gas-film relative velocity
 We_f = rho_g*|U_g - U_f|_t^2*delta/sigma.

 Detachment is activated when either dimensionless number exceeds its
 user-specified critical value. The amount of film mass transferred to
 the Lagrangian cloud is determined from the combined threshold
 exceedance using a capillary-timescale-based rate law.

 For gravity-only detachment, the injected droplet diameter follows the
 capillary-length relation used by OpenFOAM BrunDrippingInjection,
 d = dCoeff*sqrt(sigma/(rho_l*|g|)), with dCoeff = 3.3 by default.

 For shear-only detachment, the Sauter mean diameter is estimated using
 the empirical Sattelmayer prefilming-airblast correlation
 d32 = 0.67*sigma^0.75/U_g^1.57.

 If both mechanisms are active, the present model uses an SMD-inspired
 weighted harmonic interpolation between the gravity and shear diameter
 estimates. The weighting and mass-detachment rate law are model-specific
 closures introduced in this project.

References
 Brun, P. T., Damiano, A., Rieu, P., Balestra, G., & Gallaire, F. (2015).
 Rayleigh-Taylor instability under an inclined plane.
 Physics of Fluids, 27(8), 084107. doi:10.1063/1.4927857

 Lefebvre, A. H. (1988).
 Atomization and Sprays. CRC Press.
 Sattelmayer, T. (1985).
 Zum Einfluss der ausgebildeten, turbulenten Luft-Fluessigkeitsfilm-
 Stroemung auf den Filmzerfall und die Tropfenbildung am Austritt von
 Spalten geringer Hoehe. Dissertation, Universitaet Karlsruhe.

 Sattelmayer, T., & Wittig, S. (1986).
 Internal Flow Effects in Prefilming Airblast Atomizers: Mechanisms of
 Atomization and Droplet Spectra. Journal of Engineering for Gas Turbines
 and Power, 108(3), 465-472. doi:10.1115/1.3239931


SourceFiles
    hybridFilmDetachmentInjection.C

\*---------------------------------------------------------------------------*/

#ifndef hybridFilmDetachmentInjection_H
#define hybridFilmDetachmentInjection_H

#include "injectionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace surfaceFilmModels
{

/*---------------------------------------------------------------------------*\
                    Class hybridFilmDetachmentInjection Declaration
\*---------------------------------------------------------------------------*/

class hybridFilmDetachmentInjection
:
    public injectionModel
{
    // Private member functions

        //- No copy construct
        hybridFilmDetachmentInjection(const hybridFilmDetachmentInjection&) = delete;

        //- No copy assignment
        void operator=(const hybridFilmDetachmentInjection&) = delete;


protected:

    // Protected data

        //- Crtitical Bond number, default = 1
        scalar bondNumCrit_;

        //- Crticial Weber num, default = 5
        scalar weberNumCrit_;

        //- Detachment rate coefficient, default = 0.25
        scalar cm_;

        //- Coeffcient for droplet diameter, default = 3.3
        scalar dCoeff_;

        //- Diameters of particles to inject into the dripping
        scalarList diameter_;


public:

    //- Runtime type information
    TypeName("hybridFilmDetachmentInjection");


    // Constructors

        //- Construct from surface film model
        hybridFilmDetachmentInjection
        (
            surfaceFilmRegionModel& film,
            const dictionary& dict
        );


    //- Destructor
    virtual ~hybridFilmDetachmentInjection();


    // Member Functions

        //- Correct
        virtual void correct
        (
            scalarField& availableMass,
            scalarField& massToInject,
            scalarField& diameterToInject
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
