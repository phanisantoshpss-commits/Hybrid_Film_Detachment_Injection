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
    Film Dripping mass transfer model.

    If the film thickness exceeds the critical value needed to generate one or
    more drops, the equivalent mass is removed from the film.  The critical film
    thickness is calculated from the Rayleigh-Taylor stability analysis of film
    flow on an inclined plane by Brun et.al.

    Reference:
    \verbatim
        Brun, P. T., Damiano, A., Rieu, P., Balestra, G., & Gallaire, F. (2015).
        Rayleigh-Taylor instability under an inclined plane.
        Physics of Fluids (1994-present), 27(8), 084107.
    \endverbatim

    The diameter of the drops formed are obtained from the local capillary
    length multiplied by the \c dCoeff coefficient which defaults to 3.3.

    Reference:
    \verbatim
        Lefebvre, A. (1988).
        Atomization and sprays
        (Vol. 1040, No. 2756). CRC press.
    \endverbatim

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
