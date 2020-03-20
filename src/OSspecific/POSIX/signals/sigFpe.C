/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
     \\/     M anipulation  |
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

#include "sigFpe.H"
#include "error.H"
#include "jobInfo.H"
#include "OSspecific.H"
#include "IOstreams.H"

#ifdef LINUX_GNUC
    #ifndef __USE_GNU
        #define __USE_GNU
    #endif
    #include <fenv.h>
    #include <malloc.h>
#endif

#if defined(darwin64)
#include <xmmintrin.h>
#include <sys/mman.h>
#include <unistd.h>
#endif

#include <limits>

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

struct sigaction Foam::sigFpe::oldAction_;

void Foam::sigFpe::fillNan(UList<scalar>& lst)
{
    lst = std::numeric_limits<scalar>::signaling_NaN();
}

bool Foam::sigFpe::mallocNanActive_ = false;

#if defined(darwin64)
void *(*Foam::sigFpe::oldMalloc_)(struct _malloc_zone_t *zone, size_t size)
    = NULL;

void *Foam::sigFpe::nanMalloc_(struct _malloc_zone_t *zone, size_t size)
{
    void *result;

    result = oldMalloc_(zone, size);

    if (mallocNanActive_)
    {
        UList<scalar> lst(reinterpret_cast<scalar*>(result),
                          size/sizeof(scalar));
        fillNan(lst);
    }

    return result;
}
#endif

#ifdef LINUX
extern "C"
{
    extern void* __libc_malloc(size_t size);

    // Override the GLIBC malloc to support mallocNan
    void* malloc(size_t size)
    {
        if (Foam::sigFpe::mallocNanActive_)
        {
            return Foam::sigFpe::mallocNan(size);
        }
        else
        {
            return __libc_malloc(size);
        }
    }
}

void* Foam::sigFpe::mallocNan(size_t size)
{
    // Call the low-level GLIBC malloc function
    void * result = __libc_malloc(size);

    // Initialize to signalling NaN
    UList<scalar> lst(reinterpret_cast<scalar*>(result), size/sizeof(scalar));
    sigFpe::fillNan(lst);

    return result;
}
#endif


#if defined(LINUX_GNUC) || defined(darwin64)
void Foam::sigFpe::sigHandler(int)
{
    // Reset old handling
    if (sigaction(SIGFPE, &oldAction_, nullptr) < 0)
    {
        FatalErrorInFunction
            << "Cannot reset SIGFPE trapping"
            << abort(FatalError);
    }

    // Update jobInfo file
    jobInfo_.signalEnd();

    error::printStack(Perr);

    // Throw signal (to old handler)
    raise(SIGFPE);
}
#endif


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sigFpe::sigFpe()
{
    oldAction_.sa_handler = nullptr;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::sigFpe::~sigFpe()
{
    if (env("FOAM_SIGFPE"))
    {
        #if defined(LINUX_GNUC) || defined(darwin64)
        // Reset signal
        if
        (
            oldAction_.sa_handler
         && sigaction(SIGFPE, &oldAction_, nullptr) < 0
        )
        {
            FatalErrorInFunction
                << "Cannot reset SIGFPE trapping"
                << abort(FatalError);
        }
        #endif
    }

    if (env("FOAM_SETNAN"))
    {
        #if defined(LINUX) || defined(darwin64)
        // Disable initialization to NaN
        mallocNanActive_ = false;
        #endif

        #if defined(darwin64)
        // Restoring old malloc handler
        if (oldMalloc_ != NULL) {
            malloc_zone_t *zone = malloc_default_zone();

            if (zone != NULL)
            {
                mprotect(zone, getpagesize(), PROT_READ | PROT_WRITE);
                zone->malloc = oldMalloc_;
                mprotect(zone, getpagesize(), PROT_READ);
            }
        }
        #endif
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::sigFpe::set(const bool verbose)
{
    if (oldAction_.sa_handler)
    {
        FatalErrorInFunction
            << "Cannot call sigFpe::set() more than once"
            << abort(FatalError);
    }

    if (env("FOAM_SIGFPE"))
    {
        bool supported = false;

        #if defined(LINUX_GNUC) || defined(darwin64)
        supported = true;

        #if defined(LINUX_GNUC)
        feenableexcept
        (
            FE_DIVBYZERO
          | FE_INVALID
          | FE_OVERFLOW
        );
        #endif  // LINUX_GNUC
        #if defined(darwin64)
        _mm_setcsr(_MM_MASK_MASK &~
                   (_MM_MASK_OVERFLOW|_MM_MASK_INVALID|_MM_MASK_DIV_ZERO));
        #endif

        struct sigaction newAction;
        newAction.sa_handler = sigHandler;
        newAction.sa_flags = SA_NODEFER;
        sigemptyset(&newAction.sa_mask);
        if (sigaction(SIGFPE, &newAction, &oldAction_) < 0)
        {
            FatalErrorInFunction
                << "Cannot set SIGFPE trapping"
                << abort(FatalError);
        }
        #endif

        if (verbose)
        {
            if (supported)
            {
                Info<< "sigFpe : Enabling floating point exception trapping"
                    << " (FOAM_SIGFPE)." << endl;
            }
            else
            {
                Info<< "sigFpe : Floating point exception trapping"
                    << " - not supported on this platform" << endl;
            }
        }
    }


    if (env("FOAM_SETNAN"))
    {
        #if defined(LINUX) || defined(darwin64)
        mallocNanActive_ = true;
        #endif

        #if defined(darwin64)
        malloc_zone_t *zone = malloc_default_zone();

        if (zone != NULL)
        {
            oldMalloc_ = zone->malloc;
            if
            (
                mprotect(zone, getpagesize(), PROT_READ | PROT_WRITE) == 0
            )
            {
                zone->malloc = nanMalloc_;
                if
                (
                    mprotect(zone, getpagesize(), PROT_READ) == 0
                )
                {
                    mallocNanActive_ = true;
                }
            }
        }
        #endif

        if (verbose)
        {
            if (mallocNanActive_)
            {
                Info<< "SetNaN : Initialising allocated memory to NaN"
                    << " (FOAM_SETNAN)." << endl;
            }
            else
            {
                Info<< "SetNaN : Initialise allocated memory to NaN"
                    << " - not supported on this platform" << endl;
            }
        }
    }
}


// ************************************************************************* //
