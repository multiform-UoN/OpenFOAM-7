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

#ifndef doubleFloat_H
#define doubleFloat_H

#include "label.H"
#include "products.H"

#include <cmath>

#if defined(darwin64)
inline float j0f(float x) { return float(j0(double(x)));}
inline float j1f(float x) { return float(j1(double(x)));}
inline float y0f(float x) { return float(y0(double(x)));}
inline float y1f(float x) { return float(y1(double(x)));}
inline float jnf(const int n, const float s) { return float(jn(n, double(s))); }
inline float ynf(const int n, const float s) { return float(yn(n, double(s))); }

inline long double j0l(long double x)
{
    return static_cast<long double>(j0(static_cast<double>(x)));
}

inline long double j1l(long double x)
{
    return static_cast<long double>(j1(static_cast<double>(x)));
}

inline long double y0l(long double x)
{
    return static_cast<long double>(y0(static_cast<double>(x)));
}

inline long double y1l(float x)
{
    return static_cast<long double>(y1(static_cast<double>(x)));
}

inline long double jnl(const int n, const long double s)
{
    return static_cast<long double>(jn(n, static_cast<double>(s)));
}

inline long double ynl(const int n, const long double s)
{
    return static_cast<long double>(yn(n, static_cast<double>(s)));
}

#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Cmpt>
class typeOfRank<Cmpt, 0>
{
public:

    typedef Cmpt type;
};


template<class Cmpt>
class symmTypeOfRank<Cmpt, 0>
{
public:

    typedef Cmpt type;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class T>
inline bool equal(const T& s1, const T& s2)
{
    return s1 == s2;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
