/* This file is part of the Palabos library.
 *
 * Copyright (C) 2011-2017 FlowKit Sarl
 * Route d'Oron 2
 * 1010 Lausanne, Switzerland
 * E-mail contact: contact@flowkit.com
 *
 * The most recent release of Palabos can be downloaded at
 * <http://www.palabos.org/>
 *
 * The library Palabos is free software: you can redistribute it and/or
 * modify it under the terms of the GNU Affero General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * The library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

/* Main author: Orestis Malaspinas
 */

/** \file
 * 2D specialization of externalForceTemplates functions
 */
#ifndef EXTERNAL_FORCE_TEMPLATES_2D_H
#define EXTERNAL_FORCE_TEMPLATES_2D_H

#include "core/globalDefs.h"

namespace plb {

template<typename T>
struct externalForceTemplatesImpl<T, descriptors::ForcedD2Q9Descriptor<T> >
{

/// Vlasov dynamics à terminer et tester et enlever ce qui est inutile
static void addGravityLikeForce (
                Array<T,descriptors::ForcedD2Q9Descriptor<T>::q>& f,
                T* externalScalars )
{
    static const int forceBeginsAt
        = descriptors::ForcedD2Q9Descriptor<T>::ExternalField::forceBeginsAt;
    T* force = externalScalars + forceBeginsAt;
    T& fx = force[0];
    T& fy = force[1];

    // peut engendrer des fonctions de distribution négatives : attention mais semble normal
    // a comprendre si après multiplié par DeltaT voir bonnefoy v1.52
    f[1] += (-fx + fy) * (T)(1./12.);
    f[2] += (-fx     ) * (T)(1./3.);
    f[3] += (-fx - fy) * (T)(1./12.);
    f[4] += (    - fy) * (T)(1./3.);
    f[5] += ( fx - fy) * (T)(1./12.);
    f[6] += ( fx     ) * (T)(1./3.);
    f[7] += ( fx + fy) * (T)(1./12.);
    f[8] += (      fy) * (T)(1./3.);

//    T gx=fx;
//    T gy=fy;

//    f[1] += 0.0;
//    f[2] += -2.0;
//    f[3] += 0.0;
//    f[4] += -2.0;
//    f[5] += 0.0;
//    f[6] += +2.0;
//    f[7] += 0.0;
//    f[8] += +2.0;

//    // checking
//    for(int i=0; i<9; i++)
//    {
//        if (f[i]<0)
//        {
//            std::cout<< "neg" << std::endl;
//        }
//    }
}

// implementation of naive force : pas de dependance a u, lineaire au coefficient du descripteur et aux vitesses ci
// see bonnefoy v1.2 ou v1.52 pour explications
static void addNaiveForce (
                Array<T,descriptors::ForcedD2Q9Descriptor<T>::q>& f,
                T* externalScalars )
{
    static const int forceBeginsAt
        = descriptors::ForcedD2Q9Descriptor<T>::ExternalField::forceBeginsAt;
    T* force = externalScalars + forceBeginsAt;
    T& fx = force[0];
    T& fy = force[1];

    // voir bonnefoy v1.52 page 28
    // f[1]=Wi*(force_i . ci (au carré))/cs (au carré) avec ci = sqrt(krt) . ( composantes schéma ) et avec cs = sqrt(rt)
    // bizzare semble etre fi=Wi*(Fi.ei)*k avec ei=composante ci sans sqrt(krt)
    // f[0]=0;
    f[1] += (-fx + fy) * (T)(1./12.);
    f[2] += (-fx     ) * (T)(1./3.);
    f[3] += (-fx - fy) * (T)(1./12.);
    f[4] += (    - fy) * (T)(1./3.);
    f[5] += ( fx - fy) * (T)(1./12.);
    f[6] += ( fx     ) * (T)(1./3.);
    f[7] += ( fx + fy) * (T)(1./12.);
    f[8] += (      fy) * (T)(1./3.);
}

// implementation of guo force term
// voir reférénce chapter 20 "multicomponent lattice boltzmann models for biological applications" de A. montessori et Guo 2002
static void addGuoForce (
                Array<T,descriptors::ForcedD2Q9Descriptor<T>::q>& f,
                T* externalScalars,
                Array<T,descriptors::ForcedD2Q9Descriptor<T>::d> const& u, T omega, T amplitude )
{
    static const int forceBeginsAt
        = descriptors::ForcedD2Q9Descriptor<T>::ExternalField::forceBeginsAt;
    T* force = externalScalars + forceBeginsAt; // sauf erreur ici la force semble etre en newton/m3
    T mu = amplitude*((T)1-omega/(T)2);  // Guo force amplitude

    static const T oneOver3 = (T)1/(T)3; // coefficients cs
    static const T oneOver12 = (T)1/(T)12; // coefficients cs
    static const T fourOver3 = (T)4/(T)3; // coefficients cs

    const T twoUx   = (T)2*u[0];
    const T threeUx = (T)3*u[0];

    const T twoUy   = (T)2*u[1];
    const T threeUy = (T)3*u[1];

    // fi distribution en plus depends de Force F, u et omega
    // f[0]=amplitude*Wi*(1-w/2)*(-u/cs^²).Fi avec Wi= 4/9 et cs^2=1/3=k dans cette formule ci avant
    // amplitude = 1 mais il y a un facteur qui semble etre delta_t/delta_x voir pdf bonnefoy v1.52 page 28 soit inv(sqrt(krt)) mais ici dans palabos amplitude=1
    f[0] += -fourOver3*mu*(force[0]*u[0]+force[1]*u[1]); // influence force terme sur fonction de distribution f0 centrale du d2q9

    f[1] += oneOver12*mu*(force[0]*(-(T)1+twoUx-threeUy)+force[1]*((T)1+twoUy-threeUx));

    f[2] += oneOver3*mu*(force[0]*(-(T)1+twoUx)-force[1]*u[1]);

    f[3] += oneOver12*mu*(force[0]*(-(T)1+twoUx+threeUy)+force[1]*(-(T)1+twoUy+threeUx));

    f[4] += -oneOver3*mu*(force[0]*u[0]+force[1]*((T)1-twoUy));

    f[5] += oneOver12*mu*(force[0]*((T)1+twoUx-threeUy)+force[1]*(-(T)1+twoUy-threeUx));

    f[6] += oneOver3*mu*(force[0]*((T)1+twoUx)-force[1]*u[1]);

    f[7] += oneOver12*mu*(force[0]*((T)1+twoUx+threeUy)+force[1]*((T)1+twoUy+threeUx));

    f[8] += -oneOver3*mu*(force[0]*u[0]+force[1]*(-(T)1-twoUy));

//    // checking
//    for(int i=0; i<9; i++)
//    {
//        if (f[i]<0)
//        {
//            std::cout<< "neg" << std::endl;
//        }
//    }
}

static T heForcedBGKCollision (
        Array<T,descriptors::ForcedD2Q9Descriptor<T>::q>& f, T* externalScalars,
        T rhoBar, Array<T,descriptors::ForcedD2Q9Descriptor<T>::d> const& j, T omega )
{
    static const int forceBeginsAt
            = descriptors::ForcedD2Q9Descriptor<T>::ExternalField::forceBeginsAt;
    T* force = externalScalars + forceBeginsAt;
    T invRho = descriptors::ForcedD2Q9Descriptor<T>::invRho(rhoBar);

    Array<T,descriptors::ForcedD2Q9Descriptor<T>::d> uLB(j[0]*invRho, j[1]*invRho);
    const T uSqrLB = VectorTemplateImpl<T,descriptors::ForcedD2Q9Descriptor<T>::d>::normSqr(uLB);

    dynamicsTemplatesImpl<T,descriptors::D2Q9DescriptorBase<T> >::bgk_ma2_collision(f, rhoBar, j, omega);

    for (plint iPop=0; iPop < descriptors::ForcedD2Q9Descriptor<T>::q; ++iPop) {
        T ug = (descriptors::ForcedD2Q9Descriptor<T>::c[iPop][0]-uLB[0])*force[0] +
                (descriptors::ForcedD2Q9Descriptor<T>::c[iPop][1]-uLB[1])*force[1];

        f[iPop] += descriptors::ForcedD2Q9Descriptor<T>::invCs2 * ug * dynamicsTemplatesImpl<T,descriptors::ForcedD2Q9Descriptor<T> >
            ::bgk_ma2_equilibrium( iPop, (T)1, (T)1, uLB, uSqrLB );
    }
    T uSqr = util::sqr(uLB[0] + 0.5*force[0]) +
             util::sqr(uLB[1] + 0.5*force[1]);
    return uSqr;
}

// static T heForcedBGKCollision (
//         Array<T,descriptors::ForcedD2Q9Descriptor<T>::q>& f, T* externalScalars,
//         T rhoBar, Array<T,descriptors::ForcedD2Q9Descriptor<T>::d> const& j, T omega )
// {
// 	static const int forceBeginsAt
// 			= descriptors::ForcedD2Q9Descriptor<T>::ExternalField::forceBeginsAt;
// 	T* force = externalScalars + forceBeginsAt;
// 	T invRho = descriptors::ForcedD2Q9Descriptor<T>::invRho(rhoBar);
// 	const T jSqr = VectorTemplateImpl<T,descriptors::ForcedD2Q9Descriptor<T>::d>::normSqr(j);
//
// 	Array<T,descriptors::ForcedD2Q9Descriptor<T>::d> uLB(j[0]*invRho, j[1]*invRho);
// 	const T uSqrLB = VectorTemplateImpl<T,descriptors::ForcedD2Q9Descriptor<T>::d>::normSqr(uLB);
//
// 	for (plint iPop=0; iPop < descriptors::ForcedD2Q9Descriptor<T>::q; ++iPop) {
// 		T ug = (descriptors::ForcedD2Q9Descriptor<T>::c[iPop][0]-invRho*j[0])*force[0] +
// 				(descriptors::ForcedD2Q9Descriptor<T>::c[iPop][1]-invRho*j[1])*force[1];
// 		f[iPop] *= (T)1-omega;
// 		f[iPop] += omega * dynamicsTemplatesImpl<T,descriptors::ForcedD2Q9Descriptor<T> >
// 			::bgk_ma2_equilibrium( iPop, rhoBar, invRho, j, jSqr );
//
// 		f[iPop] += descriptors::ForcedD2Q9Descriptor<T>::invCs2 * ug * dynamicsTemplatesImpl<T,descriptors::ForcedD2Q9Descriptor<T> >
// 			::bgk_ma2_equilibrium( iPop, (T)1, (T)1, uLB, uSqrLB );
// 	}
// 	T uSqr = util::sqr(uLB[0] + 0.5*force[0]) +
// 			 util::sqr(uLB[1] + 0.5*force[1]);
// 	return uSqr;
// }

};

}  // namespace plb

#endif
