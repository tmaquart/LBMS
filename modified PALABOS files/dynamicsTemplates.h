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

/** \file
 * Helper functions for the implementation of LB dynamics. This file is all
 * about efficiency. The generic template code is specialized for commonly
 * used Lattices, so that a maximum performance can be taken out of each
 * case.
 */
#ifndef DYNAMICS_TEMPLATES_H
#define DYNAMICS_TEMPLATES_H

#include "core/globalDefs.h"
#include "core/cell.h"
#include "core/util.h"
#include "latticeBoltzmann/offEquilibriumTemplates.h"

namespace plb {

template<typename T, class Descriptor> struct dynamicsTemplatesImpl;

/// This structure forwards the calls to the appropriate helper class
template<typename T, template<typename U> class Descriptor>
struct dynamicsTemplates {

static T bgk_ma2_equilibrium(plint iPop, T rhoBar, T invRho, Array<T,Descriptor<T>::d> const& j, T jSqr) {
    return dynamicsTemplatesImpl<T,typename Descriptor<T>::BaseDescriptor>
        ::bgk_ma2_equilibrium(iPop, rhoBar, invRho, j, jSqr);
}

/// Solid equilibrium à terminer et tester et enlever ce qui est inutile
static T bgk_ma2_equilibrium_solid_hermite(plint iPop, T rhoBar, T invRho, Array<T,Descriptor<T>::d> const& j, T jSqr) {
    return dynamicsTemplatesImpl<T,typename Descriptor<T>::BaseDescriptor>
        ::bgk_ma2_equilibrium_solid_hermite(iPop, rhoBar, invRho, j, jSqr);
}

/// Solid equilibrium à terminer et tester et enlever ce qui est inutile
static T bgk_ma2_equilibrium_solid(plint iPop, T rhoBar, T invRho, Array<T,Descriptor<T>::d> const& j, T jSqr) {
    return dynamicsTemplatesImpl<T,typename Descriptor<T>::BaseDescriptor>
        ::bgk_ma2_equilibrium_solid(iPop, rhoBar, invRho, j, jSqr);
}

static T complete_bgk_ma2_equilibrium(plint iPop, T rhoBar, T invRho, Array<T,Descriptor<T>::d> const& j, T jSqr) {
    return dynamicsTemplatesImpl<T,typename Descriptor<T>::BaseDescriptor>
        ::complete_bgk_ma2_equilibrium(iPop, rhoBar, invRho, j, jSqr);
}

static void complete_bgk_ma2_regularize(Cell<T,Descriptor> &cell, T rhoBar, T invRho, Array<T,Descriptor<T>::d> const& j, T jSqr,
                                     Array<T,SymmetricTensor<T,Descriptor>::n> const& piNeq, T omega = (T)1, T omegaNonPhys = (T)1 )
{
    dynamicsTemplatesImpl<T,typename Descriptor<T>::BaseDescriptor>
        ::complete_bgk_ma2_regularize(cell.getRawPopulations(), rhoBar, invRho, j, jSqr, piNeq, omega, omegaNonPhys, Descriptor<T>::d+SymmetricTensor<T,Descriptor>::n);
}

static void bgk_ma2_equilibria( T rhoBar, T invRho, Array<T,Descriptor<T>::d> const& j,
                                T jSqr, Array<T,Descriptor<T>::q>& eqPop )
{
    dynamicsTemplatesImpl<T,typename Descriptor<T>::BaseDescriptor>
        ::bgk_ma2_equilibria(rhoBar, invRho, j, jSqr, eqPop);
}

static void complete_bgk_ma2_equilibria( T rhoBar, T invRho, Array<T,Descriptor<T>::d> const& j,
                                T jSqr, Array<T,Descriptor<T>::q>& eqPop )
{
    dynamicsTemplatesImpl<T,typename Descriptor<T>::BaseDescriptor>
        ::complete_bgk_ma2_equilibria(rhoBar, invRho, j, jSqr, eqPop);
}

static T bgk_ma2_collision(Cell<T,Descriptor>& cell, T rhoBar, Array<T,Descriptor<T>::d> const& j, T omega)
{
    return dynamicsTemplatesImpl<T,typename Descriptor<T>::BaseDescriptor>
        ::bgk_ma2_collision(cell.getRawPopulations(), rhoBar, j, omega);
}

/// Solid dynamics à terminer et tester et enlever ce qui est inutile
static T bgk_ma2_collision_solid_hermite(Cell<T,Descriptor>& cell, T rhoBar, Array<T,Descriptor<T>::d> const& j, T omega, T rhoMaterial)
{
    return dynamicsTemplatesImpl<T,typename Descriptor<T>::BaseDescriptor>
        ::bgk_ma2_collision_solid_hermite(cell.getRawPopulations(), rhoBar, j, omega, rhoMaterial);
}

/// Solid dynamics à terminer et tester et enlever ce qui est inutile
static T bgk_ma2_collision_solid(Cell<T,Descriptor>& cell, T rhoBar, Array<T,Descriptor<T>::d> const& j, T omega, T rhoMaterial)
{
    return dynamicsTemplatesImpl<T,typename Descriptor<T>::BaseDescriptor>
        ::bgk_ma2_collision_solid(cell.getRawPopulations(), rhoBar, j, omega, rhoMaterial);
}

/// Vlasov dynamics à terminer et tester et enlever ce qui est inutile
static T bgk_ma2_collision_vlasov(Cell<T,Descriptor>& cell, T rhoBar, Array<T,Descriptor<T>::d> const& j, T omega)
{
    return dynamicsTemplatesImpl<T,typename Descriptor<T>::BaseDescriptor>
        ::bgk_ma2_collision_vlasov(cell.getRawPopulations(), rhoBar, j, omega); // implementation of bgk_ma2_collision_vlasov ci-après ou selon descripteur dans dynamicsTemplatesXD.h
        // les classes 2D et 3D sont adaptées à certains descripteurs
}

static T complete_bgk_ma2_collision(Cell<T,Descriptor>& cell, T rhoBar, T invRho, Array<T,Descriptor<T>::d> const& j, T omega)
{
    return dynamicsTemplatesImpl<T,typename Descriptor<T>::BaseDescriptor>
        ::complete_bgk_ma2_collision(cell.getRawPopulations(), rhoBar, invRho, j, omega);
}

static T complete_regularized_bgk_ma2_collision(Cell<T,Descriptor>& cell, T rhoBar, const Array<T,Descriptor<T>::d> &j, const Array<T,SymmetricTensor<T,Descriptor>::n> &piNeq, T omega)
{
    return dynamicsTemplatesImpl<T,typename Descriptor<T>::BaseDescriptor>
        ::complete_regularized_bgk_ma2_collision(cell.getRawPopulations(), rhoBar, j, piNeq, omega);
}

static T complete_mrt_ma2_collision(Cell<T,Descriptor>& cell, plint order, T omega, T omegaNonPhys)
{
    return dynamicsTemplatesImpl<T,typename Descriptor<T>::BaseDescriptor>
        ::complete_mrt_ma2_collision(cell.getRawPopulations(), order, omega, omegaNonPhys, Descriptor<T>::d+SymmetricTensor<T,Descriptor>::n);
}

static T complete_regularized_mrt_ma2_collision(Cell<T,Descriptor>& cell, T rhoBar, const Array<T,Descriptor<T>::d> &j, const Array<T,SymmetricTensor<T,Descriptor>::n> &piNeq, plint order, T omega, T omegaNonPhys)
{
    return dynamicsTemplatesImpl<T,typename Descriptor<T>::BaseDescriptor>
        ::complete_regularized_mrt_ma2_collision(cell.getRawPopulations(), rhoBar, j, piNeq, order, omega, omegaNonPhys, Descriptor<T>::d+SymmetricTensor<T,Descriptor>::n);
}

static T consistent_smagorinsky_complete_regularized_mrt_ma2_collision(Cell<T,Descriptor>& cell, plint order,  T cSmago, T omega, T omegaNonPhys)
{
    return dynamicsTemplatesImpl<T,typename Descriptor<T>::BaseDescriptor>
        ::consistent_smagorinsky_complete_regularized_mrt_ma2_collision(cell.getRawPopulations(), order, cSmago,  omega, omegaNonPhys, Descriptor<T>::d+SymmetricTensor<T,Descriptor>::n);
}

static T truncated_mrt_ma2_collision(Cell<T,Descriptor>& cell, T omega, T omegaNonPhys)
{
    return dynamicsTemplatesImpl<T,typename Descriptor<T>::BaseDescriptor>
        ::truncated_mrt_ma2_collision(cell.getRawPopulations(), omega, omegaNonPhys, Descriptor<T>::d+SymmetricTensor<T,Descriptor>::n);
}

static T complete_mrt_smagorinsky_ma2_collision(Cell<T,Descriptor>& cell, T cSmago, T omega, T omegaNonPhys)
{
    return dynamicsTemplatesImpl<T,typename Descriptor<T>::BaseDescriptor>
        ::complete_mrt_smagorinsky_ma2_collision(cell.getRawPopulations(), cSmago, omega, omegaNonPhys, Descriptor<T>::d+SymmetricTensor<T,Descriptor>::n);
}

static T complete_mrt_smagorinsky_ma2_ext_rhoBar_j_collision(Cell<T,Descriptor>& cell, T rhoBar, Array<T,Descriptor<T>::d> const& j, plint order, T cSmago, T omega, T omegaNonPhys)
{
    return dynamicsTemplatesImpl<T,typename Descriptor<T>::BaseDescriptor>
        ::complete_mrt_smagorinsky_ma2_ext_rhoBar_j_collision(cell.getRawPopulations(), rhoBar, j, order, cSmago, omega, omegaNonPhys, Descriptor<T>::d+SymmetricTensor<T,Descriptor>::n);
}

static T truncated_mrt_smagorinsky_ma2_collision(Cell<T,Descriptor>& cell, T cSmago, T omega, T omegaNonPhys)
{
    return dynamicsTemplatesImpl<T,typename Descriptor<T>::BaseDescriptor>
        ::truncated_mrt_smagorinsky_ma2_collision(cell.getRawPopulations(), cSmago, omega, omegaNonPhys, Descriptor<T>::d+SymmetricTensor<T,Descriptor>::n);
}

static T truncated_mrt_smagorinsky_ma2_ext_rhoBar_j_collision(Cell<T,Descriptor>& cell, T rhoBar, Array<T,Descriptor<T>::d> const& j, T cSmago, T omega, T omegaNonPhys)
{
    return dynamicsTemplatesImpl<T,typename Descriptor<T>::BaseDescriptor>
        ::truncated_mrt_smagorinsky_ma2_ext_rhoBar_j_collision(cell.getRawPopulations(), rhoBar, j, cSmago, omega, omegaNonPhys, Descriptor<T>::d+SymmetricTensor<T,Descriptor>::n);
}

static T complete_mrt_ma2_ext_rhoBar_j_collision(Cell<T,Descriptor>& cell, T rhoBar, Array<T,Descriptor<T>::d> const& j, plint order, T omega, T omegaNonPhys)
{
    return dynamicsTemplatesImpl<T,typename Descriptor<T>::BaseDescriptor>
        ::complete_mrt_ma2_ext_rhoBar_j_collision(cell.getRawPopulations(), rhoBar, j, order, omega, omegaNonPhys, Descriptor<T>::d+SymmetricTensor<T,Descriptor>::n);
}

static T truncated_mrt_ma2_ext_rhoBar_j_collision(Cell<T,Descriptor>& cell, T rhoBar, Array<T,Descriptor<T>::d> const& j, T omega, T omegaNonPhys)
{
    return dynamicsTemplatesImpl<T,typename Descriptor<T>::BaseDescriptor>
        ::truncated_mrt_ma2_ext_rhoBar_j_collision(cell.getRawPopulations(), rhoBar, j, omega, omegaNonPhys, Descriptor<T>::d+SymmetricTensor<T,Descriptor>::n);
}

static T computePsiComplete(T omega) {
    return dynamicsTemplatesImpl<T,typename Descriptor<T>::BaseDescriptor>::computePsiComplete(omega);
}

static T computePsiTruncated(T omega) {
    return dynamicsTemplatesImpl<T,typename Descriptor<T>::BaseDescriptor>::computePsiTruncated(omega);
}

static T bgk_inc_collision(Cell<T,Descriptor>& cell, T rhoBar, Array<T,Descriptor<T>::d> const& j, T omega )
{
    return dynamicsTemplatesImpl<T,typename Descriptor<T>::BaseDescriptor>
        ::bgk_inc_collision(cell.getRawPopulations(), rhoBar, j, omega);
}

static T rlb_collision(Cell<T,Descriptor>& cell, T rhoBar, T invRho, Array<T,Descriptor<T>::d> const& j,
                       Array<T,SymmetricTensor<T,Descriptor>::n> const& PiNeq, T omega )
{
    return dynamicsTemplatesImpl<T,typename Descriptor<T>::BaseDescriptor>
        ::rlb_collision(cell.getRawPopulations(), rhoBar, invRho, j, PiNeq, omega);
}

static T bgk_ma2_constRho_collision(Cell<T,Descriptor>& cell,
                                T rhoBar, Array<T,Descriptor<T>::d> const& j, T ratioRho, T omega)
{
    return dynamicsTemplatesImpl<T,typename Descriptor<T>::BaseDescriptor>
               ::bgk_ma2_constRho_collision(cell.getRawPopulations(), rhoBar, j, ratioRho, omega);
}

static T precond_bgk_ma2_equilibrium(plint iPop, T rhoBar, T invRho, Array<T,Descriptor<T>::d> const& j, T jSqr, T invGamma)
{
    return dynamicsTemplatesImpl<T,typename Descriptor<T>::BaseDescriptor>
        ::precond_bgk_ma2_equilibrium(iPop, rhoBar, invRho, j, jSqr, invGamma);
}

static T precond_bgk_ma2_collision(Cell<T,Descriptor>& cell, T rhoBar, Array<T,Descriptor<T>::d> const& j, T omega, T invGamma)
{
    return dynamicsTemplatesImpl<T,typename Descriptor<T>::BaseDescriptor>
        ::precond_bgk_ma2_collision(cell.getRawPopulations(), rhoBar, j, omega, invGamma);
}

static T bgk_ma3_equilibrium(plint iPop, T rhoBar, T invRho, Array<T,Descriptor<T>::d> const& j, T jSqr, T rhoThetaBar) {
    return dynamicsTemplatesImpl<T,typename Descriptor<T>::BaseDescriptor>
        ::bgk_ma3_equilibrium(iPop, rhoBar, invRho, j, jSqr, rhoThetaBar);
}

static T bgk_ma4_equilibrium(plint iPop, T rhoBar, T invRho, Array<T,Descriptor<T>::d> const& j, T jSqr, T rhoThetaBar) {
    return dynamicsTemplatesImpl<T,typename Descriptor<T>::BaseDescriptor>
        ::bgk_ma4_equilibrium(iPop, rhoBar, invRho, j, jSqr, rhoThetaBar);
}

static T bgk_ma3_collision(Cell<T,Descriptor>& cell, T rhoBar, Array<T,Descriptor<T>::d> const& j, T rhoThetaBar, T omega)
{
    return dynamicsTemplatesImpl<T,typename Descriptor<T>::BaseDescriptor>
        ::bgk_ma3_collision(cell.getRawPopulations(), rhoBar, j, rhoThetaBar, omega);
}

static T bgk_ma4_collision(Cell<T,Descriptor>& cell, T rhoBar, Array<T,Descriptor<T>::d> const& j, T rhoThetaBar, T omega)
{
    return dynamicsTemplatesImpl<T,typename Descriptor<T>::BaseDescriptor>
        ::bgk_ma4_collision(cell.getRawPopulations(), rhoBar, j, rhoThetaBar, omega);
}

};  // struct dynamicsTemplates

/// All helper functions are inside this structure
template<typename T, class Descriptor>
struct dynamicsTemplatesImpl {

static T bgk_ma2_equilibrium(plint iPop, T rhoBar, T invRho, Array<T,Descriptor::d> const& j, T jSqr) {
    T c_j = Descriptor::c[iPop][0]*j[0];
    for (int iD=1; iD < Descriptor::d; ++iD) {
       c_j += Descriptor::c[iPop][iD]*j[iD];
    }
    return Descriptor::t[iPop] * (
           rhoBar + Descriptor::invCs2 * c_j +
           Descriptor::invCs2/(T)2 * invRho * (
               Descriptor::invCs2 * c_j*c_j - jSqr )
       );
}

/// solid dynamics à terminer et tester et enlever ce qui est inutile
static T bgk_ma2_equilibrium_solid_hermite(plint iPop, T rhoBar, T invRho, Array<T,Descriptor::d> const& j, T jSqr) {
    T c_j = Descriptor::c[iPop][0]*j[0];
    for (int iD=1; iD < Descriptor::d; ++iD) {
       c_j += Descriptor::c[iPop][iD]*j[iD];
    }
    return Descriptor::t[iPop] * (
           rhoBar + Descriptor::invCs2 * c_j +
           Descriptor::invCs2/(T)2 * invRho * (
               Descriptor::invCs2 * c_j*c_j - jSqr )
       );
}

/// solid dynamics à terminer et tester et enlever ce qui est inutile
static T bgk_ma2_equilibrium_solid(plint iPop, T rhoBar, T invRho, Array<T,Descriptor::d> const& j, T jSqr) {
    T c_j = Descriptor::c[iPop][0]*j[0];
    for (int iD=1; iD < Descriptor::d; ++iD) {
       c_j += Descriptor::c[iPop][iD]*j[iD];
    }
    return Descriptor::t[iPop] * (
           rhoBar + Descriptor::invCs2 * c_j +
           Descriptor::invCs2/(T)2 * invRho * (
               Descriptor::invCs2 * c_j*c_j - jSqr )
       );
}

static T complete_bgk_ma2_equilibrium(plint iPop, T rhoBar, T invRho, Array<T,Descriptor::d> const& j, T jSqr) {
    PLB_ASSERT(false);
    return T();
}

static void bgk_ma2_equilibria( T rhoBar, T invRho, Array<T,Descriptor::d>
        const& j,
                                T jSqr, Array<T,Descriptor::q>& eqPop )
{
    for (int iPop=0; iPop<Descriptor::q; ++iPop) {
        eqPop[iPop] = bgk_ma2_equilibrium(iPop, rhoBar, invRho, j, jSqr);
    }
}

static void complete_bgk_ma2_equilibria( T rhoBar, T invRho, Array<T,Descriptor::d> const& j,
                                         T jSqr, Array<T,Descriptor::q>& eqPop )
{
    for (int iPop=0; iPop<Descriptor::q; ++iPop) {
        eqPop[iPop] = complete_bgk_ma2_equilibrium(iPop, rhoBar, invRho, j, jSqr);
    }
}

static void complete_bgk_ma2_regularize(Array<T,Descriptor::q> &f, T rhoBar, T invRho, Array<T,Descriptor::d> const& j, T jSqr,
                                        Array<T,SymmetricTensorImpl<T,Descriptor::d>::n> const& piNeq, T omega, T omegaNonPhys)
{
    PLB_ASSERT(false);
}

static T bgk_ma2_collision(Array<T,Descriptor::q>& f, T rhoBar, Array<T,Descriptor::d> const& j, T omega) {
    T invRho = Descriptor::invRho(rhoBar);
    const T jSqr = VectorTemplateImpl<T,Descriptor::d>::normSqr(j);
    for (plint iPop=0; iPop < Descriptor::q; ++iPop) {
        f[iPop] *= (T)1-omega;
        f[iPop] += omega * dynamicsTemplatesImpl<T,Descriptor>::bgk_ma2_equilibrium (
                                iPop, rhoBar, invRho, j, jSqr );
    }
    return jSqr*invRho*invRho;
}

/// solid dynamics à terminer et tester et enlever ce qui est inutile
static T bgk_ma2_collision_solid_hermite(Array<T,Descriptor::q>& f, T rhoBar, Array<T,Descriptor::d> const& j, T omega, T rhoMaterial) {
    T invRho = Descriptor::invRho(rhoBar);
    const T jSqr = VectorTemplateImpl<T,Descriptor::d>::normSqr(j);
    // note: only for D2Q9
    plb::Array<T,9> Whermite;
    Whermite[0]=4.0/9.0;
    Whermite[1]=1.0/36.0;
    Whermite[2]=1.0/9.0;
    Whermite[3]=1.0/36.0;
    Whermite[4]=1.0/9.0;
    Whermite[5]=1.0/36.0;
    Whermite[6]=1.0/9.0;
    Whermite[7]=1.0/36.0;
    Whermite[8]=1.0/9.0;
    for (plint iPop=0; iPop < Descriptor::q; ++iPop) {
        f[iPop] *= ((T)1.0-omega);
        f[iPop] += omega*Whermite[iPop]*(rhoMaterial-1.0);
    }
    return jSqr*invRho*invRho;
}

/// solid dynamics à terminer et tester et enlever ce qui est inutile
static T bgk_ma2_collision_solid(Array<T,Descriptor::q>& f, T rhoBar, Array<T,Descriptor::d> const& j, T omega, T rhoMaterial) {
    T invRho = Descriptor::invRho(rhoBar);
    const T jSqr = VectorTemplateImpl<T,Descriptor::d>::normSqr(j);
    for (plint iPop=0; iPop < Descriptor::q; ++iPop) {
        f[iPop] *= ((T)1.0-omega);
        f[iPop] += omega*((rhoMaterial-1.0)/9.0);
    }
    return jSqr*invRho*invRho;
}

/// Vlasov dynamics à terminer et tester et enlever ce qui est inutile
static T bgk_ma2_collision_vlasov(Array<T,Descriptor::q>& f, T rhoBar, Array<T,Descriptor::d> const& j, T omega) {
    T invRho = Descriptor::invRho(rhoBar);
    const T jSqr = VectorTemplateImpl<T,Descriptor::d>::normSqr(j);
    for (plint iPop=0; iPop < Descriptor::q; ++iPop) {
        f[iPop] *= 1.0; // zero collision
    }
    return jSqr*invRho*invRho;
}

static T complete_bgk_ma2_collision(Array<T,Descriptor::q>& f, T rhoBar, T invRho, Array<T,Descriptor::d> const& j, T omega) {
    PLB_ASSERT(false);
    return T();
}

static T complete_mrt_ma2_collision(Array<T,Descriptor::q>& f, T omega, T omegaNonPhys) {
    PLB_ASSERT(false);
    return T();
}

static T truncated_mrt_ma2_collision(Array<T,Descriptor::q>& f, T omega, T omegaNonPhys) {
    PLB_ASSERT(false);
    return T();
}

static T complete_mrt_ma2_ext_rhoBar_j_collision(Array<T,Descriptor::q>& f, T rhoBar, Array<T,Descriptor::d> const &j, T omega, T omegaNonPhys, plint iPhys) {
    PLB_ASSERT(false);
    return T();
}

static T truncated_mrt_ma2_ext_rhoBar_j_collision(Array<T,Descriptor::q>& f, T rhoBar, Array<T,Descriptor::d> const &j, T omega, T omegaNonPhys, plint iPhys) {
    PLB_ASSERT(false);
    return T();
}

static T complete_mrt_smagorinsky_ma2_collision(Array<T,Descriptor::q>& f, T cSmago, T omega, T omegaNonPhys, plint iPhys) {
    PLB_ASSERT(false);
    return T();
}

static T complete_mrt_smagorinsky_ma2_ext_rhoBar_j_collision(Array<T,Descriptor::q>& f, T rhoBar, Array<T,Descriptor::d> const& j,
                                                             T cSmago, T omega, T omegaNonPhys, plint iPhys)
{
    PLB_ASSERT(false);
    return T();
}

static T truncated_mrt_smagorinsky_ma2_collision(Array<T,Descriptor::q>& f, T cSmago, T omega, T omegaNonPhys, plint iPhys) {
    PLB_ASSERT(false);
    return T();
}

static T truncated_mrt_smagorinsky_ma2_ext_rhoBar_j_collision(Array<T,Descriptor::q>& f, T rhoBar, Array<T,Descriptor::d> const& j,
                                                             T cSmago, T omega, T omegaNonPhys, plint iPhys)
{
    PLB_ASSERT(false);
    return T();
}

static T computePsiComplete(T omega) {
    PLB_ASSERT(false);
    return T();
}

static T computePsiTruncated(T omega) {
    PLB_ASSERT(false);
    return T();
}

static T bgk_inc_collision(Array<T,Descriptor::q>& f, T rhoBar, Array<T,Descriptor::d> const& j, T omega ) {
    // In incompressible BGK, the Ma^2 term is preceeded by 1/rho0 instead of 1/rho.
    // Incompressible: rho0=1
    static const T invRho0 = (T) 1;
    const T jSqr = VectorTemplateImpl<T,Descriptor::d>::normSqr(j);
    for (plint iPop=0; iPop < Descriptor::q; ++iPop) {
        f[iPop] *= (T)1-omega;
        f[iPop] += omega * dynamicsTemplatesImpl<T,Descriptor>::bgk_ma2_equilibrium (
                                iPop, rhoBar, invRho0, j, jSqr );
    }
    return jSqr;
}

static T rlb_collision (
    Array<T,Descriptor::q>& f, T rhoBar, T invRho, Array<T,Descriptor::d> const& j,
    Array<T,SymmetricTensorImpl<T,Descriptor::d>::n> const& PiNeq, T omega )
{
    const T jSqr = VectorTemplateImpl<T,Descriptor::d>::normSqr(j);
    f[0] = dynamicsTemplatesImpl<T,Descriptor>::bgk_ma2_equilibrium(0, rhoBar, invRho, j, jSqr)
                  + ((T)1-omega) *
                        offEquilibriumTemplatesImpl<T,Descriptor>::fromPiToFneq(0, PiNeq);
    for (plint iPop=1; iPop<=Descriptor::q/2; ++iPop) {
        f[iPop] = dynamicsTemplatesImpl<T,Descriptor>::
                bgk_ma2_equilibrium(iPop, rhoBar, invRho, j, jSqr);
        f[iPop+Descriptor::q/2] = dynamicsTemplatesImpl<T,Descriptor>::
                bgk_ma2_equilibrium(iPop+Descriptor::q/2, rhoBar, invRho, j, jSqr);
        T fNeq = ((T)1-omega) *
                 offEquilibriumTemplatesImpl<T,Descriptor>::fromPiToFneq(iPop, PiNeq);
        f[iPop] += fNeq;
        f[iPop+Descriptor::q/2] += fNeq;
    }
    return jSqr*invRho*invRho;
}

static T bgk_ma2_constRho_collision(Array<T,Descriptor::q>& f, T rhoBar, Array<T,Descriptor::d> const& j, T ratioRho, T omega) {
    T invRho = Descriptor::invRho(rhoBar);
    const T jSqr = VectorTemplateImpl<T,Descriptor::d>::normSqr(j);
    for (plint iPop=0; iPop < Descriptor::q; ++iPop) {
        T feq = dynamicsTemplatesImpl<T,Descriptor>::
                    bgk_ma2_equilibrium(iPop, rhoBar, invRho, j, jSqr );
        f[iPop] =
          ratioRho*feq + Descriptor::t[iPop]*(ratioRho-(T)1) +
          ((T)1-omega)*(f[iPop]-feq);
    }
    return jSqr*invRho*invRho;
}

static T precond_bgk_ma2_equilibrium(plint iPop, T rhoBar, T invRho, Array<T,Descriptor::d> const& j, T jSqr, T invGamma) {
    T c_j = Descriptor::c[iPop][0]*j[0];
    for (int iD=1; iD < Descriptor::d; ++iD) {
       c_j += Descriptor::c[iPop][iD]*j[iD];
    }
    return Descriptor::t[iPop] * (
           rhoBar + Descriptor::invCs2 * c_j +
           invGamma * Descriptor::invCs2/(T)2 * invRho * (
               Descriptor::invCs2 * c_j*c_j - jSqr )
       );
}

static T precond_bgk_ma2_collision(Array<T,Descriptor::q>& f, T rhoBar, Array<T,Descriptor::d> const& j, T omega, T invGamma) {
    T invRho = Descriptor::invRho(rhoBar);
    const T jSqr = VectorTemplateImpl<T,Descriptor::d>::normSqr(j);
    for (plint iPop=0; iPop < Descriptor::q; ++iPop) {
        f[iPop] *= (T)1-omega;
        f[iPop] += omega * dynamicsTemplatesImpl<T,Descriptor>::precond_bgk_ma2_equilibrium (
                                iPop, rhoBar, invRho, j, jSqr, invGamma );
    }
    return jSqr*invRho*invRho;
}

static T bgk_ma3_equilibrium(plint iPop, T rhoBar, T invRho, Array<T,Descriptor::d> const& j, T jSqr, T rhoThetaBar) {
    typedef Descriptor L;
    T c_j = (T)L::c[iPop][0]*j[0];
    for (int iD=1; iD < L::d; ++iD) {
        c_j += (T)L::c[iPop][iD]*j[iD];
    }

    T c_u = c_j*invRho;
    T uSqr = jSqr*invRho*invRho;
    T thetaBar = invRho*rhoThetaBar;

    return L::t[iPop] * (
        rhoBar + L::invCs2 * c_j +
        L::invCs2/(T)2 * ( invRho * (L::invCs2 * c_j*c_j - jSqr )
                                    + rhoThetaBar * ((T)L::cNormSqr[iPop]-L::cs2*(T)L::d))
        + L::invCs2*L::invCs2*L::invCs2/(T)6 * c_j * (
            c_u*c_u-(T)3*L::cs2*uSqr+(T)3*L::cs2*thetaBar*(L::cNormSqr[iPop]-L::cs2*((T)L::d+(T)2))
        )
    );
}

static T bgk_ma4_equilibrium(plint iPop, T rhoBar, T invRho, Array<T,Descriptor::d> const& j, T jSqr, T rhoThetaBar) {
    typedef Descriptor L;
    T c_j = (T)L::c[iPop][0]*j[0];
    for (int iD=1; iD < L::d; ++iD) {
        c_j += (T)L::c[iPop][iD]*j[iD];
    }

    T c_u = c_j*invRho;
    T uSqr = jSqr*invRho*invRho;
    T thetaBar = invRho*rhoThetaBar;

    return L::t[iPop] * (
        rhoBar + L::invCs2 * c_j +
        L::invCs2/(T)2 * ( invRho * (L::invCs2 * c_j*c_j - jSqr )
                                    + rhoThetaBar * ((T)L::cNormSqr[iPop]-L::cs2*(T)L::d))
        + L::invCs2*L::invCs2*L::invCs2/(T)6 * c_j * (
            c_u*c_u-(T)3*L::cs2*uSqr+(T)3*L::cs2*thetaBar*(L::cNormSqr[iPop]-L::cs2*((T)L::d+(T)2))
        )
        + L::invCs2*L::invCs2*L::invCs2*L::invCs2/(T)24 * (
            c_j*c_u*c_u*c_u - (T)6*L::cs2*uSqr*c_j*c_u + (T)3*L::cs2*L::cs2*jSqr*invRho*uSqr
            + (T)6*L::cs2*rhoThetaBar*(c_u*c_u*((T)L::cNormSqr[iPop]-L::cs2*((T)L::d+(T)4))+L::cs2*uSqr*(L::cs2*((T)L::d+(T)2)-(T)L::cNormSqr[iPop]))
            + (T)3*L::cs2*L::cs2*rhoThetaBar*rhoThetaBar*invRho*(
                (T)L::cNormSqr[iPop]*(T)L::cNormSqr[iPop]-(T)2*L::cs2*(L::d+2)*L::cNormSqr[iPop]+L::cs2*L::cs2*(T)L::d*((T)L::d+(T)2)
            )
        )
    );
}

static T bgk_ma3_collision(Array<T,Descriptor::q>& f, T rhoBar, Array<T,Descriptor::d> const& j, T rhoThetaBar, T omega) {
    T invRho = Descriptor::invRho(rhoBar);
    const T jSqr = VectorTemplateImpl<T,Descriptor::d>::normSqr(j);
    for (plint iPop=0; iPop < Descriptor::q; ++iPop) {
        f[iPop] *= (T)1-omega;
        f[iPop] += omega * dynamicsTemplatesImpl<T,Descriptor>::bgk_ma3_equilibrium (
            iPop, rhoBar, invRho, j, jSqr, rhoThetaBar);
    }
    return jSqr*invRho*invRho;
}

static T bgk_ma4_collision(Array<T,Descriptor::q>& f, T rhoBar, Array<T,Descriptor::d> const& j, T rhoThetaBar, T omega) {
    T invRho = Descriptor::invRho(rhoBar);
    const T jSqr = VectorTemplateImpl<T,Descriptor::d>::normSqr(j);
    for (plint iPop=0; iPop < Descriptor::q; ++iPop) {
        f[iPop] *= (T)1-omega;
        f[iPop] += omega * dynamicsTemplatesImpl<T,Descriptor>::bgk_ma4_equilibrium (
            iPop, rhoBar, invRho, j, jSqr, rhoThetaBar);
    }
    return jSqr*invRho*invRho;
}

};  // struct dynamicsTemplatesImpl

}  // namespace plb

#include "latticeBoltzmann/dynamicsTemplates2D.h"
#include "latticeBoltzmann/dynamicsTemplates3D.h"

#endif  // DYNAMICS_TEMPLATES_H

