///Created By T.MAQUART

/*
Constitutive laws for solids.
*/

/* This file is part of the SLD library.
 *
 * This library is mainly dedicated to research. This library uses PALABOS
 * library which has her own license. No warranty is given in terms of computational
 * efficiency and reliability of produced results. Moreover, no guarantees are given
 * concerning programmed interfaces and classes especially on the respect of C++ good practices.
 *
 * Only sequential mode is supported. SLD library has to be improved for parallelized runs
 * according to the PALABOS programming strategy.
 *
 * E-mail contact: tristan.maquart@emse.fr, tristan.maquart@hotmail.fr
 * E-mail contact: ronoel@ethz.ch
 * E-mail contact: navarro@emse.fr
 *
 * Copyright (C) <2020> <Tristan MAQUART, Romain NOEL, Laurent NAVARRO>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Affero General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

// T.MAQUART Ecole nationale supérieure des mines de Saint-Etienne : tristan.maquart@emse.fr, tristan.maquart@hotmail.fr
// R.NOEL Ecole polytechnique fédérale de Zurich : ronoel@ethz.ch
// L.NAVARRO Ecole nationale supérieure des mines de Saint-Etienne : navarro@emse.fr

#ifndef SOLIDS_CONSTITUTIVE_LAWS_H
#define SOLIDS_CONSTITUTIVE_LAWS_H

#include "palabos2D.h"
#include "palabos2D.hh"

#include <math.h>

// using namespaces in an other namespace is not recommended due to potential conflicts

namespace sld {

/// constitutive law for isotropic solids: hooke's law with lame coefficients
template <typename T>
class SolidConstitutiveLaw {
public:
    SolidConstitutiveLaw(T moduleYoung_, T poissonRatio_, T rho_, T DeltaX_, T DeltaT_, bool quiet_=true)
    {
        moduleYoung=moduleYoung_;
        poissonRatio=poissonRatio_;
        rho=rho_;
        DeltaX=DeltaX_;
        DeltaT=DeltaT_;
        quiet=quiet_;
        moduleYoungLatticeUnits=moduleYoung*DeltaX*DeltaT*DeltaT; // N/m2=kg/(m*s2)
    }
    plb::Array<T,2> GetIsotropicLameCoefficientsLatticeUnits() const
    {
        // lambda
        T lambda = (moduleYoungLatticeUnits*poissonRatio)/((1.0+poissonRatio)*(1.0-2.0*poissonRatio));

        // mu
        T mu = moduleYoungLatticeUnits/(2.0*(1.0+poissonRatio));

        // lame coefficients
        plb::Array<T,2> lameCoefficientsLatticeUnits(lambda,mu);

        return lameCoefficientsLatticeUnits;
    }
    plb::Array<T,2> GetIsotropicLameCoefficients() const
    {
        // lambda
        T lambda = (moduleYoung*poissonRatio)/((1.0+poissonRatio)*(1.0-2.0*poissonRatio));

        // mu
        T mu = moduleYoung/(2.0*(1.0+poissonRatio));

        // lame coefficients
        plb::Array<T,2> lameCoefficients(lambda,mu);

        return lameCoefficients;
    }
    T GetYoungModulus() const
    {
        return moduleYoung;
    }
    T GetYoungModulusLatticeUnits() const
    {
        return moduleYoungLatticeUnits;
    }
    T GetRho() const
    {
        return rho;
    }
    ~SolidConstitutiveLaw()
    {
        // nothing to be done if object is new (a dynamic entity) and if you won't delete other dynamic types, just calling delete(object)
        // if it is a static object, it will be destroyed after run when return 0 is called
        if (!quiet)
        {
           plb::pcout << "Destroying: ~SolidConstitutiveLaw" << std::endl;
        }
    }
private:
    T moduleYoung;
    T poissonRatio;
    T rho;
    T DeltaX;
    T DeltaT;
    T moduleYoungLatticeUnits;
    bool quiet;
};

/// constitutive law for isotropic solids: primary volume waves with lame coefficients
template <typename T>
class SolidPrimaryVolumeWaves {
public:
    SolidPrimaryVolumeWaves(T moduleYoung_, T poissonRatio_, T rho_, T DeltaT_, bool quiet_=true)
    {
        moduleYoung=moduleYoung_;
        poissonRatio=poissonRatio_;
        rho=rho_;
        DeltaT=DeltaT_;
        quiet=quiet_;

        // mu
        T mu=moduleYoung/(2.0*(1.0+poissonRatio));

        // bulk modulus
        bulkModulus=moduleYoung/(3.0*(1.0-2.0*poissonRatio));

        // wave speed
        primaryWaveVelocity=sqrt((bulkModulus+1.333333333*mu)/rho);

        // DeltaX
        DeltaX=primaryWaveVelocity*DeltaT;
    }
    T GetDeltaT() const
    {
        return DeltaT;
    }
    T GetDeltaX() const
    {
        return DeltaX;
    }
    T GetBulkModulus() const
    {
        return bulkModulus;
    }
    T GetPrimaryWaveVelocity() const
    {
        return primaryWaveVelocity;
    }
    ~SolidPrimaryVolumeWaves()
    {
        // nothing to be done if object is new (a dynamic entity) and if you won't delete other dynamic types, just calling delete(object)
        // if it is a static object, it will be destroyed after run when return 0 is called
        if (!quiet)
        {
           plb::pcout << "Destroying: ~SolidPrimaryVolumeWaves" << std::endl;
        }
    }
private:
    // constructor attributes
    T moduleYoung;
    T poissonRatio;
    T rho;
    T DeltaT;
    bool quiet;

    // computed attributes
    T DeltaX;
    T bulkModulus;
    T primaryWaveVelocity;
};

}  // namespace sld

#endif  // SOLIDS_CONSTITUTIVE_LAWS_H

