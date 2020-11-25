///Created By T.MAQUART

/*
Solids schemes stability.
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

#ifndef SOLIDS_STABILITY_CHECKER_H
#define SOLIDS_STABILITY_CHECKER_H

#include "palabos2D.h"
#include "palabos2D.hh"

// using namespaces in an other namespace is not recommended due to potential conflicts

namespace sld {

/// stability checker for cfl
template<typename T, template<typename U> class Descriptor>
class StabilityCheckerForCFL {
public:
    StabilityCheckerForCFL(bool exceptionHandle_, bool quiet_=true)
    {
        exceptionHandle=exceptionHandle_;
        quiet=quiet_;
    }
    void CheckCFLAndRaiseException(plb::MultiScalarField2D<T> initialVelocityNorm, plb::MultiBlockLattice2D<T,Descriptor>& lattice, T DeltaX, T DeltaT) const
    {
        bool stabilityCFL = true;
        T maxCFL = 0.3; // be careful, for incompressible flows there is a factor between mach number and cfl related to the speed of sound
        T maxVelocity = plb::computeMax(initialVelocityNorm, lattice.getBoundingBox());
        T maxVelocityPhysicalUnits = (maxVelocity*DeltaX)/DeltaT;

        // stability
        if (((maxVelocityPhysicalUnits*DeltaT)/(DeltaX)) > maxCFL)
        {
            stabilityCFL=false;
        }

        if (stabilityCFL==false)
        {
            std::cout << "Error: CFL number is to large: " << ((maxVelocityPhysicalUnits*DeltaT)/(DeltaX)) << std::endl;
            if (exceptionHandle==true)
            {
                throw new std::string("Error: CFL number is to large");
            }
        }

        // display stability
        std::cout << "CFL number is: " << ((maxVelocityPhysicalUnits*DeltaT)/(DeltaX)) << std::endl;
    }
    ~StabilityCheckerForCFL<T,Descriptor>()
    {
        // nothing to be done if object is new (a dynamic entity) and if you won't delete other dynamic types, just calling delete(object)
        // if it is a static object, it will be destroyed after run when return 0 is called
        if (!quiet)
        {
            plb::pcout << "Destroying: ~StabilityCheckerForCFL" << std::endl;
        }
    }
private:
    bool exceptionHandle;
    bool quiet;
};

/// stability checker with force
template<typename T, template<typename U> class Descriptor>
class StabilityCheckerByForce {
public:
    StabilityCheckerByForce(bool exceptionHandle_, bool quiet_=true)
    {
        exceptionHandle=exceptionHandle_;
        quiet=quiet_;
    }
    void AnalyseInitialStabilityAndRaiseException(T*** divStressStart, T rho, int Nx, int Ny, T DeltaX, T DeltaT, bool gaussHermite) const
    {
        // initial boolean
        bool stabilityForce = true;

        // equilibrium type and max factor
        T minMassDistributionWorstCase;
        T maxFactorD2Q9SchemeWorstCase;
        if (gaussHermite==true)
        {
            minMassDistributionWorstCase= ((rho-1.0)/36.0);
            maxFactorD2Q9SchemeWorstCase= (T)(1./6.); // due to (fx+fy)/12, like 2*fx/12
        }
        else
        {
            minMassDistributionWorstCase= ((rho-1.0)/9.0);
            maxFactorD2Q9SchemeWorstCase= (T)(1./3.);
        }

        // maximum stress divergence in lattice units
        T maxDivStressStart = findMaxValueOfThreeDimensionalArray(divStressStart, Nx, Ny);
        // T maxDivStressStartPhysicalUnits = maxDivStressStart*(DeltaX/DeltaT);
        T maxDivStressStartLatticeUnits = maxDivStressStart*((DeltaT*DeltaT)/DeltaX);

        // stability
        if (maxFactorD2Q9SchemeWorstCase*maxDivStressStartLatticeUnits > minMassDistributionWorstCase)
        {
            stabilityForce=false;
        }

        if (stabilityForce==false)
        {
            std::cout << "Error: parameters could yield negative distribution functions, decrease DeltaT" << std::endl;
            if (exceptionHandle==true)
            {
                throw new std::string("Error: parameters could yield negative distribution functions, decrease DeltaT");
            }
        }

        // display stability
        std::cout << "Ratio of force and mass distribution (%): " << (maxFactorD2Q9SchemeWorstCase*maxDivStressStartLatticeUnits*100.0)/minMassDistributionWorstCase << std::endl;
    }
    void AnalyseStabilityAndRaiseExceptionDuringIterations(plb::MultiBlockLattice2D<T,Descriptor>& lattice) const
    {
        bool stabilityForce=true;
        T eps=1.0*pow(10.0,-12.0);

        for(int t=0;t<lattice.getNy();t++)
        {
            for(int h=0;h<lattice.getNx();h++)
            {
                for (int q=1;q<Descriptor<T>::q; q++)
                {
                    if (lattice.get(h,t)[q]<(0.0-eps))
                    {
                        stabilityForce=false;
                        std::cout << "f(" << h << "," << t << ")" << "[" << q <<"] = " << lattice.get(h,t)[q] << std::endl;
                        break;
                    }
                }
                if (stabilityForce==false)
                {
                    break;
                }
            }
            if (stabilityForce==false)
            {
                break;
            }
        }

        if (stabilityForce==false)
        {
            std::cout << "Error: negative distribution function" << std::endl;
            if (exceptionHandle==true)
            {
                throw new std::string("Error: negative distribution function");
            }
        }
    }
    ~StabilityCheckerByForce<T,Descriptor>()
    {
        // nothing to be done if object is new (a dynamic entity) and if you won't delete other dynamic types, just calling delete(object)
        // if it is a static object, it will be destroyed after run when return 0 is called
        if (!quiet)
        {
            plb::pcout << "Destroying: ~StabilityCheckerByForce" << std::endl;
        }
    }
private:
    bool exceptionHandle;
    bool quiet;

    T findMaxValueOfThreeDimensionalArray(T*** divStressStart, int Nx, int Ny) const
    {
        T maxValue=0.0;

        for (int i=0; i<2; i++)
        {
            for (int j=0; j<Nx; j++)
            {
                for (int k=0; k<Ny; k++)
                {
                    if (std::abs(divStressStart[i][j][k]) > maxValue)
                    {
                        maxValue=std::abs(divStressStart[i][j][k]);
                    }
                }
            }
        }

        return maxValue;
    }
};

}  // namespace sld

#endif  // SOLIDS_STABILITY_CHECKER_H

