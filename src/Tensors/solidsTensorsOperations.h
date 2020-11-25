///Created By T.MAQUART

/*
Solids tensors operations.
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

#ifndef SOLIDS_TENSORS_OPERATIONS_H
#define SOLIDS_TENSORS_OPERATIONS_H

#include "palabos2D.h"
#include "palabos2D.hh"

#include <math.h>

// using namespaces in an other namespace is not recommended due to potential conflicts

namespace sld {

/// compute von mises stress
template<typename T, template<typename U> class Descriptor>
class ComputeVonMisesStress{
public:
    /// matching dimensions is your own problematic in order to avoid segmentation fault errors
    ComputeVonMisesStress(T*** Stress_, bool quiet_=true)
    {
        Stress=Stress_;
        quiet=quiet_;
    }
    /// multi dimensional array from constructor must match with size of the multi block lattice
    T*** ComputeVonMisesStressFromStressWithPhysicalConversion(plb::MultiBlockLattice2D<T,Descriptor>& lattice, T DeltaX, T DeltaT) const
    {
        // physical factor
        T phyFactor=DeltaX/DeltaT;

        // von mises initialization
        T*** vonmises=new T**[1];
        vonmises[0] = new T*[lattice.getNx()];
        for(int i = 0; i<lattice.getNx(); i++)
        {
            vonmises[0][i] = new T[lattice.getNy()];
        }

        // computing von mises criterion
        for (int i = 0; i<lattice.getNx(); i++)
        {
            for (int j = 0; j<lattice.getNy(); j++)
            {
               vonmises[0][i][j]=sqrt(pow(Stress[0][i][j]*phyFactor,2) - Stress[0][i][j]*phyFactor*Stress[1][i][j]*phyFactor + pow(Stress[1][i][j]*phyFactor,2) + 3.0*(pow(Stress[2][i][j]*phyFactor,2)));
            }
        }

        return vonmises;
    }
    /// multi dimensional array from constructor must match with size of the multi block lattice
    T*** ComputeVonMisesStressFromStress(plb::MultiBlockLattice2D<T,Descriptor>& lattice, T physicalFactorForSpeed) const
    {
        // von mises initialization
        T*** vonmises=new T**[1];
        vonmises[0] = new T*[lattice.getNx()];
        for(int i = 0; i<lattice.getNx(); i++)
        {
            vonmises[0][i] = new T[lattice.getNy()];
        }

        // computing von mises criterion
        for (int i = 0; i<lattice.getNx(); i++)
        {
            for (int j = 0; j<lattice.getNy(); j++)
            {
               vonmises[0][i][j]=sqrt(pow(Stress[0][i][j]*physicalFactorForSpeed,2) - Stress[0][i][j]*physicalFactorForSpeed*Stress[1][i][j]*physicalFactorForSpeed + pow(Stress[1][i][j]*physicalFactorForSpeed,2) + 3.0*(pow(Stress[2][i][j]*physicalFactorForSpeed,2)));
            }
        }

        return vonmises;
    }
    ~ComputeVonMisesStress()
    {
        // nothing to be done if object is new (a dynamic entity) and if you won't delete other dynamic types, just calling delete(object)
        // if it is a static object, it will be destroyed after run when return 0 is called
        if (!quiet)
        {
           plb::pcout << "Destroying: ~ComputeVonMisesStress" << std::endl;
        }
    }
private:
    T*** Stress;
    bool quiet;
};

/// compute speed with a backward difference scheme
template<typename T, template<typename U> class Descriptor, int nDim>
class ComputeSpeed{
public:
    /// matching dimensions is your own problematic in order to avoid segmentation fault errors
    ComputeSpeed(T*** backwarduStart_, T*** currentuStart_, bool quiet_=true)
    {
        backwarduStart=backwarduStart_;
        currentuStart=currentuStart_;
        quiet=quiet_;
    }
    /// multi dimensional array from constructor must match with size of the multi block lattice
    T*** ComputeSpeedFromDisplacement(plb::MultiBlockLattice2D<T,Descriptor>& lattice, T DeltaT) const
    {
        // speed initialization
        T*** speed=new T**[nDim];
        for(int i = 0; i<nDim; i++)
        {
            speed[i] = new T*[lattice.getNx()];
            for(int j = 0; j<lattice.getNx(); j++)
            {
                speed[i][j] = new T[lattice.getNy()];
            }
        }

        // computing speed
        for (int i = 0; i<lattice.getNx(); i++)
        {
            for (int j = 0; j<lattice.getNy(); j++)
            {
                for (int k = 0; k<nDim; k++)
                {
                    speed[k][i][j]=(currentuStart[k][i][j]-backwarduStart[k][i][j])/DeltaT;
                }
            }
        }

        return speed;
    }
    ~ComputeSpeed()
    {
        // nothing to be done if object is new (a dynamic entity) and if you won't delete other dynamic types, just calling delete(object)
        // if it is a static object, it will be destroyed after run when return 0 is called
        if (!quiet)
        {
           plb::pcout << "Destroying: ~ComputeSpeed" << std::endl;
        }
    }
private:
    T*** backwarduStart;
    T*** currentuStart;
    bool quiet;
};

/// compute acceleration with a backward difference scheme
template<typename T, template<typename U> class Descriptor, int nDim>
class ComputeAcceleration{
public:
    /// matching dimensions is your own problematic in order to avoid segmentation fault errors
    ComputeAcceleration(T*** backwardSpeed_, T*** currentSpeed_, bool quiet_=true)
    {
        backwardSpeed=backwardSpeed_;
        currentSpeed=currentSpeed_;
        quiet=quiet_;
    }
    /// multi dimensional array from constructor must match with size of the multi block lattice
    T*** ComputeAccelerationFromSpeed(plb::MultiBlockLattice2D<T,Descriptor>& lattice, T DeltaT) const
    {
        // acceleration initialization
        T*** acceleration=new T**[nDim];
        for(int i = 0; i<nDim; i++)
        {
            acceleration[i] = new T*[lattice.getNx()];
            for(int j = 0; j<lattice.getNx(); j++)
            {
                acceleration[i][j] = new T[lattice.getNy()];
            }
        }

        // computing acceleration
        for (int i = 0; i<lattice.getNx(); i++)
        {
            for (int j = 0; j<lattice.getNy(); j++)
            {
                for (int k = 0; k<nDim; k++)
                {
                    acceleration[k][i][j]=(currentSpeed[k][i][j]-backwardSpeed[k][i][j])/DeltaT;
                }
            }
        }

        return acceleration;
    }
    ~ComputeAcceleration()
    {
        // nothing to be done if object is new (a dynamic entity) and if you won't delete other dynamic types, just calling delete(object)
        // if it is a static object, it will be destroyed after run when return 0 is called
        if (!quiet)
        {
           plb::pcout << "Destroying: ~ComputeAcceleration" << std::endl;
        }
    }
private:
    T*** backwardSpeed;
    T*** currentSpeed;
    bool quiet;
};

/// compute delta of tensors
template<typename T, template<typename U> class Descriptor, int nDim>
class ComputeDeltaTensors{
public:
    /// matching dimensions is your own problematic in order to avoid segmentation fault errors
    ComputeDeltaTensors(T*** tensorOne_, T*** tensorTwo_, bool quiet_=true)
    {
        tensorOne=tensorOne_;
        tensorTwo=tensorTwo_;
        quiet=quiet_;
    }
    /// multi dimensional array from constructor must match with size of the multi block lattice
    T*** ComputeDeltaTensorsFromTwoTensors(plb::MultiBlockLattice2D<T,Descriptor>& lattice) const
    {
        // delta initialization
        T*** delta=new T**[nDim];
        for(int i = 0; i<nDim; i++)
        {
            delta[i] = new T*[lattice.getNx()];
            for(int j = 0; j<lattice.getNx(); j++)
            {
                delta[i][j] = new T[lattice.getNy()];
            }
        }

        // computing delta tensor
        for (int i = 0; i<lattice.getNx(); i++)
        {
            for (int j = 0; j<lattice.getNy(); j++)
            {
                for (int k = 0; k<nDim; k++)
                {
                    delta[k][i][j]=(tensorOne[k][i][j]-tensorTwo[k][i][j]);
                }
            }
        }

        return delta;
    }
    ~ComputeDeltaTensors()
    {
        // nothing to be done if object is new (a dynamic entity) and if you won't delete other dynamic types, just calling delete(object)
        // if it is a static object, it will be destroyed after run when return 0 is called
        if (!quiet)
        {
           plb::pcout << "Destroying: ~ComputeDeltaTensors" << std::endl;
        }
    }
private:
    T*** tensorOne;
    T*** tensorTwo;
    bool quiet;
};

}  // namespace sld

#endif  // SOLIDS_TENSORS_OPERATIONS_H

