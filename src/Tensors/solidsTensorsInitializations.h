///Created By T.MAQUART

/*
Solids tensors initializations.
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

#ifndef SOLIDS_TENSORS_INITIALIZATIONS_H
#define SOLIDS_TENSORS_INITIALIZATIONS_H

#include "palabos2D.h"
#include "palabos2D.hh"

// using namespaces in an other namespace is not recommended due to potential conflicts

namespace sld {

/// tensors initialization with n-dimensional pointers
template<typename T, template<typename U> class Descriptor>
class SolidsTensorsInitializer {
public:
    SolidsTensorsInitializer(bool quiet_=true)
    {
         quiet=quiet_;
    }
    void InitializePointers(plb::MultiBlockLattice2D<T,Descriptor>& lattice){
        // uStart
        uStart=new T**[2];
        for(int i=0; i<2; i++)
        {
            uStart[i] = new T*[lattice.getNx()];
            for(int j = 0; j<lattice.getNx(); j++)
            {
                uStart[i][j] = new T[lattice.getNy()];
            }
        }
        // uBackwardStart
        uBackwardStart=new T**[2];
        for(int i=0; i<2; i++)
        {
            uBackwardStart[i] = new T*[lattice.getNx()];
            for(int j = 0; j<lattice.getNx(); j++)
            {
                uBackwardStart[i][j] = new T[lattice.getNy()];
            }
        }
        // uSpeedBackwardStart
        uSpeedBackwardStart=new T**[2];
        for(int i=0; i<2; i++)
        {
            uSpeedBackwardStart[i] = new T*[lattice.getNx()];
            for(int j = 0; j<lattice.getNx(); j++)
            {
                uSpeedBackwardStart[i][j] = new T[lattice.getNy()];
            }
        }
        // strainStart
        strainStart=new T**[3];
        for(int i=0; i<3; i++)
        {
            strainStart[i] = new T*[lattice.getNx()];
            for(int j = 0; j<lattice.getNx(); j++)
            {
                strainStart[i][j] = new T[lattice.getNy()];
            }
        }
        // stressStart
        stressStart=new T**[3];
        for(int i=0; i<3; i++)
        {
            stressStart[i] = new T*[lattice.getNx()];
            for(int j = 0; j<lattice.getNx(); j++)
            {
                stressStart[i][j] = new T[lattice.getNy()];
            }
        }
        // divStressStart
        divStressStart=new T**[2];
        for(int i=0; i<2; i++)
        {
            divStressStart[i] = new T*[lattice.getNx()];
            for(int j = 0; j<lattice.getNx(); j++)
            {
                divStressStart[i][j] = new T[lattice.getNy()];
            }
        }
        // speedStart
        speedStart=new T**[2];
        for(int i=0; i<2; i++)
        {
            speedStart[i] = new T*[lattice.getNx()];
            for(int j = 0; j<lattice.getNx(); j++)
            {
                speedStart[i][j] = new T[lattice.getNy()];
            }
        }
        // accelerationStart
        accelerationStart=new T**[2];
        for(int i=0; i<2; i++)
        {
            accelerationStart[i] = new T*[lattice.getNx()];
            for(int j = 0; j<lattice.getNx(); j++)
            {
                accelerationStart[i][j] = new T[lattice.getNy()];
            }
        }
        // uDeltaStart
        uDeltaStart=new T**[2];
        for(int i=0; i<2; i++)
        {
            uDeltaStart[i] = new T*[lattice.getNx()];
            for(int j = 0; j<lattice.getNx(); j++)
            {
                uDeltaStart[i][j] = new T[lattice.getNy()];
            }
        }
    }
    /// they must have the same size in order to avoid segmentation fault errors
    void CopyDisplacementSpeedAccelerationTensors(T*** firstTensorToBeCopied, T*** secondTensorAsReference, plb::MultiBlockLattice2D<T,Descriptor>& lattice)
    {
        for(int i = 0; i<lattice.getNx(); i++)
        {
            for(int j = 0; j<lattice.getNy(); j++)
            {
                for (int k = 0; k<2; k++)
                {
                    firstTensorToBeCopied[k][i][j] =  secondTensorAsReference[k][i][j];
                }
            }
        }
    }
    T*** GetuStart() const
    {
         return uStart;
    }
    T*** GetuBackwardStart() const
    {
         return uBackwardStart;
    }
    T*** GetuSpeedBackwardStart() const
    {
         return uSpeedBackwardStart;
    }
    T*** GetstrainStart() const
    {
         return strainStart;
    }
    T*** GetstressStart() const
    {
         return stressStart;
    }
    T*** GetdivStressStart() const
    {
         return divStressStart;
    }
    T*** GetspeedStart() const
    {
         return speedStart;
    }
    T*** GetaccelerationStart() const
    {
         return accelerationStart;
    }
    T*** GetuDeltaStart() const
    {
         return uDeltaStart;
    }
    ~SolidsTensorsInitializer()
    {
        // nothing to be done if object is new (a dynamic entity) and if you won't delete other dynamic types, just calling delete(object)
        // if it is a static object, it will be destroyed after run when return 0 is called
        if (!quiet)
        {
           plb::pcout << "Destroying: ~SolidsTensorsInitializer" << std::endl;
        }
    }
private:
    T*** uStart;
    T*** uBackwardStart;
    T*** uSpeedBackwardStart;
    T*** strainStart;
    T*** stressStart;
    T*** divStressStart;
    T*** speedStart;
    T*** accelerationStart;
    T*** uDeltaStart;
    bool quiet;
};

}  // namespace sld

#endif  // SOLIDS_TENSORS_INITIALIZATIONS_H

