///Created By T.MAQUART

/*
Solids streaming specific modifications.
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

#ifndef SOLIDS_STREAMING_MODIFICATIONS_H
#define SOLIDS_STREAMING_MODIFICATIONS_H

#include "palabos2D.h"
#include "palabos2D.hh"

// using namespaces in an other namespace is not recommended due to potential conflicts

namespace sld {

/// modifications of a distribution function
template<typename T, template<typename U> class Descriptor>
class StreamingModificationsDistributionFunction {
public:
    StreamingModificationsDistributionFunction(bool quiet_=true)
    {
        quiet=quiet_;
    }
    void ModifyCellDistributionFunction(int indexOfDistributionFunction, int cellX, int cellY, plb::MultiBlockLattice2D<T,Descriptor>& lattice, T distributionValue) const
    {
        // get and modify distribution function
        lattice.get(cellX,cellY)[indexOfDistributionFunction]=distributionValue;
    }
    ~StreamingModificationsDistributionFunction<T,Descriptor>()
    {
        // nothing to be done if object is new (a dynamic entity) and if you won't delete other dynamic types, just calling delete(object)
        // if it is a static object, it will be destroyed after run when return 0 is called
        if (!quiet)
        {
            plb::pcout << "Destroying: ~StreamingModificationsDistributionFunction" << std::endl;
        }
    }
private:
    bool quiet;
};

}  // namespace sld

#endif  // SOLIDS_STREAMING_MODIFICATIONS_H

