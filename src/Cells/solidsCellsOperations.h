///Created By T.MAQUART

/*
Solids cells operations.
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

#ifndef SOLIDS_CELLS_OPERATIONS_H
#define SOLIDS_CELLS_OPERATIONS_H

#include "palabos2D.h"
#include "palabos2D.hh"

#include <solids2D.h> // solids headers

// using namespaces in an other namespace is not recommended due to potential conflicts

namespace sld {

/// compute displacement on the cell from speed by time integration
template <typename T, template <typename U> class Descriptor>
class ComputeCellDisplacement {
public:
    ComputeCellDisplacement(plb::Array<T,2> uStartCell_, bool quiet_=true);
    void ComputePhysicalDisplacementFromStart(plb::Cell<T, Descriptor> cell, T DeltaT, T DeltaX);
    T GetuStartCellX() const;
    T GetuStartCellY() const;
    ~ComputeCellDisplacement();
private:
    plb::Array<T,2> uStartCell;
    bool quiet;
};

/// compute stress and strain on the cell
template <typename T, template <typename U> class Descriptor>
class ComputeCellStrainStress {
public:
    ComputeCellStrainStress(plb::Array<T,3> strainStartCell_, plb::Array<T,3> stressStartCell_, plb::Array<bool,4> cellFirstVicinity_, plb::Array<bool,4> cellSecondVicinity_, int GlobalNx_, int GlobalNy_, bool PeriodicNx_, bool PeriodicNy_, bool quiet_=true);
    void ComputeStrainFromStart(T*** uStart, int globalX, int globalY, T deltaX, bool useCentralSchemeOnly, bool use0X4CentralScheme, bool largeDeformations=false);
    void ComputeStressFromStart(sld::SolidConstitutiveLaw<T>* SolidLaw, int globalX, int globalY, T deltaX, bool nonLinearLaw=false);
    T GetstrainStartCell11() const;
    T GetstrainStartCell22() const;
    T GetstrainStartCell12() const;
    T GetstressStartCell11() const;
    T GetstressStartCell22() const;
    T GetstressStartCell12() const;
    ~ComputeCellStrainStress();
private:
    plb::Array<T,3> strainStartCell;
    plb::Array<T,3> stressStartCell;
    plb::Array<bool,4> cellFirstVicinity;
    plb::Array<bool,4> cellSecondVicinity;
    int GlobalNx;
    int GlobalNy;
    bool PeriodicNx;
    bool PeriodicNy;
    bool quiet;
};

/// compute stress divergence on the cell
/// if useCentralSchemeOnly=true it can yield errors on stress divergence because we need a 2 neighborhood exclusion due to double differentiations
template <typename T, template <typename U> class Descriptor>
class ComputeCellStressDivergence {
public:
    ComputeCellStressDivergence(plb::Array<T,2> divStressStartCell_, plb::Array<bool,4> cellFirstVicinity_, plb::Array<bool,4> cellSecondVicinity_, int GlobalNx_, int GlobalNy_, bool PeriodicNx_, bool PeriodicNy_, bool quiet_=true);
    void ComputeDivergenceOfStressFromStart(T*** stressStart, int globalX, int globalY, T deltaX, bool useCentralSchemeOnly, bool use0X4CentralScheme);
    T GetdivStressStartCellX() const;
    T GetdivStressStartCellY() const;
    ~ComputeCellStressDivergence();
private:
    plb::Array<T,2> divStressStartCell;
    plb::Array<bool,4> cellFirstVicinity;
    plb::Array<bool,4> cellSecondVicinity;
    int GlobalNx;
    int GlobalNy;
    bool PeriodicNx;
    bool PeriodicNy;
    bool quiet;
};

}  // namespace sld

#endif  // SOLIDS_CELLS_OPERATIONS_H

