///Created By T.MAQUART

/*
Include headers files.
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

#include <DataProcessors/solidsDataProcessors.h> // solids data processors
#include <ConstitutiveLaws/solidsConstitutiveLaws.h> // solids constitutive laws
#include <Initialization/solidsInitialization.h> // solids distribution functions initialization
#include <BoundaryConditions/solidsBoundaryConditions.h> // solids boundary conditions
#include <Streaming/solidsStreamingModifications.h> // solids specific streaming modifications
#include <Cells/solidsCellsOperations.h> // solids cells operations
#include <Tensors/solidsTensorsOperations.h> // solids tensors operations
#include <Tensors/solidsTensorsInitializations.h> // solids tensors initializations
#include <Fields/solidsFields.h> // solids fields
#include <Stability/solidsStabilityChecker.h> // solids stability
