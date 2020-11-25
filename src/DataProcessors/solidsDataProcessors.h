///Created By T.MAQUART

/*
Solids processors for lattice boltzmann iteration steps.
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

#ifndef SOLIDS_DATA_PROCESSORS_H
#define SOLIDS_DATA_PROCESSORS_H

#include "palabos2D.h"
#include "palabos2D.hh"

#include <solids2D.h> // solids headers

// using namespaces in an other namespace is not recommended due to potential conflicts

namespace sld {

/// data processor: data processing functional for speed time integration
template<typename T, template<typename U> class Descriptor>
class ComputeDisplacementFromSpeedProcessor : public plb::BoxProcessingFunctional2D_L<T,Descriptor>{
public:
    ComputeDisplacementFromSpeedProcessor(T ***uStart_, T DeltaT_, T DeltaX_, T Nx_, T Ny_, bool debug_)
    {
        uStart=uStart_;
        DeltaT=DeltaT_;
        DeltaX=DeltaX_;
        Nx=Nx_;
        Ny=Ny_;
        debug=debug_;
    }
    virtual void process(plb::Box2D domain, plb::BlockLattice2D<T,Descriptor>& lattice) {

        // access the position of the atomic-block inside the multi-block
        plb::Dot2D relativePosition = lattice.getLocation();

        for (plb::plint iX=domain.x0; iX<=domain.x1; ++iX)
        {
            for (plb::plint iY=domain.y0; iY<=domain.y1; ++iY)
            {
                // cell
                plb::Cell<T,Descriptor> cell = lattice.get(iX,iY); // don't need a reference

                // convert local coordinates to global ones
                int globalX = iX + relativePosition.x;
                int globalY = iY + relativePosition.y;

                // note: be careful can yield errors because of hardcoded numbers like globalX==10
//                // check test and debug
//                if (globalX==1 && globalY==1 && debug)
//                {
//                   plb::pcout << "Displacement in lattice units of cell before time integration 1,1: " << uStart[0][globalX][globalY] << " " <<uStart[1][globalX][globalY] << std::endl;
//                }
//                if (globalX==10 && globalY==10 && debug)
//                {
//                   plb::pcout << "Displacement in lattice units of cell before time integration 10,10: " << uStart[0][globalX][globalY] << " " <<uStart[1][globalX][globalY] << std::endl;
//                }
//                if (globalX==(Nx-2) && globalY==(Ny-2) && debug)
//                {
//                   plb::pcout << "Displacement in lattice units of cell before time integration " << (Nx-2) << "," << (Ny-2) << ": " << uStart[0][globalX][globalY] << " " <<uStart[1][globalX][globalY] << std::endl;
//                }
//                if (globalX==(Nx-1) && globalY==(Ny-1) && debug)
//                {
//                   plb::pcout << "Displacement in lattice units of cell before time integration " << (Nx-1) << "," << (Ny-1) << ": " << uStart[0][globalX][globalY] << " " <<uStart[1][globalX][globalY] << std::endl;
//                }

                // uStart of the corresponding cell
                plb::Array<T,2> uStartCell(uStart[0][globalX][globalY], uStart[1][globalX][globalY]); // need global position of the cell

                // dynamic object construction
                ComputeCellDisplacement<T, Descriptor>* Displacement = new ComputeCellDisplacement<T, Descriptor>(uStartCell);

                // compute displacement for the cell with time integration
                Displacement->ComputePhysicalDisplacementFromStart(cell, DeltaT, DeltaX);

                // displacement update
                uStart[0][globalX][globalY]=Displacement->GetuStartCellX(); // need global position of the cell
                uStart[1][globalX][globalY]=Displacement->GetuStartCellY(); // need global position of the cell

                // note: be careful can yield errors because of hardcoded numbers like globalX==10
//                // check test and debug
//                if (globalX==1 && globalY==1 && debug)
//                {
//                   plb::pcout << "Displacement in lattice units of cell after time integration 1,1: " << uStart[0][globalX][globalY]<< " " <<uStart[1][globalX][globalY]<< std::endl;
//                }
//                if (globalX==10 && globalY==10 && debug)
//                {
//                   plb::pcout << "Displacement in lattice units of cell after time integration 10,10: " << uStart[0][globalX][globalY]<< " " <<uStart[1][globalX][globalY]<< std::endl;
//                }
//                if (globalX==(Nx-2) && globalY==(Ny-2) && debug)
//                {
//                   plb::pcout << "Displacement in lattice units of cell after time integration " << (Nx-2) << "," << (Ny-2) << ": " << uStart[0][globalX][globalY] << " " <<uStart[1][globalX][globalY] << std::endl;
//                }
//                if (globalX==(Nx-1) && globalY==(Ny-1) && debug)
//                {
//                   plb::pcout << "Displacement in lattice units of cell after time integration " << (Nx-1) << "," << (Ny-1) << ": " << uStart[0][globalX][globalY] << " " <<uStart[1][globalX][globalY] << std::endl;
//                }

                // delete dynamic object
                delete(Displacement); // si l'objet n'est pas dynamique il sera détruit dans la boucle dans son domaine de validité et le destructeur sera appelé aussi
            }
        }
    }
    virtual void getTypeOfModification(std::vector<plb::modif::ModifT>& modified) const {
        modified[0] = plb::modif::nothing; // à comprendre à quoi cela sert voir 16.3.2 palabos user guide
        // indicates what kind of cell content was modified and must be updated in a multi-block structure
    }
    virtual ComputeDisplacementFromSpeedProcessor<T,Descriptor>* clone() const
    {
        return new ComputeDisplacementFromSpeedProcessor<T,Descriptor>(*this);
    }
    // le destructeur semble hérité mais quand est appelé le delete ?
private:
    T ***uStart;
    T DeltaT;
    T DeltaX;
    T Nx;
    T Ny;
    bool debug;
};

/// data processor: data processing functional to compute strain and stress
template<typename T, template<typename U> class Descriptor>
class ComputeStrainAndStressProcessor : public plb::BoxProcessingFunctional2D_L<T,Descriptor>{
public:
    ComputeStrainAndStressProcessor(T*** uStart_, T*** strainStart_, T*** stressStart_, T DeltaX_, sld::SolidConstitutiveLaw<T>* SolidLaw_, T Nx_, T Ny_, bool PeriodicNx_, bool PeriodicNy_, bool useCentralSchemeOnly_, bool use0X4CentralScheme_, bool debug_)
    {
        uStart=uStart_;
        strainStart=strainStart_;
        stressStart=stressStart_;
        DeltaX=DeltaX_;
        SolidLaw=SolidLaw_;
        Nx=Nx_;
        Ny=Ny_;
        PeriodicNx=PeriodicNx_;
        PeriodicNy=PeriodicNy_;
        useCentralSchemeOnly=useCentralSchemeOnly_;
        use0X4CentralScheme=use0X4CentralScheme_;
        debug=debug_;
    }
    virtual void process(plb::Box2D domain, plb::BlockLattice2D<T,Descriptor>& lattice) {

        // access the position of the atomic-block inside the multi-block
        plb::Dot2D relativePosition = lattice.getLocation();

        for (plb::plint iX=domain.x0; iX<=domain.x1; ++iX)
        {
            for (plb::plint iY=domain.y0; iY<=domain.y1; ++iY)
            {
                // convert local coordinates to global ones
                int globalX = iX + relativePosition.x;
                int globalY = iY + relativePosition.y;

                // determine cell vicinity: please see palabos user guide 16.3, we can access to the vicinity with respect to the communication envelope: but in this case it is already placed in global coordinates
                plb::Array<bool,4> cellFirstVicinity(false, false, false, false);
                if (globalX==0 || globalX==(Nx-1)){
                    if (globalX==0){
                        cellFirstVicinity[0]=true;
                    }
                    else
                    {
                        cellFirstVicinity[1]=true;
                    }
                }
                if (globalY==0 || globalY==(Ny-1)){
                    if (globalY==0){
                        cellFirstVicinity[2]=true;
                    }
                    else
                    {
                        cellFirstVicinity[3]=true;
                    }
                }

                // determine cell second vicinity
                plb::Array<bool,4> cellSecondVicinity(false, false, false, false);
                if (globalX==1 || globalX==(Nx-2)){
                    if (globalX==1){
                        cellSecondVicinity[0]=true;
                    }
                    else
                    {
                        cellSecondVicinity[1]=true;
                    }
                }
                if (globalY==1 || globalY==(Ny-2)){
                    if (globalY==1){
                        cellSecondVicinity[2]=true;
                    }
                    else
                    {
                        cellSecondVicinity[3]=true;
                    }
                }

                // strainStart of the corresponding cell
                plb::Array<T,3> strainStartCell(strainStart[0][globalX][globalY], strainStart[1][globalX][globalY], strainStart[2][globalX][globalY]); // need global position of the cell

                // stressStart of the corresponding cell
                plb::Array<T,3> stressStartCell(stressStart[0][globalX][globalY], stressStart[1][globalX][globalY], stressStart[2][globalX][globalY]); // need global position of the cell

                // dynamic object construction
                ComputeCellStrainStress<T,Descriptor>* StrainStress = new ComputeCellStrainStress<T,Descriptor>(strainStartCell, stressStartCell, cellFirstVicinity, cellSecondVicinity, Nx, Ny, PeriodicNx, PeriodicNy);

                // compute strain for the cell by space derivation and speed time integration
                StrainStress->ComputeStrainFromStart(uStart, globalX, globalY, DeltaX, useCentralSchemeOnly, use0X4CentralScheme);

                // compute stress for the cell
                StrainStress->ComputeStressFromStart(SolidLaw, globalX, globalY, DeltaX);

                // strain update
                strainStart[0][globalX][globalY]=StrainStress->GetstrainStartCell11(); // need global position of the cell
                strainStart[1][globalX][globalY]=StrainStress->GetstrainStartCell22(); // need global position of the cell
                strainStart[2][globalX][globalY]=StrainStress->GetstrainStartCell12(); // need global position of the cell

                // stress update
                stressStart[0][globalX][globalY]=StrainStress->GetstressStartCell11(); // need global position of the cell
                stressStart[1][globalX][globalY]=StrainStress->GetstressStartCell22(); // need global position of the cell
                stressStart[2][globalX][globalY]=StrainStress->GetstressStartCell12(); // need global position of the cell

                // delete dynamic object
                delete(StrainStress); // si l'objet n'est pas dynamique il sera détruit dans la boucle dans son domaine de validité et le destructeur sera appelé aussi
            }
        }
    }
    virtual void getTypeOfModification(std::vector<plb::modif::ModifT>& modified) const {
        modified[0] = plb::modif::nothing; // à comprendre à quoi cela sert voir 16.3.2 palabos user guide
        // indicates what kind of cell content was modified and must be updated in a multi-block structure
    }
    virtual ComputeStrainAndStressProcessor<T,Descriptor>* clone() const
    {
        return new ComputeStrainAndStressProcessor<T,Descriptor>(*this);
    }
    // le destructeur semble hérité mais quand est appelé le delete ?
private:
    T*** uStart;
    T*** strainStart;
    T*** stressStart;
    T DeltaX;
    SolidConstitutiveLaw<T>* SolidLaw;
    T Nx;
    T Ny;
    bool PeriodicNx;
    bool PeriodicNy;
    bool useCentralSchemeOnly;
    bool use0X4CentralScheme;
    bool debug;
};

/// data processor: data processing functional for force application which is space time dependent
template<typename T, template<typename U> class Descriptor>
class ComputeStressDivergenceAndApplySpaceTimeDependentExternalForceProcessor : public plb::BoxProcessingFunctional2D_L<T,Descriptor>{
public:
    ComputeStressDivergenceAndApplySpaceTimeDependentExternalForceProcessor(T ***stressStart_, T ***divStressStart_, T DeltaX_, T DeltaT_, T Nx_, T Ny_, bool PeriodicNx_, bool PeriodicNy_, bool useCentralSchemeOnly_, bool use0X4CentralScheme_, T convergenceFactor_, bool debug_)
    {
        stressStart=stressStart_;
        divStressStart=divStressStart_;
        DeltaX=DeltaX_;
        DeltaT=DeltaT_;
        Nx=Nx_;
        Ny=Ny_;
        PeriodicNx=PeriodicNx_;
        PeriodicNy=PeriodicNy_;
        useCentralSchemeOnly=useCentralSchemeOnly_;
        use0X4CentralScheme=use0X4CentralScheme_;
        convergenceFactor=convergenceFactor_;
        debug=debug_;
    }
    virtual void process(plb::Box2D domain, plb::BlockLattice2D<T,Descriptor>& lattice) {

        // initialization
        T forcex=0.0;
        T forcey=0.0;

        // access the position of the atomic-block inside the multi-block
        plb::Dot2D relativePosition = lattice.getLocation();

        for (plb::plint iX=domain.x0; iX<=domain.x1; ++iX)
        {
            for (plb::plint iY=domain.y0; iY<=domain.y1; ++iY)
            {
                // cell
                plb::Cell<T,Descriptor>& cell = lattice.get(iX,iY); // need a reference

                // convert local coordinates to global ones
                int globalX = iX + relativePosition.x;
                int globalY = iY + relativePosition.y;

                // determine cell vicinity: please see palabos user guide 16.3, we can access to the vicinity with respect to the communication envelope: but in this case it is already placed in global coordinates
                plb::Array<bool,4> cellFirstVicinity(false, false, false, false);
                if (globalX==0 || globalX==(Nx-1)){
                    if (globalX==0){
                        cellFirstVicinity[0]=true;
                    }
                    else
                    {
                        cellFirstVicinity[1]=true;
                    }
                }
                if (globalY==0 || globalY==(Ny-1)){
                    if (globalY==0){
                        cellFirstVicinity[2]=true;
                    }
                    else
                    {
                        cellFirstVicinity[3]=true;
                    }
                }

                // determine cell second vicinity
                plb::Array<bool,4> cellSecondVicinity(false, false, false, false);
                if (globalX==1 || globalX==(Nx-2)){
                    if (globalX==1){
                        cellSecondVicinity[0]=true;
                    }
                    else
                    {
                        cellSecondVicinity[1]=true;
                    }
                }
                if (globalY==1 || globalY==(Ny-2)){
                    if (globalY==1){
                        cellSecondVicinity[2]=true;
                    }
                    else
                    {
                        cellSecondVicinity[3]=true;
                    }
                }

                // divStressStart of the corresponding cell
                plb::Array<T,2> divStressStartCell(divStressStart[0][globalX][globalY], divStressStart[1][globalX][globalY]); // need global position of the cell

                // dynamic object construction
                ComputeCellStressDivergence<T,Descriptor>* StressDivergence = new ComputeCellStressDivergence<T,Descriptor>(divStressStartCell, cellFirstVicinity, cellSecondVicinity, Nx, Ny, PeriodicNx, PeriodicNy);

                // compute divergence of stress
                StressDivergence->ComputeDivergenceOfStressFromStart(stressStart, globalX, globalY, DeltaX, useCentralSchemeOnly, use0X4CentralScheme);

                // stress as force
                forcex=StressDivergence->GetdivStressStartCellX();
                forcey=StressDivergence->GetdivStressStartCellY();

                // stress divergence update
                divStressStart[0][globalX][globalY]=StressDivergence->GetdivStressStartCellX(); // need global position of the cell
                divStressStart[1][globalX][globalY]=StressDivergence->GetdivStressStartCellY(); // need global position of the cell

                // non-uniform gravity force: space and time dependent with internal stress
                T latticeUnitsForcex = convergenceFactor*forcex*((DeltaT*DeltaT)/DeltaX); // in lattice units
                T latticeUnitsForcey = convergenceFactor*forcey*((DeltaT*DeltaT)/DeltaX); // in lattice units
                plb::Array<T,2> forcefield(latticeUnitsForcex, latticeUnitsForcey);
                cell.setExternalField(Descriptor<T>::ExternalField::forceBeginsAt, Descriptor<T>::ExternalField::sizeOfForce, &forcefield[0]);

                // delete dynamic object
                delete(StressDivergence); // si l'objet n'est pas dynamique il sera détruit dans la boucle dans son domaine de validité et le destructeur sera appelé aussi
            }
        }
    }
    virtual void getTypeOfModification(std::vector<plb::modif::ModifT>& modified) const {
        modified[0] = plb::modif::allVariables; // à comprendre à quoi cela sert voir 16.3.2 palabos user guide
        // indicates what kind of cell content was modified and must be updated in a multi-block structure
    }
    virtual ComputeStressDivergenceAndApplySpaceTimeDependentExternalForceProcessor<T,Descriptor>* clone() const
    {
        return new ComputeStressDivergenceAndApplySpaceTimeDependentExternalForceProcessor<T,Descriptor>(*this);
    }
    // le destructeur semble hérité mais quand est appelé le delete ?
private:
    T ***stressStart;
    T ***divStressStart;
    T DeltaX;
    T DeltaT;
    T Nx;
    T Ny;
    bool PeriodicNx;
    bool PeriodicNy;
    bool useCentralSchemeOnly;
    bool use0X4CentralScheme;
    T convergenceFactor;
    bool debug;
};

/// data processor: data processing functional for clearing external forces
template<typename T, template<typename U> class Descriptor>
class ClearExternalForcesProcessor : public plb::BoxProcessingFunctional2D_L<T,Descriptor>{
public:
    ClearExternalForcesProcessor(int direction_, bool maxNxNy_, int includedMin_, int excludedMax_, int Nx_, int Ny_, bool clearOnlyCellsOfBoundaryConditions_, bool debug_)
    {
        direction=direction_;
        maxNxNy=maxNxNy_;
        includedMin=includedMin_;
        excludedMax=excludedMax_;
        Nx=Nx_;
        Ny=Ny_;
        clearOnlyCellsOfBoundaryConditions=clearOnlyCellsOfBoundaryConditions_;
        debug=debug_;
    }
    virtual void process(plb::Box2D domain, plb::BlockLattice2D<T,Descriptor>& lattice) {

        // access the position of the atomic-block inside the multi-block
        plb::Dot2D relativePosition = lattice.getLocation();

        for (plb::plint iX=domain.x0; iX<=domain.x1; ++iX)
        {
            for (plb::plint iY=domain.y0; iY<=domain.y1; ++iY)
            {
                // cell
                plb::Cell<T,Descriptor>& cell = lattice.get(iX,iY); // need a reference

                // zero force field
                plb::Array<T,2> forcefieldZero(0.0, 0.0);

                // convert local coordinates to global ones
                int globalX = iX + relativePosition.x;
                int globalY = iY + relativePosition.y;

                if (clearOnlyCellsOfBoundaryConditions==true)
                {
                    switch (direction)
                    {
                        case 0:
                            if (maxNxNy)
                            {
                                if (globalX==(Nx-1))
                                {
                                    if ((globalY>=includedMin) && (globalY<excludedMax))
                                    {
                                        cell.setExternalField(Descriptor<T>::ExternalField::forceBeginsAt, Descriptor<T>::ExternalField::sizeOfForce, &forcefieldZero[0]);
                                    }
                                }
                            }
                            else
                            {
                                if (globalX==0)
                                {
                                    if ((globalY>=includedMin) && (globalY<excludedMax))
                                    {
                                        cell.setExternalField(Descriptor<T>::ExternalField::forceBeginsAt, Descriptor<T>::ExternalField::sizeOfForce, &forcefieldZero[0]);
                                    }
                                }
                            }
                        case 1:
                            if (maxNxNy)
                            {
                                if (globalY==(Ny-1))
                                {
                                    if ((globalX>=includedMin) && (globalX<excludedMax))
                                    {
                                        cell.setExternalField(Descriptor<T>::ExternalField::forceBeginsAt, Descriptor<T>::ExternalField::sizeOfForce, &forcefieldZero[0]);
                                    }
                                }
                            }
                            else
                            {
                                if (globalY==0)
                                {
                                    if ((globalX>=includedMin) && (globalX<excludedMax))
                                    {
                                        cell.setExternalField(Descriptor<T>::ExternalField::forceBeginsAt, Descriptor<T>::ExternalField::sizeOfForce, &forcefieldZero[0]);
                                    }
                                }
                            }
                    }
                }
                else
                {
                    cell.setExternalField(Descriptor<T>::ExternalField::forceBeginsAt, Descriptor<T>::ExternalField::sizeOfForce, &forcefieldZero[0]);
                }
            }
        }
    }
    virtual void getTypeOfModification(std::vector<plb::modif::ModifT>& modified) const {
        modified[0] = plb::modif::allVariables; // à comprendre à quoi cela sert voir 16.3.2 palabos user guide
        // indicates what kind of cell content was modified and must be updated in a multi-block structure
    }
    virtual ClearExternalForcesProcessor<T,Descriptor>* clone() const
    {
        return new ClearExternalForcesProcessor<T,Descriptor>(*this);
    }
    // le destructeur semble hérité mais quand est appelé le delete ?
private:
    int direction;
    bool maxNxNy;
    int includedMin;
    int excludedMax;
    int Nx;
    int Ny;
    bool clearOnlyCellsOfBoundaryConditions;
    bool debug;
};

/// data processor: data processing functional for computing specific equilibrium near boundary conditions of solids
template<typename T, template<typename U> class Descriptor>
class EquilibriumOfBoundariesWithConcatenatedRelaxationParameter : public plb::BoxProcessingFunctional2D_L<T,Descriptor>{
public:
    EquilibriumOfBoundariesWithConcatenatedRelaxationParameter(T omegaBoundaryConditions_, T rho_, bool gaussHermite_, int direction_, bool maxNxNy_, int includedMin_, int excludedMax_, int Nx_, int Ny_, bool onlyEquilibrateCellsOfBoundaryConditions_, bool debug_)
    {
        omegaBoundaryConditions=omegaBoundaryConditions_;
        rho=rho_;
        gaussHermite=gaussHermite_;
        direction=direction_;
        maxNxNy=maxNxNy_;
        includedMin=includedMin_;
        excludedMax=excludedMax_;
        Nx=Nx_;
        Ny=Ny_;
        onlyEquilibrateCellsOfBoundaryConditions=onlyEquilibrateCellsOfBoundaryConditions_;
        debug=debug_;
    }
    virtual void process(plb::Box2D domain, plb::BlockLattice2D<T,Descriptor>& lattice) {

        // access the position of the atomic-block inside the multi-block
        plb::Dot2D relativePosition = lattice.getLocation();

        for (plb::plint iX=domain.x0; iX<=domain.x1; ++iX)
        {
            for (plb::plint iY=domain.y0; iY<=domain.y1; ++iY)
            {
                // cell
                plb::Cell<T,Descriptor>& cell = lattice.get(iX,iY); // need a reference

                // convert local coordinates to global ones
                int globalX = iX + relativePosition.x;
                int globalY = iY + relativePosition.y;

                if (onlyEquilibrateCellsOfBoundaryConditions==true)
                {
                    switch (direction)
                    {
                        case 0:
                            if (maxNxNy)
                            {
                                if (globalX==(Nx-1))
                                {
                                    if ((globalY>=includedMin) && (globalY<excludedMax))
                                    {
                                        if (gaussHermite==true)
                                        {
                                            T rhoBar=rho-1.0;
                                            plb::Array<T,9> populationAtEquilibrium;
                                            populationAtEquilibrium[0]=cell[0]-omegaBoundaryConditions*cell[0] +omegaBoundaryConditions*(4.0/9.0)*rhoBar;
                                            populationAtEquilibrium[1]=cell[1]-omegaBoundaryConditions*cell[1] +omegaBoundaryConditions*(1.0/36.0)*rhoBar;
                                            populationAtEquilibrium[2]=cell[2]-omegaBoundaryConditions*cell[2] +omegaBoundaryConditions*(1.0/9.0)*rhoBar;
                                            populationAtEquilibrium[3]=cell[3]-omegaBoundaryConditions*cell[3] +omegaBoundaryConditions*(1.0/36.0)*rhoBar;
                                            populationAtEquilibrium[4]=cell[4]-omegaBoundaryConditions*cell[4] +omegaBoundaryConditions*(1.0/9.0)*rhoBar;
                                            populationAtEquilibrium[5]=cell[5]-omegaBoundaryConditions*cell[5] +omegaBoundaryConditions*(1.0/36.0)*rhoBar;
                                            populationAtEquilibrium[6]=cell[6]-omegaBoundaryConditions*cell[6] +omegaBoundaryConditions*(1.0/9.0)*rhoBar;
                                            populationAtEquilibrium[7]=cell[7]-omegaBoundaryConditions*cell[7] +omegaBoundaryConditions*(1.0/36.0)*rhoBar;
                                            populationAtEquilibrium[8]=cell[8]-omegaBoundaryConditions*cell[8] +omegaBoundaryConditions*(1.0/9.0)*rhoBar;
                                            cell.setPopulations(populationAtEquilibrium);
                                        }
                                        else
                                        {
                                            T rhoBar=(rho-1.0)/9.0;
                                            plb::Array<T,9> populationAtEquilibrium;
                                            populationAtEquilibrium[0]=cell[0]-omegaBoundaryConditions*cell[0] +omegaBoundaryConditions*rhoBar;
                                            populationAtEquilibrium[1]=cell[1]-omegaBoundaryConditions*cell[1] +omegaBoundaryConditions*rhoBar;
                                            populationAtEquilibrium[2]=cell[2]-omegaBoundaryConditions*cell[2] +omegaBoundaryConditions*rhoBar;
                                            populationAtEquilibrium[3]=cell[3]-omegaBoundaryConditions*cell[3] +omegaBoundaryConditions*rhoBar;
                                            populationAtEquilibrium[4]=cell[4]-omegaBoundaryConditions*cell[4] +omegaBoundaryConditions*rhoBar;
                                            populationAtEquilibrium[5]=cell[5]-omegaBoundaryConditions*cell[5] +omegaBoundaryConditions*rhoBar;
                                            populationAtEquilibrium[6]=cell[6]-omegaBoundaryConditions*cell[6] +omegaBoundaryConditions*rhoBar;
                                            populationAtEquilibrium[7]=cell[7]-omegaBoundaryConditions*cell[7] +omegaBoundaryConditions*rhoBar;
                                            populationAtEquilibrium[8]=cell[8]-omegaBoundaryConditions*cell[8] +omegaBoundaryConditions*rhoBar;
                                            cell.setPopulations(populationAtEquilibrium);
                                        }
                                    }
                                }
                            }
                            else
                            {
                                if (globalX==0)
                                {
                                    if ((globalY>=includedMin) && (globalY<excludedMax))
                                    {
                                        if (gaussHermite==true)
                                        {
                                            T rhoBar=rho-1.0;
                                            plb::Array<T,9> populationAtEquilibrium;
                                            populationAtEquilibrium[0]=cell[0]-omegaBoundaryConditions*cell[0] +omegaBoundaryConditions*(4.0/9.0)*rhoBar;
                                            populationAtEquilibrium[1]=cell[1]-omegaBoundaryConditions*cell[1] +omegaBoundaryConditions*(1.0/36.0)*rhoBar;
                                            populationAtEquilibrium[2]=cell[2]-omegaBoundaryConditions*cell[2] +omegaBoundaryConditions*(1.0/9.0)*rhoBar;
                                            populationAtEquilibrium[3]=cell[3]-omegaBoundaryConditions*cell[3] +omegaBoundaryConditions*(1.0/36.0)*rhoBar;
                                            populationAtEquilibrium[4]=cell[4]-omegaBoundaryConditions*cell[4] +omegaBoundaryConditions*(1.0/9.0)*rhoBar;
                                            populationAtEquilibrium[5]=cell[5]-omegaBoundaryConditions*cell[5] +omegaBoundaryConditions*(1.0/36.0)*rhoBar;
                                            populationAtEquilibrium[6]=cell[6]-omegaBoundaryConditions*cell[6] +omegaBoundaryConditions*(1.0/9.0)*rhoBar;
                                            populationAtEquilibrium[7]=cell[7]-omegaBoundaryConditions*cell[7] +omegaBoundaryConditions*(1.0/36.0)*rhoBar;
                                            populationAtEquilibrium[8]=cell[8]-omegaBoundaryConditions*cell[8] +omegaBoundaryConditions*(1.0/9.0)*rhoBar;
                                            cell.setPopulations(populationAtEquilibrium);
                                        }
                                        else
                                        {
                                            T rhoBar=(rho-1.0)/9.0;
                                            plb::Array<T,9> populationAtEquilibrium;
                                            populationAtEquilibrium[0]=cell[0]-omegaBoundaryConditions*cell[0] +omegaBoundaryConditions*rhoBar;
                                            populationAtEquilibrium[1]=cell[1]-omegaBoundaryConditions*cell[1] +omegaBoundaryConditions*rhoBar;
                                            populationAtEquilibrium[2]=cell[2]-omegaBoundaryConditions*cell[2] +omegaBoundaryConditions*rhoBar;
                                            populationAtEquilibrium[3]=cell[3]-omegaBoundaryConditions*cell[3] +omegaBoundaryConditions*rhoBar;
                                            populationAtEquilibrium[4]=cell[4]-omegaBoundaryConditions*cell[4] +omegaBoundaryConditions*rhoBar;
                                            populationAtEquilibrium[5]=cell[5]-omegaBoundaryConditions*cell[5] +omegaBoundaryConditions*rhoBar;
                                            populationAtEquilibrium[6]=cell[6]-omegaBoundaryConditions*cell[6] +omegaBoundaryConditions*rhoBar;
                                            populationAtEquilibrium[7]=cell[7]-omegaBoundaryConditions*cell[7] +omegaBoundaryConditions*rhoBar;
                                            populationAtEquilibrium[8]=cell[8]-omegaBoundaryConditions*cell[8] +omegaBoundaryConditions*rhoBar;
                                            cell.setPopulations(populationAtEquilibrium);
                                        }
                                    }
                                }
                            }
                        case 1:
                            if (maxNxNy)
                            {
                                if (globalY==(Ny-1))
                                {
                                    if ((globalX>=includedMin) && (globalX<excludedMax))
                                    {
                                        if (gaussHermite==true)
                                        {
                                            T rhoBar=rho-1.0;
                                            plb::Array<T,9> populationAtEquilibrium;
                                            populationAtEquilibrium[0]=cell[0]-omegaBoundaryConditions*cell[0] +omegaBoundaryConditions*(4.0/9.0)*rhoBar;
                                            populationAtEquilibrium[1]=cell[1]-omegaBoundaryConditions*cell[1] +omegaBoundaryConditions*(1.0/36.0)*rhoBar;
                                            populationAtEquilibrium[2]=cell[2]-omegaBoundaryConditions*cell[2] +omegaBoundaryConditions*(1.0/9.0)*rhoBar;
                                            populationAtEquilibrium[3]=cell[3]-omegaBoundaryConditions*cell[3] +omegaBoundaryConditions*(1.0/36.0)*rhoBar;
                                            populationAtEquilibrium[4]=cell[4]-omegaBoundaryConditions*cell[4] +omegaBoundaryConditions*(1.0/9.0)*rhoBar;
                                            populationAtEquilibrium[5]=cell[5]-omegaBoundaryConditions*cell[5] +omegaBoundaryConditions*(1.0/36.0)*rhoBar;
                                            populationAtEquilibrium[6]=cell[6]-omegaBoundaryConditions*cell[6] +omegaBoundaryConditions*(1.0/9.0)*rhoBar;
                                            populationAtEquilibrium[7]=cell[7]-omegaBoundaryConditions*cell[7] +omegaBoundaryConditions*(1.0/36.0)*rhoBar;
                                            populationAtEquilibrium[8]=cell[8]-omegaBoundaryConditions*cell[8] +omegaBoundaryConditions*(1.0/9.0)*rhoBar;
                                            cell.setPopulations(populationAtEquilibrium);
                                        }
                                        else
                                        {
                                            T rhoBar=(rho-1.0)/9.0;
                                            plb::Array<T,9> populationAtEquilibrium;
                                            populationAtEquilibrium[0]=cell[0]-omegaBoundaryConditions*cell[0] +omegaBoundaryConditions*rhoBar;
                                            populationAtEquilibrium[1]=cell[1]-omegaBoundaryConditions*cell[1] +omegaBoundaryConditions*rhoBar;
                                            populationAtEquilibrium[2]=cell[2]-omegaBoundaryConditions*cell[2] +omegaBoundaryConditions*rhoBar;
                                            populationAtEquilibrium[3]=cell[3]-omegaBoundaryConditions*cell[3] +omegaBoundaryConditions*rhoBar;
                                            populationAtEquilibrium[4]=cell[4]-omegaBoundaryConditions*cell[4] +omegaBoundaryConditions*rhoBar;
                                            populationAtEquilibrium[5]=cell[5]-omegaBoundaryConditions*cell[5] +omegaBoundaryConditions*rhoBar;
                                            populationAtEquilibrium[6]=cell[6]-omegaBoundaryConditions*cell[6] +omegaBoundaryConditions*rhoBar;
                                            populationAtEquilibrium[7]=cell[7]-omegaBoundaryConditions*cell[7] +omegaBoundaryConditions*rhoBar;
                                            populationAtEquilibrium[8]=cell[8]-omegaBoundaryConditions*cell[8] +omegaBoundaryConditions*rhoBar;
                                            cell.setPopulations(populationAtEquilibrium);
                                        }
                                    }
                                }
                            }
                            else
                            {
                                if (globalY==0)
                                {
                                    if ((globalX>=includedMin) && (globalX<excludedMax))
                                    {
                                        if (gaussHermite==true)
                                        {
                                            T rhoBar=rho-1.0;
                                            plb::Array<T,9> populationAtEquilibrium;
                                            populationAtEquilibrium[0]=cell[0]-omegaBoundaryConditions*cell[0] +omegaBoundaryConditions*(4.0/9.0)*rhoBar;
                                            populationAtEquilibrium[1]=cell[1]-omegaBoundaryConditions*cell[1] +omegaBoundaryConditions*(1.0/36.0)*rhoBar;
                                            populationAtEquilibrium[2]=cell[2]-omegaBoundaryConditions*cell[2] +omegaBoundaryConditions*(1.0/9.0)*rhoBar;
                                            populationAtEquilibrium[3]=cell[3]-omegaBoundaryConditions*cell[3] +omegaBoundaryConditions*(1.0/36.0)*rhoBar;
                                            populationAtEquilibrium[4]=cell[4]-omegaBoundaryConditions*cell[4] +omegaBoundaryConditions*(1.0/9.0)*rhoBar;
                                            populationAtEquilibrium[5]=cell[5]-omegaBoundaryConditions*cell[5] +omegaBoundaryConditions*(1.0/36.0)*rhoBar;
                                            populationAtEquilibrium[6]=cell[6]-omegaBoundaryConditions*cell[6] +omegaBoundaryConditions*(1.0/9.0)*rhoBar;
                                            populationAtEquilibrium[7]=cell[7]-omegaBoundaryConditions*cell[7] +omegaBoundaryConditions*(1.0/36.0)*rhoBar;
                                            populationAtEquilibrium[8]=cell[8]-omegaBoundaryConditions*cell[8] +omegaBoundaryConditions*(1.0/9.0)*rhoBar;
                                            cell.setPopulations(populationAtEquilibrium);
                                        }
                                        else
                                        {
                                            T rhoBar=(rho-1.0)/9.0;
                                            plb::Array<T,9> populationAtEquilibrium;
                                            populationAtEquilibrium[0]=cell[0]-omegaBoundaryConditions*cell[0] +omegaBoundaryConditions*rhoBar;
                                            populationAtEquilibrium[1]=cell[1]-omegaBoundaryConditions*cell[1] +omegaBoundaryConditions*rhoBar;
                                            populationAtEquilibrium[2]=cell[2]-omegaBoundaryConditions*cell[2] +omegaBoundaryConditions*rhoBar;
                                            populationAtEquilibrium[3]=cell[3]-omegaBoundaryConditions*cell[3] +omegaBoundaryConditions*rhoBar;
                                            populationAtEquilibrium[4]=cell[4]-omegaBoundaryConditions*cell[4] +omegaBoundaryConditions*rhoBar;
                                            populationAtEquilibrium[5]=cell[5]-omegaBoundaryConditions*cell[5] +omegaBoundaryConditions*rhoBar;
                                            populationAtEquilibrium[6]=cell[6]-omegaBoundaryConditions*cell[6] +omegaBoundaryConditions*rhoBar;
                                            populationAtEquilibrium[7]=cell[7]-omegaBoundaryConditions*cell[7] +omegaBoundaryConditions*rhoBar;
                                            populationAtEquilibrium[8]=cell[8]-omegaBoundaryConditions*cell[8] +omegaBoundaryConditions*rhoBar;
                                            cell.setPopulations(populationAtEquilibrium);
                                        }
                                    }
                                }
                            }
                    }
                }
                else if (globalX==0 || globalX==(Nx-1) || globalY==0 || globalY==(Ny-1))
                {
                    if (gaussHermite==true)
                    {
                        T rhoBar=rho-1.0;
                        plb::Array<T,9> populationAtEquilibrium;
                        populationAtEquilibrium[0]=cell[0]-omegaBoundaryConditions*cell[0] +omegaBoundaryConditions*(4.0/9.0)*rhoBar;
                        populationAtEquilibrium[1]=cell[1]-omegaBoundaryConditions*cell[1] +omegaBoundaryConditions*(1.0/36.0)*rhoBar;
                        populationAtEquilibrium[2]=cell[2]-omegaBoundaryConditions*cell[2] +omegaBoundaryConditions*(1.0/9.0)*rhoBar;
                        populationAtEquilibrium[3]=cell[3]-omegaBoundaryConditions*cell[3] +omegaBoundaryConditions*(1.0/36.0)*rhoBar;
                        populationAtEquilibrium[4]=cell[4]-omegaBoundaryConditions*cell[4] +omegaBoundaryConditions*(1.0/9.0)*rhoBar;
                        populationAtEquilibrium[5]=cell[5]-omegaBoundaryConditions*cell[5] +omegaBoundaryConditions*(1.0/36.0)*rhoBar;
                        populationAtEquilibrium[6]=cell[6]-omegaBoundaryConditions*cell[6] +omegaBoundaryConditions*(1.0/9.0)*rhoBar;
                        populationAtEquilibrium[7]=cell[7]-omegaBoundaryConditions*cell[7] +omegaBoundaryConditions*(1.0/36.0)*rhoBar;
                        populationAtEquilibrium[8]=cell[8]-omegaBoundaryConditions*cell[8] +omegaBoundaryConditions*(1.0/9.0)*rhoBar;
                        cell.setPopulations(populationAtEquilibrium);
                    }
                    else
                    {
                        T rhoBar=(rho-1.0)/9.0;
                        plb::Array<T,9> populationAtEquilibrium;
                        populationAtEquilibrium[0]=cell[0]-omegaBoundaryConditions*cell[0] +omegaBoundaryConditions*rhoBar;
                        populationAtEquilibrium[1]=cell[1]-omegaBoundaryConditions*cell[1] +omegaBoundaryConditions*rhoBar;
                        populationAtEquilibrium[2]=cell[2]-omegaBoundaryConditions*cell[2] +omegaBoundaryConditions*rhoBar;
                        populationAtEquilibrium[3]=cell[3]-omegaBoundaryConditions*cell[3] +omegaBoundaryConditions*rhoBar;
                        populationAtEquilibrium[4]=cell[4]-omegaBoundaryConditions*cell[4] +omegaBoundaryConditions*rhoBar;
                        populationAtEquilibrium[5]=cell[5]-omegaBoundaryConditions*cell[5] +omegaBoundaryConditions*rhoBar;
                        populationAtEquilibrium[6]=cell[6]-omegaBoundaryConditions*cell[6] +omegaBoundaryConditions*rhoBar;
                        populationAtEquilibrium[7]=cell[7]-omegaBoundaryConditions*cell[7] +omegaBoundaryConditions*rhoBar;
                        populationAtEquilibrium[8]=cell[8]-omegaBoundaryConditions*cell[8] +omegaBoundaryConditions*rhoBar;
                        cell.setPopulations(populationAtEquilibrium);
                    }
                }
            }
        }
    }
    virtual void getTypeOfModification(std::vector<plb::modif::ModifT>& modified) const {
        modified[0] = plb::modif::allVariables; // à comprendre à quoi cela sert voir 16.3.2 palabos user guide
        // indicates what kind of cell content was modified and must be updated in a multi-block structure
    }
    virtual EquilibriumOfBoundariesWithConcatenatedRelaxationParameter<T,Descriptor>* clone() const
    {
        return new EquilibriumOfBoundariesWithConcatenatedRelaxationParameter<T,Descriptor>(*this);
    }
    // le destructeur semble hérité mais quand est appelé le delete ?
private:
    T omegaBoundaryConditions;
    T rho;
    bool gaussHermite;
    int direction;
    bool maxNxNy;
    int includedMin;
    int excludedMax;
    int Nx;
    int Ny;
    bool onlyEquilibrateCellsOfBoundaryConditions;
    bool debug;
};

/// data processor: data processing functional for dynamic force application which is space time dependent
/// in test and development
template<typename T, template<typename U> class Descriptor>
class ComputeDynamicEquilibriumForceAndApplySpaceTimeDependentExternalForceProcessor : public plb::BoxProcessingFunctional2D_L<T,Descriptor>{
public:
    ComputeDynamicEquilibriumForceAndApplySpaceTimeDependentExternalForceProcessor(int studyDynamicType_, T ***uStart_, T ***speedStart_, T ***accelerationStart_, T Nx_, T Ny_, T DeltaX_, T DeltaT_, sld::SolidConstitutiveLaw<T>* SolidLaw_, T alphaRayleigh_, T betaRayleigh_, bool debug_)
    {
        studyDynamicType=studyDynamicType_;
        uStart=uStart_;
        speedStart=speedStart_;
        accelerationStart=accelerationStart_;
        Nx=Nx_;
        Ny=Ny_;
        DeltaX=DeltaX_;
        DeltaT=DeltaT_;
        SolidLaw=SolidLaw_;
        alphaRayleigh=alphaRayleigh_;
        betaRayleigh=betaRayleigh_;
        debug=debug_;
    }
    virtual void process(plb::Box2D domain, plb::BlockLattice2D<T,Descriptor>& lattice) {

        // access the position of the atomic-block inside the multi-block
        plb::Dot2D relativePosition = lattice.getLocation();

        for (plb::plint iX=domain.x0; iX<=domain.x1; ++iX)
        {
            for (plb::plint iY=domain.y0; iY<=domain.y1; ++iY)
            {
                // initialization
                T forcex=0.0;
                T forcey=0.0;

                // cell
                plb::Cell<T,Descriptor>& cell = lattice.get(iX,iY); // need a reference

                // convert local coordinates to global ones
                int globalX = iX + relativePosition.x;
                int globalY = iY + relativePosition.y;

                // note: unités et dimensions à revoir

                // study dynamic type
                if (studyDynamicType==dynamicType::quasiStatic)
                {
                    // compute stiffness contribution like ES/L
                    // note: without poisson ratio for the moment
                    forcex+=((SolidLaw->GetYoungModulus()*DeltaX)/DeltaX)*uStart[0][globalX][globalY];
                    forcey+=((SolidLaw->GetYoungModulus()*DeltaX)/DeltaX)*uStart[1][globalX][globalY];

                }
                else if (studyDynamicType==dynamicType::semiDynamic)
                {
                    // compute stiffness contribution like ES/L
                    // note: without poisson ratio for the moment
                    forcex+=((SolidLaw->GetYoungModulus()*DeltaX)/DeltaX)*uStart[0][globalX][globalY];
                    forcey+=((SolidLaw->GetYoungModulus()*DeltaX)/DeltaX)*uStart[1][globalX][globalY];

                    // compute damping contribution with rayleigh coefficients
                    forcex+=((alphaRayleigh*(SolidLaw->GetRho()*DeltaX*DeltaX))+(betaRayleigh*(SolidLaw->GetRho()*DeltaX*DeltaX)))*speedStart[0][globalX][globalY];
                    forcey+=((alphaRayleigh*(SolidLaw->GetRho()*DeltaX*DeltaX))+(betaRayleigh*(SolidLaw->GetRho()*DeltaX*DeltaX)))*speedStart[1][globalX][globalY];
                }
                else if (studyDynamicType==dynamicType::fullDynamic)
                {
                    // compute stiffness contribution like ES/L
                    // note: without poisson ratio for the moment
                    forcex+=((SolidLaw->GetYoungModulus()*DeltaX)/DeltaX)*uStart[0][globalX][globalY];
                    forcey+=((SolidLaw->GetYoungModulus()*DeltaX)/DeltaX)*uStart[1][globalX][globalY];

                    // compute damping contribution with rayleigh coefficients
                    forcex+=((alphaRayleigh*(SolidLaw->GetRho()*DeltaX*DeltaX))+(betaRayleigh*(SolidLaw->GetRho()*DeltaX*DeltaX)))*speedStart[0][globalX][globalY];
                    forcey+=((alphaRayleigh*(SolidLaw->GetRho()*DeltaX*DeltaX))+(betaRayleigh*(SolidLaw->GetRho()*DeltaX*DeltaX)))*speedStart[1][globalX][globalY];

                    // compute inertial contribution
                    forcex+=(SolidLaw->GetRho()*DeltaX*DeltaX)*accelerationStart[0][globalX][globalY];
                    forcey+=(SolidLaw->GetRho()*DeltaX*DeltaX)*accelerationStart[1][globalX][globalY];
                }
                else
                {
                    throw new std::string("Error: wrong study dynamic type");
                }

                // non-uniform gravity force: space and time dependent
                T latticeUnitsForcex = forcex*((DeltaT*DeltaT)/DeltaX); // in lattice units
                T latticeUnitsForcey = forcey*((DeltaT*DeltaT)/DeltaX); // in lattice units
                plb::Array<T,2> forcefield(latticeUnitsForcex, latticeUnitsForcey);
                cell.setExternalField(Descriptor<T>::ExternalField::forceBeginsAt, Descriptor<T>::ExternalField::sizeOfForce, &forcefield[0]);
            }
        }
    }
    virtual void getTypeOfModification(std::vector<plb::modif::ModifT>& modified) const {
        modified[0] = plb::modif::allVariables; // à comprendre à quoi cela sert voir 16.3.2 palabos user guide
        // indicates what kind of cell content was modified and must be updated in a multi-block structure
    }
    virtual ComputeDynamicEquilibriumForceAndApplySpaceTimeDependentExternalForceProcessor<T,Descriptor>* clone() const
    {
        return new ComputeDynamicEquilibriumForceAndApplySpaceTimeDependentExternalForceProcessor<T,Descriptor>(*this);
    }
    // le destructeur semble hérité mais quand est appelé le delete ?
private:
    enum dynamicType{quasiStatic=1, semiDynamic=2, fullDynamic=3};
    int studyDynamicType;
    T ***uStart;
    T ***speedStart;
    T ***accelerationStart;
    T Nx;
    T Ny;
    T DeltaX;
    T DeltaT;
    sld::SolidConstitutiveLaw<T>* SolidLaw;
    T alphaRayleigh;
    T betaRayleigh;
    bool debug;
};

}  // namespace sld

#endif  // SOLIDS_DATA_PROCESSORS_H

