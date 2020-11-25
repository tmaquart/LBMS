///Created By T.MAQUART

/*
Solids cells operations implementation file.
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

#ifndef SOLIDS_CELLS_OPERATIONS_CPP
#define SOLIDS_CELLS_OPERATIONS_CPP

#include "palabos2D.h"
#include "palabos2D.hh"

#include <solids2D.h> // solids headers

// using namespaces in an other namespace is not recommended due to potential conflicts

using namespace sld;

/// compute displacement on the cell from speed by time integration
template <typename T, template <typename U> class Descriptor>
ComputeCellDisplacement<T,Descriptor>::ComputeCellDisplacement(plb::Array<T,2> uStartCell_, bool quiet_)
{
    uStartCell=uStartCell_;
    quiet=quiet_;
}

template <typename T, template <typename U> class Descriptor>
void ComputeCellDisplacement<T,Descriptor>::ComputePhysicalDisplacementFromStart(plb::Cell<T, Descriptor> cell, T DeltaT, T DeltaX)
{
    plb::Array<T,2> vCell(0.0, 0.0);
    cell.computeVelocity(vCell); // in lattice units

    uStartCell[0]+=vCell[0]*((DeltaX*DeltaT)/DeltaT); // in physical units
    uStartCell[1]+=vCell[1]*((DeltaX*DeltaT)/DeltaT); // in physical units
}

template <typename T, template <typename U> class Descriptor>
T ComputeCellDisplacement<T,Descriptor>::GetuStartCellX() const
{
    return uStartCell[0];
}

template <typename T, template <typename U> class Descriptor>
T ComputeCellDisplacement<T,Descriptor>::GetuStartCellY() const
{
    return uStartCell[1];
}

template <typename T, template <typename U> class Descriptor>
ComputeCellDisplacement<T,Descriptor>::~ComputeCellDisplacement()
{
    // nothing to be done if object is new (a dynamic entity) and if you won't delete other dynamic types, just calling delete(object)
    // if it is a static object, it will be destroyed after run when return 0 is called
    if (!quiet)
    {
        plb::pcout << "Destroying: ~ComputeCellDisplacement" << std::endl;
    }
}

/// compute stress and strain on the cell
template <typename T, template <typename U> class Descriptor>
ComputeCellStrainStress<T,Descriptor>::ComputeCellStrainStress(plb::Array<T,3> strainStartCell_, plb::Array<T,3> stressStartCell_, plb::Array<bool,4> cellFirstVicinity_, plb::Array<bool,4> cellSecondVicinity_, int GlobalNx_, int GlobalNy_, bool PeriodicNx_, bool PeriodicNy_, bool quiet_)
{
    strainStartCell=strainStartCell_;
    stressStartCell=stressStartCell_;
    cellFirstVicinity=cellFirstVicinity_;
    cellSecondVicinity=cellSecondVicinity_;
    GlobalNx=GlobalNx_;
    GlobalNy=GlobalNy_;
    PeriodicNx=PeriodicNx_;
    PeriodicNy=PeriodicNy_;
    quiet=quiet_;
}

template <typename T, template <typename U> class Descriptor>
void ComputeCellStrainStress<T,Descriptor>::ComputeStrainFromStart(T*** uStart, int globalX, int globalY, T DeltaX, bool useCentralSchemeOnly, bool use0X4CentralScheme, bool largeDeformations)
{
    if (largeDeformations==false)
    {
        // handle exceptions
        if (cellFirstVicinity[0]==true && cellFirstVicinity[1]==true)
        {
            throw new std::string("Error: mismatching cell vicinity");
        }
        if (cellFirstVicinity[2]==true && cellFirstVicinity[3]==true)
        {
            throw new std::string("Error: mismatching cell vicinity");
        }

        // if no vicinity problems
        if (cellFirstVicinity[0]==false && cellFirstVicinity[1]==false && cellFirstVicinity[2]==false && cellFirstVicinity[3]==false)
        {
            // central difference schemes
            if (cellSecondVicinity[0]==false && cellSecondVicinity[1]==false && cellSecondVicinity[2]==false && cellSecondVicinity[3]==false && use0X4CentralScheme==true)
            {
                // highly accurate central difference scheme + O(DeltaX^4)
                // integrating dx all deformations yield the maximum displacement
                strainStartCell[0]=(8.0*uStart[0][globalX+1][globalY]-8.0*uStart[0][globalX-1][globalY]-uStart[0][globalX+2][globalY]+uStart[0][globalX-2][globalY])/(12*DeltaX);
                strainStartCell[1]=(8.0*uStart[1][globalX][globalY+1]-8.0*uStart[1][globalX][globalY-1]-uStart[1][globalX][globalY+2]+uStart[1][globalX][globalY-2])/(12*DeltaX);
                T strainStartCell2First=(8.0*uStart[0][globalX][globalY+1]-8.0*uStart[0][globalX][globalY-1]-uStart[0][globalX][globalY+2]+uStart[0][globalX][globalY-2])/(12*DeltaX);
                T strainStartCell2Second=(8.0*uStart[1][globalX+1][globalY]-8.0*uStart[1][globalX-1][globalY]-uStart[1][globalX+2][globalY]+uStart[1][globalX-2][globalY])/(12*DeltaX);
                strainStartCell[2]=0.5*(strainStartCell2First+strainStartCell2Second);
            }
            else
            {
                // central difference scheme + O(DeltaX^2)
                // integrating dx all deformations yield the maximum displacement
                strainStartCell[0]=(uStart[0][globalX+1][globalY]-uStart[0][globalX-1][globalY])/(2*DeltaX);
                strainStartCell[1]=(uStart[1][globalX][globalY+1]-uStart[1][globalX][globalY-1])/(2*DeltaX);
                strainStartCell[2]=0.5*((uStart[0][globalX][globalY+1]-uStart[0][globalX][globalY-1])/(2*DeltaX) + (uStart[1][globalX+1][globalY]-uStart[1][globalX-1][globalY])/(2*DeltaX));
            }
            return;
        }
        else if (useCentralSchemeOnly==false) // depending on boundary conditions
        {
            // Nx=0
            if (cellFirstVicinity[0]==true)
            {
                if (cellFirstVicinity[2]==false && cellFirstVicinity[3]==false)
                {
                    if (PeriodicNx || PeriodicNy)
                    {
                        if (PeriodicNx && PeriodicNy)
                        {
                            int Xup = globalX+1;
                            int Xdwn = GlobalNx-1;
                            int Yup = globalY+1;
                            int Ydwn = globalY-1;
                            strainStartCell[0]=(uStart[0][Xup][globalY]-uStart[0][Xdwn][globalY])/(2*DeltaX);
                            strainStartCell[1]=(uStart[1][globalX][Yup]-uStart[1][globalX][Ydwn])/(2*DeltaX);
                            strainStartCell[2]=0.5*((uStart[0][globalX][Yup]-uStart[0][globalX][Ydwn])/(2*DeltaX) + (uStart[1][Xup][globalY]-uStart[1][Xdwn][globalY])/(2*DeltaX));
                            return;
                        }
                        else if (PeriodicNx)
                        {
                            int Xup = globalX+1;
                            int Xdwn = GlobalNx-1;
                            int Yup = globalY+1;
                            int Ydwn = globalY-1;
                            strainStartCell[0]=(uStart[0][Xup][globalY]-uStart[0][Xdwn][globalY])/(2*DeltaX);
                            strainStartCell[1]=(uStart[1][globalX][Yup]-uStart[1][globalX][Ydwn])/(2*DeltaX);
                            strainStartCell[2]=0.5*((uStart[0][globalX][Yup]-uStart[0][globalX][Ydwn])/(2*DeltaX) + (uStart[1][Xup][globalY]-uStart[1][Xdwn][globalY])/(2*DeltaX));
                            return;
                        }
                        else
                        {
                            int Yup = globalY+1;
                            int Ydwn = globalY-1;
                            strainStartCell[0]=(uStart[0][globalX+1][globalY]-uStart[0][globalX][globalY])/(DeltaX); // forward Nx scheme
                            strainStartCell[1]=(uStart[1][globalX][Yup]-uStart[1][globalX][Ydwn])/(2*DeltaX);
                            strainStartCell[2]=0.5*((uStart[0][globalX][Yup]-uStart[0][globalX][Ydwn])/(2*DeltaX) + (uStart[1][globalX+1][globalY]-uStart[1][globalX][globalY])/(DeltaX)); // forward Nx scheme
                            return;
                        }
                    }
                    else
                    {
                        // strainStartCell[0]=(uStart[0][globalX+1][globalY]-uStart[0][globalX][globalY])/(DeltaX); // forward Nx scheme
                        strainStartCell[0]=(4.0*uStart[0][globalX+1][globalY]-3.0*uStart[0][globalX][globalY]-uStart[0][globalX+2][globalY])/(2*DeltaX); // three-point forward Nx scheme O(DeltaX^2)
                        strainStartCell[1]=(uStart[1][globalX][globalY+1]-uStart[1][globalX][globalY-1])/(2*DeltaX); // central scheme
                        // strainStartCell[2]=0.5*((uStart[0][globalX][globalY+1]-uStart[0][globalX][globalY-1])/(2*DeltaX) + (uStart[1][globalX+1][globalY]-uStart[1][globalX][globalY])/(DeltaX)); // forward Nx scheme
                        strainStartCell[2]=0.5*((uStart[0][globalX][globalY+1]-uStart[0][globalX][globalY-1])/(2*DeltaX) + (4.0*uStart[1][globalX+1][globalY]-3.0*uStart[1][globalX][globalY]-uStart[1][globalX+2][globalY])/(2*DeltaX)); // three-point forward Nx scheme O(DeltaX^2)
                        return;
                    }
                }
                else if (cellFirstVicinity[2]==true && cellFirstVicinity[3]==false)
                {
                    if (PeriodicNx || PeriodicNy)
                    {
                        if (PeriodicNx && PeriodicNy)
                        {
                            int Xup = globalX+1;
                            int Xdwn = GlobalNx-1;
                            int Yup = globalY+1;
                            int Ydwn = GlobalNy-1;
                            strainStartCell[0]=(uStart[0][Xup][globalY]-uStart[0][Xdwn][globalY])/(2*DeltaX);
                            strainStartCell[1]=(uStart[1][globalX][Yup]-uStart[1][globalX][Ydwn])/(2*DeltaX);
                            strainStartCell[2]=0.5*((uStart[0][globalX][Yup]-uStart[0][globalX][Ydwn])/(2*DeltaX) + (uStart[1][Xup][globalY]-uStart[1][Xdwn][globalY])/(2*DeltaX));
                            return;
                        }
                        else if (PeriodicNx)
                        {
                            int Xup = globalX+1;
                            int Xdwn = GlobalNx-1;
                            strainStartCell[0]=(uStart[0][Xup][globalY]-uStart[0][Xdwn][globalY])/(2*DeltaX);
                            strainStartCell[1]=(uStart[1][globalX][globalY+1]-uStart[1][globalX][globalY])/(DeltaX); // forward Ny scheme
                            strainStartCell[2]=0.5*((uStart[0][globalX][globalY+1]-uStart[0][globalX][globalY])/(DeltaX) + (uStart[1][Xup][globalY]-uStart[1][Xdwn][globalY])/(2*DeltaX)); // forward Ny scheme
                            return;
                        }
                        else
                        {
                            int Yup = globalY+1;
                            int Ydwn = GlobalNy-1;
                            strainStartCell[0]=(uStart[0][globalX+1][globalY]-uStart[0][globalX][globalY])/(DeltaX); // forward Nx scheme
                            strainStartCell[1]=(uStart[1][globalX][Yup]-uStart[1][globalX][Ydwn])/(2*DeltaX);
                            strainStartCell[2]=0.5*((uStart[0][globalX][Yup]-uStart[0][globalX][Ydwn])/(2*DeltaX) + (uStart[1][globalX+1][globalY]-uStart[1][globalX][globalY])/(DeltaX)); // forward Nx scheme
                            return;
                        }
                    }
                    else
                    {
                        // strainStartCell[0]=(uStart[0][globalX+1][globalY]-uStart[0][globalX][globalY])/(DeltaX); // forward Nx scheme
                        strainStartCell[0]=(4.0*uStart[0][globalX+1][globalY]-3.0*uStart[0][globalX][globalY]-uStart[0][globalX+2][globalY])/(2*DeltaX); // three-point forward Nx scheme O(DeltaX^2)
                        // strainStartCell[1]=(uStart[1][globalX][globalY+1]-uStart[1][globalX][globalY])/(DeltaX); // forward Ny scheme
                        strainStartCell[1]=(4.0*uStart[1][globalX][globalY+1]-3.0*uStart[1][globalX][globalY]-uStart[1][globalX][globalY+2])/(2*DeltaX); // three-point forward Ny scheme O(DeltaX^2)
                        // strainStartCell[2]=0.5*((uStart[0][globalX][globalY+1]-uStart[0][globalX][globalY])/(DeltaX) + (uStart[1][globalX+1][globalY]-uStart[1][globalX][globalY])/(DeltaX)); // forward Nx Ny scheme
                        strainStartCell[2]=0.5*((4.0*uStart[0][globalX][globalY+1]-3.0*uStart[0][globalX][globalY]-uStart[0][globalX][globalY+2])/(2*DeltaX) + (4.0*uStart[1][globalX+1][globalY]-3.0*uStart[1][globalX][globalY]-uStart[1][globalX+2][globalY])/(2*DeltaX)); // three-point forward Nx Ny scheme O(DeltaX^2)
                        return;
                    }
                }
                else if (cellFirstVicinity[2]==false && cellFirstVicinity[3]==true)
                {
                    if (PeriodicNx || PeriodicNy)
                    {
                        if (PeriodicNx && PeriodicNy)
                        {
                            int Xup = globalX+1;
                            int Xdwn = GlobalNx-1;
                            int Yup = 0;
                            int Ydwn = globalY-1;
                            strainStartCell[0]=(uStart[0][Xup][globalY]-uStart[0][Xdwn][globalY])/(2*DeltaX);
                            strainStartCell[1]=(uStart[1][globalX][Yup]-uStart[1][globalX][Ydwn])/(2*DeltaX);
                            strainStartCell[2]=0.5*((uStart[0][globalX][Yup]-uStart[0][globalX][Ydwn])/(2*DeltaX) + (uStart[1][Xup][globalY]-uStart[1][Xdwn][globalY])/(2*DeltaX));
                            return;
                        }
                        else if (PeriodicNx)
                        {
                            int Xup = globalX+1;
                            int Xdwn = GlobalNx-1;
                            strainStartCell[0]=(uStart[0][Xup][globalY]-uStart[0][Xdwn][globalY])/(2*DeltaX);
                            strainStartCell[1]=(uStart[1][globalX][globalY]-uStart[1][globalX][globalY-1])/(DeltaX); // backward Ny scheme
                            strainStartCell[2]=0.5*((uStart[0][globalX][globalY]-uStart[0][globalX][globalY-1])/(DeltaX) + (uStart[1][Xup][globalY]-uStart[1][Xdwn][globalY])/(2*DeltaX)); // backward Ny scheme
                            return;
                        }
                        else
                        {
                            int Yup = 0;
                            int Ydwn = globalY-1;
                            strainStartCell[0]=(uStart[0][globalX+1][globalY]-uStart[0][globalX][globalY])/(DeltaX); // forward Nx scheme
                            strainStartCell[1]=(uStart[1][globalX][Yup]-uStart[1][globalX][Ydwn])/(2*DeltaX);
                            strainStartCell[2]=0.5*((uStart[0][globalX][Yup]-uStart[0][globalX][Ydwn])/(2*DeltaX) + (uStart[1][globalX+1][globalY]-uStart[1][globalX][globalY])/(DeltaX));  // forward Nx scheme
                            return;
                        }
                    }
                    else
                    {
                        // strainStartCell[0]=(uStart[0][globalX+1][globalY]-uStart[0][globalX][globalY])/(DeltaX); // forward Nx scheme
                        strainStartCell[0]=(4.0*uStart[0][globalX+1][globalY]-3.0*uStart[0][globalX][globalY]-uStart[0][globalX+2][globalY])/(2*DeltaX); // three-point forward Nx scheme O(DeltaX^2)
                        // strainStartCell[1]=(uStart[1][globalX][globalY]-uStart[1][globalX][globalY-1])/(DeltaX); // backward Ny scheme
                        strainStartCell[1]=(3.0*uStart[1][globalX][globalY]-4.0*uStart[1][globalX][globalY-1]+uStart[1][globalX][globalY-2])/(2*DeltaX); // three-point backward Ny scheme O(DeltaX^2)
                        // strainStartCell[2]=0.5*((uStart[0][globalX][globalY]-uStart[0][globalX][globalY-1])/(DeltaX) + (uStart[1][globalX+1][globalY]-uStart[1][globalX][globalY])/(DeltaX)); // forward Nx backward Ny scheme
                        strainStartCell[2]=0.5*((3.0*uStart[0][globalX][globalY]-4.0*uStart[0][globalX][globalY-1]+uStart[0][globalX][globalY-2])/(2*DeltaX) + (4.0*uStart[1][globalX+1][globalY]-3.0*uStart[1][globalX][globalY]-uStart[1][globalX+2][globalY])/(2*DeltaX)); // three-point forward Nx backward Ny scheme O(DeltaX^2)
                        return;
                    }
                }
            }
            // Nx-1
            if (cellFirstVicinity[1]==true)
            {
                if (cellFirstVicinity[2]==false && cellFirstVicinity[3]==false)
                {
                    if (PeriodicNx || PeriodicNy)
                    {
                        if (PeriodicNx && PeriodicNy)
                        {
                            // not implemented
                            throw new std::string("Error: not implemented");
                        }
                        else if (PeriodicNx)
                        {
                            int Xup = 0;
                            int Xdwn = globalX-1;
                            int Yup = globalY+1;
                            int Ydwn = globalY-1;
                            strainStartCell[0]=(uStart[0][Xup][globalY]-uStart[0][Xdwn][globalY])/(2*DeltaX);
                            strainStartCell[1]=(uStart[1][globalX][Yup]-uStart[1][globalX][Ydwn])/(2*DeltaX);
                            strainStartCell[2]=0.5*((uStart[0][globalX][Yup]-uStart[0][globalX][Ydwn])/(2*DeltaX) + (uStart[1][Xup][globalY]-uStart[1][Xdwn][globalY])/(2*DeltaX));
                            return;
                        }
                        else
                        {
                            // not implemented
                            throw new std::string("Error: not implemented");
                        }
                    }
                    else
                    {
                        // strainStartCell[0]=(uStart[0][globalX][globalY]-uStart[0][globalX-1][globalY])/(DeltaX); // backward Nx scheme
                        strainStartCell[0]=(3.0*uStart[0][globalX][globalY]-4.0*uStart[0][globalX-1][globalY]+uStart[0][globalX-2][globalY])/(2*DeltaX); // three-point backward Nx scheme O(DeltaX^2)
                        strainStartCell[1]=(uStart[1][globalX][globalY+1]-uStart[1][globalX][globalY-1])/(2*DeltaX); // central scheme
                        // strainStartCell[2]=0.5*((uStart[0][globalX][globalY+1]-uStart[0][globalX][globalY-1])/(2*DeltaX) + (uStart[1][globalX][globalY]-uStart[1][globalX-1][globalY])/(DeltaX)); // backward Nx scheme
                        strainStartCell[2]=0.5*((uStart[0][globalX][globalY+1]-uStart[0][globalX][globalY-1])/(2*DeltaX) + (3.0*uStart[1][globalX][globalY]-4.0*uStart[1][globalX-1][globalY]+uStart[1][globalX-2][globalY])/(2*DeltaX)); // three-point backward Nx scheme O(DeltaX^2)
                        return;
                    }
                }
                else if (cellFirstVicinity[2]==true && cellFirstVicinity[3]==false)
                {
                    if (PeriodicNx || PeriodicNy)
                    {
                        if (PeriodicNx && PeriodicNy)
                        {
                            // not implemented
                            throw new std::string("Error: not implemented");
                        }
                        else if (PeriodicNx)
                        {
                            int Xup = 0;
                            int Xdwn = globalX-1;
                            strainStartCell[0]=(uStart[0][Xup][globalY]-uStart[0][Xdwn][globalY])/(2*DeltaX);
                            strainStartCell[1]=(uStart[1][globalX][globalY+1]-uStart[1][globalX][globalY])/(DeltaX); // forward Ny scheme
                            strainStartCell[2]=0.5*((uStart[0][globalX][globalY+1]-uStart[0][globalX][globalY])/(DeltaX) + (uStart[1][Xup][globalY]-uStart[1][Xdwn][globalY])/(2*DeltaX)); // forward Ny scheme
                            return;
                        }
                        else
                        {
                            // not implemented
                            throw new std::string("Error: not implemented");
                        }
                    }
                    else
                    {
                        // strainStartCell[0]=(uStart[0][globalX][globalY]-uStart[0][globalX-1][globalY])/(DeltaX); // backward Nx scheme
                        strainStartCell[0]=(3.0*uStart[0][globalX][globalY]-4.0*uStart[0][globalX-1][globalY]+uStart[0][globalX-2][globalY])/(2*DeltaX); // three-point backward Nx scheme O(DeltaX^2)
                        // strainStartCell[1]=(uStart[1][globalX][globalY+1]-uStart[1][globalX][globalY])/(DeltaX); // forward Ny scheme
                        strainStartCell[1]=(4.0*uStart[1][globalX][globalY+1]-3.0*uStart[1][globalX][globalY]-uStart[1][globalX][globalY+2])/(2*DeltaX); // three-point forward Ny scheme O(DeltaX^2)
                        // strainStartCell[2]=0.5*((uStart[0][globalX][globalY+1]-uStart[0][globalX][globalY])/(DeltaX) + (uStart[1][globalX][globalY]-uStart[1][globalX-1][globalY])/(DeltaX)); // backward Nx forward Ny scheme
                        strainStartCell[2]=0.5*((4.0*uStart[0][globalX][globalY+1]-3.0*uStart[0][globalX][globalY]-uStart[0][globalX][globalY+2])/(2*DeltaX) + (3.0*uStart[1][globalX][globalY]-4.0*uStart[1][globalX-1][globalY]+uStart[1][globalX-2][globalY])/(2*DeltaX)); // three-point backward Nx forward Ny scheme O(DeltaX^2)
                        return;
                    }
                }
                else if (cellFirstVicinity[2]==false && cellFirstVicinity[3]==true)
                {
                    if (PeriodicNx || PeriodicNy)
                    {
                        if (PeriodicNx && PeriodicNy)
                        {
                            // not implemented
                            throw new std::string("Error: not implemented");
                        }
                        else if (PeriodicNx)
                        {
                            int Xup = 0;
                            int Xdwn = globalX-1;
                            strainStartCell[0]=(uStart[0][Xup][globalY]-uStart[0][Xdwn][globalY])/(2*DeltaX);
                            strainStartCell[1]=(uStart[1][globalX][globalY]-uStart[1][globalX][globalY-1])/(DeltaX); // backward Ny scheme
                            strainStartCell[2]=0.5*((uStart[0][globalX][globalY]-uStart[0][globalX][globalY-1])/(DeltaX) + (uStart[1][Xup][globalY]-uStart[1][Xdwn][globalY])/(2*DeltaX)); // backward Ny scheme
                            return;
                        }
                        else
                        {
                            // not implemented
                            throw new std::string("Error: not implemented");
                        }
                    }
                    else
                    {
                        // strainStartCell[0]=(uStart[0][globalX][globalY]-uStart[0][globalX-1][globalY])/(DeltaX); // backward Nx scheme
                        strainStartCell[0]=(3.0*uStart[0][globalX][globalY]-4.0*uStart[0][globalX-1][globalY]+uStart[0][globalX-2][globalY])/(2*DeltaX); // three-point backward Nx scheme O(DeltaX^2)
                        // strainStartCell[1]=(uStart[1][globalX][globalY]-uStart[1][globalX][globalY-1])/(DeltaX); // backward Ny scheme
                        strainStartCell[1]=(3.0*uStart[1][globalX][globalY]-4.0*uStart[1][globalX][globalY-1]+uStart[1][globalX][globalY-2])/(2*DeltaX); // three-point backward Ny scheme O(DeltaX^2)
                        // strainStartCell[2]=0.5*((uStart[0][globalX][globalY]-uStart[0][globalX][globalY-1])/(DeltaX) + (uStart[1][globalX][globalY]-uStart[1][globalX-1][globalY])/(DeltaX)); // backward Nx Ny scheme
                        strainStartCell[2]=0.5*((3.0*uStart[0][globalX][globalY]-4.0*uStart[0][globalX][globalY-1]+uStart[0][globalX][globalY-2])/(2*DeltaX) + (3.0*uStart[1][globalX][globalY]-4.0*uStart[1][globalX-1][globalY]+uStart[1][globalX-2][globalY])/(2*DeltaX)); // three-point backward Nx Ny scheme O(DeltaX^2)
                        return;
                    }
                }
            }
            // Ny=0
            if (cellFirstVicinity[2]==true)
            {
                if (PeriodicNx || PeriodicNy)
                {
                    if (PeriodicNx && PeriodicNy)
                    {
                        // not implemented
                        throw new std::string("Error: not implemented");
                    }
                    else if (PeriodicNx)
                    {
                        int Xup = globalX+1;
                        int Xdwn = globalX-1;
                        strainStartCell[0]=(uStart[0][Xup][globalY]-uStart[0][Xdwn][globalY])/(2*DeltaX);
                        strainStartCell[1]=(uStart[1][globalX][globalY+1]-uStart[1][globalX][globalY])/(DeltaX); // forward Ny scheme
                        strainStartCell[2]=0.5*((uStart[0][globalX][globalY+1]-uStart[0][globalX][globalY])/(DeltaX) + (uStart[1][Xup][globalY]-uStart[1][Xdwn][globalY])/(2*DeltaX)); // forward Ny scheme
                        return;
                    }
                    else
                    {
                        // not implemented
                        throw new std::string("Error: not implemented");
                    }
                }
                else
                {
                    strainStartCell[0]=(uStart[0][globalX+1][globalY]-uStart[0][globalX-1][globalY])/(2*DeltaX); // central scheme
                    // strainStartCell[1]=(uStart[1][globalX][globalY+1]-uStart[1][globalX][globalY])/(DeltaX); // forward Ny scheme
                    strainStartCell[1]=(4.0*uStart[1][globalX][globalY+1]-3.0*uStart[1][globalX][globalY]-uStart[1][globalX][globalY+2])/(2*DeltaX); // three-point forward Ny scheme O(DeltaX^2)
                    // strainStartCell[2]=0.5*((uStart[0][globalX][globalY+1]-uStart[0][globalX][globalY])/(DeltaX) + (uStart[1][globalX+1][globalY]-uStart[1][globalX-1][globalY])/(2*DeltaX)); // forward Ny scheme
                    strainStartCell[2]=0.5*((4.0*uStart[0][globalX][globalY+1]-3.0*uStart[0][globalX][globalY]-uStart[0][globalX][globalY+2])/(2*DeltaX) + (uStart[1][globalX+1][globalY]-uStart[1][globalX-1][globalY])/(2*DeltaX)); // three-point forward Ny scheme O(DeltaX^2)
                    return;
                }
            }
            // Ny-1
            if (cellFirstVicinity[3]==true)
            {
                if (PeriodicNx || PeriodicNy)
                {
                    if (PeriodicNx && PeriodicNy)
                    {
                        // not implemented
                        throw new std::string("Error: not implemented");
                    }
                    else if (PeriodicNx)
                    {
                        int Xup = globalX+1;
                        int Xdwn = globalX-1;
                        strainStartCell[0]=(uStart[0][Xup][globalY]-uStart[0][Xdwn][globalY])/(2*DeltaX);
                        strainStartCell[1]=(uStart[1][globalX][globalY]-uStart[1][globalX][globalY-1])/(DeltaX); // backward Ny scheme
                        strainStartCell[2]=0.5*((uStart[0][globalX][globalY]-uStart[0][globalX][globalY-1])/(DeltaX) + (uStart[1][Xup][globalY]-uStart[1][Xdwn][globalY])/(2*DeltaX)); // backward Ny scheme
                        return;
                    }
                    else
                    {
                        // not implemented
                        throw new std::string("Error: not implemented");
                    }
                }
                else
                {
                    strainStartCell[0]=(uStart[0][globalX+1][globalY]-uStart[0][globalX-1][globalY])/(2*DeltaX); // central scheme
                    // strainStartCell[1]=(uStart[1][globalX][globalY]-uStart[1][globalX][globalY-1])/(DeltaX); // backward Ny scheme
                    strainStartCell[1]=(3.0*uStart[1][globalX][globalY]-4.0*uStart[1][globalX][globalY-1]+uStart[1][globalX][globalY-2])/(2*DeltaX); // three-point backward Ny scheme O(DeltaX^2)
                    // strainStartCell[2]=0.5*((uStart[0][globalX][globalY]-uStart[0][globalX][globalY-1])/(DeltaX) + (uStart[1][globalX+1][globalY]-uStart[1][globalX-1][globalY])/(2*DeltaX)); // backward Ny scheme
                    strainStartCell[2]=0.5*((3.0*uStart[0][globalX][globalY]-4.0*uStart[0][globalX][globalY-1]+uStart[0][globalX][globalY-2])/(2*DeltaX) + (uStart[1][globalX+1][globalY]-uStart[1][globalX-1][globalY])/(2*DeltaX)); // three-point backward Ny scheme O(DeltaX^2)
                    return;
                }
            }
        }
    }
    else
    {
        // not implemented for large deformations
        throw new std::string("Error: not implemented: need physical displacement because of non-linear behavior");
    }
}

template <typename T, template <typename U> class Descriptor>
void ComputeCellStrainStress<T,Descriptor>::ComputeStressFromStart(SolidConstitutiveLaw<T>* SolidLaw, int globalX, int globalY, T DeltaX, bool nonLinearLaw)
{
    if (nonLinearLaw==false)
    {
        // get lame coefficients from class
        plb::Array<T,2> LocalLameCoefficients = SolidLaw->GetIsotropicLameCoefficients();

        // local static variables
        T MatB[2][2] = {{strainStartCell[0],strainStartCell[2]},{strainStartCell[2],strainStartCell[1]}};
        T MatC[2][2] = {{0.0,0.0},{0.0,0.0}};
        T MatI[2][2] = {{1.0,0.0},{0.0,1.0}};
        T traceStrain = MatB[0][0] + MatB[1][1];

        // matrix scalar multiplication
        for (int i = 0; i < 2; ++i)
        {
            for (int j = 0; j < 2; ++j)
            {
                MatC[i][j] = LocalLameCoefficients[0]*traceStrain*MatI[i][j] + 2.0*LocalLameCoefficients[1]*MatB[i][j];
            }
        }

        // injection
        stressStartCell[0]=MatC[0][0];
        stressStartCell[1]=MatC[1][1];
        stressStartCell[2]=MatC[0][1];
    }
    else
    {
        // not implemented for non linear law
        throw new std::string("Error: not implemented: need physical deformation because of non-linear law");
    }
}

template <typename T, template <typename U> class Descriptor>
T ComputeCellStrainStress<T,Descriptor>::GetstrainStartCell11() const
{
    return strainStartCell[0];
}

template <typename T, template <typename U> class Descriptor>
T ComputeCellStrainStress<T,Descriptor>::GetstrainStartCell22() const
{
    return strainStartCell[1];
}

template <typename T, template <typename U> class Descriptor>
T ComputeCellStrainStress<T,Descriptor>::GetstrainStartCell12() const
{
    return strainStartCell[2];
}

template <typename T, template <typename U> class Descriptor>
T ComputeCellStrainStress<T,Descriptor>::GetstressStartCell11() const
{
    return stressStartCell[0];
}

template <typename T, template <typename U> class Descriptor>
T ComputeCellStrainStress<T,Descriptor>::GetstressStartCell22() const
{
    return stressStartCell[1];
}

template <typename T, template <typename U> class Descriptor>
T ComputeCellStrainStress<T,Descriptor>::GetstressStartCell12() const
{
    return stressStartCell[2];
}

template <typename T, template <typename U> class Descriptor>
ComputeCellStrainStress<T,Descriptor>::~ComputeCellStrainStress()
{
    // nothing to be done if object is new (a dynamic entity) and if you won't delete other dynamic types, just calling delete(object)
    // if it is a static object, it will be destroyed after run when return 0 is called
    if (!quiet)
    {
        plb::pcout << "Destroying: ~ComputeCellStrainStress" << std::endl;
    }
}

/// compute stress divergence on the cell
/// if useCentralSchemeOnly==true it can yield errors on stress divergence because we need a 2 neighborhood exclusion due to double differentiations
template <typename T, template <typename U> class Descriptor>
ComputeCellStressDivergence<T,Descriptor>::ComputeCellStressDivergence(plb::Array<T,2> divStressStartCell_, plb::Array<bool,4> cellFirstVicinity_, plb::Array<bool,4> cellSecondVicinity_, int GlobalNx_, int GlobalNy_, bool PeriodicNx_, bool PeriodicNy_, bool quiet_)
{
    divStressStartCell=divStressStartCell_;
    cellFirstVicinity=cellFirstVicinity_;
    cellSecondVicinity=cellSecondVicinity_;
    GlobalNx=GlobalNx_;
    GlobalNy=GlobalNy_;
    PeriodicNx=PeriodicNx_;
    PeriodicNy=PeriodicNy_;
    quiet=quiet_;
}

template <typename T, template <typename U> class Descriptor>
void ComputeCellStressDivergence<T,Descriptor>::ComputeDivergenceOfStressFromStart(T*** stressStart, int globalX, int globalY, T DeltaX, bool useCentralSchemeOnly, bool use0X4CentralScheme)
{
    // handle exceptions
    if (cellFirstVicinity[0]==true && cellFirstVicinity[1]==true)
    {
        throw new std::string("Error: mismatching cell vicinity");
    }
    if (cellFirstVicinity[2]==true && cellFirstVicinity[3]==true)
    {
        throw new std::string("Error: mismatching cell vicinity");
    }

    // if no vicinity problems
    if (cellFirstVicinity[0]==false && cellFirstVicinity[1]==false && cellFirstVicinity[2]==false && cellFirstVicinity[3]==false)
    {
        // central difference schemes
        if (cellSecondVicinity[0]==false && cellSecondVicinity[1]==false && cellSecondVicinity[2]==false && cellSecondVicinity[3]==false && use0X4CentralScheme==true)
        {
            // highly accurate central difference scheme + O(DeltaX^4)
            T divStressStartCell0First=(8.0*stressStart[0][globalX+1][globalY]-8.0*stressStart[0][globalX-1][globalY]-stressStart[0][globalX+2][globalY]+stressStart[0][globalX-2][globalY])/(12*DeltaX);
            T divStressStartCell0Second=(8.0*stressStart[2][globalX][globalY+1]-8.0*stressStart[2][globalX][globalY-1]-stressStart[2][globalX][globalY+2]+stressStart[2][globalX][globalY-2])/(12*DeltaX);
            T divStressStartCell1First=(8.0*stressStart[2][globalX+1][globalY]-8.0*stressStart[2][globalX-1][globalY]-stressStart[2][globalX+2][globalY]+stressStart[2][globalX-2][globalY])/(12*DeltaX);
            T divStressStartCell1Second=(8.0*stressStart[1][globalX][globalY+1]-8.0*stressStart[1][globalX][globalY-1]-stressStart[1][globalX][globalY+2]+stressStart[1][globalX][globalY-2])/(12*DeltaX);
            divStressStartCell[0]=divStressStartCell0First+divStressStartCell0Second;
            divStressStartCell[1]=divStressStartCell1First+divStressStartCell1Second;
        }
        else
        {
            // central difference scheme + O(DeltaX^2)
            divStressStartCell[0]=(stressStart[0][globalX+1][globalY]-stressStart[0][globalX-1][globalY])/(2*DeltaX) + (stressStart[2][globalX][globalY+1]-stressStart[2][globalX][globalY-1])/(2*DeltaX);
            divStressStartCell[1]=(stressStart[2][globalX+1][globalY]-stressStart[2][globalX-1][globalY])/(2*DeltaX) + (stressStart[1][globalX][globalY+1]-stressStart[1][globalX][globalY-1])/(2*DeltaX);
        }
        return;
    }
    else if (useCentralSchemeOnly==false) // depending on boundary conditions
    {
        // Nx=0
        if (cellFirstVicinity[0]==true)
        {
            if (cellFirstVicinity[2]==false && cellFirstVicinity[3]==false)
            {
                if (PeriodicNx || PeriodicNy)
                {
                    if (PeriodicNx && PeriodicNy)
                    {
                        // not implemented
                        throw new std::string("Error: not implemented");
                    }
                    else if (PeriodicNx)
                    {
                        int Xup = globalX+1;
                        int Xdwn = GlobalNx-1;
                        int Yup = globalY+1;
                        int Ydwn = globalY-1;
                        divStressStartCell[0]=(stressStart[0][Xup][globalY]-stressStart[0][Xdwn][globalY])/(2*DeltaX) + (stressStart[2][globalX][Yup]-stressStart[2][globalX][Ydwn])/(2*DeltaX);
                        divStressStartCell[1]=(stressStart[2][Xup][globalY]-stressStart[2][Xdwn][globalY])/(2*DeltaX) + (stressStart[1][globalX][Yup]-stressStart[1][globalX][Ydwn])/(2*DeltaX);
                        return;
                    }
                    else
                    {
                        // not implemented
                        throw new std::string("Error: not implemented");
                    }
                }
                else
                {
                    // divStressStartCell[0]=(stressStart[0][globalX+1][globalY]-stressStart[0][globalX][globalY])/(DeltaX) + (stressStart[2][globalX][globalY+1]-stressStart[2][globalX][globalY-1])/(2*DeltaX); // forward Nx
                    // divStressStartCell[1]=(stressStart[2][globalX+1][globalY]-stressStart[2][globalX][globalY])/(DeltaX) + (stressStart[1][globalX][globalY+1]-stressStart[1][globalX][globalY-1])/(2*DeltaX); // forward Nx
                    divStressStartCell[0]=(4.0*stressStart[0][globalX+1][globalY]-3.0*stressStart[0][globalX][globalY]-stressStart[0][globalX+2][globalY])/(2*DeltaX) + (stressStart[2][globalX][globalY+1]-stressStart[2][globalX][globalY-1])/(2*DeltaX); // three-point forward Nx scheme O(DeltaX^2)
                    divStressStartCell[1]=(4.0*stressStart[2][globalX+1][globalY]-3.0*stressStart[2][globalX][globalY]-stressStart[2][globalX+2][globalY])/(2*DeltaX) + (stressStart[1][globalX][globalY+1]-stressStart[1][globalX][globalY-1])/(2*DeltaX); // three-point forward Nx scheme O(DeltaX^2)
                    return;
                }
            }
            else if (cellFirstVicinity[2]==true && cellFirstVicinity[3]==false)
            {
                if (PeriodicNx || PeriodicNy)
                {
                    if (PeriodicNx && PeriodicNy)
                    {
                        // not implemented
                        throw new std::string("Error: not implemented");
                    }
                    else if (PeriodicNx)
                    {
                        int Xup = globalX+1;
                        int Xdwn = GlobalNx-1;
                        divStressStartCell[0]=(stressStart[0][Xup][globalY]-stressStart[0][Xdwn][globalY])/(2*DeltaX) + (stressStart[2][globalX][globalY+1]-stressStart[2][globalX][globalY])/(DeltaX); // forward Ny
                        divStressStartCell[1]=(stressStart[2][Xup][globalY]-stressStart[2][Xdwn][globalY])/(2*DeltaX) + (stressStart[1][globalX][globalY+1]-stressStart[1][globalX][globalY])/(DeltaX); // forward Ny
                        return;
                    }
                    else
                    {
                        // not implemented
                        throw new std::string("Error: not implemented");
                    }
                }
                else
                {
                    // divStressStartCell[0]=(stressStart[0][globalX+1][globalY]-stressStart[0][globalX][globalY])/(DeltaX) + (stressStart[2][globalX][globalY+1]-stressStart[2][globalX][globalY])/(DeltaX); // forward Nx Ny
                    // divStressStartCell[1]=(stressStart[2][globalX+1][globalY]-stressStart[2][globalX][globalY])/(DeltaX) + (stressStart[1][globalX][globalY+1]-stressStart[1][globalX][globalY])/(DeltaX); // forward Nx Ny
                    divStressStartCell[0]=(4.0*stressStart[0][globalX+1][globalY]-3.0*stressStart[0][globalX][globalY]-stressStart[0][globalX+2][globalY])/(2*DeltaX) + (4.0*stressStart[2][globalX][globalY+1]-3.0*stressStart[2][globalX][globalY]-stressStart[2][globalX][globalY+2])/(2*DeltaX); // three-point forward Nx Ny scheme O(DeltaX^2)
                    divStressStartCell[1]=(4.0*stressStart[2][globalX+1][globalY]-3.0*stressStart[2][globalX][globalY]-stressStart[2][globalX+2][globalY])/(2*DeltaX) + (4.0*stressStart[1][globalX][globalY+1]-3.0*stressStart[1][globalX][globalY]-stressStart[1][globalX][globalY+2])/(2*DeltaX); // three-point forward Nx Ny scheme O(DeltaX^2)
                    return;
                }
            }
            else if (cellFirstVicinity[2]==false && cellFirstVicinity[3]==true)
            {
                if (PeriodicNx || PeriodicNy)
                {
                    if (PeriodicNx && PeriodicNy)
                    {
                        // not implemented
                        throw new std::string("Error: not implemented");
                    }
                    else if (PeriodicNx)
                    {
                        int Xup = globalX+1;
                        int Xdwn = GlobalNx-1;
                        divStressStartCell[0]=(stressStart[0][Xup][globalY]-stressStart[0][Xdwn][globalY])/(2*DeltaX) + (stressStart[2][globalX][globalY]-stressStart[2][globalX][globalY-1])/(DeltaX); // backward Ny
                        divStressStartCell[1]=(stressStart[2][Xup][globalY]-stressStart[2][Xdwn][globalY])/(2*DeltaX) + (stressStart[1][globalX][globalY]-stressStart[1][globalX][globalY-1])/(DeltaX); // backward Ny
                        return;
                    }
                    else
                    {
                        // not implemented
                        throw new std::string("Error: not implemented");
                    }
                }
                else
                {
                    // divStressStartCell[0]=(stressStart[0][globalX+1][globalY]-stressStart[0][globalX][globalY])/(DeltaX) + (stressStart[2][globalX][globalY]-stressStart[2][globalX][globalY-1])/(DeltaX); // forward Nx backward Ny
                    // divStressStartCell[1]=(stressStart[2][globalX+1][globalY]-stressStart[2][globalX][globalY])/(DeltaX) + (stressStart[1][globalX][globalY]-stressStart[1][globalX][globalY-1])/(DeltaX); // forward Nx backward Ny
                    divStressStartCell[0]=(4.0*stressStart[0][globalX+1][globalY]-3.0*stressStart[0][globalX][globalY]-stressStart[0][globalX+2][globalY])/(2*DeltaX) + (3.0*stressStart[2][globalX][globalY]-4.0*stressStart[2][globalX][globalY-1]+stressStart[2][globalX][globalY-2])/(2*DeltaX); // three-point forward Nx backward Ny scheme O(DeltaX^2)
                    divStressStartCell[1]=(4.0*stressStart[2][globalX+1][globalY]-3.0*stressStart[2][globalX][globalY]-stressStart[2][globalX+2][globalY])/(2*DeltaX) + (3.0*stressStart[1][globalX][globalY]-4.0*stressStart[1][globalX][globalY-1]+stressStart[1][globalX][globalY-2])/(2*DeltaX); // three-point forward Nx backward Ny scheme O(DeltaX^2)
                    return;
                }
            }
        }
        // Nx-1
        if (cellFirstVicinity[1]==true)
        {
            if (cellFirstVicinity[2]==false && cellFirstVicinity[3]==false)
            {
                if (PeriodicNx || PeriodicNy)
                {
                    if (PeriodicNx && PeriodicNy)
                    {
                        // not implemented
                        throw new std::string("Error: not implemented");
                    }
                    else if (PeriodicNx)
                    {
                        int Xup = 0;
                        int Xdwn = globalX-1;
                        int Yup = globalY+1;
                        int Ydwn = globalY-1;
                        divStressStartCell[0]=(stressStart[0][Xup][globalY]-stressStart[0][Xdwn][globalY])/(2*DeltaX) + (stressStart[2][globalX][Yup]-stressStart[2][globalX][Ydwn])/(2*DeltaX);
                        divStressStartCell[1]=(stressStart[2][Xup][globalY]-stressStart[2][Xdwn][globalY])/(2*DeltaX) + (stressStart[1][globalX][Yup]-stressStart[1][globalX][Ydwn])/(2*DeltaX);
                        return;
                    }
                    else
                    {
                        // not implemented
                        throw new std::string("Error: not implemented");
                    }
                }
                else
                {
                    // divStressStartCell[0]=(stressStart[0][globalX][globalY]-stressStart[0][globalX-1][globalY])/(DeltaX) + (stressStart[2][globalX][globalY+1]-stressStart[2][globalX][globalY-1])/(2*DeltaX); // backward Nx
                    // divStressStartCell[1]=(stressStart[2][globalX][globalY]-stressStart[2][globalX-1][globalY])/(DeltaX) + (stressStart[1][globalX][globalY+1]-stressStart[1][globalX][globalY-1])/(2*DeltaX); // backward Nx
                    divStressStartCell[0]=(3.0*stressStart[0][globalX][globalY]-4.0*stressStart[0][globalX-1][globalY]+stressStart[0][globalX-2][globalY])/(2*DeltaX) + (stressStart[2][globalX][globalY+1]-stressStart[2][globalX][globalY-1])/(2*DeltaX); // three-point backward Nx scheme O(DeltaX^2)
                    divStressStartCell[1]=(3.0*stressStart[2][globalX][globalY]-4.0*stressStart[2][globalX-1][globalY]+stressStart[2][globalX-2][globalY])/(2*DeltaX) + (stressStart[1][globalX][globalY+1]-stressStart[1][globalX][globalY-1])/(2*DeltaX); // three-point backward Nx scheme O(DeltaX^2)
                    return;
                }
            }
            else if (cellFirstVicinity[2]==true && cellFirstVicinity[3]==false)
            {
                if (PeriodicNx || PeriodicNy)
                {
                    if (PeriodicNx && PeriodicNy)
                    {
                        // not implemented
                        throw new std::string("Error: not implemented");
                    }
                    else if (PeriodicNx)
                    {
                        int Xup = 0;
                        int Xdwn = globalX-1;
                        divStressStartCell[0]=(stressStart[0][Xup][globalY]-stressStart[0][Xdwn][globalY])/(2*DeltaX) + (stressStart[2][globalX][globalY+1]-stressStart[2][globalX][globalY])/(DeltaX); // forward Ny
                        divStressStartCell[1]=(stressStart[2][Xup][globalY]-stressStart[2][Xdwn][globalY])/(2*DeltaX) + (stressStart[1][globalX][globalY+1]-stressStart[1][globalX][globalY])/(DeltaX); // forward Ny
                        return;
                    }
                    else
                    {
                        // not implemented
                        throw new std::string("Error: not implemented");
                    }
                }
                else
                {
                    // divStressStartCell[0]=(stressStart[0][globalX][globalY]-stressStart[0][globalX-1][globalY])/(DeltaX) + (stressStart[2][globalX][globalY+1]-stressStart[2][globalX][globalY])/(DeltaX); // backward Nx forward Ny
                    // divStressStartCell[1]=(stressStart[2][globalX][globalY]-stressStart[2][globalX-1][globalY])/(DeltaX) + (stressStart[1][globalX][globalY+1]-stressStart[1][globalX][globalY])/(DeltaX); // backward Nx forward Ny
                    divStressStartCell[0]=(3.0*stressStart[0][globalX][globalY]-4.0*stressStart[0][globalX-1][globalY]+stressStart[0][globalX-2][globalY])/(2*DeltaX) + (4.0*stressStart[2][globalX][globalY+1]-3.0*stressStart[2][globalX][globalY]-stressStart[2][globalX][globalY+2])/(2*DeltaX); // three-point backward Nx forward Ny scheme O(DeltaX^2)
                    divStressStartCell[1]=(3.0*stressStart[2][globalX][globalY]-4.0*stressStart[2][globalX-1][globalY]+stressStart[2][globalX-2][globalY])/(2*DeltaX) + (4.0*stressStart[1][globalX][globalY+1]-3.0*stressStart[1][globalX][globalY]-stressStart[1][globalX][globalY+2])/(2*DeltaX); // three-point backward Nx forward Ny scheme O(DeltaX^2)
                    return;
                }
            }
            else if (cellFirstVicinity[2]==false && cellFirstVicinity[3]==true)
            {
                if (PeriodicNx || PeriodicNy)
                {
                    if (PeriodicNx && PeriodicNy)
                    {
                        // not implemented
                        throw new std::string("Error: not implemented");
                    }
                    else if (PeriodicNx)
                    {
                        int Xup = 0;
                        int Xdwn = globalX-1;
                        divStressStartCell[0]=(stressStart[0][Xup][globalY]-stressStart[0][Xdwn][globalY])/(2*DeltaX) + (stressStart[2][globalX][globalY]-stressStart[2][globalX][globalY-1])/(DeltaX); // backward Ny
                        divStressStartCell[1]=(stressStart[2][Xup][globalY]-stressStart[2][Xdwn][globalY])/(2*DeltaX) + (stressStart[1][globalX][globalY]-stressStart[1][globalX][globalY-1])/(DeltaX); // backward Ny
                        return;
                    }
                    else
                    {
                        // not implemented
                        throw new std::string("Error: not implemented");
                    }
                }
                else
                {
                    // divStressStartCell[0]=(stressStart[0][globalX][globalY]-stressStart[0][globalX-1][globalY])/(DeltaX) + (stressStart[2][globalX][globalY]-stressStart[2][globalX][globalY-1])/(DeltaX); // backward Nx Ny
                    // divStressStartCell[1]=(stressStart[2][globalX][globalY]-stressStart[2][globalX-1][globalY])/(DeltaX) + (stressStart[1][globalX][globalY]-stressStart[1][globalX][globalY-1])/(DeltaX); // backward Nx Ny
                    divStressStartCell[0]=(3.0*stressStart[0][globalX][globalY]-4.0*stressStart[0][globalX-1][globalY]+stressStart[0][globalX-2][globalY])/(2*DeltaX) + (3.0*stressStart[2][globalX][globalY]-4.0*stressStart[2][globalX][globalY-1]+stressStart[2][globalX][globalY-2])/(2*DeltaX); // three-point backward Nx Ny scheme O(DeltaX^2)
                    divStressStartCell[1]=(3.0*stressStart[2][globalX][globalY]-4.0*stressStart[2][globalX-1][globalY]+stressStart[2][globalX-2][globalY])/(2*DeltaX) + (3.0*stressStart[1][globalX][globalY]-4.0*stressStart[1][globalX][globalY-1]+stressStart[1][globalX][globalY-2])/(2*DeltaX); // three-point backward Nx Ny scheme O(DeltaX^2)
                    return;
                }
            }
        }
        // Ny=0
        if (cellFirstVicinity[2]==true)
        {
            if (PeriodicNx || PeriodicNy)
            {
                if (PeriodicNx && PeriodicNy)
                {
                    // not implemented
                    throw new std::string("Error: not implemented");
                }
                else if (PeriodicNx)
                {
                    int Xup = globalX+1;
                    int Xdwn = globalX-1;
                    divStressStartCell[0]=(stressStart[0][Xup][globalY]-stressStart[0][Xdwn][globalY])/(2*DeltaX) + (stressStart[2][globalX][globalY+1]-stressStart[2][globalX][globalY])/(DeltaX); // forward Ny
                    divStressStartCell[1]=(stressStart[2][Xup][globalY]-stressStart[2][Xdwn][globalY])/(2*DeltaX) + (stressStart[1][globalX][globalY+1]-stressStart[1][globalX][globalY])/(DeltaX); // forward Ny
                    return;
                }
                else
                {
                    // not implemented
                    throw new std::string("Error: not implemented");
                }
            }
            else
            {
                // divStressStartCell[0]=(stressStart[0][globalX+1][globalY]-stressStart[0][globalX-1][globalY])/(2*DeltaX) + (stressStart[2][globalX][globalY+1]-stressStart[2][globalX][globalY])/(DeltaX); // forward Ny
                // divStressStartCell[1]=(stressStart[2][globalX+1][globalY]-stressStart[2][globalX-1][globalY])/(2*DeltaX) + (stressStart[1][globalX][globalY+1]-stressStart[1][globalX][globalY])/(DeltaX); // forward Ny
                divStressStartCell[0]=(stressStart[0][globalX+1][globalY]-stressStart[0][globalX-1][globalY])/(2*DeltaX) + (4.0*stressStart[2][globalX][globalY+1]-3.0*stressStart[2][globalX][globalY]-stressStart[2][globalX][globalY+2])/(2*DeltaX); // three-point forward Ny scheme O(DeltaX^2)
                divStressStartCell[1]=(stressStart[2][globalX+1][globalY]-stressStart[2][globalX-1][globalY])/(2*DeltaX) + (4.0*stressStart[1][globalX][globalY+1]-3.0*stressStart[1][globalX][globalY]-stressStart[1][globalX][globalY+2])/(2*DeltaX); // three-point forward Ny scheme O(DeltaX^2)
                return;
            }
        }
        // Ny-1
        if (cellFirstVicinity[3]==true)
        {
            if (PeriodicNx || PeriodicNy)
            {
                if (PeriodicNx && PeriodicNy)
                {
                    // not implemented
                    throw new std::string("Error: not implemented");
                }
                else if (PeriodicNx)
                {
                    int Xup = globalX+1;
                    int Xdwn = globalX-1;
                    divStressStartCell[0]=(stressStart[0][Xup][globalY]-stressStart[0][Xdwn][globalY])/(2*DeltaX) + (stressStart[2][globalX][globalY]-stressStart[2][globalX][globalY-1])/(DeltaX); // backward Ny
                    divStressStartCell[1]=(stressStart[2][Xup][globalY]-stressStart[2][Xdwn][globalY])/(2*DeltaX) + (stressStart[1][globalX][globalY]-stressStart[1][globalX][globalY-1])/(DeltaX); // backward Ny
                    return;
                }
                else
                {
                    // not implemented
                    throw new std::string("Error: not implemented");
                }
            }
            else
            {
                // divStressStartCell[0]=(stressStart[0][globalX+1][globalY]-stressStart[0][globalX-1][globalY])/(2*DeltaX) + (stressStart[2][globalX][globalY]-stressStart[2][globalX][globalY-1])/(DeltaX); // backward Ny
                // divStressStartCell[1]=(stressStart[2][globalX+1][globalY]-stressStart[2][globalX-1][globalY])/(2*DeltaX) + (stressStart[1][globalX][globalY]-stressStart[1][globalX][globalY-1])/(DeltaX); // backward Ny
                divStressStartCell[0]=(stressStart[0][globalX+1][globalY]-stressStart[0][globalX-1][globalY])/(2*DeltaX) + (3.0*stressStart[2][globalX][globalY]-4.0*stressStart[2][globalX][globalY-1]+stressStart[2][globalX][globalY-2])/(2*DeltaX); // three-point backward Ny scheme O(DeltaX^2)
                divStressStartCell[1]=(stressStart[2][globalX+1][globalY]-stressStart[2][globalX-1][globalY])/(2*DeltaX) + (3.0*stressStart[1][globalX][globalY]-4.0*stressStart[1][globalX][globalY-1]+stressStart[1][globalX][globalY-2])/(2*DeltaX); // three-point backward Ny scheme O(DeltaX^2)
                return;
            }
        }
    }
}

template <typename T, template <typename U> class Descriptor>
T ComputeCellStressDivergence<T,Descriptor>::GetdivStressStartCellX() const
{
    return divStressStartCell[0];
}

template <typename T, template <typename U> class Descriptor>
T ComputeCellStressDivergence<T,Descriptor>::GetdivStressStartCellY() const
{
    return divStressStartCell[1];
}

template <typename T, template <typename U> class Descriptor>
ComputeCellStressDivergence<T,Descriptor>::~ComputeCellStressDivergence()
{
    // nothing to be done if object is new (a dynamic entity) and if you won't delete other dynamic types, just calling delete(object)
    // if it is a static object, it will be destroyed after run when return 0 is called
    if (!quiet)
    {
        plb::pcout << "Destroying: ~ComputeCellStressDivergence" << std::endl;
    }
}

#endif // SOLIDS_CELLS_OPERATIONS_CPP

