///Created By T.MAQUART

/*
Solids boundary conditions.
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

#ifndef SOLIDS_BOUNDARY_CONDITIONS_H
#define SOLIDS_BOUNDARY_CONDITIONS_H

#include "palabos2D.h"
#include "palabos2D.hh"

#include <math.h> // round

// using namespaces in an other namespace is not recommended due to potential conflicts

namespace sld {

/// data processor: data processing functional for boundary velocity condition
template<typename T, template<typename U> class Descriptor>
class BoundaryVelocityConditionProcessor : public plb::BoxProcessingFunctional2D_L<T,Descriptor>{
public:
    BoundaryVelocityConditionProcessor(T speedlatticeUnits_, T rho_, int direction_, bool maxNxNy_, int includedMin_, int excludedMax_, int Nx_, int Ny_, bool gaussHermite_, bool triangularOption_)
    {
        speedlatticeUnits=speedlatticeUnits_;
        rho=rho_;
        direction=direction_;
        maxNxNy=maxNxNy_;
        includedMin=includedMin_;
        excludedMax=excludedMax_;
        Nx=Nx_;
        Ny=Ny_;
        gaussHermite=gaussHermite_;
        triangularOption=triangularOption_;
    }
    virtual void process(plb::Box2D domain, plb::BlockLattice2D<T,Descriptor>& lattice)
    {
        // access the position of the atomic-block inside the multi-block
        plb::Dot2D relativePosition = lattice.getLocation();

        for (plb::plint iX=domain.x0; iX<=domain.x1; ++iX)
        {
            for (plb::plint iY=domain.y0; iY<=domain.y1; ++iY)
            {
                // cell
                plb::Cell<T,Descriptor>& cell = lattice.get(iX,iY); // need a reference to modify cell attributes

                // convert local coordinates to global ones
                plb::plint globalX = iX + relativePosition.x;
                plb::plint globalY = iY + relativePosition.y;

                // deal with values
                T rhoBar;
                plb::Array<T,9> Whermite;
                if (gaussHermite==true)
                {
                    rhoBar = (rho-1.0);
                    Whermite[0]=4.0/9.0;
                    Whermite[1]=1.0/36.0;
                    Whermite[2]=1.0/9.0;
                    Whermite[3]=1.0/36.0;
                    Whermite[4]=1.0/9.0;
                    Whermite[5]=1.0/36.0;
                    Whermite[6]=1.0/9.0;
                    Whermite[7]=1.0/36.0;
                    Whermite[8]=1.0/9.0;
                }
                else
                {
                    rhoBar = (rho-1.0)/9.0;
                    Whermite[0]=1.0;
                    Whermite[1]=1.0;
                    Whermite[2]=1.0;
                    Whermite[3]=1.0;
                    Whermite[4]=1.0;
                    Whermite[5]=1.0;
                    Whermite[6]=1.0;
                    Whermite[7]=1.0;
                    Whermite[8]=1.0;
                }

                // deal with boundary conditions
                // note: mass conservative boundary conditions
                switch (direction)
                {
                    case 0:
                        if (maxNxNy)
                        {
                            if (globalX==(Nx-1))
                            {
                                if ((globalY>=includedMin) && (globalY<excludedMax))
                                {
                                    if (triangularOption==true)
                                    {
                                        int maxTriangle=std::round((((excludedMax-1.0)-includedMin)/2.0) + includedMin);
                                        T triangularValue;
                                        if (globalY<=maxTriangle)
                                        {
                                            triangularValue = ((globalY-includedMin)/((((excludedMax-1.0)-includedMin)/2.0)))*speedlatticeUnits;
                                            T eps=maxTriangle*pow(10.0,-12.0);
                                            if (triangularValue>(0.0-eps) && triangularValue<(0.0+eps))
                                            {
                                                triangularValue = ((1.0)/((((excludedMax-1.0)-includedMin)/2.0)))*speedlatticeUnits;
                                            }
                                        }
                                        else
                                        {
                                            triangularValue = (std::abs((globalY-(excludedMax-1.0)))/((((excludedMax-1.0)-includedMin)/2.0)))*speedlatticeUnits;
                                            T eps=maxTriangle*pow(10.0,-12.0);
                                            if (triangularValue>(0.0-eps) && triangularValue<(0.0+eps))
                                            {
                                                triangularValue = (std::abs((1.0))/((((excludedMax-1.0)-includedMin)/2.0)))*speedlatticeUnits;
                                            }
                                        }
                                        cell[1]=rhoBar*Whermite[1]+triangularValue*0.70;
                                        cell[2]=rhoBar*Whermite[2]+triangularValue;
                                        cell[3]=rhoBar*Whermite[3]+triangularValue*0.70;
                                        cell[4]=rhoBar*Whermite[4];
                                        cell[5]=rhoBar*Whermite[5]-triangularValue*0.70;
                                        cell[6]=rhoBar*Whermite[6]-triangularValue;
                                        cell[7]=rhoBar*Whermite[7]-triangularValue*0.70;
                                        cell[8]=rhoBar*Whermite[8];
                                    }
                                    else
                                    {
                                        cell[1]=rhoBar*Whermite[1]+speedlatticeUnits*0.70;
                                        cell[2]=rhoBar*Whermite[2]+speedlatticeUnits;
                                        cell[3]=rhoBar*Whermite[3]+speedlatticeUnits*0.70;
                                        cell[4]=rhoBar*Whermite[4];
                                        cell[5]=rhoBar*Whermite[5]-speedlatticeUnits*0.70;
                                        cell[6]=rhoBar*Whermite[6]-speedlatticeUnits;
                                        cell[7]=rhoBar*Whermite[7]-speedlatticeUnits*0.70;
                                        cell[8]=rhoBar*Whermite[8];
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
                                    if (triangularOption==true)
                                    {
                                        int maxTriangle=std::round((((excludedMax-1.0)-includedMin)/2.0) + includedMin);
                                        T triangularValue;
                                        if (globalY<=maxTriangle)
                                        {
                                            triangularValue = ((globalY-includedMin)/((((excludedMax-1.0)-includedMin)/2.0)))*speedlatticeUnits;
                                            T eps=maxTriangle*pow(10.0,-12.0);
                                            if (triangularValue>(0.0-eps) && triangularValue<(0.0+eps))
                                            {
                                                triangularValue = ((1.0)/((((excludedMax-1.0)-includedMin)/2.0)))*speedlatticeUnits;
                                            }
                                        }
                                        else
                                        {
                                            triangularValue = (std::abs((globalY-(excludedMax-1.0)))/((((excludedMax-1.0)-includedMin)/2.0)))*speedlatticeUnits;
                                            T eps=maxTriangle*pow(10.0,-12.0);
                                            if (triangularValue>(0.0-eps) && triangularValue<(0.0+eps))
                                            {
                                                triangularValue = (std::abs((1.0))/((((excludedMax-1.0)-includedMin)/2.0)))*speedlatticeUnits;
                                            }
                                        }
                                        cell[1]=rhoBar*Whermite[1]-triangularValue*0.70;
                                        cell[2]=rhoBar*Whermite[2]-triangularValue;
                                        cell[3]=rhoBar*Whermite[3]-triangularValue*0.70;
                                        cell[4]=rhoBar*Whermite[4];
                                        cell[5]=rhoBar*Whermite[5]+triangularValue*0.70;
                                        cell[6]=rhoBar*Whermite[6]+triangularValue;
                                        cell[7]=rhoBar*Whermite[7]+triangularValue*0.70;
                                        cell[8]=rhoBar*Whermite[8];
                                    }
                                    else
                                    {
                                        cell[1]=rhoBar*Whermite[1]-speedlatticeUnits*0.70;
                                        cell[2]=rhoBar*Whermite[2]-speedlatticeUnits;
                                        cell[3]=rhoBar*Whermite[3]-speedlatticeUnits*0.70;
                                        cell[4]=rhoBar*Whermite[4];
                                        cell[5]=rhoBar*Whermite[5]+speedlatticeUnits*0.70;
                                        cell[6]=rhoBar*Whermite[6]+speedlatticeUnits;
                                        cell[7]=rhoBar*Whermite[7]+speedlatticeUnits*0.70;
                                        cell[8]=rhoBar*Whermite[8];
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
                                    if (triangularOption==true)
                                    {
                                        int maxTriangle=std::round((((excludedMax-1.0)-includedMin)/2.0) + includedMin);
                                        T triangularValue;
                                        if (globalX<=maxTriangle)
                                        {
                                            triangularValue = ((globalX-includedMin)/((((excludedMax-1.0)-includedMin)/2.0)))*speedlatticeUnits;
                                            T eps=maxTriangle*pow(10.0,-12.0);
                                            if (triangularValue>(0.0-eps) && triangularValue<(0.0+eps))
                                            {
                                                triangularValue = ((1.0)/((((excludedMax-1.0)-includedMin)/2.0)))*speedlatticeUnits;
                                            }
                                        }
                                        else
                                        {
                                            triangularValue = (std::abs((globalX-(excludedMax-1.0)))/((((excludedMax-1.0)-includedMin)/2.0)))*speedlatticeUnits;
                                            T eps=maxTriangle*pow(10.0,-12.0);
                                            if (triangularValue>(0.0-eps) && triangularValue<(0.0+eps))
                                            {
                                                triangularValue = (std::abs((1.0))/((((excludedMax-1.0)-includedMin)/2.0)))*speedlatticeUnits;
                                            }
                                        }
                                        cell[1]=rhoBar*Whermite[1]-triangularValue*0.70;
                                        cell[2]=rhoBar*Whermite[2];
                                        cell[3]=rhoBar*Whermite[3]+triangularValue*0.70;
                                        cell[4]=rhoBar*Whermite[4]+triangularValue;
                                        cell[5]=rhoBar*Whermite[5]+triangularValue*0.70;
                                        cell[6]=rhoBar*Whermite[6];
                                        cell[7]=rhoBar*Whermite[7]-triangularValue*0.70;
                                        cell[8]=rhoBar*Whermite[8]-triangularValue;
                                    }
                                    else
                                    {
                                        cell[1]=rhoBar*Whermite[1]-speedlatticeUnits*0.70;
                                        cell[2]=rhoBar*Whermite[2];
                                        cell[3]=rhoBar*Whermite[3]+speedlatticeUnits*0.70;
                                        cell[4]=rhoBar*Whermite[4]+speedlatticeUnits;
                                        cell[5]=rhoBar*Whermite[5]+speedlatticeUnits*0.70;
                                        cell[6]=rhoBar*Whermite[6];
                                        cell[7]=rhoBar*Whermite[7]-speedlatticeUnits*0.70;
                                        cell[8]=rhoBar*Whermite[8]-speedlatticeUnits;
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
                                    if (triangularOption==true)
                                    {
                                        int maxTriangle=std::round((((excludedMax-1.0)-includedMin)/2.0) + includedMin);
                                        T triangularValue;
                                        if (globalX<=maxTriangle)
                                        {
                                            triangularValue = ((globalX-includedMin)/((((excludedMax-1.0)-includedMin)/2.0)))*speedlatticeUnits;
                                            T eps=maxTriangle*pow(10.0,-12.0);
                                            if (triangularValue>(0.0-eps) && triangularValue<(0.0+eps))
                                            {
                                                triangularValue = ((1.0)/((((excludedMax-1.0)-includedMin)/2.0)))*speedlatticeUnits;
                                            }
                                        }
                                        else
                                        {
                                            triangularValue = (std::abs((globalX-(excludedMax-1.0)))/((((excludedMax-1.0)-includedMin)/2.0)))*speedlatticeUnits;
                                            T eps=maxTriangle*pow(10.0,-12.0);
                                            if (triangularValue>(0.0-eps) && triangularValue<(0.0+eps))
                                            {
                                                triangularValue = (std::abs((1.0))/((((excludedMax-1.0)-includedMin)/2.0)))*speedlatticeUnits;
                                            }
                                        }
                                        cell[1]=rhoBar*Whermite[1]+triangularValue*0.70;
                                        cell[2]=rhoBar*Whermite[2];
                                        cell[3]=rhoBar*Whermite[3]-triangularValue*0.70;
                                        cell[4]=rhoBar*Whermite[4]-triangularValue;
                                        cell[5]=rhoBar*Whermite[5]-triangularValue*0.70;
                                        cell[6]=rhoBar*Whermite[6];
                                        cell[7]=rhoBar*Whermite[7]+triangularValue*0.70;
                                        cell[8]=rhoBar*Whermite[8]+triangularValue;
                                    }
                                    else
                                    {
                                        cell[1]=rhoBar*Whermite[1]+speedlatticeUnits*0.70;
                                        cell[2]=rhoBar*Whermite[2];
                                        cell[3]=rhoBar*Whermite[3]-speedlatticeUnits*0.70;
                                        cell[4]=rhoBar*Whermite[4]-speedlatticeUnits;
                                        cell[5]=rhoBar*Whermite[5]-speedlatticeUnits*0.70;
                                        cell[6]=rhoBar*Whermite[6];
                                        cell[7]=rhoBar*Whermite[7]+speedlatticeUnits*0.70;
                                        cell[8]=rhoBar*Whermite[8]+speedlatticeUnits;
                                    }
                                }
                            }
                        }
                }
            }
        }
    }
    virtual void getModificationPattern(std::vector<bool>& isWritten) const {
        isWritten[0] = true;
    }
    virtual void getTypeOfModification(std::vector<plb::modif::ModifT>& modified) const {
        modified[0] = plb::modif::allVariables; // à comprendre à quoi cela sert voir 16.3.2 palabos user guide
        // indicates what kind of cell content was modified and must be updated in a multi-block structure
    }
    virtual BoundaryVelocityConditionProcessor<T,Descriptor>* clone() const
    {
        return new BoundaryVelocityConditionProcessor<T,Descriptor>(*this);
    }
    // le destructeur semble hérité mais quand est appelé le delete ?
private:
    T speedlatticeUnits;
    T rho;
    int direction;
    bool maxNxNy;
    int includedMin;
    int excludedMax;
    int Nx;
    int Ny;
    bool gaussHermite;
    bool triangularOption;
};

/// boundary velocity condition without cell indexation
template<typename T, template<typename U> class Descriptor>
class BoundaryVelocityCondition {
public:
    BoundaryVelocityCondition(T uspeedf_, T rho_, bool gaussHermite_, bool quiet_=true)
    {
        uspeedf=uspeedf_;
        rho=rho_;
        gaussHermite=gaussHermite_;
        quiet=quiet_;
    }
    void ApplyVelocityBoundaryCondition(plb::MultiBlockLattice2D<T,Descriptor>& lattice, int direction, bool maxNxNy, int includedMin, int excludedMax, bool triangularOption) const
    {
        T rhoBar;
        plb::Array<T,9> Whermite;
        if (gaussHermite==true)
        {
            rhoBar = (rho-1.0);
            Whermite[0]=4.0/9.0;
            Whermite[1]=1.0/36.0;
            Whermite[2]=1.0/9.0;
            Whermite[3]=1.0/36.0;
            Whermite[4]=1.0/9.0;
            Whermite[5]=1.0/36.0;
            Whermite[6]=1.0/9.0;
            Whermite[7]=1.0/36.0;
            Whermite[8]=1.0/9.0;
        }
        else
        {
            rhoBar = (rho-1.0)/9.0;
            Whermite[0]=1.0;
            Whermite[1]=1.0;
            Whermite[2]=1.0;
            Whermite[3]=1.0;
            Whermite[4]=1.0;
            Whermite[5]=1.0;
            Whermite[6]=1.0;
            Whermite[7]=1.0;
            Whermite[8]=1.0;
        }

        switch (direction)
        {
            case 0:
                if (maxNxNy)
                {
                    // velocity boundary condition in negative x direction
                    // note: mass conservative boundary conditions
                    for(int h=includedMin;h<excludedMax;h++)
                    {
                        if (triangularOption==true)
                        {
                            int maxTriangle=std::round((((excludedMax-1.0)-includedMin)/2.0) + includedMin);
                            T triangularValue;
                            if (h<=maxTriangle)
                            {
                                triangularValue = ((h-includedMin)/((((excludedMax-1.0)-includedMin)/2.0)))*uspeedf;
                                T eps=maxTriangle*pow(10.0,-12.0);
                                if (triangularValue>(0.0-eps) && triangularValue<(0.0+eps))
                                {
                                    triangularValue = ((1.0)/((((excludedMax-1.0)-includedMin)/2.0)))*uspeedf;
                                }
                            }
                            else
                            {
                                triangularValue = (std::abs((h-(excludedMax-1.0)))/((((excludedMax-1.0)-includedMin)/2.0)))*uspeedf;
                                T eps=maxTriangle*pow(10.0,-12.0);
                                if (triangularValue>(0.0-eps) && triangularValue<(0.0+eps))
                                {
                                    triangularValue = (std::abs((1.0))/((((excludedMax-1.0)-includedMin)/2.0)))*uspeedf;
                                }
                            }
                            lattice.get(lattice.getNx()-1,h)[1]=rhoBar*Whermite[1]+triangularValue*0.70;
                            lattice.get(lattice.getNx()-1,h)[2]=rhoBar*Whermite[2]+triangularValue;
                            lattice.get(lattice.getNx()-1,h)[3]=rhoBar*Whermite[3]+triangularValue*0.70;
                            lattice.get(lattice.getNx()-1,h)[4]=rhoBar*Whermite[4];
                            lattice.get(lattice.getNx()-1,h)[5]=rhoBar*Whermite[5]-triangularValue*0.70;
                            lattice.get(lattice.getNx()-1,h)[6]=rhoBar*Whermite[6]-triangularValue;
                            lattice.get(lattice.getNx()-1,h)[7]=rhoBar*Whermite[7]-triangularValue*0.70;
                            lattice.get(lattice.getNx()-1,h)[8]=rhoBar*Whermite[8];
                        }
                        else
                        {
                            lattice.get(lattice.getNx()-1,h)[1]=rhoBar*Whermite[1]+uspeedf*0.70;
                            lattice.get(lattice.getNx()-1,h)[2]=rhoBar*Whermite[2]+uspeedf;
                            lattice.get(lattice.getNx()-1,h)[3]=rhoBar*Whermite[3]+uspeedf*0.70;
                            lattice.get(lattice.getNx()-1,h)[4]=rhoBar*Whermite[4];
                            lattice.get(lattice.getNx()-1,h)[5]=rhoBar*Whermite[5]-uspeedf*0.70;
                            lattice.get(lattice.getNx()-1,h)[6]=rhoBar*Whermite[6]-uspeedf;
                            lattice.get(lattice.getNx()-1,h)[7]=rhoBar*Whermite[7]-uspeedf*0.70;
                            lattice.get(lattice.getNx()-1,h)[8]=rhoBar*Whermite[8];
                        }
                    }
                }
                else
                {
                    // velocity boundary condition in positive x direction
                    // note: mass conservative boundary conditions
                    for(int h=includedMin;h<excludedMax;h++)
                    {
                        if (triangularOption==true)
                        {
                            int maxTriangle=std::round((((excludedMax-1.0)-includedMin)/2.0) + includedMin);
                            T triangularValue;
                            if (h<=maxTriangle)
                            {
                                triangularValue = ((h-includedMin)/((((excludedMax-1.0)-includedMin)/2.0)))*uspeedf;
                                T eps=maxTriangle*pow(10.0,-12.0);
                                if (triangularValue>(0.0-eps) && triangularValue<(0.0+eps))
                                {
                                    triangularValue = ((1.0)/((((excludedMax-1.0)-includedMin)/2.0)))*uspeedf;
                                }
                            }
                            else
                            {
                                triangularValue = (std::abs((h-(excludedMax-1.0)))/((((excludedMax-1.0)-includedMin)/2.0)))*uspeedf;
                                T eps=maxTriangle*pow(10.0,-12.0);
                                if (triangularValue>(0.0-eps) && triangularValue<(0.0+eps))
                                {
                                    triangularValue = (std::abs((1.0))/((((excludedMax-1.0)-includedMin)/2.0)))*uspeedf;
                                }
                            }
                            lattice.get(0,h)[1]=rhoBar*Whermite[1]-triangularValue*0.70;
                            lattice.get(0,h)[2]=rhoBar*Whermite[2]-triangularValue;
                            lattice.get(0,h)[3]=rhoBar*Whermite[3]-triangularValue*0.70;
                            lattice.get(0,h)[4]=rhoBar*Whermite[4];
                            lattice.get(0,h)[5]=rhoBar*Whermite[5]+triangularValue*0.70;
                            lattice.get(0,h)[6]=rhoBar*Whermite[6]+triangularValue;
                            lattice.get(0,h)[7]=rhoBar*Whermite[7]+triangularValue*0.70;
                            lattice.get(0,h)[8]=rhoBar*Whermite[8];
                        }
                        else
                        {
                            lattice.get(0,h)[1]=rhoBar*Whermite[1]-uspeedf*0.70;
                            lattice.get(0,h)[2]=rhoBar*Whermite[2]-uspeedf;
                            lattice.get(0,h)[3]=rhoBar*Whermite[3]-uspeedf*0.70;
                            lattice.get(0,h)[4]=rhoBar*Whermite[4];
                            lattice.get(0,h)[5]=rhoBar*Whermite[5]+uspeedf*0.70;
                            lattice.get(0,h)[6]=rhoBar*Whermite[6]+uspeedf;
                            lattice.get(0,h)[7]=rhoBar*Whermite[7]+uspeedf*0.70;
                            lattice.get(0,h)[8]=rhoBar*Whermite[8];
                        }
                    }
                }
            case 1:
                if (maxNxNy)
                {
                    // velocity boundary condition in negative y direction
                    // note: mass conservative boundary conditions
                    for(int h=includedMin;h<excludedMax;h++)
                    {
                        if (triangularOption==true)
                        {
                            int maxTriangle=std::round((((excludedMax-1.0)-includedMin)/2.0) + includedMin);
                            T triangularValue;
                            if (h<=maxTriangle)
                            {
                                triangularValue = ((h-includedMin)/((((excludedMax-1.0)-includedMin)/2.0)))*uspeedf;
                                T eps=maxTriangle*pow(10.0,-12.0);
                                if (triangularValue>(0.0-eps) && triangularValue<(0.0+eps))
                                {
                                    triangularValue = ((1.0)/((((excludedMax-1.0)-includedMin)/2.0)))*uspeedf;
                                }
                            }
                            else
                            {
                                triangularValue = (std::abs((h-(excludedMax-1.0)))/((((excludedMax-1.0)-includedMin)/2.0)))*uspeedf;
                                T eps=maxTriangle*pow(10.0,-12.0);
                                if (triangularValue>(0.0-eps) && triangularValue<(0.0+eps))
                                {
                                    triangularValue = (std::abs((1.0))/((((excludedMax-1.0)-includedMin)/2.0)))*uspeedf;
                                }
                            }
                            lattice.get(h,lattice.getNy()-1)[1]=rhoBar*Whermite[1]-triangularValue*0.70;
                            lattice.get(h,lattice.getNy()-1)[2]=rhoBar*Whermite[2];
                            lattice.get(h,lattice.getNy()-1)[3]=rhoBar*Whermite[3]+triangularValue*0.70;
                            lattice.get(h,lattice.getNy()-1)[4]=rhoBar*Whermite[4]+triangularValue;
                            lattice.get(h,lattice.getNy()-1)[5]=rhoBar*Whermite[5]+triangularValue*0.70;
                            lattice.get(h,lattice.getNy()-1)[6]=rhoBar*Whermite[6];
                            lattice.get(h,lattice.getNy()-1)[7]=rhoBar*Whermite[7]-triangularValue*0.70;
                            lattice.get(h,lattice.getNy()-1)[8]=rhoBar*Whermite[8]-triangularValue;
                        }
                        else
                        {
                            lattice.get(h,lattice.getNy()-1)[1]=rhoBar*Whermite[1]-uspeedf*0.70;
                            lattice.get(h,lattice.getNy()-1)[2]=rhoBar*Whermite[2];
                            lattice.get(h,lattice.getNy()-1)[3]=rhoBar*Whermite[3]+uspeedf*0.70;
                            lattice.get(h,lattice.getNy()-1)[4]=rhoBar*Whermite[4]+uspeedf;
                            lattice.get(h,lattice.getNy()-1)[5]=rhoBar*Whermite[5]+uspeedf*0.70;
                            lattice.get(h,lattice.getNy()-1)[6]=rhoBar*Whermite[6];
                            lattice.get(h,lattice.getNy()-1)[7]=rhoBar*Whermite[7]-uspeedf*0.70;
                            lattice.get(h,lattice.getNy()-1)[8]=rhoBar*Whermite[8]-uspeedf;
                        }
                    }
                }
                else
                {
                    // velocity boundary condition in positive y direction
                    // note: mass conservative boundary conditions
                    for(int h=includedMin;h<excludedMax;h++)
                    {
                        if (triangularOption==true)
                        {
                            int maxTriangle=std::round((((excludedMax-1.0)-includedMin)/2.0) + includedMin);
                            T triangularValue;
                            if (h<=maxTriangle)
                            {
                                triangularValue = ((h-includedMin)/((((excludedMax-1.0)-includedMin)/2.0)))*uspeedf;
                                T eps=maxTriangle*pow(10.0,-12.0);
                                if (triangularValue>(0.0-eps) && triangularValue<(0.0+eps))
                                {
                                    triangularValue = ((1.0)/((((excludedMax-1.0)-includedMin)/2.0)))*uspeedf;
                                }
                            }
                            else
                            {
                                triangularValue = (std::abs((h-(excludedMax-1.0)))/((((excludedMax-1.0)-includedMin)/2.0)))*uspeedf;
                                T eps=maxTriangle*pow(10.0,-12.0);
                                if (triangularValue>(0.0-eps) && triangularValue<(0.0+eps))
                                {
                                    triangularValue = (std::abs((1.0))/((((excludedMax-1.0)-includedMin)/2.0)))*uspeedf;
                                }
                            }
                            lattice.get(h,0)[1]=rhoBar*Whermite[1]+triangularValue*0.70;
                            lattice.get(h,0)[2]=rhoBar*Whermite[2];
                            lattice.get(h,0)[3]=rhoBar*Whermite[3]-triangularValue*0.70;
                            lattice.get(h,0)[4]=rhoBar*Whermite[4]-triangularValue;
                            lattice.get(h,0)[5]=rhoBar*Whermite[5]-triangularValue*0.70;
                            lattice.get(h,0)[6]=rhoBar*Whermite[6];
                            lattice.get(h,0)[7]=rhoBar*Whermite[7]+triangularValue*0.70;
                            lattice.get(h,0)[8]=rhoBar*Whermite[8]+triangularValue;
                        }
                        else
                        {
                            lattice.get(h,0)[1]=rhoBar*Whermite[1]+uspeedf*0.70;
                            lattice.get(h,0)[2]=rhoBar*Whermite[2];
                            lattice.get(h,0)[3]=rhoBar*Whermite[3]-uspeedf*0.70;
                            lattice.get(h,0)[4]=rhoBar*Whermite[4]-uspeedf;
                            lattice.get(h,0)[5]=rhoBar*Whermite[5]-uspeedf*0.70;
                            lattice.get(h,0)[6]=rhoBar*Whermite[6];
                            lattice.get(h,0)[7]=rhoBar*Whermite[7]+uspeedf*0.70;
                            lattice.get(h,0)[8]=rhoBar*Whermite[8]+uspeedf;
                        }
                    }
                }
        }
    }
    ~BoundaryVelocityCondition<T,Descriptor>()
    {
        // nothing to be done if object is new (a dynamic entity) and if you won't delete other dynamic types, just calling delete(object)
        // if it is a static object, it will be destroyed after run when return 0 is called
        if (!quiet)
        {
            plb::pcout << "Destroying: ~BoundaryVelocityCondition" << std::endl;
        }
    }
private:
    T uspeedf;
    T rho;
    bool gaussHermite;
    bool quiet;
};

/// boundary force condition without cell indexation
/// only for distributed mass on all fis
// note: need to be validated: force arrays and other tests are mandatory
template<typename T, template<typename U> class Descriptor>
class BoundaryForceCondition {
public:
    BoundaryForceCondition(plb::Array<T,2> forcef_, T rho_, bool quiet_=true)
    {
        forcef=forcef_;
        rho=rho_;
        quiet=quiet_;
    }
    void ApplyForceBoundaryCondition(plb::MultiBlockLattice2D<T,Descriptor>& lattice, int direction, bool maxNxNy, int includedMin, int excludedMax) const
    {
        // dealing with force
        plb::Array<T,8> force;
        force[0] = (-forcef[0]+forcef[1])*(T)(1./12.);
        force[1] = (-forcef[0])*(T)(1./3.);
        force[2] = (-forcef[0]-forcef[1])*(T)(1./12.);
        force[3] = (-forcef[1])*(T)(1./3.);
        force[4] = (forcef[0]-forcef[1])*(T)(1./12.);
        force[5] = (forcef[0])*(T)(1./3.);
        force[6] = (forcef[0]+forcef[1])*(T)(1./12.);
        force[7] = (forcef[1])*(T)(1./3.);

        switch (direction)
        {
            case 0:
                if (maxNxNy)
                {
                    // force boundary condition in negative x direction
                    // note: mass conservative boundary conditions
                    for(int h=includedMin;h<excludedMax;h++)
                    {
                        for (int q=1; q<Descriptor<T>::q; q++)
                        {
                            lattice.get(lattice.getNx()-1,h)[q]=((rho-1.0)/9.0)+force[q-1];
                        }
                    }
                }
                else
                {
                    // force boundary condition in positive x direction
                    // note: mass conservative boundary conditions
                    for(int h=includedMin;h<excludedMax;h++)
                    {
                        for (int q=1; q<Descriptor<T>::q; q++)
                        {
                            lattice.get(0,h)[q]=((rho-1.0)/9.0)+force[q-1];
                        }
                    }
                }
            case 1:
                if (maxNxNy)
                {
                    // force boundary condition in negative y direction
                    // note: mass conservative boundary conditions
                    for(int h=includedMin;h<excludedMax;h++)
                    {
                        for (int q=1; q<Descriptor<T>::q; q++)
                        {
                            lattice.get(h,lattice.getNy()-1)[q]=((rho-1.0)/9.0)+force[q-1];
                        }
                    }
                }
                else
                {
                    // force boundary condition in positive y direction
                    // note: mass conservative boundary conditions
                    for(int h=includedMin;h<excludedMax;h++)
                    {
                        for (int q=1; q<Descriptor<T>::q; q++)
                        {
                            lattice.get(h,0)[q]=((rho-1.0)/9.0)+force[q-1];
                        }
                    }
                }
        }
    }
    ~BoundaryForceCondition<T,Descriptor>()
    {
        // nothing to be done if object is new (a dynamic entity) and if you won't delete other dynamic types, just calling delete(object)
        // if it is a static object, it will be destroyed after run when return 0 is called
        if (!quiet)
        {
            plb::pcout << "Destroying: ~BoundaryForceCondition" << std::endl;
        }
    }
private:
    plb::Array<T,2> forcef;
    T rho;
    bool quiet;
};

/// data processor: data processing functional for boundary body force condition
template<typename T, template<typename U> class Descriptor>
class BoundaryBodyForceConditionProcessor : public plb::BoxProcessingFunctional2D_L<T,Descriptor>{
public:
    BoundaryBodyForceConditionProcessor(plb::Array<T,2> physicalBodyForce_, T DeltaX_, T DeltaT_, int direction_, bool maxNxNy_, int includedMin_, int excludedMax_, int Nx_, int Ny_)
    {
        physicalBodyForce=physicalBodyForce_;
        DeltaX=DeltaX_;
        DeltaT=DeltaT_;
        direction=direction_;
        maxNxNy=maxNxNy_;
        includedMin=includedMin_;
        excludedMax=excludedMax_;
        Nx=Nx_;
        Ny=Ny_;
    }
    virtual void process(plb::Box2D domain, plb::BlockLattice2D<T,Descriptor>& lattice)
    {
        // access the position of the atomic-block inside the multi-block
        plb::Dot2D relativePosition = lattice.getLocation();

        for (plb::plint iX=domain.x0; iX<=domain.x1; ++iX)
        {
            for (plb::plint iY=domain.y0; iY<=domain.y1; ++iY)
            {
                // cell
                plb::Cell<T,Descriptor>& cell = lattice.get(iX,iY); // need a reference to modify cell attributes

                // convert local coordinates to global ones
                plb::plint globalX = iX + relativePosition.x;
                plb::plint globalY = iY + relativePosition.y;

                // deal with boundary conditions
                // note: mass conservative boundary conditions
                switch (direction)
                {
                    case 0:
                        if (maxNxNy)
                        {
                            if (globalX==(Nx-1))
                            {
                                if ((globalY>=includedMin) && (globalY<excludedMax))
                                {
                                    T latticeUnitsForcex = physicalBodyForce[0]*((DeltaT*DeltaT)/DeltaX); // in lattice units
                                    T latticeUnitsForcey = physicalBodyForce[1]*((DeltaT*DeltaT)/DeltaX); // in lattice units
                                    plb::Array<T,2> forcefield(latticeUnitsForcex, latticeUnitsForcey);
                                    cell.setExternalField(Descriptor<T>::ExternalField::forceBeginsAt, Descriptor<T>::ExternalField::sizeOfForce, &forcefield[0]);
                                }
                            }
                        }
                        else
                        {
                            if (globalX==0)
                            {
                                if ((globalY>=includedMin) && (globalY<excludedMax))
                                {
                                    T latticeUnitsForcex = physicalBodyForce[0]*((DeltaT*DeltaT)/DeltaX); // in lattice units
                                    T latticeUnitsForcey = physicalBodyForce[1]*((DeltaT*DeltaT)/DeltaX); // in lattice units
                                    plb::Array<T,2> forcefield(latticeUnitsForcex, latticeUnitsForcey);
                                    cell.setExternalField(Descriptor<T>::ExternalField::forceBeginsAt, Descriptor<T>::ExternalField::sizeOfForce, &forcefield[0]);
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
                                    T latticeUnitsForcex = physicalBodyForce[0]*((DeltaT*DeltaT)/DeltaX); // in lattice units
                                    T latticeUnitsForcey = physicalBodyForce[1]*((DeltaT*DeltaT)/DeltaX); // in lattice units
                                    plb::Array<T,2> forcefield(latticeUnitsForcex, latticeUnitsForcey);
                                    cell.setExternalField(Descriptor<T>::ExternalField::forceBeginsAt, Descriptor<T>::ExternalField::sizeOfForce, &forcefield[0]);
                                }
                            }
                        }
                        else
                        {
                            if (globalY==0)
                            {
                                if ((globalX>=includedMin) && (globalX<excludedMax))
                                {
                                    T latticeUnitsForcex = physicalBodyForce[0]*((DeltaT*DeltaT)/DeltaX); // in lattice units
                                    T latticeUnitsForcey = physicalBodyForce[1]*((DeltaT*DeltaT)/DeltaX); // in lattice units
                                    plb::Array<T,2> forcefield(latticeUnitsForcex, latticeUnitsForcey);
                                    cell.setExternalField(Descriptor<T>::ExternalField::forceBeginsAt, Descriptor<T>::ExternalField::sizeOfForce, &forcefield[0]);
                                }
                            }
                        }
                }
            }
        }
    }
    virtual void getModificationPattern(std::vector<bool>& isWritten) const {
        isWritten[0] = true;
    }
    virtual void getTypeOfModification(std::vector<plb::modif::ModifT>& modified) const {
        modified[0] = plb::modif::allVariables; // à comprendre à quoi cela sert voir 16.3.2 palabos user guide
        // indicates what kind of cell content was modified and must be updated in a multi-block structure
    }
    virtual BoundaryBodyForceConditionProcessor<T,Descriptor>* clone() const
    {
        return new BoundaryBodyForceConditionProcessor<T,Descriptor>(*this);
    }
    // le destructeur semble hérité mais quand est appelé le delete ?
private:
    plb::Array<T,2> physicalBodyForce;
    T DeltaX;
    T DeltaT;
    int direction;
    bool maxNxNy;
    int includedMin;
    int excludedMax;
    int Nx;
    int Ny;
};

}  // namespace sld

#endif  // SOLIDS_BOUNDARY_CONDITIONS_H

