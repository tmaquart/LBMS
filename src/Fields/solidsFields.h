///Created By T.MAQUART

/*
Solids fields for vtk rendering or other purposes.
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

#ifndef SOLIDS_FIELDS_H
#define SOLIDS_FIELDS_H

#include "palabos2D.h"
#include "palabos2D.hh"

// using namespaces in an other namespace is not recommended due to potential conflicts

namespace sld {

/// compute 2D multi scalar field from a 3 dimensional pointer array
template<typename T, template<typename U> class Descriptor>
class ComputeMultiScalarField2DFrom3DimensionalPointerArray {
public:
    /// matching dimensions is your own problematic in order to avoid segmentation fault errors
    ComputeMultiScalarField2DFrom3DimensionalPointerArray(T*** multiscalarfield_, int component_=0, bool quiet_=true)
    {
        multiscalarfield=multiscalarfield_;
        component=component_;
        quiet=quiet_;
    }
    /// multi scalar field from constructor must match with size of the multi block lattice
    plb::MultiScalarField2D<T>* GetMultiScalarField(plb::MultiBlockLattice2D<T,Descriptor>& lattice) const
    {
        // creating palabos object
        plb::MultiScalarField2D<T>* plbmultiscalarfield = new plb::MultiScalarField2D<T>(lattice.getNx(), lattice.getNy());

        // filling scalar field
        for (int i=0; i<lattice.getNx(); i++)
        {
            for (int j=0; j<lattice.getNy(); j++)
            {
                plbmultiscalarfield->get(i, j)=multiscalarfield[component][i][j];
            }
        }
        return plbmultiscalarfield;
    }
    ~ComputeMultiScalarField2DFrom3DimensionalPointerArray()
    {
        // nothing to be done if object is new (a dynamic entity) and if you won't delete other dynamic types, just calling delete(object)
        // if it is a static object, it will be destroyed after run when return 0 is called
        if (!quiet)
        {
           plb::pcout << "Destroying: ~ComputeMultiScalarField2DFrom3DimensionalPointerArray" << std::endl;
        }
    }
private:
    T*** multiscalarfield;
    int component;
    bool quiet;
};

/// compute 2D multi tensor field from a 3 dimensional pointer array
template<typename T, template<typename U> class Descriptor, int nDim>
class ComputeMultiTensorField2DFrom3DimensionalPointerArray {
public:
    /// matching dimensions is your own problematic in order to avoid segmentation fault errors
    ComputeMultiTensorField2DFrom3DimensionalPointerArray(T*** multitensorfield_, bool quiet_=true)
    {
        multitensorfield=multitensorfield_;
        quiet=quiet_;
    }
    /// multi tensor field from constructor must match with size of the multi block lattice
    plb::MultiTensorField2D<T, nDim>* GetMultiTensorField(plb::MultiBlockLattice2D<T,Descriptor>& lattice) const
    {
        // creating palabos object
        plb::MultiTensorField2D<T, nDim>* plbmultitensorfield = new plb::MultiTensorField2D<T, nDim>(lattice.getNx(), lattice.getNy());

        // filling tensor field
        for (int i=0; i<lattice.getNx(); i++)
        {
            for (int j=0; j<lattice.getNy(); j++)
            {
                for (int k=0; k<nDim; k++)
                {
                    plbmultitensorfield->get(i, j)[k]=multitensorfield[k][i][j];
                }
            }
        }
        return plbmultitensorfield;
    }
    ~ComputeMultiTensorField2DFrom3DimensionalPointerArray()
    {
        // nothing to be done if object is new (a dynamic entity) and if you won't delete other dynamic types, just calling delete(object)
        // if it is a static object, it will be destroyed after run when return 0 is called
        if (!quiet)
        {
           plb::pcout << "Destroying: ~ComputeMultiTensorField2DFrom3DimensionalPointerArray" << std::endl;
        }
    }
private:
    T*** multitensorfield;
    bool quiet;
};

/// processing functional: parallelized version for building multi scalar field from 3 dimensional pointer
template<typename T, template<typename U> class Descriptor>
class ParallelizedComputeMultiScalarField2DFrom3DimensionalPointerArray : public plb::BoxProcessingFunctional2D_L<T,Descriptor>{
public:
    ParallelizedComputeMultiScalarField2DFrom3DimensionalPointerArray(T*** multiscalarfield_, plb::MultiScalarField2D<T>* plbmultiscalarfield_, int component_=0)
    {
        multiscalarfield=multiscalarfield_;
        plbmultiscalarfield=plbmultiscalarfield_;
        component=component_;
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

                // filling scalar field
                plbmultiscalarfield->get(globalX, globalY)=multiscalarfield[component][globalX][globalY];
            }
        }
    }
    virtual void getTypeOfModification(std::vector<plb::modif::ModifT>& modified) const {
        modified[0] = plb::modif::nothing; // à comprendre à quoi cela sert voir 16.3.2 palabos user guide
        // indicates what kind of cell content was modified and must be updated in a multi-block structure
    }
    virtual ParallelizedComputeMultiScalarField2DFrom3DimensionalPointerArray<T,Descriptor>* clone() const
    {
        return new ParallelizedComputeMultiScalarField2DFrom3DimensionalPointerArray<T,Descriptor>(*this);
    }
    // le destructeur semble hérité mais quand est appelé le delete ?
private:
    T*** multiscalarfield;
    plb::MultiScalarField2D<T>* plbmultiscalarfield;
    int component;
};

/// processing functional: parallelized version for building multi tensor field from 3 dimensional pointer
template<typename T, template<typename U> class Descriptor, int nDim>
class ParallelizedComputeMultiTensorField2DFrom3DimensionalPointerArray : public plb::BoxProcessingFunctional2D_L<T,Descriptor>{
public:
    ParallelizedComputeMultiTensorField2DFrom3DimensionalPointerArray(T*** multitensorfield_, plb::MultiTensorField2D<T,nDim>* plbmultitensorfield_)
    {
        multitensorfield=multitensorfield_;
        plbmultitensorfield=plbmultitensorfield_;
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

                // filling tensor field
                for (int k=0; k<nDim; k++)
                {
                    plbmultitensorfield->get(globalX, globalY)[k]=multitensorfield[k][globalX][globalY];
                }
            }
        }
    }
    virtual void getTypeOfModification(std::vector<plb::modif::ModifT>& modified) const {
        modified[0] = plb::modif::nothing; // à comprendre à quoi cela sert voir 16.3.2 palabos user guide
        // indicates what kind of cell content was modified and must be updated in a multi-block structure
    }
    virtual ParallelizedComputeMultiTensorField2DFrom3DimensionalPointerArray<T,Descriptor,nDim>* clone() const
    {
        return new ParallelizedComputeMultiTensorField2DFrom3DimensionalPointerArray<T,Descriptor,nDim>(*this);
    }
    // le destructeur semble hérité mais quand est appelé le delete ?
private:
    T*** multitensorfield;
    plb::MultiTensorField2D<T,nDim>* plbmultitensorfield;
};

/// processing functional: parallelized version for building multi tensor field from 3 dimensional pointer for VTK rendering
template<typename T, template<typename U> class Descriptor, int nDim>
class ParallelizedComputeMultiTensorField2DFrom3DimensionalPointerArrayVTK : public plb::BoxProcessingFunctional2D_L<T,Descriptor>{
public:
    ParallelizedComputeMultiTensorField2DFrom3DimensionalPointerArrayVTK(T*** multitensorfield_, plb::MultiTensorField2D<T,nDim>* plbmultitensorfield_)
    {
        multitensorfield=multitensorfield_;
        plbmultitensorfield=plbmultitensorfield_;
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

                // filling tensor field
                for (int k=0; k<(nDim-1); k++)
                {
                    plbmultitensorfield->get(globalX, globalY)[k]=multitensorfield[k][globalX][globalY];
                }

                // filling third or last fake dimension
                plbmultitensorfield->get(globalX, globalY)[nDim-1]=0.0;
            }
        }
    }
    virtual void getTypeOfModification(std::vector<plb::modif::ModifT>& modified) const {
        modified[0] = plb::modif::nothing; // à comprendre à quoi cela sert voir 16.3.2 palabos user guide
        // indicates what kind of cell content was modified and must be updated in a multi-block structure
    }
    virtual ParallelizedComputeMultiTensorField2DFrom3DimensionalPointerArrayVTK<T,Descriptor,nDim>* clone() const
    {
        return new ParallelizedComputeMultiTensorField2DFrom3DimensionalPointerArrayVTK<T,Descriptor,nDim>(*this);
    }
    // le destructeur semble hérité mais quand est appelé le delete ?
private:
    T*** multitensorfield;
    plb::MultiTensorField2D<T,nDim>* plbmultitensorfield;
};

}  // namespace sld

#endif  // SOLIDS_FIELDS_H

