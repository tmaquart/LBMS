/////Created By T.MAQUART
//
///*
//Module as header.
//*/
//
///* This file is part of the SLD library.
// *
// * This library is mainly dedicated to research. This library uses PALABOS
// * library which has her own license. No warranty is given in terms of computational
// * efficiency and reliability of produced results. Moreover, no guarantees are given
// * concerning programmed interfaces and classes especially on the respect of C++ good practices.
// *
// * Only sequential mode is supported. SLD library has to be improved for parallelized runs
// * according to the PALABOS programming strategy.
// *
// * E-mail contact: tristan.maquart@emse.fr, tristan.maquart@hotmail.fr
// * E-mail contact: ronoel@ethz.ch
// * E-mail contact: navarro@emse.fr
// *
// * Copyright (C) <2020> <Tristan MAQUART, Romain NOEL, Laurent NAVARRO>
// *
// * This program is free software: you can redistribute it and/or modify
// * it under the terms of the GNU Affero General Public License as published by
// * the Free Software Foundation, either version 3 of the License, or
// * (at your option) any later version.
// *
// * This program is distributed in the hope that it will be useful,
// * but WITHOUT ANY WARRANTY; without even the implied warranty of
// * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// * GNU Affero General Public License for more details.
// *
// * You should have received a copy of the GNU Affero General Public License
// * along with this program. If not, see <http://www.gnu.org/licenses/>.
//*/
//
//// T.MAQUART Ecole nationale supérieure des mines de Saint-Etienne : tristan.maquart@emse.fr, tristan.maquart@hotmail.fr
//// R.NOEL Ecole polytechnique fédérale de Zurich : ronoel@ethz.ch
//// L.NAVARRO Ecole nationale supérieure des mines de Saint-Etienne : navarro@emse.fr
//
//#ifndef STATIC_EQUILIBRIUM_ISO_H
//#define STATIC_EQUILIBRIUM_ISO_H
//
//#include "palabos2D.h"
//#include "palabos2D.hh"
//
//// using namespaces in an other namespace is not recommended due to potential conflicts
//
//namespace sld {
//
//// boundary velocity to zero
//template<typename T>
//class UnSetVelocityByZero {
//public:
//    UnSetVelocityByZero(int nx_, int ny_)
//    {
//        nx=nx_;
//        ny=ny_;
//    }
//    void operator()(plb::plint iX, plb::plint iY, plb::Array<T,2>& u) const {
//        u[0] = 0.0;
//        u[1] = 0.0;
//    }
//private:
//    int nx;
//    int ny;
//};
//
//// boundary velocity
//template<typename T>
//class CustomizedBoundaryVelocityConditions {
//public:
//    CustomizedBoundaryVelocityConditions(int nx_, int ny_, T physicalSpeed_, T physicalFactorForSpeed_, T DeltaX_, T DeltaT_, int includedMinXOrY_, int excludedMaxXOrY_)
//    {
//        nx=nx_;
//        ny=ny_;
//        physicalSpeed=physicalSpeed_;
//        physicalFactorForSpeed=physicalFactorForSpeed_;
//        DeltaX=DeltaX_;
//        DeltaT=DeltaT_;
//        includedMinXOrY=includedMinXOrY_;
//        excludedMaxXOrY=excludedMaxXOrY_;
//    }
//    void operator()(plb::plint iX, plb::plint iY, plb::Array<T,2>& u) const {
//
//        // y max loading
//        if (iY==(ny-1))
//        {
//            if (iX<includedMinXOrY || iX>=excludedMaxXOrY)
//            {
//                u[0] = 0.0;
//                u[1] = 0.0;
//            }
//            else
//            {
//                // rectangle load
//                u[0] = 0.0;
//                u[1] = physicalSpeed*(DeltaT/DeltaX)*physicalFactorForSpeed; // in lattice units
//
////                // triangle load
////                u[0] = 0.0;
////                T stepV =(physicalSpeed*(DeltaT/DeltaX)*physicalFactorForSpeed)/(excludedMaxXOrY-includedMinXOrY);
////                u[1] = stepV*(iX-(includedMinXOrY-1)); // in lattice units
//            }
//        }
//        else
//        {
//            u[0] = 0.0;
//            u[1] = 0.0;
//        }
//    }
//private:
//    int nx;
//    int ny;
//    T physicalSpeed;
//    T physicalFactorForSpeed;
//    T DeltaX;
//    T DeltaT;
//    int includedMinXOrY;
//    int excludedMaxXOrY;
//};
//
//template <typename T, template <typename U> class Descriptor>
//void UnSetBoundaryVelocity(plb::MultiBlockLattice2D<T,Descriptor>& lattice, plb::OnLatticeBoundaryCondition2D<T,Descriptor>& boundaryCondition, int includedMinXOrY, int excludedMaxXOrY)
//{
//    // change value of the velocity
//    boundaryCondition.setVelocityConditionOnBlockBoundaries(lattice, plb::Box2D(1, lattice.getNx()-2, lattice.getNy()-1, lattice.getNy()-1));
//    boundaryCondition.setVelocityConditionOnBlockBoundaries(lattice, plb::Box2D(1, lattice.getNx()-2, 0, 0));
//    boundaryCondition.setVelocityConditionOnBlockBoundaries(lattice, plb::Box2D(0, 0, 0, lattice.getNy()-1));
//    boundaryCondition.setVelocityConditionOnBlockBoundaries(lattice, plb::Box2D(lattice.getNx()-1, lattice.getNx()-1, 0, lattice.getNy()-1));
//    plb::setBoundaryVelocity(lattice, lattice.getBoundingBox(), UnSetVelocityByZero<T>(lattice.getNx(), lattice.getNy()));
//}
//
//template <typename T, template <typename U> class Descriptor>
//void Setup(plb::MultiBlockLattice2D<T,Descriptor>& lattice, plb::OnLatticeBoundaryCondition2D<T,Descriptor>& boundaryCondition, T rho, int includedMinXOrY, int excludedMaxXOrY, bool gaussHermite, T physicalFactorForSpeed, T DeltaX, T DeltaT)
//{
//    // solid equilibrium with a data processor
//    plb::applyProcessingFunctional(new sld::InitDistributionFunctionsProcessor<T,Descriptor>(rho, gaussHermite), lattice.getBoundingBox(), lattice);
//
//    // palabos boundary conditions
//    T physicalSpeed=-480.0; // physical speed
//    boundaryCondition.setVelocityConditionOnBlockBoundaries(lattice, plb::Box2D(1, lattice.getNx()-2, lattice.getNy()-1, lattice.getNy()-1));
//    boundaryCondition.setVelocityConditionOnBlockBoundaries(lattice, plb::Box2D(1, lattice.getNx()-2, 0, 0));
//    boundaryCondition.setVelocityConditionOnBlockBoundaries(lattice, plb::Box2D(0, 0, 0, lattice.getNy()-1));
//    boundaryCondition.setVelocityConditionOnBlockBoundaries(lattice, plb::Box2D(lattice.getNx()-1, lattice.getNx()-1, 0, lattice.getNy()-1));
//    plb::setBoundaryVelocity(lattice, lattice.getBoundingBox(), CustomizedBoundaryVelocityConditions<T>(lattice.getNx(), lattice.getNy(), physicalSpeed, physicalFactorForSpeed, DeltaX, DeltaT, includedMinXOrY, excludedMaxXOrY));
//
//    lattice.initialize(); // first internal communications
//}
//
//template <typename T>
//void deleteThreeDimensionalArray(T*** arrayToDelete, int aDim, int bDim)
//{
//    for(int i=0; i<aDim; i++)
//    {
//        for(int j=0; j<bDim; j++)
//        {
//            delete[] arrayToDelete[i][j];
//        }
//        delete[] arrayToDelete[i];
//    }
//    delete[] arrayToDelete;
//}
//
//template <typename T, template <typename U> class Descriptor>
//T*** computeBackwardSpeed(plb::MultiBlockLattice2D<T,Descriptor>& lattice, T*** uStart, T*** backwarduStart, T DeltaT)
//{
//    // time
//    T dt = DeltaT;
//
//    // computing speed
//    sld::ComputeSpeed<T,Descriptor,2>* speed = new ComputeSpeed<T,Descriptor,2>(backwarduStart, uStart);
//    T*** speedComputed=speed->ComputeSpeedFromDisplacement(lattice, dt);
//
//    // delete sld namespace dynamic objects
//    delete(speed);
//
//    return speedComputed;
//}
//
//}  // namespace sld
//
//#endif  // STATIC_EQUILIBRIUM_ISO_H
//
