/////Created By T.MAQUART
//
///*
//vlasov simulations.
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
//#include "palabos2D.h"
//#include "palabos2D.hh"
//#include <vector>
//#include <cmath>
//#include <iostream>
//#include <fstream>
//#include <iomanip>
//
//using namespace plb;
//using namespace plb::descriptors;
//using namespace std;
//
//typedef double T;
//#define DESCRIPTOR descriptors::ForcedD2Q9Descriptor // with external force or not
//
//T ConstantVelocityX(T latticeU) {
//    return latticeU;
//}
//
//T ConstantDensity(T rhoLattice) {
//    return rhoLattice;
//}
//
//template<typename T>
//class ConstantVelocityAndDensityX {
//public:
//    ConstantVelocityAndDensityX(T LatticeU_, T rhoLattice_)
//    {
//     LatticeU=LatticeU_;
//     rhoLattice=rhoLattice_;
//    }
//    void operator()(plint iX, plint iY, T& rho, Array<T,2>& u) const {
//        rho = ConstantDensity(rhoLattice);
//        u[0] = ConstantVelocityX(LatticeU);
//        u[1] = T();
//    }
//private:
//    T LatticeU;
//    T rhoLattice;
//};
//
//template<typename T>
//class CylinderShapeDomain2D : public plb::DomainFunctional2D {
//public:
//    CylinderShapeDomain2D(plb::plint cx_, plb::plint cy_, plb::plint radius)
//        : cx(cx_),
//          cy(cy_),
//          radiusSqr(plb::util::sqr(radius))
//    { }
//    virtual bool operator() (plb::plint iX, plb::plint iY) const {
//        return plb::util::sqr(iX-cx) + plb::util::sqr(iY-cy) <= radiusSqr;
//    }
//    virtual CylinderShapeDomain2D<T>* clone() const {
//        return new CylinderShapeDomain2D<T>(*this);
//    }
//private:
//    plb::plint cx;
//    plb::plint cy;
//    plb::plint radiusSqr;
//};
//
//void Setup( MultiBlockLattice2D<T,DESCRIPTOR>& lattice, int nx, int ny, T latticeU, T rhoLattice, int resolution, T diameter, T DeltaT, T DeltaX)
//{
//    // initial equilibrium
//    initializeAtEquilibrium(lattice, lattice.getBoundingBox(), ConstantVelocityAndDensityX<T>(latticeU, rhoLattice));
//
//    // geometry
//    plint cx=nx/2;
//    plint cy=ny/2;
//    plint radius=int((diameter/2)*resolution);
//
//    // bounceback dynamics
//    defineDynamics(lattice, lattice.getBoundingBox(), new CylinderShapeDomain2D<T>(cx,cy,radius), new plb::BounceBack<T,DESCRIPTOR>);
//
//    // external gravity force
//    T gravity=9.81*(DeltaT*DeltaT/DeltaX);
//
//    // gravity
//    setExternalVector(lattice, lattice.getBoundingBox(), DESCRIPTOR<T>::ExternalField::forceBeginsAt, Array<T,2>(0.0,-rhoLattice*gravity));
//
//    lattice.initialize();
//}
//
//void writeGif(MultiBlockLattice2D<T,DESCRIPTOR>& lattice, plint iter)
//{
//    ImageWriter<T> imageWriter("leeloo");
//    imageWriter.writeScaledGif(createFileName("u", iter, 6),
//                               *computeVelocityNorm(lattice) );
//}
//
//void writeVTK(MultiBlockLattice2D<T,DESCRIPTOR>& lattice, T DeltaX, T DeltaT, plint iter)
//{
//    T dx = DeltaX;
//    T dt = DeltaT;
//    VtkImageOutput2D<T> vtkOut(createFileName("vtk", iter, 6), dx);
//    vtkOut.writeData<float>(*computeDensity(lattice), "Density", 1);
//    vtkOut.writeData<float>(*computeVelocityNorm(lattice), "velocityNorm", dx/dt);
//    vtkOut.writeData<2,float>(*computeVelocity(lattice), "velocity", dx/dt);
//}
//
//int main(int argc, char* argv[]) {
//    plbInit(&argc, &argv);
//    global::directories().setOutputDir("./tmp/");
//
//    T PhysicalU=0.2;
//    T cylinderDiameter=0.2;
//    T rho=1000; // density
//    const T DeltaT=0.0001; // physical step
//    T lx=1;
//    T ly=1;
//    const T resolution=500;
//    plint nx=lx*resolution;
//    plint ny=ly*resolution;
//    T DeltaX=lx/nx;
//    // note: omega must be non-zero to avoid arithmetics bug with palabos library: to be investigated
//    T omega=0.0001; // omega has no influence on vlasov dynamics
//    const T maxT=0.05; // physical time
//    const plint savestatIter=10;
//    const plint statIter=10;
//
//    // displaying
//    pcout << "resolution: " << resolution << endl;
//    pcout << "DeltaT: " << DeltaT << endl;
//    pcout << "DeltaX: " << DeltaX << endl;
//    pcout << "omega: " << omega << endl;
//
//    // without external force term
//    // MultiBlockLattice2D<T, DESCRIPTOR> lattice(nx, ny, new plb::VlasovDynamics<T,DESCRIPTOR>(omega));
//
//    // with external force term: if external gravity force is uncommented
//    // note: simple force term
//    // MultiBlockLattice2D<T, DESCRIPTOR> lattice(nx, ny, new plb::GravityLikeExternalForceVlasovDynamics<T,DESCRIPTOR>(omega));
//
//    // with external force term: if external gravity force is uncommented
//    // note: guo force term
//    MultiBlockLattice2D<T, DESCRIPTOR> lattice(nx, ny, new plb::GuoExternalForceVlasovDynamics<T,DESCRIPTOR>(omega));
//
//    // periodicity
//    lattice.periodicity().toggleAll(true);
//
//    // setup
//    Setup(lattice, nx, ny, PhysicalU*(DeltaT/DeltaX), rho, resolution, cylinderDiameter, DeltaT, DeltaX);
//
//    // main loop
//    for (plint iT=0; iT*DeltaT<maxT; ++iT) {
//
//        if (iT%savestatIter==0) {
//            pcout << "Saving Gif ..." << endl;
//            writeGif(lattice, iT);
//        }
//
//        if (iT%savestatIter==0 && iT>0) {
//            pcout << "Saving VTK file ..." << endl;
//            writeVTK(lattice, DeltaX, DeltaT, iT);
//        }
//
//        if (iT%savestatIter==0) {
//            pcout << "step " << iT
//                  << "; t=" << iT*DeltaT;
//        }
//
//        lattice.collideAndStream();
//
//        if (iT%statIter==0) {
//            pcout << "; av energy ="
//                  << setprecision(10) << getStoredAverageEnergy<T>(lattice)
//                  << "; av rho ="
//                  << getStoredAverageDensity<T>(lattice) << endl;
//        }
//    }
//}
