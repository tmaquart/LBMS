/////Created By T.MAQUART
//
///*
//Solids simulations.
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
//#include <solids2D.h> // solids headers
//#include <solids2D.cpp> // solids implementations
//
//#include "soliddynamics.h" // module header for this file
//
//using namespace plb;
//using namespace plb::descriptors;
//using namespace std;
//using namespace sld; // namespace for solids built on the top of palabos: use in sequential mode only for the moment
//
//typedef double T;
//#define DESCRIPTOR ForcedD2Q9Descriptor // for external force: mandatory for vlasov or solid dynamics
//
//void writeVTKResultsForSolids(plb::MultiBlockLattice2D<T,DESCRIPTOR>& lattice, T*** uStart, T*** strainStart, T*** stressStart, T*** divStressStart, T*** uBackwardStart, T*** uSpeedBackwardStart, plb::plint iter, T DeltaX, T DeltaT, T InvPhysicalFactorForSpeed)
//{
//    // time and space
//    T dx = DeltaX;
//    T dt = DeltaT;
//
//    // computing von mises stress
//    sld::ComputeVonMisesStress<T,DESCRIPTOR>* vonMisesStress = new ComputeVonMisesStress<T,DESCRIPTOR>(stressStart);
//    T*** vonMisesStressComputed=vonMisesStress->ComputeVonMisesStressFromStress(lattice, InvPhysicalFactorForSpeed);
//
//    // computing speed
//    sld::ComputeSpeed<T,DESCRIPTOR,2>* speed = new ComputeSpeed<T,DESCRIPTOR,2>(uBackwardStart, uStart);
//    T*** speedComputed=speed->ComputeSpeedFromDisplacement(lattice, dt);
//
//    // computing acceleration
//    sld::ComputeAcceleration<T,DESCRIPTOR,2>* acceleration = new ComputeAcceleration<T,DESCRIPTOR,2>(uSpeedBackwardStart, speedComputed);
//    T*** accelerationComputed=acceleration->ComputeAccelerationFromSpeed(lattice, dt);
//
//    // initialization
//    plb::MultiScalarField2D<T>* plbvonMisesStressMultiScalarField = new plb::MultiScalarField2D<T>(lattice.getNx(), lattice.getNy());
//    plb::MultiTensorField2D<T,2>* plbuMultiTensorField = new plb::MultiTensorField2D<T,2>(lattice.getNx(), lattice.getNy());
//    plb::MultiTensorField2D<T,3>* plbstrainMultiTensorField = new plb::MultiTensorField2D<T,3>(lattice.getNx(), lattice.getNy());
//    plb::MultiTensorField2D<T,3>* plbstressMultiTensorField  = new plb::MultiTensorField2D<T,3>(lattice.getNx(), lattice.getNy());
//    plb::MultiTensorField2D<T,2>* plbdivStressMultiTensorField = new plb::MultiTensorField2D<T,2>(lattice.getNx(), lattice.getNy());
//    plb::MultiTensorField2D<T,3>* plbu3DMultiTensorField  = new plb::MultiTensorField2D<T,3>(lattice.getNx(), lattice.getNy());
//    plb::MultiTensorField2D<T,2>* plbspeedMultiTensorField  = new plb::MultiTensorField2D<T,2>(lattice.getNx(), lattice.getNy());
//    plb::MultiTensorField2D<T,2>* plbaccelerationMultiTensorField  = new plb::MultiTensorField2D<T,2>(lattice.getNx(), lattice.getNy());
//
//    // apply functionals from namespace sld
//    plb::applyProcessingFunctional(new sld::ParallelizedComputeMultiScalarField2DFrom3DimensionalPointerArray<T,DESCRIPTOR>(vonMisesStressComputed, plbvonMisesStressMultiScalarField), lattice.getBoundingBox(), lattice);
//    plb::applyProcessingFunctional(new sld::ParallelizedComputeMultiTensorField2DFrom3DimensionalPointerArray<T,DESCRIPTOR,2>(uStart, plbuMultiTensorField), lattice.getBoundingBox(), lattice);
//    plb::applyProcessingFunctional(new sld::ParallelizedComputeMultiTensorField2DFrom3DimensionalPointerArray<T,DESCRIPTOR,3>(strainStart, plbstrainMultiTensorField), lattice.getBoundingBox(), lattice);
//    plb::applyProcessingFunctional(new sld::ParallelizedComputeMultiTensorField2DFrom3DimensionalPointerArray<T,DESCRIPTOR,3>(stressStart, plbstressMultiTensorField), lattice.getBoundingBox(), lattice);
//    plb::applyProcessingFunctional(new sld::ParallelizedComputeMultiTensorField2DFrom3DimensionalPointerArray<T,DESCRIPTOR,2>(divStressStart, plbdivStressMultiTensorField), lattice.getBoundingBox(), lattice);
//    plb::applyProcessingFunctional(new sld::ParallelizedComputeMultiTensorField2DFrom3DimensionalPointerArrayVTK<T,DESCRIPTOR,3>(uStart, plbu3DMultiTensorField), lattice.getBoundingBox(), lattice);
//    plb::applyProcessingFunctional(new sld::ParallelizedComputeMultiTensorField2DFrom3DimensionalPointerArray<T,DESCRIPTOR,2>(speedComputed, plbspeedMultiTensorField), lattice.getBoundingBox(), lattice);
//    plb::applyProcessingFunctional(new sld::ParallelizedComputeMultiTensorField2DFrom3DimensionalPointerArray<T,DESCRIPTOR,2>(accelerationComputed, plbaccelerationMultiTensorField), lattice.getBoundingBox(), lattice);
//
//    // classic palabos vtk output
//    plb::VtkImageOutput2D<T> vtkOut(plb::createFileName("sld", iter, 8), dx);
//    vtkOut.writeData(*plbvonMisesStressMultiScalarField, "Von Mises Stress (N/m2) (physical units)", 1.0); // (NON-LINEAR)
//    vtkOut.writeData<2,float>(*plbuMultiTensorField, "Displacement (m) (physical units)", InvPhysicalFactorForSpeed); // (LINEAR)
//    vtkOut.writeData<3,float>(*plbstrainMultiTensorField, "Strain (m/m) % (physical units)", InvPhysicalFactorForSpeed); // (LINEAR)
//    vtkOut.writeData<3,float>(*plbstressMultiTensorField , "Stress (N/m2) (physical units)", InvPhysicalFactorForSpeed); // (LINEAR)
//    vtkOut.writeData<2,float>(*plbdivStressMultiTensorField, "Stress Divergence (N/m3) (physical units)", InvPhysicalFactorForSpeed); // (LINEAR)
//    vtkOut.writeData<2,float>(*plbdivStressMultiTensorField, "Stress Divergence (lattice units without speed factor)", (dt*dt)/dx); // (LINEAR)
//    vtkOut.writeData<2,float>(*plbdivStressMultiTensorField, "Stress Divergence (lattice units)", InvPhysicalFactorForSpeed*(dt*dt)/dx); // (LINEAR)
//    vtkOut.writeData<3,float>(*plbu3DMultiTensorField, "Displacement (m) (physical units) 3D Warping", InvPhysicalFactorForSpeed); // (LINEAR)
//    vtkOut.writeData<2,float>(*plbspeedMultiTensorField, "Speed (m/s) (physical units)", InvPhysicalFactorForSpeed); // (LINEAR)
//    vtkOut.writeData<2,float>(*plbaccelerationMultiTensorField, "Acceleration (m/s2) (physical units)", InvPhysicalFactorForSpeed); // (LINEAR)
//
//    // delete plb namespace dynamic objects
//    delete(plbvonMisesStressMultiScalarField);
//    delete(plbuMultiTensorField);
//    delete(plbu3DMultiTensorField);
//    delete(plbspeedMultiTensorField);
//    delete(plbaccelerationMultiTensorField);
//    delete(plbstrainMultiTensorField);
//    delete(plbstressMultiTensorField);
//    delete(plbdivStressMultiTensorField);
//
//    // delete sld namespace dynamic objects
//    delete(vonMisesStress);
//    delete(speed);
//    delete(acceleration);
//
//    // free pointers from three dimensional arrays
//    sld::deleteThreeDimensionalArray(vonMisesStressComputed, 1, lattice.getNx());
//    sld::deleteThreeDimensionalArray(speedComputed, 2, lattice.getNx());
//    sld::deleteThreeDimensionalArray(accelerationComputed, 2, lattice.getNx());
//}
//
//int main(int argc, char* argv[]) {
//    plbInit(&argc, &argv);
//
//    global::directories().setOutputDir("./tmp/");
//
///// INPUTS ///
//
//    //////////////////////// main inputs ////////////////////////
//    const T DeltaT=0.000001; // physical step
//    const T maxT=0.02; // max study physical time in seconds
//    int nodesPerMeter=166; // number of nodes per meter
//    T lx=1.0; // physical dimensions x
//    T ly=0.2; // physical dimensions y
//
//    //////////////////////// output and checking specifications ////////////////////////
//    const plint savestatIter=100; // stats and save stepping
//    const plint statIter=100; // other stats stepping
//    const plint stabilityIter=100; // for stability checking
//
//    //////////////////////// factor ////////////////////////
//    const T physicalFactorForSpeed=0.001; // factor for speed, in case of negative distribution functions
//
//    //////////////////////// relaxation parameter ////////////////////////
//    // note: has no influence on vlasov dynamics and solid dynamics due to body force computation, but near boundaries it is significant
//    // note: omega must be non-zero to avoid arithmetics bug with palabos library: to be investigated
//    const T omega=0.001; // solid relaxation parameter
//
//    //////////////////////// finite difference scheme options ////////////////////////
//    bool use0X4CentralScheme=true; // 0(DeltaX^4) accurate finite difference scheme, it will take a little bit more time than the standard 0(DeltaX^2) scheme
//
//    //////////////////////// material properties ////////////////////////
//    T rho=7860; // density kg/m3
//    T convergenceForceFactor=0.00005; // intensity of the force
//    T moduleYoung=210.0*pow(10.0, 9.0); // N/m2
//    T poissonRatio=0.30; // without units
//
//    //////////////////////// boundary conditions values and options ////////////////////////
//    int includedMinXOrY=25; // min condition application
//    int excludedMaxXOrY=141; // max condition application
//    int direction=1; // direction of the load (must be consistent with palabos boundary conditions in .h module)
//    int maxNxNy=true; // sign of the load (must be consistent with palabos boundary conditions in .h module)
//    bool onlyEquilibrateCellsOfBoundaryConditions=false; // equilibrate only cells which have velocity boundary conditions
//    T omegaBoundaryConditions=1.0; // relaxation parameter for boundaries in addition of omega (concatenated relaxation parameters)
//    bool relaxationOfBoundaries=true; // must be true if you want to equilibrate boundaries
//    T maxIterationBoundaryCondition=1000; // after this iteration number, boundary condition is replaced by an other condition value
//
//    //////////////////////// solid equilibrium type ////////////////////////
//    bool gaussHermite=true; // true for gauss hermite distribution and related dynamic, false for equilibrated dynamic
//
///// END INPUTS ///
//
///// OTHER RELATED INPUTS AND DECLARATIONS ///
//
//    //////////////////////// computed inputs ////////////////////////
//    T DeltaX=1.0/nodesPerMeter; // DeltaX
//    plint nx = lx/DeltaX; // number of nodes x
//    plint ny = ly/DeltaX; // number of nodes y
//
//    //////////////////////// constitutive law: isotropic law with a symmetric stress tensor ////////////////////////
//    sld::SolidConstitutiveLaw<T>* SolidLaw = new sld::SolidConstitutiveLaw<T>(moduleYoung, poissonRatio, rho, DeltaX, DeltaT);
//
//    //////////////////////// displaying ////////////////////////
//    pcout << "DeltaT: " << DeltaT << " s" << endl;
//    pcout << "DeltaX: " << DeltaX << " m" << endl;
//    pcout << "nx: " << nx << " number of nodes x" << endl;
//    pcout << "ny: " << ny << " number of nodes y" << endl;
//    pcout << "lx: " << lx << " lx" << endl;
//    pcout << "ly: " << ly << " ly" << endl;
//    pcout << "Factor of gravity: " << (DeltaT*DeltaT)/DeltaX << " s2/m" << endl;
//    pcout << "Lambda: " << SolidLaw->GetIsotropicLameCoefficients()[0] << endl;
//    pcout << "Mu: " << SolidLaw->GetIsotropicLameCoefficients()[1] << endl;
//
//    //////////////////////// stability checker for solids ////////////////////////
//    sld::StabilityCheckerByForce<T,DESCRIPTOR>* StabilityCheckerForce = new sld::StabilityCheckerByForce<T,DESCRIPTOR>(true, true);
//    sld::StabilityCheckerForCFL<T,DESCRIPTOR>* StabilityCheckerCFL  = new sld::StabilityCheckerForCFL<T,DESCRIPTOR>(true, true);
//
//    //////////////////////// solid lattice boltzmann dynamic type ////////////////////////
//    // note: for external force we need a specific solid dynamic
//    // note: dynamic with hermite distribution and guo external force
//    MultiBlockLattice2D<T,DESCRIPTOR> lattice(nx, ny, new GuoExternalForceSolidHermiteDynamics<T, DESCRIPTOR>(omega, rho));
//
//    // note: for external force we need a specific solid dynamic
//    // note: dynamic with hermite distribution and simple external force
//    // MultiBlockLattice2D<T,DESCRIPTOR> lattice(nx, ny, new GravityLikeExternalForceSolidHermiteDynamics<T, DESCRIPTOR>(omega, rho));
//
//    // note: for external force we need a specific solid dynamic
//    // note: gaussHermite must be false
//    // note: dynamic with equilibrated distribution and simple external force
//    // MultiBlockLattice2D<T,DESCRIPTOR> lattice(nx, ny, new GravityLikeExternalForceSolidDynamics<T, DESCRIPTOR>(omega, rho));
//
//    // note: for external force we need a specific solid dynamic
//    // note: vlasov dynamics without collisions: omega is not used
//    // note: simple external force
//    // note: gaussHermite option sets the initial distribution
//    // MultiBlockLattice2D<T,DESCRIPTOR> lattice(nx, ny, new GravityLikeExternalForceVlasovDynamics<T, DESCRIPTOR>(omega));
//
//    // note: for external force we need a specific solid dynamic
//    // note: vlasov dynamics without collisions: omega is not used
//    // note: guo external force
//    // note: gaussHermite option sets the initial distribution
//    // MultiBlockLattice2D<T,DESCRIPTOR> lattice(nx, ny, new GuoExternalForceVlasovDynamics<T, DESCRIPTOR>(omega));
//
//    //////////////////////// set perodicity ////////////////////////
//    lattice.periodicity().toggle(0, false); // default bounceback
//    lattice.periodicity().toggle(1, false); // default bounceback
//
//    //////////////////////// palabos boundary condition ////////////////////////
//    OnLatticeBoundaryCondition2D<T,DESCRIPTOR>* boundaryCondition = createLocalBoundaryCondition2D<T,DESCRIPTOR>();
//
//    //////////////////////// setup ////////////////////////
//    Setup(lattice, *boundaryCondition, rho, includedMinXOrY, excludedMaxXOrY, gaussHermite, physicalFactorForSpeed, DeltaX, DeltaT);
//
//    //////////////////////// initialization of solids tensors ////////////////////////
//    sld::SolidsTensorsInitializer<T,DESCRIPTOR>* SolidsTensors =  new sld::SolidsTensorsInitializer<T,DESCRIPTOR>();
//    SolidsTensors->InitializePointers(lattice);
//    T*** uStart=SolidsTensors->GetuStart();
//    T*** uBackwardStart=SolidsTensors->GetuBackwardStart();
//    T*** uSpeedBackwardStart=SolidsTensors->GetuSpeedBackwardStart();
//    T*** strainStart=SolidsTensors->GetstrainStart();
//    T*** stressStart=SolidsTensors->GetstressStart();
//    T*** divStressStart=SolidsTensors->GetdivStressStart();
//
//    //////////////////////// data processors integration ////////////////////////
//    // note: data processor: with a specific processor level to be executed manually: need to be integrated to be executed at each time step
//    plb::integrateProcessingFunctional(new sld::ComputeDisplacementFromSpeedProcessor<T,DESCRIPTOR>(uStart, DeltaT, DeltaX, lattice.getNx(), lattice.getNy(), false), lattice.getBoundingBox(), lattice, -3);
//    // note: data processor: with a specific processor level to be executed manually: need to be integrated to be executed at each time step
//    plb::integrateProcessingFunctional(new sld::ComputeStrainAndStressProcessor<T,DESCRIPTOR>(uStart, strainStart, stressStart, DeltaX, SolidLaw, lattice.getNx(), lattice.getNy(), false, false, false, use0X4CentralScheme, false), lattice.getBoundingBox(), lattice, -2);
//    // note: data processor: with a specific processor level to be executed manually: need to be integrated to be executed at each time step
//    plb::integrateProcessingFunctional(new sld::ComputeStressDivergenceAndApplySpaceTimeDependentExternalForceProcessor<T,DESCRIPTOR>(stressStart, divStressStart, DeltaX, DeltaT, lattice.getNx(), lattice.getNy(), false, false, false, use0X4CentralScheme, convergenceForceFactor, false), lattice.getBoundingBox(), lattice, -1);
//    // note: data processor: with a specific processor level to be executed manually: need to be integrated to be executed at each time step
//    plb::integrateProcessingFunctional(new sld::EquilibriumOfBoundariesWithConcatenatedRelaxationParameter<T,DESCRIPTOR>(omegaBoundaryConditions, rho, gaussHermite, direction, maxNxNy, includedMinXOrY, excludedMaxXOrY, lattice.getNx(), lattice.getNy(), onlyEquilibrateCellsOfBoundaryConditions, false), lattice.getBoundingBox(), lattice, -4);
//
//    //////////////////////// data processors execution at iT=0 ////////////////////////
//    // note: compute initial displacement by time integration iT=0, one step earlier for displacement to integrate initial velocity condition
//    lattice.executeInternalProcessors(-3); // execute the data processor at a specific level
//
//    // note: compute strain and stress
//    lattice.executeInternalProcessors(-2); // execute the data processor at a specific level
//
//    //////////////////////// initial CFL checker: need a velocity boundary condition ////////////////////////
//    StabilityCheckerCFL->CheckCFLAndRaiseException(*computeVelocityNorm(lattice), lattice, DeltaX, DeltaT);
//
///// END OTHER RELATED INPUTS AND DECLARATIONS ///
//
///// SOLID LATTICE BOLTZMANN LOOP ///
//
//    // loop
//    for (plint iT=0; iT*DeltaT<maxT; ++iT) {
//
//        // stepping
//        if (iT%savestatIter==0) {
//            pcout << endl;
//            pcout << "step " << iT
//                  << "; t=" << iT*DeltaT << endl;
//        }
//
//        // vtk saving for solids
//        if (iT%savestatIter==0) {
//            pcout << "Saving solids VTK file ..." << endl;
//            writeVTKResultsForSolids(lattice, uStart, strainStart, stressStart, divStressStart, uBackwardStart, uSpeedBackwardStart, iT, DeltaX, DeltaT, 1.0/physicalFactorForSpeed);
//        }
//
//        // for acceleration and speed rendering purposes
//        uSpeedBackwardStart=computeBackwardSpeed(lattice, uStart, uBackwardStart, DeltaT);
//        SolidsTensors->CopyDisplacementSpeedAccelerationTensors(uBackwardStart, uStart, lattice);
//
//        // apply space time dependent force: at this time related functions of cells are modified (f+force is then considered for velocity): it may change velocity at this line
//        lattice.executeInternalProcessors(-1); // execute the data processor at a specific level
//
//        // relaxation of boundaries
//        if (relaxationOfBoundaries==true)
//        {
//            // relaxation of boundaries
//            lattice.executeInternalProcessors(-4);
//        }
//
//        // collide: external forces are already integrated with externalfield
//        lattice.collide(); // apply external forces into distribution functions
//
//        // stream
//        lattice.stream();
//
//        // palabos boundary conditions
//        if (iT==maxIterationBoundaryCondition)
//        {
//            UnSetBoundaryVelocity(lattice, *boundaryCondition, includedMinXOrY, excludedMaxXOrY);
//        }
//
//        // execute customized internal processor, need boundary velocity condition or force at initial condition for level -3
//        lattice.executeInternalProcessors(-3); // execute the data processor at a specific level
//        lattice.executeInternalProcessors(-2); // execute the data processor at a specific level
//
//        // stability checker
//        if (iT%stabilityIter==0 && iT>0) {
//            StabilityCheckerForce->AnalyseStabilityAndRaiseExceptionDuringIterations(lattice);
//            StabilityCheckerCFL->CheckCFLAndRaiseException(*computeVelocityNorm(lattice), lattice, DeltaX, DeltaT);
//        }
//
//        // initial stability checker
//        if (iT==0)
//        {
//            StabilityCheckerForce->AnalyseInitialStabilityAndRaiseException(divStressStart, rho, lattice.getNx(), lattice.getNy(), DeltaX, DeltaT, gaussHermite);
//        }
//
//        // useful stats
//        if (iT%statIter==0) {
//            pcout << "; av energy ="
//                  << setprecision(17) << getStoredAverageEnergy<T>(lattice)
//                  << "; av rho ="
//                  << setprecision(17) << getStoredAverageDensity<T>(lattice) << endl;
//        }
//    }
//
///// END SOLID LATTICE BOLTZMANN LOOP ///
//
///// FREE MEMORY ///
//
//    // free palabos boundary condition
//    delete boundaryCondition;
//
//    // free pointers from three dimensional arrays
//    sld::deleteThreeDimensionalArray(uStart, 2, lattice.getNx());
//    sld::deleteThreeDimensionalArray(uBackwardStart, 2, lattice.getNx());
//    sld::deleteThreeDimensionalArray(uSpeedBackwardStart, 2, lattice.getNx());
//    sld::deleteThreeDimensionalArray(strainStart, 3, lattice.getNx());
//    sld::deleteThreeDimensionalArray(stressStart, 3, lattice.getNx());
//    sld::deleteThreeDimensionalArray(divStressStart, 2, lattice.getNx());
//
//    // free dynamic types from namespace sld
//    delete(StabilityCheckerForce); // need to delete explicitly because of dynamic type
//    delete(StabilityCheckerCFL); // need to delete explicitly because of dynamic type
//    delete(SolidsTensors); // need to delete explicitly because of dynamic type
//    delete(SolidLaw); // need to delete explicitly because of dynamic type
//    return 0; // calling destructors for static objects
//
///// END FREE MEMORY ///
//
//}
