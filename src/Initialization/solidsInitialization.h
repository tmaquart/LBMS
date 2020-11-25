///Created By T.MAQUART

/*
Initialization of distribution functions specialized for solids.
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

#ifndef SOLIDS_INITIALIZATION_H
#define SOLIDS_INITIALIZATION_H

#include "palabos2D.h"
#include "palabos2D.hh"

// using namespaces in an other namespace is not recommended due to potential conflicts

namespace sld {

/// initialization of distribution functions
template<typename T, template<typename U> class Descriptor>
class InitDistributionFunctionsProcessor : public plb::BoxProcessingFunctional2D_L<T,Descriptor>{
public:
    InitDistributionFunctionsProcessor(T rho_, bool gaussHermite_)
    {
        rho=rho_;
        gaussHermite=gaussHermite_;
    }
    virtual void process(plb::Box2D domain, plb::BlockLattice2D<T,Descriptor>& lattice)
    {
        for (plb::plint iX=domain.x0; iX<=domain.x1; ++iX)
        {
            for (plb::plint iY=domain.y0; iY<=domain.y1; ++iY)
            {
                // cell
                plb::Cell<T,Descriptor>& cell = lattice.get(iX,iY); // need a reference to modify cell attributes

                if (gaussHermite==true)
                {
                    // population with gauss hermite weights
                    T rhoBar=(rho-(T)1.0);
                    plb::Array<T,9> populationAtEquilibrium;
                    populationAtEquilibrium[0]=(4.0/9.0)*rhoBar;
                    populationAtEquilibrium[1]=(1.0/36.0)*rhoBar;
                    populationAtEquilibrium[2]=(1.0/9.0)*rhoBar;
                    populationAtEquilibrium[3]=(1.0/36.0)*rhoBar;
                    populationAtEquilibrium[4]=(1.0/9.0)*rhoBar;
                    populationAtEquilibrium[5]=(1.0/36.0)*rhoBar;
                    populationAtEquilibrium[6]=(1.0/9.0)*rhoBar;
                    populationAtEquilibrium[7]=(1.0/36.0)*rhoBar;
                    populationAtEquilibrium[8]=(1.0/9.0)*rhoBar;
                    cell.setPopulations(populationAtEquilibrium);
                }
                else
                {
                    // population like spring
                    T popRhoBar=((rho-(T)1.0)/(T)9.0);
                    plb::Array<T,9> populationAtEquilibrium;
                    for (int i=0; i<9; i++)
                    {
                        populationAtEquilibrium[i]=popRhoBar;
                    }
                    cell.setPopulations(populationAtEquilibrium);
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
    virtual InitDistributionFunctionsProcessor<T,Descriptor>* clone() const
    {
        return new InitDistributionFunctionsProcessor<T,Descriptor>(*this);
    }
    // le destructeur semble hérité mais quand est appelé le delete ?
private:
    T rho;
    bool gaussHermite;
};

/// simple initialization of distribution functions
/// deprecated due to other data processors and not for hermite dynamics
template<typename T, template<typename U> class Descriptor>
class InitializeLatticeVlasov {
public:
    InitializeLatticeVlasov(T rho_, bool quiet_=true)
    {
        rho=rho_;
        quiet=quiet_;
    }
    void InitializeSolidAtEquilibrium(plb::MultiBlockLattice2D<T,Descriptor>& lattice, bool f0Density) const
    {
        // be careful: it can modify imposed velocity conditions and computed density
        if (f0Density==false)
        {
            // density imposed on all fis and no velocity everywhere
            T rhof=(rho-1.0)/9.0; // due to rhobar
            for(int t=0;t<lattice.getNy();t++)
            {
                for(int h=0;h<lattice.getNx();h++)
                {
                    lattice.get(h,t)[0]=rhof;
                    for (int q=1;q<Descriptor<T>::q; q++)
                    {
                        lattice.get(h,t)[q]=rhof;
                    }
                }
            }
        }
        else{
            // density imposed on f0 and no velocity everywhere: can yield errors with negative distribution functions
            for(int t=0;t<lattice.getNy();t++)
            {
                for(int h=0;h<lattice.getNx();h++)
                {
                    lattice.get(h,t)[0]=rho-1.0; // due to rhobar
                    for (int q=1;q<Descriptor<T>::q; q++)
                    {
                        lattice.get(h,t)[q]=0.0; // set all other fi to zero, if lattice is a multiblock then iX and iY are global coordinates
                    }
                }
            }
        }
    }
    ~InitializeLatticeVlasov<T,Descriptor>()
    {
        // nothing to be done if object is new (a dynamic entity) and if you won't delete other dynamic types, just calling delete(object)
        // if it is a static object, it will be destroyed after run when return 0 is called
        if (!quiet)
        {
            plb::pcout << "Destroying: ~InitializeLatticeVlasov" << std::endl;
        }
    }
private:
    T rho;
    bool quiet;
};

}  // namespace sld

#endif  // SOLIDS_INITIALIZATION_H

