/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2014 Stanford University and the Authors.           *
 * Authors: Peter Eastman                                                     *
 * Contributors:                                                              *
 *                                                                            *
 * Permission is hereby granted, free of charge, to any person obtaining a    *
 * copy of this software and associated documentation files (the "Software"), *
 * to deal in the Software without restriction, including without limitation  *
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,   *
 * and/or sell copies of the Software, and to permit persons to whom the      *
 * Software is furnished to do so, subject to the following conditions:       *
 *                                                                            *
 * The above copyright notice and this permission notice shall be included in *
 * all copies or substantial portions of the Software.                        *
 *                                                                            *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR *
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,   *
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL    *
 * THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,    *
 * DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR      *
 * OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE  *
 * USE OR OTHER DEALINGS IN THE SOFTWARE.                                     *
 * -------------------------------------------------------------------------- */

#include "CudaFeedInKernels.h"
#include "CudaFeedInKernelSources.h"
#include "openmm/internal/ContextImpl.h"
#include "openmm/cuda/CudaForceInfo.h"
#include <math.h>

using namespace FeedInPlugin;
using namespace OpenMM;
using namespace std;

class CudaFeedInForceInfo : public CudaForceInfo {
public:
    CudaFeedInForceInfo(const FeedInForce& force) : force(force) {
    }

private:
    const FeedInForce& force;
};

CudaCalcFeedInForceKernel::CudaCalcFeedInForceKernel(std::string name, const OpenMM::Platform& platform, OpenMM::CudaContext& cu, const OpenMM::System& system) :
        CalcFeedInForceKernel(name, platform), hasInitializedKernel(false), cu(cu), system(system) {}

void CudaCalcFeedInForceKernel::copyForceToContext(ContextImpl& context, const std::vector<std::vector<double> > &forces_vec) {
    cu.setAsCurrent();
    if (forces_vec.size() != context.getSystem().getNumParticles())
        throw OpenMMException("copyForceToContext: number of particles != number of forces provided.");
    vector<double> tmp_vec;
    for (int i=0; i<forces_vec.size(); ++i){
        if (forces_vec[i].size() != 3)
            throw OpenMMException("execute: number of force components is wrong (paritcle index=" + to_string(i) + "). Expecting 3, got " + to_string(forces_vec[i].size()) + ".");
        for (int j=0; j<3; ++j){
            tmp_vec.push_back(forces_vec[i][j]);
        }
    }
    forces_in->upload(tmp_vec);
}

CudaCalcFeedInForceKernel::~CudaCalcFeedInForceKernel() {
}

void CudaCalcFeedInForceKernel::initialize(const System& system, const FeedInForce& force) {
    cu.setAsCurrent();

    forces_in = CudaArray::create<double>(cu, 3*system.getNumParticles(), "forces_in");
    // Initialize with 0.0
    forces_in->upload(vector<double>(3*system.getNumParticles(), 0.0));

    // Create kernel
    std::map<std::string, std::string> replacements;
    std::map<std::string, std::string> defines;
    defines["NUM_ATOMS"] = cu.intToString(system.getNumParticles());
    defines["PADDED_NUM_ATOMS"] = cu.intToString(cu.getPaddedNumAtoms());
    CUmodule module = cu.createModule(cu.replaceStrings(CudaFeedInKernelSources::FeedInForce, replacements), defines);
    computeForceKernel = cu.getKernel(module, "computeForceKernel");
    cu.addForce(new CudaFeedInForceInfo(force));
}

double CudaCalcFeedInForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    void* args_force[] = {
        &cu.getForce().getDevicePointer(),
        &forces_in->getDevicePointer(),
        &cu.getEnergyBuffer().getDevicePointer()
    };
    cu.executeKernel(computeForceKernel, args_force, context.getSystem().getNumParticles());
    return 0.0;
}

