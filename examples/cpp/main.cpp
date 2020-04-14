// -----------------------------------------------------------------------------
//           OpenMM(tm) HelloArgon example in C++ (June 2009)
// -----------------------------------------------------------------------------
// This program demonstrates a simple molecular simulation using the OpenMM
// API for GPU-accelerated molecular dynamics simulation. The primary goal is
// to make sure you can compile, link, and run with OpenMM and view the output.
// The example is available in C++, C, and Fortran 95.
//
// The system modeled here is a small number of argon atoms in a vacuum.
// A multi-frame PDB file is written to stdout which  can be read by VMD or
// other visualization tool to produce an animation of the resulting trajectory.
// -----------------------------------------------------------------------------

#include "OpenMM.h"
#include "FeedInForce.h"
#include <iostream>

// Forward declaration of routine for printing one frame of the
// trajectory, defined later in this source file.
void writePdbFrame(int frameNum, const OpenMM::State&);

void simulate()
{
    // Load any shared libraries containing CPU/GPU implementations.
    // If OpenMM is installed into a default directory, then it can be done with:
    //  OpenMM::Platform::loadPluginsFromDirectory( OpenMM::Platform::getDefaultPluginsDirectory() );
    // If non-standard location is used, explicit path to OpenMM/lib/plugins directory should be provided.
    OpenMM::Platform::loadPluginsFromDirectory("/home/yevhen/programs/openmm/lib/plugins");

    // Create a System with nonbonded forces.
    OpenMM::System system;
    // Create CustomNonbondedForce. We use it to implement repulsive potential between particles.
    // Only keep repulsive part of L-J potential
    OpenMM::CustomNonbondedForce* nonbond = new OpenMM::CustomNonbondedForce("(sigma/r)^12; sigma=0.5*(sigma1+sigma2)"); // http://docs.openmm.org/latest/userguide/theory.html#custom-forces
                                                                                                                         // http://docs.openmm.org/latest/userguide/theory.html#writing-custom-expressions
    nonbond->addPerParticleParameter("sigma");
    // Here we don't use cutoff for repulsive potential, but adding cutoff might speed things up a bit as interaction
    // strength decays very fast.
    nonbond->setNonbondedMethod(OpenMM::CustomNonbondedForce::NonbondedMethod::NoCutoff);
    // We can add cutoffs if needed.
    /*
    nonbond->setNonbondedMethod(OpenMM::CustomNonbondedForce::NonbondedMethod::CutoffNonPeriodic);
    float cutoff_distance = 2.0;
    nonbond->setCutoffDistance(cutoff_distance);
    */
    // We can also add swithching function to make the interaction smoothly go to zero over at swithching distance.
    /*
    nonbond->setUseSwitchingFunction(true);
    float switching_distance = 1.95;          // Switching distance has to be less than cutoff_distance
    nonbond->setSwitchingDistance(switching_distance);
    */
    system.addForce(nonbond);
    /* FeedInForce is an OpenMM Force object that serves as the interface to feed forces into the simulation. It stores
     * additional forces that we want to introduce internally and has a method
     * .updateForceInContext(OpenMM::Context context, vector< vector<double> > in_forces). in_forces should be a vector of
     * 3d vectors - one 3d vector per particle and each 3d vector should contain X,Y,Z components of the force.
     * updateForceInContext() will copy forces from in_forces into internal variable that will keep them
     * until updateForceInContext() is called again. Every simulation step OpenMM will add forces stored in FeedInForce to
     * the simulation. You don't need to call updateForceInContext() every step if the forces didn't change - OpenMM will
     * use the forces from the last call of updateForceInContext().
    */
    FeedInPlugin::FeedInForce* in_force = new FeedInPlugin::FeedInForce();
    system.addForce(in_force);

    // Create three atoms.
    std::vector<OpenMM::Vec3> initPosInNm(3);
    for (int a = 0; a < 3; ++a)
    {
        initPosInNm[a] = OpenMM::Vec3(0.5*a, 0.0, 0.0); // positions, nm
        // Populate the system with particles.
        // system.addParticle(mass in atomic mass units)
        system.addParticle(100.0);
        // Add particles to CustomNonbondedForce. Here we define repulsion strength sigma for each particle. If particles
        // have different value of sigma, combination rule sigma=0.5*(sigma1+sigma2) will be applied as defined in
        // CustomNonbondedForce force expression.
        float sigma = 1.0;
        nonbond->addParticle(std::vector<double>(1 , sigma));
    }

    float temperature   = 298.0;    // Temperature (K)
    float frictionCoeff = 1.0;      // Friction coefficient (1/ps)
    float step_size     = 0.002;    // Time step (ps)
    OpenMM::BrownianIntegrator integrator(temperature, frictionCoeff, step_size);

    OpenMM::Platform& platform = OpenMM::Platform::getPlatformByName("CUDA");
    OpenMM::Context context(system, integrator, platform);
    printf( "REMARK  Using OpenMM platform %s\n",
        context.getPlatform().getName().c_str() );

    // Set starting positions of the atoms. Leave time and velocity at zero.
    context.setPositions(initPosInNm);

    std::vector< std::vector<double> > forces_vec; // Will use this vector to feed additional forces into the simulation
    for (int i=0; i<system.getNumParticles(); ++i){
        std::vector<double> force_on_particle;
        force_on_particle.push_back(1.0);          // X component of the force
        force_on_particle.push_back(2.0);          // Y component of the force
        force_on_particle.push_back(3.0);          // Z component of the force
        forces_vec.push_back(force_on_particle);
    }

    // Simulate.
    for (int frameNum=1; frameNum<5;++frameNum) {
        // Feed additional force into the simulation.
        in_force->updateForceInContext(context, forces_vec); // Here we actually add force to the system.
                                                             // We will initialize forces_vec with forces produced by substrate concentration field and update forces when needed.
        integrator.step(1); // Advance simulation for 1 step. If we don't want to update external forces every step, we can do integrator.step(n)
        // Output current state information.
        OpenMM::State state = context.getState(OpenMM::State::Forces+OpenMM::State::Positions+OpenMM::State::Energy);
        std::cout << "Forces: " << std::endl;
        for (int i=0; i<state.getForces().size(); ++i){
            std::cout << state.getForces()[i] << std::endl;
        }
        std::cout << "Energy " << state.getPotentialEnergy() << std::endl;
        writePdbFrame(frameNum, state); // output coordinates
    }
}

int main()
{
    try {
        simulate();
        return 0; // success!
    }
    // Catch and report usage and runtime errors detected by OpenMM and fail.
    catch(const std::exception& e) {
        printf("EXCEPTION: %s\n", e.what());
        return 1; // failure!
    }
}

// Handy homebrew PDB writer for quick-and-dirty trajectory output.
void writePdbFrame(int frameNum, const OpenMM::State& state)
{
    // Reference atomic positions in the OpenMM State.
    const std::vector<OpenMM::Vec3>& posInNm = state.getPositions();

    // Use PDB MODEL cards to number trajectory frames
    printf("MODEL     %d\n", frameNum); // start of frame
    for (int a = 0; a < (int)posInNm.size(); ++a)
    {
        printf("ATOM  %5d  AR   AR     1    ", a+1); // atom number
        printf("%8.3f%8.3f%8.3f  1.00  0.00\n",      // coordinates
            // "*10" converts nanometers to Angstroms
            posInNm[a][0]*10, posInNm[a][1]*10, posInNm[a][2]*10);
    }
    printf("ENDMDL\n"); // end of frame
}
