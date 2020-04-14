#####################################################################################
# This is an example usage of the interface to feed forces into OpenMM with Python. #
#####################################################################################

from simtk.unit import *
from simtk.openmm import *
from simtk.openmm.app import *
import numpy
from FeedInplugin import FeedInForce


def caclSubstrateForces(positions):
    # Do some magic here that will return substrate field-related force acting on each particle
    forces = []
    for i in range(3):
        force_on_particle = [0.1, 0.2, 0.3]
        forces.append(force_on_particle)
    return forces


# Here, we use a pdb just to create topology latter on.
# Alternatively one can create topology object from scratch.
pdb = PDBFile('test_input.pdb')

# Create an OpenMM System object
system = System()

# Create CustomNonbondedForce. We use it to implement repulsive potential between particles.
# Only keep repulsive part of L-J potential
nonbond = CustomNonbondedForce("(sigma/r)^12; sigma=0.5*(sigma1+sigma2)")
nonbond.addPerParticleParameter("sigma")
# Here we don't use cutoff for repulsive potential, but adding cutoff might speed things up a bit as interaction
# strength decays very fast.
nonbond.setNonbondedMethod(CustomNonbondedForce.NoCutoff)
# Add force to the System
system.addForce(nonbond)

# FeedInForce is an OpenMM Force object that serves as the interface to feed forces into the simulation. It stores
# additional forces that we want to introduce internally and has a method
# .updateForceInContext(OpenMM::Context context, vector< vector<double> > in_forces). in_forces should be a vector of
# 3d vectors - one 3d vector per particle and each 3d vector should contain X,Y,Z components of the force.
# updateForceInContext() will copy forces from in_forces into internal variable that will keep them
# until updateForceInContext() is called again. Every simulation step OpenMM will add forces stored in FeedInForce to
# the simulation. You don't need to call updateForceInContext() every step if the forces didn't change - OpenMM will
# use the forces from the last call of updateForceInContext().
in_force = FeedInForce()
# Add in_force to the system
system.addForce(in_force)

num_particles = len(pdb.getPositions())
initPosInNm = []
for i in range(num_particles):
    initPosInNm.append([0.5 * i, 0.0, 0.0])
    # Populate the system with particles.
    # system.addParticle(mass in atomic mass units)
    system.addParticle(100.0)  # 100.0 is a particle mass
    # Add particles to CustomNonbondedForce. Here we define repulsion strength sigma for each particle. If particles
    # have different value of sigma, combination rule sigma=0.5*(sigma1+sigma2) will be applied as defined in
    # CustomNonbondedForce force expression.
    sigma = 1.0
    nonbond.addParticle([sigma])

# Simulation parameters
temperature = 298.0  # K
frictionCoeff = 1.0  # 1/ps
step_size = 0.002  # ps

# Create integrator
integrator = BrownianIntegrator(temperature, frictionCoeff, step_size)
# Create platform
platform = Platform.getPlatformByName('CUDA')
# Create simulation
simulation = Simulation(pdb.topology, system, integrator, platform)
print('REMARK  Using OpenMM platform %s' % simulation.context.getPlatform().getName())
# Set initial positions for the particles
simulation.context.setPositions(initPosInNm)

# Simulate
# We can add dcd or other reporters to the simulation to get desired output with
# simulation.reporters.append(DCDReporter(file_name, reportInterval, append=False, enforcePeriodicBox=None))
for num_inter in range(1, 3):
    # forces_vec will contain external forces that we want to feed into OpenMM.
    # Get current positions of the particles for concentration field/forces calculation
    forces_vec = caclSubstrateForces(simulation.context.getState(getPositions=True).getPositions())
    # Feed external forces into OpenMM
    in_force.updateForceInContext(simulation.context, forces_vec)
    # Advance simulation for 1 steps
    integrator.step(1)
    state = simulation.context.getState(getEnergy=True, getForces=True)
    print(state.getForces())
