%module FeedInplugin

%import(module="simtk.openmm") "swig/OpenMMSwigHeaders.i"
%include "swig/typemaps.i"

/*
 * The following lines are needed to handle std::vector.
 * Similar lines may be needed for vectors of vectors or
 * for other STL types like maps.
 */

%include "std_vector.i"
namespace std {
  %template(vectord) vector<double>;
  %template(VecVecdouble) vector< vector<double> >;
  %template(vectori) vector<int>;
};

%{
#include "FeedInForce.h"
#include "OpenMM.h"
#include "OpenMMAmoeba.h"
#include "OpenMMDrude.h"
#include "openmm/RPMDIntegrator.h"
#include "openmm/RPMDMonteCarloBarostat.h"
%}

%pythoncode %{
import simtk.openmm as mm
import simtk.unit as unit
%}

/*
 * Convert C++ exceptions to Python exceptions.
*/
%exception {
    try {
        $action
    } catch (std::exception &e) {
        PyErr_SetString(PyExc_Exception, const_cast<char*>(e.what()));
        return NULL;
    }
}


namespace FeedInPlugin {

class FeedInForce : public OpenMM::Force {
public:
    FeedInForce();
    
    void updateForceInContext(OpenMM::Context& context, std::vector<std::vector<double> > const &inForces);

    /*
     * Add methods for casting a Force to an FeedInForce.
    */
    %extend {
        static FeedInPlugin::FeedInForce& cast(OpenMM::Force& force) {
            return dynamic_cast<FeedInPlugin::FeedInForce&>(force);
        }

        static bool isinstance(OpenMM::Force& force) {
            return (dynamic_cast<FeedInPlugin::FeedInForce*>(&force) != NULL);
        }
    }
};

}
