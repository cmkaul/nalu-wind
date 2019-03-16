
#ifndef ABLDAMPINGALGORITHM_H
#define ABLDAMPINGALGORITHM_H

#include "NaluParsing.h"
#include "FieldTypeDef.h"

#include "stk_mesh/base/Selector.hpp"

#include <string>
#include <vector>
#include <stdexcept>
#include <sstream>
#include <unordered_set>

namespace sierra {
namespace nalu {

class Realm;

/**
 * \brief ABL Damping Source terms for Momentum and Temperature equations
 *
 * This class parses the user inputs and provides a common interface to the
 * momentum and temperature ABL damping source term implementations within Nalu.
 * The ABL forcing capability is turned on by the presence of a sub-section
 * titled `abl_damping` within the Realm section of the Nalu input file.
 *
 * ```
 *   abl_forcing:
 *     momentum:
 *       relaxation_factor: 1.0
 *       thickness: 500.0
 *
 *     temperature:
 *       relaxation_factor: 1.0
 *       thickness: 500.0
 * ```
 *
 * There are two optional sub-sections in `abl_damping`: "momentum" and
 * "temperature".
 */
class ABLDampingAlgorithm
{
public:
  template <typename T>
  using Array2D = std::vector<std::vector<T>>;

  /**
   * Types of ABL damping available
   */
  enum ABLDampingTypes {
    OFF = 0,        //!< No ABL Damping applied
    ON = 1,         //!< Use Damping
    NUM_ABL_DAMPING_TYPES //!< Guard
  };

  ABLDampingAlgorithm(Realm&, const YAML::Node&);

  ~ABLDampingAlgorithm();

  //! Parse input file for user options and initialize
  void load(const YAML::Node&);

  //! Initialize ABL damping (steps after mesh creation)
  void initialize();

  //! Execute field transfers, compute planar averaging, and determine source
  //! terms at desired levels.
  void execute();

  //! Evaluate the ABL damping source contribution at a node
  //! CK: I think this will need to be modified
  //! CK: I think vector here is the 3 components, we will need 3*nz
  void eval_momentum_source(
    const double,        //!< Height of the node from terrain
    std::vector<double>& //!< Source vector to be populated
    );

  //! Evaluate the ABL forcing source contribution (temperature)
  //! CK: Will need to modify this (ie double to a vector)
  void eval_temperature_source(
    const double, //!< Height of the node from terrain
    double&       //!< Temperature source term to be populated
    );

  inline bool momentumDampingOn() { return (momSrcType_ != OFF); }

  inline bool temperatureDampingOn() { return (tempSrcType_ != OFF); }

  inline bool ablDampingOn()
  {
    return (momentumDampingOn() || temperatureDampingOn());
  }

private:
  ABLDampingAlgorithm();
  ABLDampingAlgorithm(const ABLDampingAlgorithm&);

  //! Utility function to parse momentum forcing options from input file.
  void load_momentum_info(const YAML::Node&);

  //! Helper method to parse temperature forcing options from input file.
  void load_temperature_info(const YAML::Node&);

  //! Create 2-D interpolation lookup tables from YAML data structures
  //! CK: I don't think I will need this
  void create_interp_arrays(
    const std::vector<double>::size_type,
    const Array2D<double>&,
    std::vector<double>&,
    Array2D<double>&);

  //! Compute mean velocity and estimate source term for a given timestep
  //! CK: this should be something I need
  void compute_momentum_sources();

  //! Compute average planar temperature and estimate source term
  //! CK: this should be something I need
  void compute_temperature_sources();

  //! Reference to Realm
  Realm& realm_;

  //! Momentum Forcing Source Type
  ABLDampingTypes momSrcType_;

  //! Temperature Forcing Source Type
  ABLDampingTypes tempSrcType_;

  //! Relaxation factor for momentum sources
  //! CK: make this a vector (was double) to allow height-dependence
  double alphaMomentumMax_;
  std::vector<double> alphaMomentum_;

  //! Relaxation factor for temperature sources
  //! CK: make this a vector (was a double) to allow height-dependence
  double alphaTemperatureMax_;
  std::vector<double> alphaTemperature_;

  //! Thickness of the momentum damping layer 
  double momThickness_; // 

  //! Thickness of the temperature damping layer
  double tempThickness_; // Array of shape [num_Theights]

  //! Planar average velocity calculated on the surface [num_UHeights, 3]
  //! CK: This will need to be changed to all z levels within momThickness
  // or simply at all z levels
  //! Or in other words num_UHeights=nz
  //! CK: should enforce UmeanCalc[:,2]=0 ?
  Array2D<double> UmeanCalc_;

  //! Planar average density calculated on the surface [num_UHeights]
  //! CK: this seems unnecessary for Boussinesq sys, but may be forward looking
  std::vector<double> rhoMeanCalc_;

  //! U source as a function of height [3,num_UHeights]
  Array2D<double> USource_;

  //! Planar average temperature calculated on the surface [num_THeights]
  std::vector<double> TmeanCalc_;

  //! T source as a function of height [num_THeights]
  std::vector<double> TSource_;

  //! Write frequency for source term output
  //! CK: I am not sure this is needed. I believe this is related to .dat files
  int outputFreq_{10};

  //! Format string specifier indicating the file name for output. The
  //! specification takes one `%s` specifier that is used to populate Ux, Uy,
  //! Uz, T. Default is "abl_sources_%s.dat"
  std::string outFileFmt_{"abl_%s_damping.dat"};
};
}
}

#endif /* ABLDAMPINGALGORITHM_H */
