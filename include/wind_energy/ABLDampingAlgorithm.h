
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
 * momentum and temperature ABL forcing source term implementations within Nalu.
 * The ABL forcing capability is turned on by the presence of a sub-section
 * titled `abl_damping` within the Realm section of the Nalu input file.
 *
 * ```
 *   abl_damping:
 *     momentum:
 *       type: given_profile
 *       relaxation_factor: 1.0
 *       heights: [80.0]
 *       velocity_x:
 *         - [0.0, 10.0]                 # [Time0, vxH0, ... , vxHN]
 *         - [100000.0, 10.0]            # [TimeN, vxH0, ... , vxHN]
 *
 *       velocity_y:
 *         - [0.0, 0.0]
 *         - [10000.0, 0.0]
 *
 *       velocity_z:
 *         - [0.0, 0.0]
 *         - [10000.0, 0.0]
 *
 *     temperature:
 *       type: given_profile
 *       relaxation_factor: 1.0
 *       heights: [80.0]
 *       temperature:
 *         - [0.0, 300.0]
 *         - [10000.0, 300.0]
 * ```
 *
 * There are two optional sub-sections in `abl_forcing`: "momentum" and
 * "temperature".
 */
class ABLDampingAlgorithm
{
public:
  template <typename T>
  using Array2D = std::vector<std::vector<T>>;

  /**
   * Types of ABL forcing available
   */
  enum ABLDampingTypes {
    OFF = 0,              //!< No ABL damping applied
    GIVEN_PROFILE = 1,     //!< Target profile provided by user
    MEAN_PROFILE = 2,      //!< Target profile is the current mean profile (WIP)
    NUM_ABL_DAMPING_TYPES //!< Guard
  };

  //! Smoothly increasing damping coefficient profile for Momentum
  double gammaMomentum_;
  std::vector<double> dampingCoeffMomentum;

  //! Smoothly increasing damping coefficient profile for Temperature
  double gammaTemperature_;
  std::vector<double> dampingCoeffTemperature;
 
  //! Base of the damping layer for momentum
  double minDampingHeightMomentum;

  //! Base of the damping layer for temperature
  double minDampingHeightTemperature;

  //! Target velocity profile on the ABL vertical grid [nHeights, 3]
  Array2D<double> UDamp;

  //! Target temperature profile on the ABL vertical grid [nHeights]
  std::vector<double> TDamp;

  ABLDampingAlgorithm(Realm&, const YAML::Node&);

  ~ABLDampingAlgorithm();

  //! Parse input file for user options and initialize
  void load(const YAML::Node&);

  //! Initialize ABL damping (steps after mesh creation)
  void initialize();

  //! Execute field transfers, compute planar averaging, and determine source
  //! terms at desired levels.
  //! CK: TODO: Keep this methods but need to update descriptive comment
  void execute();

  /* CK: Tentatively removing
  //! Evaluate the ABL forcing source contribution at a node
  void eval_momentum_source(
    const double,        //!< Height of the node from terrain
    std::vector<double>& //!< Source vector to be populated
    );

  //! Evaluate the ABL forcing source contribution (temperature)
  void eval_temperature_source(
    const double, //!< Height of the node from terrain
    double&       //!< Temperature source term to be populated
    );
   */

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
  void create_interp_arrays(
    const std::vector<double>::size_type,
    const Array2D<double>&,
    std::vector<double>&,
    Array2D<double>&);

  //! Compute mean velocity and estimate source term for a given timestep
  void compute_momentum_target_profile();

  //! Compute average planar temperature and estimate source term
  void compute_temperature_target_profile();

  //! Reference to Realm
  Realm& realm_;

  //! Momentum Forcing Source Type
  ABLDampingTypes momSrcType_;

  //! Temperature Forcing Source Type
  ABLDampingTypes tempSrcType_;

  //! Maximum relaxation factor for momentum sources
  double gammaMomentum_;


  //! Maximum relaxation factor for temperature sources
  double gammaTemperature_;

  //! Heights where velocity information is provided
  std::vector<double> velHeights_; // Array of shape [num_Uheights]

  //! Heights where temperature information is provided
  std::vector<double> tempHeights_; // Array of shape [num_Theights]

  //! Times where velocity information is available
  std::vector<double> velXTimes_; // Arrays of shape [num_Utimes]
  std::vector<double> velYTimes_;
  std::vector<double> velZTimes_;

  //! Times where temperature information is available
  std::vector<double> tempTimes_; // Array of shape [num_Ttimes]

  // The following arrays are shaped [num_UHeights, num_Utimes]
  Array2D<double> velX_;
  Array2D<double> velY_;
  Array2D<double> velZ_;
  // The temperature array is shaped [num_Theights, num_Ttimes]
  Array2D<double> temp_;

  
  /* CK: staging for removal
  //! Planar average density calculated on the surface [num_UHeights]
  std::vector<double> rhoMeanCalc_;

  //! T source as a function of height [num_THeights]
  std::vector<double> TSource_;
  */

  /* CK: TODO: Not sure on whether to keep the output or not
  //! Write frequency for source term output
  int outputFreq_{10};

  //! Format string specifier indicating the file name for output. The
  //! specification takes one `%s` specifier that is used to populate Ux, Uy,
  //! Uz, T. Default is "abl_sources_%s.dat"
  std::string outFileFmt_{"abl_%s_sources.dat"};
  */
};
}
}

#endif /* ABLDAMPINGALGORITHM_H */
