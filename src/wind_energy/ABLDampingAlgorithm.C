
#include "wind_energy/ABLDampingAlgorithm.h"
#include "Realm.h"
#include "NaluEnv.h"
#include "xfer/Transfer.h"
#include "xfer/Transfers.h"
#include "utils/LinearInterpolation.h"
#include "wind_energy/BdyLayerStatistics.h"
//CK I think this is needed
#include "wind_energy/BdyHeightAlgorithm.h"


// stk_mesh/base/fem
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/GetBuckets.hpp>
#include <stk_mesh/base/Selector.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Part.hpp>

#include <stk_io/IossBridge.hpp>

#include <stk_util/parallel/ParallelReduce.hpp>

#include <boost/format.hpp>
#include <fstream>
#include <iostream>
#include <iomanip>

namespace sierra {
namespace nalu {

ABLDampingAlgorithm::ABLDampingAlgorithm(Realm& realm, const YAML::Node& node)
  : minDampingHeightMomentum(-999.),
    minDampingHeightTemperature(-999.),
    realm_(realm),
    momSrcType_(ABLDampingAlgorithm::OFF),
    tempSrcType_(ABLDampingAlgorithm::OFF),
    gammaMomentum_(0.001),
    gammaTemperature_(0.001),
    velHeights_(0),
    tempHeights_(0),
    velXTimes_(0),
    velYTimes_(0),
    velZTimes_(0),
    tempTimes_(0),
    velX_(0),
    velY_(0),
    velZ_(0),
    temp_(0)
{
  if (realm_.bdyLayerStats_ == nullptr)
    throw std::runtime_error("ABL Damping requires ABL Boundary Layer statistics");
  load(node);
}

ABLDampingAlgorithm::~ABLDampingAlgorithm()
{}

void
ABLDampingAlgorithm::load(const YAML::Node& node)
{
  /* CK: Output, staging for removal
  get_if_present(node, "output_frequency", outputFreq_, outputFreq_);
  get_if_present(node, "output_format", outFileFmt_, outFileFmt_);
  */
  NaluEnv::self().naluOutputP0()<< "ABLDampingAlgorithm::load." << std::endl;
 

  if (node["momentum"])
    load_momentum_info(node["momentum"]);

  if (node["temperature"])
    load_temperature_info(node["temperature"]);
}

void
ABLDampingAlgorithm::load_momentum_info(const YAML::Node& node)
{
  std::string mom_type = node["type"].as<std::string>();
  NaluEnv::self().naluOutputP0()<< "ABLDampingAlgorithm::load_momentum_info." << std::endl;
  if (mom_type == "given_profile") {
    momSrcType_ = ABLDampingAlgorithm::GIVEN_PROFILE;
  } else if (mom_type == "mean_profile") {
    momSrcType_ = ABLDampingAlgorithm::MEAN_PROFILE;
    throw std::runtime_error(
      "ABLDampingAlgorithm: mean_profile type specification for momentum. "
      "This type is not ready yet.");
  } else {
    throw std::runtime_error(
      "ABLDampingAlgorithm: Invalid type specification for momentum. ");
  }
  get_if_present(node, "relaxation_factor", gammaMomentum_, gammaMomentum_);
  get_if_present(node, "base_height", minDampingHeightMomentum, 
                   minDampingHeightMomentum);
  // auto nHeights = velHeights_.size();
  if (momSrcType_ == GIVEN_PROFILE){
    get_required<std::vector<double>>(node, "heights", velHeights_);
    auto nHeights = velHeights_.size();
    // Load momentum source time histories in temporary data structures, test
    // consistency of data with heights and then recast them into 2-D lookup
    // arrays.
    Array2D<double> vxtmp, vytmp, vztmp;
    get_required<Array2D<double>>(node, "velocity_x", vxtmp);
    get_required<Array2D<double>>(node, "velocity_y", vytmp);
    get_required<Array2D<double>>(node, "velocity_z", vztmp);

    create_interp_arrays(nHeights, vxtmp, velXTimes_, velX_);
    create_interp_arrays(nHeights, vytmp, velYTimes_, velY_);
    create_interp_arrays(nHeights, vztmp, velZTimes_, velZ_);
  } else if (momSrcType_ == MEAN_PROFILE){
    //! Do something else?
  }
}

void
ABLDampingAlgorithm::load_temperature_info(const YAML::Node& node)
{
  std::string temp_type = node["type"].as<std::string>();

  NaluEnv::self().naluOutputP0()<< "ABLDampingAlgorithm::load._temperature_info" << std::endl;
  if (temp_type == "given_profile") {
    tempSrcType_ = ABLDampingAlgorithm::GIVEN_PROFILE;
  } else if (temp_type == "mean_profile") {
    tempSrcType_ = ABLDampingAlgorithm::MEAN_PROFILE;
    throw std::runtime_error(
      "ABLDampingAlgorithm: mean_profile type specification for temperature. "
      "This type is not ready yet.");
  } else {
    throw std::runtime_error(
      "ABLDampingAlgorithm: Invalid type specification for temperature. ");
  }
  get_if_present(node, "relaxation_factor", gammaTemperature_, gammaTemperature_);
  get_if_present(node, "base_height", minDampingHeightTemperature, minDampingHeightTemperature);
  if (tempSrcType_ == GIVEN_PROFILE){  
    get_required<std::vector<double>>(node, "heights", tempHeights_);
    auto nHeights = tempHeights_.size();
    // Load temperature source time histories, check data consistency and create
    // interpolation lookup tables.
    Array2D<double> temp;
    get_required<Array2D<double>>(node, "temperature", temp);
    
    create_interp_arrays(nHeights, temp, tempTimes_, temp_);

  }else if (tempSrcType_ == MEAN_PROFILE){
    //! Do something else?

  }

}

void
ABLDampingAlgorithm::create_interp_arrays(
  const std::vector<double>::size_type nHeights,
  const Array2D<double>& inpArr,
  std::vector<double>& outTimes,
  Array2D<double>& outValues)
{
  /* The input vector is shaped [nTimes, nHeights+1]. We transform it to two
   * arrays:
   *    time[nTimes] = inp[nTimes,0], and
   *    value[nHeights, nTimes] -> swap rows/cols from input
   */
    NaluEnv::self().naluOutputP0() << "ABLDampingAlgorithm::create_interp_arrays" << std::endl;
  // Check that all timesteps contain values for all the heights
  for (auto vx : inpArr) {
    ThrowAssert((vx.size() == (nHeights + 1)));
  }

  auto nTimes = inpArr.size();
  outTimes.resize(nTimes);
  outValues.resize(nHeights);
  for (std::vector<double>::size_type i = 0; i < nHeights; i++) {
    outValues[i].resize(nTimes);
  }
  for (std::vector<double>::size_type i = 0; i < nTimes; i++) {
    outTimes[i] = inpArr[i][0];
    for (std::vector<double>::size_type j = 0; j < nHeights; j++) {
      outValues[j][i] = inpArr[i][j + 1];
    }
  }
}

void
ABLDampingAlgorithm::initialize()
{ 
  NaluEnv::self().naluOutputP0()<< "ABLDampingAlgorithm::initialize" << std::endl;
  auto*  bdyLayerStats = realm_.bdyLayerStats_;
  //! CK not sure about this! 
  const std::vector<double>& ablHeights = bdyLayerStats->abl_heights();
  const int nAblHeights = bdyLayerStats->abl_num_levels();
  const int ndim = realm_.spatialDimension_;
  const double zTop = ablHeights[nAblHeights-1];
  const double zdMom = zTop - minDampingHeightMomentum;
  const double zdTemp = zTop - minDampingHeightTemperature;
  const double halfPi = 0.5 * std::acos(-1.0);
  double sinArg;
  
  

  if (momSrcType_ != OFF) {
    NaluEnv::self().naluOutputP0()
      << "ABL Damping active for Momentum Equations\n"<< std::endl;
    dampingCoeffMomentum.resize(nAblHeights);
    UDamp.resize(nAblHeights);

    for (int i = 0; i < nAblHeights; i++) {
      UDamp[i].resize(ndim);
      if (ablHeights[i]< minDampingHeightMomentum){
        dampingCoeffMomentum[i] = 0.0;
      }else{
        sinArg = halfPi * (1.0 - (zTop - ablHeights[i]) / zdMom);
        dampingCoeffMomentum[i] = gammaMomentum_ * std::pow(std::sin(sinArg),2);
      }
    }
  }

  if (tempSrcType_ != OFF) {
    NaluEnv::self().naluOutputP0()
      << "ABL Damping active for Temperature Equation\n"
      << std::endl;
    dampingCoeffTemperature.resize(nAblHeights);
    TDamp.resize(nAblHeights);
    for (int i = 0; i < nAblHeights; i++) {
      if (ablHeights[i]< minDampingHeightMomentum){
        dampingCoeffTemperature[i] = 0.0;
      }else{
        sinArg = halfPi * (1.0 - (zTop - ablHeights[i]) / zdTemp);
        dampingCoeffTemperature[i] = gammaTemperature_ * std::pow(std::sin(sinArg),2);
      }
    }
  }
  /*
  // Prepare output files to dump sources when computed during precursor phase
  if (( NaluEnv::self().parallel_rank() == 0 ) &&
      ( momSrcType_ == COMPUTED )) {
    std::string uxname((boost::format(outFileFmt_)%"Ux").str());
    std::string uyname((boost::format(outFileFmt_)%"Uy").str());
    std::string uzname((boost::format(outFileFmt_)%"Uz").str());
    std::fstream uxFile, uyFile, uzFile;
    uxFile.open(uxname.c_str(), std::fstream::out);
    uyFile.open(uyname.c_str(), std::fstream::out);
    uzFile.open(uzname.c_str(), std::fstream::out);

    uxFile << "# Time, " ;
    uyFile << "# Time, " ;
    uzFile << "# Time, " ;
    for (size_t ih = 0; ih < velHeights_.size(); ih++) {
      uxFile << velHeights_[ih] << " ";
      uyFile << velHeights_[ih] << " ";
      uzFile << velHeights_[ih] << " ";
    }
    uxFile << std::endl ;
    uyFile << std::endl ;
    uzFile << std::endl ;
    uxFile.close();
    uyFile.close();
    uzFile.close();
  } */
}

void
ABLDampingAlgorithm::execute()
{  
  NaluEnv::self().naluOutputP0()<< "ABLDampingAlgorithm::execute" << std::endl;
  if (momentumDampingOn())
    compute_momentum_target_profile();

  if (temperatureDampingOn())
    compute_temperature_target_profile();
}

void
ABLDampingAlgorithm::compute_momentum_target_profile()

{
  const double currTime = realm_.get_current_time();
  auto* bdyLayerStats = realm_.bdyLayerStats_;
  //! CK not sure about this! 
  const std::vector<double>& ablHeights = bdyLayerStats->abl_heights();
  const int nAblHeights = bdyLayerStats->abl_num_levels();
  const int nInputHeights = velHeights_.size();
  std::vector<double> timeInterpVelX;
  std::vector<double> timeInterpVelY; 
  std::vector<double> timeInterpVelZ; 

  timeInterpVelX.resize(nInputHeights);
  timeInterpVelY.resize(nInputHeights);
  timeInterpVelZ.resize(nInputHeights);
  

  if (momSrcType_ == MEAN_PROFILE) {
    for (int ih=0; ih < nInputHeights; ih++) {
      bdyLayerStats->velocity(ablHeights[ih], UDamp[ih].data());
    }
  }else if( momSrcType_ == GIVEN_PROFILE){
      //! First interp in time
      for (int ih = 0; ih < nInputHeights; ih++) {
        // Interpolate the velocities from the table to the current time
        utils::linear_interp(velXTimes_, velX_[ih], currTime, timeInterpVelX[ih]);
        utils::linear_interp(velYTimes_, velY_[ih], currTime, timeInterpVelY[ih]);
        utils::linear_interp(velZTimes_, velZ_[ih], currTime, timeInterpVelZ[ih]);
      }
      //! Second interp in height
      for (int ih=0; ih<nAblHeights; ih++) {
        utils::linear_interp(velHeights_, timeInterpVelX, ablHeights[ih], UDamp[ih][0]);
        utils::linear_interp(velHeights_, timeInterpVelY, ablHeights[ih], UDamp[ih][1]);
        utils::linear_interp(velHeights_, timeInterpVelZ, ablHeights[ih], UDamp[ih][2]);
      }

  }

 

  /*
  const int tcount = realm_.get_time_step_count();
  if (( NaluEnv::self().parallel_rank() == 0 ) &&
      ( momSrcType_ == COMPUTED ) &&
      ( tcount % outputFreq_ == 0)) {
    std::string uxname((boost::format(outFileFmt_)%"Ux").str());
    std::string uyname((boost::format(outFileFmt_)%"Uy").str());
    std::string uzname((boost::format(outFileFmt_)%"Uz").str());
    std::fstream uxFile, uyFile, uzFile;
    uxFile.open(uxname.c_str(), std::fstream::app);
    uyFile.open(uyname.c_str(), std::fstream::app);
    uzFile.open(uzname.c_str(), std::fstream::app);

    uxFile << std::setw(12) << currTime << " ";
    uyFile << std::setw(12) << currTime << " ";
    uzFile << std::setw(12) << currTime << " ";
    for (size_t ih = 0; ih < velHeights_.size(); ih++) {
      uxFile << std::setprecision(6)
             << std::setw(15)
             << USource_[0][ih] << " ";
      uyFile << std::setprecision(6)
             << std::setw(15)
             << USource_[1][ih] << " ";
      uzFile << std::setprecision(6)
             << std::setw(15)
             << USource_[2][ih] << " ";
    }
    uxFile << std::endl;
    uyFile << std::endl;
    uzFile << std::endl;
    uxFile.close();
    uyFile.close();
    uzFile.close();
  } */
}


void
ABLDampingAlgorithm::compute_temperature_target_profile()
{
  const double currTime = realm_.get_current_time();
  auto* bdyLayerStats = realm_.bdyLayerStats_;
  //! CK not sure about this! 
  const std::vector<double>& ablHeights = bdyLayerStats->abl_heights();
  const int nAblHeights = bdyLayerStats->abl_num_levels();
  const int nInputHeights = tempHeights_.size();
  std::vector<double> timeInterpTemp;

  timeInterpTemp.resize(nInputHeights);

  if (tempSrcType_ == MEAN_PROFILE) {
    for (size_t ih=0; ih < tempHeights_.size(); ih++) {
      bdyLayerStats->temperature(tempHeights_[ih], &TDamp[ih]);
    }
  }else if (tempSrcType_ == GIVEN_PROFILE){
    //! First interp in time
    for (int ih = 0; ih < nInputHeights; ih++) {
      utils::linear_interp(tempTimes_, temp_[ih], currTime, timeInterpTemp[ih]);
    }
    //! Second interp in height
    for (int ih =0; ih < nAblHeights; ih++){
      utils::linear_interp(tempHeights_, timeInterpTemp, ablHeights[ih], TDamp[ih]);

    }
  }




/*
  const int tcount = realm_.get_time_step_count();
  if (( NaluEnv::self().parallel_rank() == 0 ) &&
      ( tempSrcType_ == COMPUTED ) &&
      ( tcount % outputFreq_ == 0)) {
    std::string fname((boost::format(outFileFmt_)%"T").str());
    std::fstream tFile;
    tFile.open(fname.c_str(), std::fstream::app);

    tFile << currTime << " ";
    for (size_t ih = 0; ih < tempHeights_.size(); ih++) {
      tFile << std::setprecision(6)
            << std::setw(15)
            << TSource_[ih] << " ";
    }
    tFile << std::endl;
    tFile.close();
  }
  */
}


} // namespace nalu
} // namespace sierra
