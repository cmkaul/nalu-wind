
#include "wind_energy/ABLDampingAlgorithm.h"
#include "Realm.h"
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
  : realm_(realm),
    momSrcType_(ABLDampingAlgorithm::OFF),
    tempSrcType_(ABLDampingAlgorithm::OFF),
    alphaMomentumMax_(1.0),
    alphaTemperatureMax_(1.0);
    alphaMomentum_(1.0),
    alphaTemperature_(1.0),
    momThickness_(0),
    tempThickness_(0),
    UmeanCalc_(0),
    rhoMeanCalc_(0),
    USource_(0),
    TmeanCalc_(0),
    TSource_(0)
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
  get_if_present(node, "output_frequency", outputFreq_, outputFreq_);
  get_if_present(node, "output_format", outFileFmt_, outFileFmt_);

  if (node["momentum"])
    momSrcType_ = ABLDampingAlgorithm::ON
    load_momentum_info(node["momentum"]);

  if (node["temperature"])
    tempSrcType = ABLDampingAlgorithm::OFF
    load_temperature_info(node["temperature"]);
}

void
ABLForcingAlgorithm::load_momentum_info(const YAML::Node& node)
{
  get_if_present(node, "relaxation_factor", alphaMomentumMax_, alphaMomentumMax_);
  get_required<std::vector<double>>(node, "thickness", momThickness_);
  // auto nHeights = velHeights_.size();
}

void
ABLForcingAlgorithm::load_temperature_info(const YAML::Node& node)
{
  get_if_present(node, "relaxation_factor", alphaTemperatureMax_, alphaTemperatureMax_);
  get_required<std::vector<double>>(node, "thickness", tempThickness_);
  //(CK)get_required<std::vector<double>>(node, "heights", tempHeights_);
  //(CK)auto nHeights = tempHeights_.size();



  //(CK)TSource_.resize(nHeights);
  //(CK)if (tempSrcType_ == COMPUTED)
  //(CK)  TmeanCalc_.resize(nHeights);
}

void
ABLForcingAlgorithm::initialize()
{ 
  auto* bdyLayerStats = realm_.bdyLayerStats_;
  //! CK not sure about this! (right way to get heights_ array from bdyLayerStats?
  // Is resulting heights_ the values or a reference? 
  auto heights_ = bdyLayerStats->heights_
   
  if (momSrcType_ != OFF) {
    NaluEnv::self().naluOutputP0()
      << "ABL Forcing active for Momentum Equations\n"
      << "\t Number of planes: " << velHeights_.size()
      << "\n\t Number of time steps: " << velXTimes_.size() << std::endl;
  }

  const int ndim = realm_.spatialDimension_;
  if (momSrcType_ == COMPUTED) {
    UmeanCalc_.resize(nHeights);
    rhoMeanCalc_.resize(nHeights);
    
    for (size_t i = 0; i < nHeights; i++) {
      UmeanCalc_[i].resize(ndim);
    }
  }
  
  USource_.resize(ndim);
  for (int i = 0; i < ndim; i++) {
    USource_[i].resize(nHeights);
  }




  if (tempSrcType_ != OFF) {
    NaluEnv::self().naluOutputP0()
      << "ABL Forcing active for Temperature Equation\n"
      << "\t Number of planes: " << tempHeights_.size()
      << "\n\t Number of time steps: " << tempTimes_.size() << std::endl
      << std::endl;
  }

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
  }
}

void
ABLForcingAlgorithm::execute()
{
  if (momentumForcingOn())
    compute_momentum_sources();

  if (temperatureForcingOn())
    compute_temperature_sources();
}

void
ABLForcingAlgorithm::compute_momentum_sources()
{
  const double dt = realm_.get_time_step();
  const double currTime = realm_.get_current_time();

  if (momSrcType_ == COMPUTED) {
    auto* bdyLayerStats = realm_.bdyLayerStats_;
    for (size_t ih=0; ih < velHeights_.size(); ih++) {
      bdyLayerStats->velocity(velHeights_[ih], UmeanCalc_[ih].data());
      bdyLayerStats->density(velHeights_[ih], &rhoMeanCalc_[ih]);
    }
    for (size_t ih = 0; ih < velHeights_.size(); ih++) {
      double xval, yval;

      // Interpolate the velocities from the table to the current time
      utils::linear_interp(velXTimes_, velX_[ih], currTime, xval);
      utils::linear_interp(velYTimes_, velY_[ih], currTime, yval);

      // Compute the momentum source
      // Momentum source in the x direction
      USource_[0][ih] = rhoMeanCalc_[ih] * (alphaMomentum_ / dt) *
                          (xval - UmeanCalc_[ih][0]);
      // Momentum source in the y direction
      USource_[1][ih] = rhoMeanCalc_[ih] * (alphaMomentum_ / dt) *
                          (yval - UmeanCalc_[ih][1]);

      // No momentum source in z-direction
      USource_[2][ih] = 0.0;

    }

  } else {
    for (size_t ih = 0; ih < velHeights_.size(); ih++) {
      utils::linear_interp(velXTimes_, velX_[ih], currTime, USource_[0][ih]);
      utils::linear_interp(velYTimes_, velY_[ih], currTime, USource_[1][ih]);
      utils::linear_interp(velZTimes_, velZ_[ih], currTime, USource_[2][ih]);
    }
  }

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
  }
}

void
ABLForcingAlgorithm::compute_temperature_sources()
{
  const double dt = realm_.get_time_step();
  const double currTime = realm_.get_current_time();

  if (tempSrcType_ == COMPUTED) {
    auto* bdyLayerStats = realm_.bdyLayerStats_;
    for (size_t ih=0; ih < tempHeights_.size(); ih++) {
      bdyLayerStats->temperature(tempHeights_[ih], &TmeanCalc_[ih]);
    }
    for (size_t ih = 0; ih < tempHeights_.size(); ih++) {
      double tval;
      utils::linear_interp(tempTimes_, temp_[ih], currTime, tval);
      TSource_[ih] = (alphaTemperature_ / dt) * (tval - TmeanCalc_[ih]);
    }
  } else {
    for (size_t ih = 0; ih < tempHeights_.size(); ih++) {
      utils::linear_interp(tempTimes_, temp_[ih], currTime, TSource_[ih]);
    }
  }

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
}

void
ABLForcingAlgorithm::eval_momentum_source(
  const double zp, std::vector<double>& momSrc)
{
  const int nDim = realm_.spatialDimension_;
  if (velHeights_.size() == 1) {
    // Constant source term throughout the domain
    for (int i = 0; i < nDim; i++) {
      momSrc[i] = USource_[i][0];
    }
  } else {
    // Linearly interpolate source term within the planes, maintain constant
    // source term above and below the heights provided
    for (int i = 0; i < nDim; i++) {
      utils::linear_interp(velHeights_, USource_[i], zp, momSrc[i]);
    }
  }
}

void
ABLForcingAlgorithm::eval_temperature_source(const double zp, double& tempSrc)
{
  if (tempHeights_.size() == 1) {
    tempSrc = TSource_[0];
  } else {
    utils::linear_interp(tempHeights_, TSource_, zp, tempSrc);
  }
}

} // namespace nalu
} // namespace sierra
