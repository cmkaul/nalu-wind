
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
ABLDampingAlgorithm::initialize()
{ 
  auto* bdyLayerStats = realm_.bdyLayerStats_;
  //! CK not sure about this! (right way to get heights_ array from bdyLayerStats?
  // Is resulting heights_ the values or a reference? 
  auto heights_ = bdyLayerStats->heights_;
  const size_t nHeights= heights_.size();
   
  if (momSrcType_ != OFF) {
    NaluEnv::self().naluOutputP0()
      << "ABL Damping active for Momentum Equations\n"<< std::endl;
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
      << "ABL Damping active for Temperature Equation\n"
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
ABLDampingAlgorithm::execute()
{
  if (momentumForcingOn())
    compute_momentum_sources();

  if (temperatureForcingOn())
    compute_temperature_sources();
}

void
ABLDampingAlgorithm::compute_momentum_sources()
{
  const double dt = realm_.get_time_step();
  const int nDim = realm_.spatialDimension_;
  const size_t nHeights = heights_.size();
  auto& meta = real.meta_data();
  auto& bulk = realm.bulk_data();
  stk::mesh::Selector sel = meta.locally_owned_part()
    & stk::mesh::selectUnion(fluidParts)
    & !(realm_.get_inactive_selector());
  const auto bkts = bulk.get_buckets(stk::topology::NORE_RANK, sel);  

  
  VectorFieldType* velocity = meta.get_field<VectorFieldType>(
    stk::topology::NODE_RANK, "velocity")
  


  // CK: At various points, it is assumed that bdyLayerStats have been initialized  
  auto* bdyLayerStats = realm_.bdyLayerStats_;
  // Gather the mean profile
  for (size_t ih=0; ih < nHeights; ih++) {
    //CK: not 100% sure this if is a good idea
    if (heights_[ih] > heights_[nHeights-1] - momThickess){
      bdyLayerStats->velocity(heights_[ih], UmeanCalc_[ih].data());
      bdyLayerStats->density(heights_[ih], &rhoMeanCalc_[ih]);
    }
  }
  // CK Get the local source terms--note that the Usource array might not be useful!
  // CK Should I just add the  source terms on here???
  for (auto b: bkts) {
    for (size_t in=0; in < b->size(); in++){
      auto node = (*b)[in];
      int ih = *stk::mesh::field_data(*heightIndex_, node);
      double* vel = stk::mesh::field_data(*velocity, node);
      if (heights_[ih] > heights_[nHeights-1] - momThickess){
        for (int d=0; d < nDim_; d++){
           momSrc[d] = rhoMeanCalc_[ih] * (alphaMomentum_/dt) *
                         (vel[d] -UmeanCalc_[ih][d]);
         }
      }
      else{
      for (int d=0; d < nDim_; d++){
            momSrc[d] =0.0;
       }
     }
    }
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
  const int ih, std, std::vector<double>& momSrc)
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
