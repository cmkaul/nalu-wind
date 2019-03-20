
#include "MomentumABLDampSrcNodeSuppAlg.h"
#include "Realm.h"
#include "wind_energy/ABLDampingAlgorithm.h"

// stk_mesh/base/fem
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Field.hpp>

namespace sierra {
namespace nalu {

MomentumABLDampSrcNodeSuppAlg::MomentumABLDampSrcNodeSuppAlg(
  Realm& realm, ABLDampingAlgorithm* abldamp)
  : SupplementalAlgorithm(realm),
    ablSrc_(ablsrc),
    nDim_(realm_.meta_data().spatial_dimension())
{
  // save off fields
  stk::mesh::MetaData & meta = realm_.meta_data();
  VectorFieldType* coords_ = meta.get_field<VectorFieldType>(
    stk::topology::NODE_RANK, "coordinates");
  ScalarFieldType* density = meta.get_field<ScalarFieldType>(
    stk::topology::NODE_RANK, "density");
  densityNP1_ = &(density->field_of_state(stk::mesh::StateNP1));
  VectorFieldType* velocity_ = meta.get_field<VectorFieldType>(
    stk::topology::NODE_RANK, "velocity");
  velocityNP1_ = &(velocity__>field_of_state(stk::mesh::StateNP1));
  ScalarFieldType* dualNodalVolume_ = meta.get_field<ScalarFieldType>(
    stk::topology::NODE_RANK, "dual_nodal_volume");
  ScalarIntFieldType* heightIndex_ = meta.get_field<ScalarIntFieldType>(
    stk::topology::NODE_RANK, "bdy_layer_height_index_field")
  const int nDim_ = meta.spatial_dimension(); 
  
  }

void
MomentumABLDampSrcNodeSuppAlg::node_execute(
  double*  /* lhs */, double* rhs, stk::mesh::Entity node)
{
  const double dualVol = *stk::mesh::field_data(*dualNodalVolume_, node);
  const double* pt = stk::mesh::field_data(*coords_, node);
  const double rhoNP1 = *stk::mesh::field_data(*densityNP1_, node);
  const double* vel = *stk::mesh::field_data(*velocityNP1_, node);
  const int ih = stk::mesh::field_data(*heightIndex_, node);
  std::vector<double> momSrc(nDim_);
  
  //CK: tenteative here
  const double dampCoeff = abldamp->alphaMomentum[ih]
  // CK: TODO: Ensure that mean velocity is up to date

  const double* velMean = abldamp->UmeanCalc_[ih].data()

 
  ///// Evaluate the momentum source
  ////ablSrc_->eval_momentum_source(ih, vel, momSrc);

  //CK Add a check for the height being high enough or ensure by dampCoef??
  // Add the momentum source into the RHS
  for (int i = 0; i < nDim_; i++) {
    // CK: IMPORTANT! TODO: Check on density being used
    // TODO: Check on sign
    rhs[i] += dualVol * rhoNP1* dampCoeff * (velMean[d]-vel[d]);
  }
}

} // namespace nalu
} // namespace sierra
