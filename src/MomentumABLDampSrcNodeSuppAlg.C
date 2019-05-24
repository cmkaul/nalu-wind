
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
    ablDamp_(abldamp),
    nDim_(realm.meta_data().spatial_dimension())
{
  // Save some fields
  // get the realm meta data
  stk::mesh::MetaData& meta = realm.meta_data();
  // Get the coordinates 
  coords_ = meta.get_field<VectorFieldType>(
    stk::topology::NODE_RANK, "coordinates");
  // Cet the density
  ScalarFieldType* density_ = meta.get_field<ScalarFieldType>(
    stk::topology::NODE_RANK, "density");
   // Set access to the correct update state
  densityNP1_ = &(density_->field_of_state(stk::mesh::StateNP1));
  // Get the velocity field
  VectorFieldType* velocity_ = meta.get_field<VectorFieldType>(
    stk::topology::NODE_RANK, "velocity");
  // Set access to the correct  update state
  velocityNP1_ = &(velocity_->field_of_state(stk::mesh::StateNP1));
  // Get the node volume for weighting the source term 
  dualNodalVolume_ = meta.get_field<ScalarFieldType>(
    stk::topology::NODE_RANK, "dual_nodal_volume");
  /* Get the index in the unique height array
     created in bdyLayerStats->setup */
  //CK: TODO: I can't figure out when that is executed relative to this!
  heightIndex_ = meta.get_field<ScalarIntFieldType>(
    stk::topology::NODE_RANK, "bdy_layer_height_index_field");

  
  }

void
MomentumABLDampSrcNodeSuppAlg::node_execute(
  double*  /* lhs */, double* rhs, stk::mesh::Entity node)
{
  //Access nodal values
  const double* pt = stk::mesh::field_data(*coords_, node);
  const double dualVol = *stk::mesh::field_data(*dualNodalVolume_, node);
  const double rhoNP1 = *stk::mesh::field_data(*densityNP1_, node);
  const double* vel = stk::mesh::field_data(*velocityNP1_, node);
  const int ih = *stk::mesh::field_data(*heightIndex_, node);
  
  //CK: tentative here
  //Getting some values needed from the damping algorithm
  const double dampHeight = ablDamp_->minDampingHeightMomentum;
  const double dampCoeff = ablDamp_->dampingCoeffMomentum[ih];
  const double* dampVel = ablDamp_->UDamp[ih].data();

  for (int i = 0; i < 2; i++) {
    ThrowAssertMsg(std::isfinite(rhs[i]), "Inf or NAN rhs before damp");
    rhs[i] += dualVol * rhoNP1* dampCoeff * (dampVel[i]-vel[i]);
    ThrowAssertMsg(std::isfinite(rhs[i]), "Inf or NAN rhs after damp");
      
  }
  ThrowAssertMsg(std::isfinite(rhs[2]), "Inf or NAN rhs before damp");
  rhs[2] -= dualVol * rhoNP1 * (dampVel[2]);
  ThrowAssertMsg(std::isfinite(rhs[2]), "Inf or NAN rhs after damp");

  // If below the minimum damping heihght, return without doing anything
  /*
  if (pt[nDim_-1] < dampHeight){

  }else{  
    // Add the momentum source into the RHS
    for (int i = 0; i < nDim_; i++) {
      ThrowAssertMsg(std::isfinite(rhs[i]), "Inf or NAN rhs before damp");
      
      rhs[i] += dualVol * rhoNP1* dampCoeff * (dampVel[i]-vel[i]);
      ThrowAssertMsg(std::isfinite(rhs[i]), "Inf or NAN rhs after damp");
      
    }
    
  }
  */


}

} // namespace nalu
} // namespace sierra
