/******************************************************************************/
/* This software is released under the license detailed in the file, LICENSE, */
/* located in the top-level Nalu directory structure.                         */
/******************************************************************************/

#include "EnthalpyABLDampSrcNodeSuppAlg.h"
#include "Realm.h"
#include "wind_energy/ABLDampingAlgorithm.h"

// stk_mesh/base/fem
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Field.hpp>

namespace sierra {
namespace nalu {

EnthalpyABLDampSrcNodeSuppAlg::EnthalpyABLDampSrcNodeSuppAlg(
  Realm& realm, ABLDampingAlgorithm* abldamp)
  : SupplementalAlgorithm(realm),
    ablDamp_(abldamp),
    nDim_(realm.meta_data().spatial_dimension())
{
  // Save some fields
  // ! get the realm meta data
  stk::mesh::MetaData& meta = realm.meta_data();
  //! Get the coordinates  
  VectorFieldType* coords_ = meta.get_field<VectorFieldType>(
    stk::topology::NODE_RANK, "coordinates");
  //! Get the density
  ScalarFieldType* density = realm.meta_data().get_field<ScalarFieldType>(
    stk::topology::NODE_RANK, "density");
  //! Set access to the correct update state
  densityNP1_ = &(density->field_of_state(stk::mesh::StateNP1));
  // Get the node volume for weighting the source term 
  dualNodalVolume_ = meta.get_field<ScalarFieldType>(
    stk::topology::NODE_RANK, "dual_nodal_volume");

  // Get the temperature field
  ScalarFieldType* temperature_ = meta.get_field<ScalarFieldType>(
    stk::topology::NODE_RANK, "dual_nodal_volume");
  //set acces to the correct update state
  temperatureNP1_ = &(temperature_->field_of_state(stk::mesh::StateNP1));
  // Get the specific heat
  specificHeat_ = meta.get_field<ScalarFieldType>(
    stk::topology::NODE_RANK, "specific_heat");
  /* Get the index in the unique height array
     created in bdyLayerStats->setup */
  //CK: TODO: I can't figure out when that is executed relative to this!
  heightIndex_ = meta.get_field<ScalarIntFieldType>(
    stk::topology::NODE_RANK, "bdy_layer_height_index_field");
}

void
EnthalpyABLDampSrcNodeSuppAlg::node_execute(
  double* lhs, double* rhs, stk::mesh::Entity node)
{
  const double* pt = stk::mesh::field_data(*coords_, node);
  const double dualVol = *stk::mesh::field_data(*dualNodalVolume_, node);
  const double rhoNP1 = *stk::mesh::field_data(*densityNP1_, node);
  const double temp = *stk::mesh::field_data(*temperatureNP1_, node);
  const int ih = stk::mesh::field_data(*heightIndex_, node);
  const double spHeat = *stk::mesh::field_data(*specificHeat_, node);
  
  //CK: tentative here
  //Getting some values needed from the damping algorithm
  const double dampHeight = ablDamp_->minDampingHeightTemperature
  const double dampCoeff = ablDamp_->dampingCoeffTemperature[ih]
  const double dampTemp = ablDamp_->TDamp[ih]
  // If below the minimum damping heihght, return without doing anything
  if (pt[nDim_-1] < dampHeight) {


  }else{  
    rhs[0] += dualVol * rhoNP1 * spHeat * dampCoeff * (dampTemp-temp);
    lhs[0] += 0.0;

  }


}

} // nalu
} // sierra
