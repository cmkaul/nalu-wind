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
  Realm& realm, ABLForcingAlgorithm* abldamp)
  : SupplementalAlgorithm(realm),
    ablDamp_(abldamp),
    nDim_(realm.meta_data().spatial_dimension())
{
  // Save some fields
  // ! CK : Left off here
  ScalarFieldType* density = realm.meta_data().get_field<ScalarFieldType>(
    stk::topology::NODE_RANK, "density");
  density_ = &(density->field_of_state(stk::mesh::StateNP1));



  coords_(realm.meta_data().get_field<VectorFieldType>(
      stk::topology::NODE_RANK, "coordinates")),
    dualNodalVolume_(realm.meta_data().get_field<ScalarFieldType>(
      stk::topology::NODE_RANK, "dual_nodal_volume")),
    density_(NULL),
    specificHeat_(realm.meta_data().get_field<ScalarFieldType>(
      stk::topology::NODE_RANK, "specific_heat")),
}

void
EnthalpyABLSrcNodeSuppAlg::node_execute(
  double* lhs, double* rhs, stk::mesh::Entity node)
{
  const double dualVol = *stk::mesh::field_data(*dualNodalVolume_, node);
  const double* pt = stk::mesh::field_data(*coords_, node);
  const double rhoNP1 = *stk::mesh::field_data(*density_, node);
  const double spHeat = *stk::mesh::field_data(*specificHeat_, node);
  double tempSrc;

  ablSrc_->eval_temperature_source(pt[nDim_ - 1], tempSrc);

  rhs[0] += dualVol * rhoNP1 * spHeat * tempSrc;
  lhs[0] += 0.0;
}

} // nalu
} // sierra
