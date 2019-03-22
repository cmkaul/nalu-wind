/******************************************************************************/
/* This software is released under the license detailed in the file, LICENSE, */
/* located in the top-level Nalu directory structure.                         */
/******************************************************************************/

#ifndef ENTHALPYABLDAMPSRCNODESUPPALG_H
#define ENTHALPYABLDAMPSRCNODESUPPALG_H

#include "SupplementalAlgorithm.h"
#include "FieldTypeDef.h"

#include <stk_mesh/base/Entity.hpp>

namespace sierra {
namespace nalu {

// Forward declarations
class Realm;
class ABLDampingAlgorithm;

class EnthalpyABLDampSrcNodeSuppAlg : public SupplementalAlgorithm
{
public:
  EnthalpyABLDampSrcNodeSuppAlg(Realm&, ABLDampingAlgorithm*);

  virtual ~EnthalpyABLDampSrcNodeSuppAlg() {}

  virtual void setup() {}

  virtual void node_execute(double*, double*, stk::mesh::Entity);

private:
  EnthalpyABLDampSrcNodeSuppAlg();
  EnthalpyABLDampSrcNodeSuppAlg(const EnthalpyABLDampSrcNodeSuppAlg&);

  //! Pointer to ABL Damping algorithm object
  ABLDampingAlgorithm* ablDamp_;

  //! Pointer to mesh coordinates
  VectorFieldType* coords_;

  //! Pointer to the density
  ScalarFieldType* densityNP1_

  //! Pointer to dual volume of the mesh
  ScalarFieldType* dualNodalVolume_;

  //! Pointer to the density of the mesh
  ScalarFieldType* density_;

  //! Pointer to specific heat
  ScalarFieldType* specificHeat_;

  //! Pointer to the temperature
  ScalarFieldType* temperatureNP1

  //! Pointer to the ABL height index field
  ScalarIntFieldType* heightIndex_

  const int nDim_;
};

} // namespace nalu
} // namespace sierra

#endif /* ENTHALPYABLDAMPSRCNODESUPPALG_H */
