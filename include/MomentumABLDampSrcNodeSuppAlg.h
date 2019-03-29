/******************************************************************************/
/* This software is released under the license detailed in the file, LICENSE, */
/* located in the top-level Nalu directory structure.                         */
/******************************************************************************/

#ifndef MOMENTUMABLDAMPSRCNODESUPPALG_H
#define MOMENTUMABLDAMPSRCNODESUPPALG_H

#include <SupplementalAlgorithm.h>
#include <FieldTypeDef.h>

#include <stk_mesh/base/Entity.hpp>

namespace sierra {
namespace nalu {

class Realm;
class ABLDampingAlgorithm;

class MomentumABLDampSrcNodeSuppAlg : public SupplementalAlgorithm
{
public:
  MomentumABLDampSrcNodeSuppAlg(Realm&, ABLDampingAlgorithm*);

  virtual ~MomentumABLDampSrcNodeSuppAlg() {}

  virtual void setup() {}

  virtual void node_execute(double*, double*, stk::mesh::Entity);

private:
  MomentumABLDampSrcNodeSuppAlg();
  MomentumABLDampSrcNodeSuppAlg(const MomentumABLDampSrcNodeSuppAlg&);

  //! Pointer to ABL Forcing Algorithm object
  ABLDampingAlgorithm* ablDamp_;

  //! Pointer to the mesh coordinates
  VectorFieldType* coords_;

  //! Pointer to the density at state np1
  ScalarFieldType* densityNP1_;
  
  //! Pointer to the velocity at state np1
  VectorFieldType* velocityNP1_;

  //! Pointer to the dual volume of the mesh
  ScalarFieldType* dualNodalVolume_;
  
  //! Pointer to the height index array created in bdyLayerStats
  ScalarIntFieldType* heightIndex_;

  //! Spatial dimension of the computational mesh
  const int nDim_;




  
};
}
}

#endif /* MOMENTUMABLDAMPSRCNODESUPPALG_H */
