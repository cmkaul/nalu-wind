/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#ifndef AssembleFaceScalarFluxBCSolverAlgorithm_h
#define AssembleFaceScalarFluxBCSolverAlgorithm_h

#include<SolverAlgorithm.h>
#include<FieldTypeDef.h>

namespace stk {
namespace mesh {
class Part;
}
}

namespace sierra{
namespace nalu{

class LinearSystem;
class Realm;

class AssembleFaceScalarFluxBCSolverAlgorithm : public SolverAlgorithm
{
public:

  AssembleFaceScalarFluxBCSolverAlgorithm(
    Realm &realm,
    stk::mesh::Part *part,
    EquationSystem *eqSystem,
    GenericFieldType *bcScalarQ);
  virtual ~AssembleFaceScalarFluxBCSolverAlgorithm() {}
  virtual void initialize_connectivity();
  virtual void execute();

private:

  GenericFieldType *bcScalarQ_;
  GenericFieldType *exposedAreaVec_;
};

}
}


#endif /* ASSEMBLESCALARELEMDIFFBCSOLVERALGORITHM_H_ */
