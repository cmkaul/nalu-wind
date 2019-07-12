/*------------------------------------------------------------------------*/
/*  Copyright 2014 Sandia Corporation.                                    */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/


#ifndef ComputeABLWallFluxesAlgorithm_h
#define ComputeABLWallFluxesAlgorithm_h

#include <Algorithm.h>
#include <NaluParsing.h>
#include <FieldTypeDef.h>
#include <ABLProfileFunction.h>

// stk
#include <stk_mesh/base/Part.hpp>

#include <vector>
#include <iostream>

namespace sierra{
namespace nalu{

class Realm;

class ComputeABLWallFluxesAlgorithm : public Algorithm
{
public:
  template <typename T>
  using ListArray = std::vector<std::vector<T>>;

  // Constructor
  ComputeABLWallFluxesAlgorithm(
    Realm &realm,
    const YAML::Node&,
    stk::mesh::Part *part,
    const bool &useShifted,
    const double &gravity,
    const double &z0,
    const double &Tref);

  // Destructor
  virtual ~ComputeABLWallFluxesAlgorithm();

  // Parse input file for user options and initialize
  void load(const YAML::Node&);

  void execute();

  void zero_nodal_fields();

  void compute_utau(
      const double &up, const double &zp,
      const double &qsurf, const ABLProfileFunction *ABLProfFun,
      double &utau);

  //! Compute the surface fluxes.
  void compute_fluxes_given_surface_temperature(
      const double tol,
      const double &up, 
      const double &Tp,
      const double Tsurface,
      const double rho,
      const double Cp,
      const double &zp,
      const ABLProfileFunction *ABLProfFun,
      double &utau, 
      double &qsurf);

  void compute_fluxes_given_surface_heating(
      const double tol,
      const double &up,
      const double &Tp,
      double &Tsurface,
      const double rho,
      const double Cp,
      const double &zp,
      const ABLProfileFunction *ABLProfFun,
      double &utau, 
      const double &qsurf);
  
  void normalize_nodal_fields();

  const bool useShifted_;
  const double z0_;
  const double Tref_;
  const double gravity_;
  const double alpha_h_;
  const double beta_m_;
  const double beta_h_;
  const double gamma_m_;
  const double gamma_h_;
  const double kappa_;
  const int maxIteration_;
  const double tolerance_;

  bool useLESSamplingHeight_{false};
  double lesModelRefHeight_{0.0};
  double lesModelRefHeightTemp_{0.0};

  VectorFieldType *velocity_;
  ScalarFieldType *temperature_;
  VectorFieldType *bcVelocity_;
  ScalarFieldType *bcHeatFlux_;
  ScalarFieldType *density_;
  ScalarFieldType *specificHeat_;
  VectorFieldType *coordinates_;
  GenericFieldType *exposedAreaVec_;
  GenericFieldType *wallFrictionVelocityBip_;
  GenericFieldType *wallHeatFluxBip_;
  GenericFieldType *wallNormalDistanceBip_;
  ScalarFieldType *assembledWallNormalDistance_;
  ScalarFieldType *assembledWallArea_;

private:
  // Break the flux/surface temperature vs. time input table into vectors
  // of each quantity and store in the following vectors.
  std::vector<double> tableTimes_;
  std::vector<double> tableFluxes_;
  std::vector<double> tableSurfaceTemperatures_;
  std::vector<double> tableWeights_;
};

} // namespace nalu
} // namespace Sierra

#endif
