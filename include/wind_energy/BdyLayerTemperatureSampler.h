/*------------------------------------------------------------------------*/
/*  Copyright 2014 National Renewable Energy Laboratory.                  */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#ifndef BDYLAYERTEMPERATURESAMPLER_H
#define BDYLAYERTEMPERATURESAMPLER_H

#include "AlgorithmDriver.h"

#include "stk_mesh/base/Entity.hpp"
#include "stk_mesh/base/Ghosting.hpp"
#include "stk_search/BoundingBox.hpp"
#include "stk_search/IdentProc.hpp"
#include "stk_search/SearchMethod.hpp"

#include "yaml-cpp/yaml.h"

#include <vector>
#include <unordered_map>
#include <memory>
#include <cstdint>

namespace stk{
namespace mesh {
class Part;
typedef std::vector<Part*> PartVector;
}
}

namespace sierra {
namespace nalu {

typedef stk::search::IdentProc<uint64_t, int> SearchKey;
typedef stk::search::Point<double> Point;
typedef stk::search::Sphere<double> Sphere;
typedef stk::search::Box<double> Box;
typedef std::pair<Point, SearchKey> BoundingPoint;
typedef std::pair<Sphere, SearchKey> BoundingSphere;
typedef std::pair<Box, SearchKey> BoundingBox;

class Realm;
class SolverAlgorithm;
struct WallUserData;

struct TargetPointInfo
{
  std::vector<double> nodalCoords_;
  std::vector<double> isoParCoords_;
  stk::mesh::Entity elem_;
  double minDistance_;

  TargetPointInfo(const int nDim)
    : nodalCoords_(nDim),
      isoParCoords_(nDim)
  {}
};

/** Sample velocity at a given height for use with LES wall models
 *
 */
class BdyLayerTemperatureSampler: public AlgorithmDriver
{
public:
  BdyLayerTemperatureSampler(
      Realm&,
      WallUserData&);

  virtual ~BdyLayerTemperatureSampler();

  virtual void pre_work();

  void set_wall_func_algorithm(SolverAlgorithm* alg)
  { wallFuncAlg_ = alg; }

  /** Calculate interpolated velocity at the desired sampling point
   */
  void get_temperature(
      const stk::mesh::Entity,
      double*);

protected:
  void determine_node_elem_mapping();

  /** Determine point field where the enclosing element must be determined.
   *
   *  This method processes the wall sidesets, applies the user-provided offset
   *  vector to determine the spatial location where the velocity must be
   *  sampled. The point field is converted into bounding spheres required by
   *  the STK::Search interface.
   */
  void determine_bounding_spheres();

  /** Populate a vector of bounding boxes around elements on a meshblock
   *
   */
  void determine_bounding_boxes();

  /** Perform search and determine the (point, element) map used to interpolate
   * velocity field upon query.
   */
  void search_and_update_map();

  /** Determine the optimal element used for interpolating velocity field
   */
  void finalize_search(std::vector<std::pair<SearchKey, SearchKey>>&);

private:
  BdyLayerTemperatureSampler() = delete;
  BdyLayerTemperatureSampler(const BdyLayerTemperatureSampler&) = delete;

  SolverAlgorithm* wallFuncAlg_{nullptr};

  //! List of mesh blocks where the target element search is performed
  std::vector<std::string> searchPartNames_;

  //! The normal offset to be applied to determine the height above ground where
  //! the velocity is sampled for the LES wall model
  std::vector<double> tempOffsetVector_{0.0, 0.0, 0.0};

  //! List of (element, mpi_rank) mapping for ghosting elements to other MPI
  //! ranks
  std::vector<stk::mesh::EntityProc> elemsToGhost_;

  //! Custom ghosting object that handles ghosting the elements to the desired
  //! MPI rank
  stk::mesh::Ghosting* samplerElemGhosting_{nullptr};

  //! Target nodes whose bounding element must be determined
  std::vector<BoundingSphere> boundingPoints_;

  //! List of candidate element bounding boxes where search is performed
  std::vector<BoundingBox> boundingBoxes_;

  //! Mapping of a sideset node to the element that is used to determine the
  //! velocity used in the LES model.
  std::map<stk::mesh::EntityId, std::unique_ptr<TargetPointInfo>> targetElemMap_;

  stk::search::SearchMethod searchMethod_{stk::search::KDTREE};

  bool doInit_{true};
};

}  // nalu
}  // sierra

#endif /* BDYLAYERTEMPERATURESAMPLER_H */
