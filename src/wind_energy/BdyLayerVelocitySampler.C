/*------------------------------------------------------------------------*/
/*  Copyright 2014 National Renewable Energy Laboratory.                  */
/*  This software is released under the license detailed                  */
/*  in the file, LICENSE, which is located in the top-level Nalu          */
/*  directory structure                                                   */
/*------------------------------------------------------------------------*/

#include "wind_energy/BdyLayerVelocitySampler.h"
#include "Realm.h"
#include "NaluParsing.h"
#include "SolverAlgorithm.h"
#include "master_element/MasterElement.h"
#include "master_element/MasterElementFactory.h"

#include "stk_mesh/base/MetaData.hpp"
#include "stk_mesh/base/BulkData.hpp"
#include "stk_mesh/base/Part.hpp"
#include "stk_mesh/base/Field.hpp"
#include "stk_search/CoarseSearch.hpp"
#include "stk_search/IdentProc.hpp"

#include <limits>

namespace sierra {
namespace nalu {

BdyLayerVelocitySampler::BdyLayerVelocitySampler(
  Realm& realm,
  WallUserData& wallUserData
) : AlgorithmDriver(realm),
    searchPartNames_(wallUserData.ablTargetPartNames_),
    offsetVector_(wallUserData.velOffsetVector_)
{
  if (offsetVector_.size() != realm_.meta_data().spatial_dimension())
    throw std::runtime_error(
      "BdyLayerVelocitySampler:: Invalid offset vector provided");
}

BdyLayerVelocitySampler::~BdyLayerVelocitySampler()
{}

void
BdyLayerVelocitySampler::pre_work()
{
  if (doInit_) {
    determine_node_elem_mapping();
    doInit_ = false;
  }
}

void
BdyLayerVelocitySampler::get_velocity(
  const stk::mesh::Entity node, double* velocity)
{
  auto& meta = realm_.meta_data();
  auto& bulk = realm_.bulk_data();
  const int nDim = meta.spatial_dimension();

  const VectorFieldType* vfield =
    meta.get_field<VectorFieldType>(stk::topology::NODE_RANK, "velocity");
  const stk::mesh::EntityId nodeID = bulk.identifier(node);

  auto it = targetElemMap_.find(nodeID);
  if (it == targetElemMap_.end())
    throw std::runtime_error(
      "BdyLayerVelocitySampler:: Invalid node requested");

  auto elem = it->second->elem_;
  auto& isoParCoords = it->second->isoParCoords_;

  const stk::topology topo = bulk.bucket(elem).topology();
  const int num_nodes = bulk.num_nodes(elem);
  const auto* elemNodes = bulk.begin_nodes(elem);
  std::vector<double> elemVelocity(num_nodes * nDim);

  for (int in = 0; in < num_nodes; in++) {
    const auto enode = elemNodes[in];
    const double* vel = stk::mesh::field_data(*vfield, enode);

    for (int d = 0; d < nDim; d++) {
      const int offset = d * num_nodes + in;
      elemVelocity[offset] = vel[d];
    }
  }

  MasterElement* meSCS =
    sierra::nalu::MasterElementRepo::get_surface_master_element(topo);
  meSCS->interpolatePoint(
    nDim, isoParCoords.data(), elemVelocity.data(), velocity);
}

void
BdyLayerVelocitySampler::determine_node_elem_mapping()
{
  auto& bulk = realm_.bulk_data();
  // Reset data structures befor processing
  boundingPoints_.clear();
  boundingBoxes_.clear();
  elemsToGhost_.clear();
  targetElemMap_.clear();

  bulk.modification_begin();
  {
    if (samplerElemGhosting_ == nullptr) {
      std::string ghostingName = "nalu_abl_les_model_ghosting";
      samplerElemGhosting_ = &(bulk.create_ghosting(ghostingName));
    }
    else {
      bulk.destroy_ghosting(*samplerElemGhosting_);
    }
  }
  bulk.modification_end();

  determine_bounding_spheres();

  determine_bounding_boxes();

  search_and_update_map();
}

void
BdyLayerVelocitySampler::determine_bounding_spheres()
{
  auto& meta = realm_.meta_data();
  auto& bulk = realm_.bulk_data();
  const int iproc = bulk.parallel_rank();
  const int nDim = meta.spatial_dimension();

  stk::mesh::Selector sel = stk::mesh::selectUnion(wallFuncAlg_->partVec_);
  const auto& bkts = bulk.get_buckets(
    stk::topology::NODE_RANK, sel);
  const VectorFieldType* coords = meta.get_field<VectorFieldType>(
    stk::topology::NODE_RANK, realm_.get_coordinates_name());

  for (auto b: bkts) {
    for (size_t in = 0; in < b->size(); in++) {
      Point searchPt;
      auto node = (*b)[in];
      auto nodeID = bulk.identifier(node);
      double* crd = stk::mesh::field_data(*coords, node);

      // Apply offset to the wall node to determine the spatial location where
      // velocity is sampled.
      for (int d=0; d < nDim; d++) {
        searchPt[d] = crd[d] + offsetVector_[d];
      }

      // Create STK search data structures and add this point to the search vector
      SearchKey key(nodeID, iproc);
      BoundingSphere bSphere(Sphere(searchPt, 1.0), key);
      boundingPoints_.push_back(bSphere);
    }
  }
}

void
BdyLayerVelocitySampler::determine_bounding_boxes()
{
  auto& meta = realm_.meta_data();
  auto& bulk = realm_.bulk_data();
  const int iproc = bulk.parallel_rank();
  const int nDim = meta.spatial_dimension();

  stk::mesh::PartVector searchParts;
  for (auto pname : searchPartNames_) {
    stk::mesh::Part* part = meta.get_part(pname);
    if (part == nullptr)
      throw std::runtime_error(
        "BdyLayerVelocitySampler:: Cannot find required part " + pname);
    searchParts.push_back(part);
  }

  stk::mesh::Selector sel =
    meta.locally_owned_part() & stk::mesh::selectUnion(searchParts);
  const auto& bkts = bulk.get_buckets(stk::topology::ELEMENT_RANK, sel);
  const VectorFieldType* coords = meta.get_field<VectorFieldType>(
    stk::topology::NODE_RANK, realm_.get_coordinates_name());

  Point minPt, maxPt;
  for (auto b: bkts) {
    for (size_t ie =0; ie < b->size(); ie++) {
      auto elem = (*b)[ie];
      const int num_nodes = bulk.num_nodes(elem);

      // Reset min max points before processing this element
      for (int d=0; d < nDim; d++) {
        minPt[d] = std::numeric_limits<double>::max();
        maxPt[d] = std::numeric_limits<double>::lowest();
      }

      const stk::mesh::Entity* nodes = bulk.begin_nodes(elem);
      for (int in=0; in < num_nodes; in++) {
        auto node = nodes[in];

        const double* crd = stk::mesh::field_data(*coords, node);

        for (int d=0; d<nDim; d++) {
          minPt[d] = std::min(minPt[d], crd[d]);
          maxPt[d] = std::max(maxPt[d], crd[d]);
        }
      }

      SearchKey key(bulk.identifier(elem), iproc);
      Box box(minPt, maxPt);
      BoundingBox bBox(box, key);
      boundingBoxes_.push_back(bBox);
    }
  }
}

void
BdyLayerVelocitySampler::search_and_update_map()
{
  auto& bulk = realm_.bulk_data();
  std::vector<std::pair<SearchKey, SearchKey>> searchKeyPair;
  const int iproc = bulk.parallel_rank();

  // Use STK search to determine the point-element pair. This coarse search
  // might return multiple candidate elements because we are using bounding
  // boxes around the element. The actual enclosing element and the unique
  // (point, element) pair will be determined in the finalize_search method.
  stk::search::coarse_search(boundingPoints_, boundingBoxes_, searchMethod_,
                             bulk.parallel(), searchKeyPair);

  // Ensure that the elements are ghosted to the appropriate MPI ranks before
  // proceeding with finalize_search.
  for (auto it: searchKeyPair) {
    const stk::mesh::EntityId elemID = it.second.id();
    const int ptProc = it.first.proc();
    const int elemProc = it.second.proc();

    if ((elemProc == iproc) && (ptProc != iproc)) {
      // The bounding element resides on different MPI rank, needs ghosting
      auto elem = bulk.get_entity(stk::topology::ELEMENT_RANK, elemID);

      if (!bulk.is_valid(elem))
        throw std::runtime_error(
          "BdyLayerVelocitySampler:: Invalid element encountered during search");

      stk::mesh::EntityProc elemProc(elem, ptProc);
      elemsToGhost_.push_back(elemProc);
    }
  }

  int numGhostElems = elemsToGhost_.size();
  int gNumGhostElems = 0;
  stk::all_reduce_sum(bulk.parallel(), &numGhostElems, &gNumGhostElems, 1);

  if (gNumGhostElems > 0) {
    bulk.modification_begin();
    bulk.change_ghosting(*samplerElemGhosting_, elemsToGhost_);
    bulk.modification_end();
  }

  // Determine unique (point, enclosing element) pairs
  finalize_search(searchKeyPair);
}

void
BdyLayerVelocitySampler::finalize_search(
  std::vector<std::pair<SearchKey, SearchKey>>& searchKeyPair)
{
  auto& meta = realm_.meta_data();
  auto& bulk = realm_.bulk_data();
  const int iproc = bulk.parallel_rank();
  const int nDim = meta.spatial_dimension();
  const VectorFieldType* coords = meta.get_field<VectorFieldType>(
    stk::topology::NODE_RANK, realm_.get_coordinates_name());

  std::vector<double> isoParCoords(nDim);
  std::vector<double> nodalCoords(nDim);
  std::vector<double> elemCoords;

  for (auto it: searchKeyPair) {
    const int ptProc = it.first.proc();
    // Skip nodes not owned by this proc
    if (iproc != ptProc) continue;

    const stk::mesh::EntityId nodeID = it.first.id();
    const stk::mesh::EntityId elemID = it.second.id();
    const auto elem = bulk.get_entity(stk::topology::ELEM_RANK, elemID);
    const stk::topology topo = bulk.bucket(elem).topology();
    const auto node = bulk.get_entity(stk::topology::NODE_RANK, nodeID);

    // Check if this node was already processed with another element before, if
    // not create new data structure to hold the point, element pair.
    double curDist = std::numeric_limits<double>::max();
    auto ptPtr = targetElemMap_.find(nodeID);
    TargetPointInfo* tgtInfo{nullptr};
    if (ptPtr == targetElemMap_.end()) {
      targetElemMap_[nodeID] = std::unique_ptr<TargetPointInfo>(
        new TargetPointInfo(nDim));
      tgtInfo = targetElemMap_[nodeID].get();
    } else {
      tgtInfo = targetElemMap_[nodeID].get();

      // Reset current minimum distance with the one stored in the data structure
      curDist = tgtInfo->minDistance_;
    }

    // Get the coordinates of the point where the velocity is to be sampled.
    const double* nodeCrd = stk::mesh::field_data(*coords, node);
    for (int d=0; d < nDim; d++) {
      nodalCoords[d] = nodeCrd[d] + offsetVector_[d];
    }

    // Populate the coordinates of the element nodes
    const int num_nodes = bulk.num_nodes(elem);
    const auto* elemNodes = bulk.begin_nodes(elem);
    elemCoords.resize(num_nodes * nDim);
    for (int in=0; in < num_nodes; in++) {
      const auto enode = elemNodes[in];
      const double* crd = stk::mesh::field_data(*coords, enode);
      for (int d=0; d < nDim; d++) {
        const int offset = d * num_nodes + in;
        elemCoords[offset] = crd[d];
      }
    }

    MasterElement* meSCS = sierra::nalu::MasterElementRepo::get_surface_master_element(topo);
    const double newDistance = meSCS->isInElement(
      elemCoords.data(), nodalCoords.data(), isoParCoords.data());

    // If the distance is less than the previous candidate update the mapping
    // and replace with new element.
    if (newDistance < curDist) {
      tgtInfo->elem_ = elem;
      tgtInfo->minDistance_ = newDistance;

      for (int d=0; d < nDim; d++) {
        tgtInfo->nodalCoords_[d] = nodalCoords[d];
        tgtInfo->isoParCoords_[d] = isoParCoords[d];
      }
    }
  }

  if (targetElemMap_.size() != boundingPoints_.size())
    throw std::runtime_error("BdyLayerVelocitySampler:: Cannot find bounding elements for all wall nodes.");
}

}  // nalu
}  // sierra
