#ifdef __GPUCOMPILE__
#include "Garfield/GPUInterface.hh"
#include "TetrahedralTreeGPU.h"
#include "GPUFunctions.h"
#else
#include "Garfield/TetrahedralTree.hh"
#include <iostream>
#endif

namespace Garfield {

#ifndef __GPUCOMPILE__
std::vector<int> TetrahedralTree::emptyBlock = {};

/**
TetrahedralTree.cc
This class stores the mesh nodes and elements in an Octree data
structure to optimize the element search operations

Author: Ali Sheharyar
Organization: Texas A&M University at Qatar
*/
TetrahedralTree::TetrahedralTree(const Vec3& origin, const Vec3& halfDimension)
    : m_origin(origin), m_halfDimension(halfDimension) {
  m_min.x = origin.x - halfDimension.x;
  m_min.y = origin.y - halfDimension.y;
  m_min.z = origin.z - halfDimension.z;
  m_max.x = origin.x + halfDimension.x;
  m_max.y = origin.y + halfDimension.y;
  m_max.z = origin.z + halfDimension.z;

  // Initially, there are no children
  for (int i = 0; i < 8; ++i) children[i] = nullptr;
}

TetrahedralTree::~TetrahedralTree() {
  // Recursively destroy octants
  for (int i = 0; i < 8; ++i) delete children[i];
}

// Check if a box overlaps with this node
bool TetrahedralTree::DoesBoxOverlap(const double bb[6]) const {
  if (m_max.x < bb[0] || m_max.y < bb[1] || m_max.z < bb[2]) return false;
  if (m_min.x > bb[3] || m_min.y > bb[4] || m_min.z > bb[5]) return false;
  return true;
}
#endif

// Determine which octant of the tree would contain 'point'
__GPULABEL__ int __TETRAHEDRALTREECLASS__::GetOctantContainingPoint(const __VEC3CLASS__& point) const
{
  int oct = 0;
  if (point.x >= m_origin.x) oct |= 4;
  if (point.y >= m_origin.y) oct |= 2;
  if (point.z >= m_origin.z) oct |= 1;
  return oct;
}

__GPULABEL__ bool __TETRAHEDRALTREECLASS__::IsLeafNode() const {
  // We are a leaf if we have no children. Since we either have none, or
  // all eight, it is sufficient to just check the first.
  return children[0] == nullptr;
}

#ifndef __GPUCOMPILE__
void TetrahedralTree::InsertMeshNode(Vec3 point, const int index) {
  // Check if it is a leaf node.
  if (!IsLeafNode()) {
    // We are at an interior node. Insert recursively into the
    // appropriate child octant.
    int octant = GetOctantContainingPoint(point);
    children[octant]->InsertMeshNode(point, index);
    return;
  }
  
  // Add the new point if the block is not full.
  if (nodes.size() < BlockCapacity) {
    nodes.push_back(std::make_pair(point, index));
    return;
  } 
  // Block is full, so we need to partition it.
  // Split the current node and create new empty trees for each child octant.
  for (int i = 0; i < 8; ++i) {
    // Compute new bounding box for this child
    Vec3 newOrigin = m_origin;
    newOrigin.x += m_halfDimension.x * (i & 4 ? .5f : -.5f);
    newOrigin.y += m_halfDimension.y * (i & 2 ? .5f : -.5f);
    newOrigin.z += m_halfDimension.z * (i & 1 ? .5f : -.5f);
    children[i] = new TetrahedralTree(newOrigin, m_halfDimension * .5f);
  }

  // Move the mesh nodes from the partitioned node (now marked as interior) to
  // its children.
  while (!nodes.empty()) {
    auto node = nodes.back();
    nodes.pop_back();
    const int oct = GetOctantContainingPoint(node.first);
    children[oct]->InsertMeshNode(node.first, node.second);
  }
  // Insert the new point in the appropriate octant.
  children[GetOctantContainingPoint(point)]->InsertMeshNode(point, index);
}

void TetrahedralTree::InsertMeshElement(const double bb[6], const int index) {
  if (IsLeafNode()) {
    // Add the element to the list of this octant.
    elements.push_back(index);
    return;
  } 
  // Check which children overlap with the element's bounding box.
  for (int i = 0; i < 8; i++) {
    if (!children[i]->DoesBoxOverlap(bb)) continue;
    children[i]->InsertMeshElement(bb, index);
  }
}
#endif

// It returns the list of tetrahedrons that intersects in a bounding box (Octree
// block) that contains the
// point passed as input.
#ifdef __GPUCOMPILE__
__device__ void TetrahedralTreeGPU::GetElementsInBlock(const Vec3GPU& point, const int *&tet_list_elems, int &num_elems) const {
    const TetrahedralTreeGPU* octreeNode = GetBlockFromPoint(point);

    if (octreeNode) {
        tet_list_elems = octreeNode->elements;
        num_elems = octreeNode->numelements;
        return;
    }

    tet_list_elems = nullptr;
    num_elems = 0;
  }
#else
const std::vector<int>& TetrahedralTree::GetElementsInBlock(const Vec3& point) const {
  const TetrahedralTree* octreeNode = GetBlockFromPoint(point);

  if (octreeNode) {
    return octreeNode->elements;
  }

  return emptyBlock;
}
#endif

// check if the point is inside the domain.
// This function is only executed at root to ensure that input point is inside
// the mesh's bounding box
// If we don't check this, the case when root is leaf node itself will return
// wrong block
__GPULABEL__ const __TETRAHEDRALTREECLASS__* __TETRAHEDRALTREECLASS__::GetBlockFromPoint(const __VEC3CLASS__& point) const
{
  if (!(m_min.x <= point.x && point.x <= m_max.x &&
        m_min.y <= point.y && point.y <= m_max.y &&
        m_min.z <= point.z && point.z <= m_max.z))
    return nullptr;

  return GetBlockFromPointHelper(point);
}

__GPULABEL__ const __TETRAHEDRALTREECLASS__* __TETRAHEDRALTREECLASS__::GetBlockFromPointHelper(
    const __VEC3CLASS__& point) const
{
  // If we're at a leaf node, it means, the point is inside this block
  if (IsLeafNode()) return this;
  // We are at the interior node, so check which child octant contains the
  // point
  int octant = GetOctantContainingPoint(point);
  return children[octant]->GetBlockFromPointHelper(point);
}
#ifndef __GPUCOMPILE__
#ifndef USEGPU
double TetrahedralTree::CreateGPUTransferObject(TetrahedralTreeGPU *&tree_gpu)
{
  tree_gpu = nullptr;
  return 0;
}
#endif
#endif
}
