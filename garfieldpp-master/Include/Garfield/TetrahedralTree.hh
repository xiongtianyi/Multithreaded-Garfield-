// Include this header if we're compiling with the GPU or this is the first time without
#if defined(__GPUCOMPILE__) || !defined(TETRAHEDRAL_TREE_H)

#if !defined(__GPUCOMPILE__) && !defined(TETRAHEDRAL_TREE_H)
#define TETRAHEDRAL_TREE_H
#endif

// undefine everything first 
#ifdef __TETRAHEDRALTREECLASS__
#undef __TETRAHEDRALTREECLASS__
#undef __VEC3CLASS__
#undef __GPULABEL__
#endif

#ifdef __GPUCOMPILE__

#define __TETRAHEDRALTREECLASS__ TetrahedralTreeGPU
#define __VEC3CLASS__ Vec3GPU
#define __GPULABEL__ __device__

#else

class TetrahedralTreeGPU;

#define __TETRAHEDRALTREECLASS__ TetrahedralTree
#define __VEC3CLASS__ Vec3
#define __GPULABEL__

#endif

#include <cstddef>
#include <vector>
#include <utility>

namespace Garfield {

// TODO: replace this class with ROOT's TVector3 class

struct __VEC3CLASS__ {
  float x = 0., y = 0., z = 0.;

  __GPULABEL__ __VEC3CLASS__() {}
  __GPULABEL__ __VEC3CLASS__(float _x, float _y, float _z) : x(_x), y(_y), z(_z) {}

  __GPULABEL__ __VEC3CLASS__ operator+(const __VEC3CLASS__& r) const {
    return __VEC3CLASS__(x + r.x, y + r.y, z + r.z);
  }

  __GPULABEL__ __VEC3CLASS__ operator-(const __VEC3CLASS__& r) const {
    return __VEC3CLASS__(x - r.x, y - r.y, z - r.z);
  }

  __GPULABEL__ __VEC3CLASS__& operator+=(const __VEC3CLASS__& r) {
    x += r.x;
    y += r.y;
    z += r.z;
    return *this;
  }

  __GPULABEL__ __VEC3CLASS__& operator-=(const __VEC3CLASS__& r) {
    x -= r.x;
    y -= r.y;
    z -= r.z;
    return *this;
  }

  __GPULABEL__ __VEC3CLASS__ operator*(float r) const { return __VEC3CLASS__(x * r, y * r, z * r); }

  __GPULABEL__ __VEC3CLASS__ operator/(float r) const { return __VEC3CLASS__(x / r, y / r, z / r); }
};

/**

\brief Helper class for searches in field maps.

This class stores the mesh nodes and elements in an Octree data
structure to optimize the element search operations

Author: Ali Sheharyar

Organization: Texas A&M University at Qatar

*/
class __TETRAHEDRALTREECLASS__ {
 public:
 #ifdef __GPUCOMPILE__
  // Constructor
  TetrahedralTreeGPU() = default;

  // Destructor
  ~TetrahedralTreeGPU() {};
 #else
  // Constructor
  TetrahedralTree(const Vec3& origin, const Vec3& halfDimension);

  /// Destructor
  ~TetrahedralTree();
#endif

  #ifndef __GPUCOMPILE__
  // Insert a mesh node (a vertex/point) to the tree
  void InsertMeshNode(Vec3 point, const int index);

  /// Insert a mesh element with given bounding box and index to the tree.
  void InsertMeshElement(const double bb[6], const int index);

  /// Create and initialise GPU Transfer class
  double CreateGPUTransferObject(TetrahedralTreeGPU *&tree_gpu);
  #endif

 private:
  static std::vector<int> emptyBlock;

  // Physical centre of this tree node.
  __VEC3CLASS__ m_origin;
  // Half the width/height/depth of this tree node. 
  __VEC3CLASS__ m_halfDimension;
  // Storing min and max points for convenience
  __VEC3CLASS__ m_min, m_max;  

  // The tree has up to eight children and can additionally store
  // a list of mesh nodes and mesh elements.
  // Pointers to child octants.
  __TETRAHEDRALTREECLASS__* children[8];  

  // Children follow a predictable pattern to make accesses simple.
  // Here, - means less than 'origin' in that dimension, + means greater than.
  // child:	0 1 2 3 4 5 6 7
  // x:     - - - - + + + +
  // y:     - - + + - - + +
  // z:     - + - + - + - +

  #ifndef __GPUCOMPILE__
  std::vector<std::pair<__VEC3CLASS__, int> > nodes;
  #endif

  #ifdef __GPUCOMPILE__
  int* elements{nullptr};
  int numelements{0};
  #else
  std::vector<int> elements;
  #endif

  static const size_t BlockCapacity = 10;

  #ifndef __GPUCOMPILE__
  // Check if the given box overlaps with this tree node.
  bool DoesBoxOverlap(const double bb[6]) const;
  #endif
  // Check if this tree node is a leaf or intermediate node.
  __GPULABEL__ bool IsLeafNode() const;


public:
  // Get all tetrahedra linked to a block corresponding to the given point
  #ifdef __GPUCOMPILE__
  __device__ void GetElementsInBlock(const Vec3GPU& point, const int *&tet_list_elems, int &num_elems) const;
  #else
  const std::vector<int>& GetElementsInBlock(const Vec3& point) const;
  #endif

private:
  __GPULABEL__ int GetOctantContainingPoint(const __VEC3CLASS__& point) const;
  // Get a block containing the input point
  __GPULABEL__ const __TETRAHEDRALTREECLASS__* GetBlockFromPoint(const __VEC3CLASS__& point) const;
  // A helper function used by the function above.
  // Called recursively on the child nodes.
  __GPULABEL__ const __TETRAHEDRALTREECLASS__* GetBlockFromPointHelper(const __VEC3CLASS__& point) const;  

#ifdef __GPUCOMPILE__
  friend class TetrahedralTree;
#endif
};

}

#endif
