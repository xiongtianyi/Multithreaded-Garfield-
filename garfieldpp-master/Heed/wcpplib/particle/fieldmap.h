#ifndef FIELDMAP_H
#define FIELDMAP_H

#include "wcpplib/geometry/vec.h"
#include <limits>

namespace Heed {

/// Retrieve electric and magnetic field.

class fieldmap {
 public:
  fieldmap() = default;
  virtual ~fieldmap() = default;
  virtual void evaluate(const point& /*pt*/, vec& efield, vec& bfield, double& mrange) const {
    efield.x = bfield.x = 0.;
    efield.y = bfield.y = 0.;
    efield.z = bfield.z = 0.;
    mrange = std::numeric_limits<double>::max();
  }
  virtual bool inside(const point& /*pt*/) { return true; }

};
}

#endif
