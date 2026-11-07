#pragma once
#include "material/MaterialData.h"

namespace phonomc {
// Expand a validated unique first-BZ grid to an equal-weight representation of
// all shortest reciprocal-translation representatives. Every original q gets
// L samples, every distinct representative L/m samples. Thus normalization by
// L*Nq implements exactly the inverse-multiplicity BZ quadrature, without
// changing the physical occupation or the uniform carrier sampling contract.
// Returns L. The expanded table is internal, not the raw input-file schema.
unsigned complete_callaway_bz(MaterialData& data);
} // namespace phonomc
