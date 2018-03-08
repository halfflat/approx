#pragma once

#include <cmath>

namespace fmaselect {

template <typename V>
V fma(V a, V b, V c) {
#ifdef NOFMA
#warning "Replacing fma() with explicit multiply then add."
    return a*b+c;

#else
    return std::fma(a, b, c);
#endif
}

} // namespace fmaselect
