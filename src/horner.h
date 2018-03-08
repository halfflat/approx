#pragma once

#include "fmaselect.h"

template <typename V>
static inline double horner(V x, V a0) {
    return a0;
}

template <typename V, typename... T>
static double horner(V x, V a0, T... tail) {
    using fmaselect::fma;
    return fma(x, horner(x, tail...), a0);
}

template <typename V>
static inline double horner1(V x, V a0) {
    return x+a0;
}

template <typename V, typename... T>
static double horner1(V x, V a0, T... tail) {
    using fmaselect::fma;
    return fma(x, horner1(x, tail...), a0);
}

