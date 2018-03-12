#pragma once

// Wrap a unary function V -> V with metadata.

#include <functional>
#include <string>

template <typename V>
struct unary_fn {
    std::function<V (V)> fn;
    std::string name;

    using value_type = V;

    template <typename Fn>
    explicit unary_fn(Fn fn): fn(fn), name("<>") {}

    template <typename Fn>
    unary_fn(Fn fn, std::string name): fn(fn), name(std::move(name)) {}

    V operator()(V x) const { return fn(x); }
};

