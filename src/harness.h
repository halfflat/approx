#pragma once

#include <iostream>
#include <iomanip>
#include <random>

#include "common_opt.h"
#include "unary_fn.h"
#include "ulp_check.h"

template <typename V>
void harness(
    std::ostream& out,
    const common_opt& opt,
    const unary_fn<V>& ref,
    const unary_fn<V>& eval,
    V lb, V ub)
{
    using std::cout;
    using std::setw;

    std::mt19937_64 R(opt.seed);
    std::uniform_real_distribution<V> u(lb, ub);

    auto gen  = [&]() { return u(R); };

    if (opt.raw) {
        out << "# lb=" << lb << "; ub=" << ub << '\n';
        out << "# x " << ref.name << " " << eval.name << '\n';
        for (std::size_t i = 0; i<opt.N; ++i) {
            double x = gen();
            out << std::setprecision(18)
                << setw(25) << x << setw(25) << ref(x) << setw(25) << eval(x) << '\n';
        }
    }
    else {
        auto result = ulp_check(opt.N, ref.fn, eval.fn, gen);

        out << ref.name << " vs " << eval.name << " on [" << lb << ", " << ub << "]\n";
        pretty_print(out, result);
        out << '\n';
    }
}
