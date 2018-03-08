#include <iostream>
#include <iomanip>
#include <random>

#include "ulp_check.h"
#include "common_opt.h"

template <typename V, typename Ref, typename Eval>
void harness(
    std::ostream& out,
    const common_opt& opt,
    const char* ref_name,  Ref ref,
    const char* eval_name, Eval eval,
    V lb, V ub)
{
    using std::cout;
    using std::setw;

    std::mt19937_64 R(opt.seed);
    std::uniform_real_distribution<V> u(lb, ub);

    auto gen  = [&]() { return u(R); };

    if (opt.raw) {
        out << "# lb=" << lb << "; ub=" << ub << '\n';
        out << "# x " << ref_name << " " << eval_name << '\n';
        for (std::size_t i = 0; i<opt.N; ++i) {
            double x = gen();
            out << std::setprecision(18)
                << setw(25) << x << setw(25) << ref(x) << setw(25) << eval(x) << '\n';
        }
    }
    else {
        auto result = ulp_check(opt.N, ref, eval, gen);

        out << ref_name << " vs " << eval_name << " on [" << lb << ", " << ub << "]\n";
        pretty_print(out, result);
        out << '\n';
    }
}

