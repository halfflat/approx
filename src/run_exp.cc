#include <cfloat>

#include "fmaselect.h"
#include "horner.h"
#include "common_opt.h"
#include "harness.h"

// exp(x) following Cephes library algorithm;
// order 6 rational poly approximation R(x)/R(-x)
// with even coeffs Q and odd coeffs P below.

constexpr double P0exp = 9.99999999999999999910E-1;
constexpr double P1exp = 3.02994407707441961300E-2;
constexpr double P2exp = 1.26177193074810590878E-4;

constexpr double Q0exp = 2.00000000000000000009E0;
constexpr double Q1exp = 2.27265548208155028766E-1;
constexpr double Q2exp = 2.52448340349684104192E-3;
constexpr double Q3exp = 3.00198505138664455042E-6;

// Cancellation-corrected ln(2) = ln2C1 + ln2C2:

constexpr double ln2C1 = 6.93145751953125E-1;
constexpr double ln2C2 = 1.42860682030941723212E-6;

// 1/ln(2):

constexpr double ln2inv = 1.4426950408889634073599;

// Min and max argument values for finite and normal
// double-precision exponential.

constexpr double exp_minarg = -708.3964185322641;
constexpr double exp_maxarg = 709.782712893384;

using fmaselect::fma;

double exp_c(double x) {
    double n = std::floor(fma(ln2inv, x, 0.5));
    double g = fma(n, -ln2C1, x);
    g = fma(n, -ln2C2, g);

    auto gg = g*g;

    auto odd = g * horner(gg, P0exp, P1exp, P2exp);
    auto even = horner(gg, Q0exp, Q1exp, Q2exp, Q3exp);

    // Compute R(g)/R(-g) - 1 = 2*g*P(g^2) / (Q(g^2)-g*P(g^2))

    auto expg = fma(2., odd/(even-odd), 1.);
    return std::scalbn(expg, (int)n);
}

#if 0
void harness(bool raw, int N, double lb, double ub, int seed=12345) {
    using std::cout;
    using std::setw;

    std::mt19937_64 R(seed);
    std::uniform_real_distribution<> u(lb, ub);

    auto gen  = [&]() { return u(R); };
    auto ref  = [](double x) { return std::exp(x); };
    auto eval = &exp_c;

    if (raw) {
        cout << "# lb=" << lb << "; ub=" << ub << '\n';
        cout << "# x exp, exp_c\n";
        for (std::size_t i = 0; i<N; ++i) {
            double x = gen();
            cout << std::setprecision(18)
                 << setw(25) << x << setw(25) << ref(x) << setw(25) << eval(x) << '\n';
        }
    }
    else {
        auto result = ulp_check(N, ref, eval, gen);

        cout << "exp vs exp_c on [" << lb << ", " << ub << "]\n";
        pretty_print(std::cout, result);
        cout << '\n';
    }
}
#endif

int main(int argc, char** argv) {
    auto opt = common_options(argc, argv);
    if (opt.help) {
        std::cout << "usage: run_exp [OPTIONS]\n" << common_option_summary;
        return 0;
    }

    auto exp = [](double x) { return std::exp(x); };

    harness<double>(std::cout, opt, "std::exp", exp, "exp_c", exp_c, -0.1, 0.1);
    harness<double>(std::cout, opt, "std::exp", exp, "exp_c", exp_c, exp_minarg, exp_maxarg);
    return 0;
}
