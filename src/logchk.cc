#include <cmath>
#include <cfloat>
#include <cstdio>
#include <random>

//using std::fma;
double fma(double a, double b, double c) { return a*b+c; }

static inline double horner(double x, double a0) {
    return a0;
}

template <typename... T>
static double horner(double x, double a0, T... tail) {
    return fma(x, horner(x, tail...), a0);
}

static inline double horner1(double x, double a0) {
    return x+a0;
}

template <typename... T>
static double horner1(double x, double a0, T... tail) {
    return fma(x, horner1(x, tail...), a0);
}

constexpr double sqrt2 = 1.41421356237309504880;

// Cancellation-corrected ln(2) = ln2C3 + ln2C4:

constexpr double ln2C3 = 0.693359375;
constexpr double ln2C4 = -2.121944400546905827679e-4;

// Polynomial coefficients (Q is also order 5,
// but monic).

constexpr double P0log = 7.70838733755885391666E0;
constexpr double P1log = 1.79368678507819816313e1;
constexpr double P2log = 1.44989225341610930846e1;
constexpr double P3log = 4.70579119878881725854e0;
constexpr double P4log = 4.97494994976747001425e-1;
constexpr double P5log = 1.01875663804580931796e-4;

constexpr double Q0log = 2.31251620126765340583E1;
constexpr double Q1log = 7.11544750618563894466e1;
constexpr double Q2log = 8.29875266912776603211e1;
constexpr double Q3log = 4.52279145837532221105e1;
constexpr double Q4log = 1.12873587189167450590e1;

double log_c(double x) {
    double g = std::logb(x);
    double u = std::scalbn(x, -(int)g);

    if (u>sqrt2) {
        u /= 2;
        g += 1;
    }

    double z = u-1;
    auto pz = horner(z, P0log, P1log, P2log, P3log, P4log, P5log);
    auto qz = horner1(z, Q0log, Q1log, Q2log, Q3log, Q4log);

    auto z2 = z*z;
    auto z3 = z2*z;

    auto r = z3*pz/qz;
    r = fma(z2, -0.5, fma(g, ln2C4, r)) + z;
    r = fma(g, ln2C3, r);

    return r;
}

void run_check(double lb, double ub, int n) {
    std::minstd_rand R(12345);
    std::uniform_real_distribution<> u(lb, ub);

    for (int i = 0; i<n; ++i) {
        double x = u(R);
        double lx = std::log(x);
        double cx = log_c(x);

        double err = std::abs(lx-cx);
        double ulp = std::scalbn(1., std::logb(lx)-52);
        int err_ulp = std::ceil(err/ulp);

        if (err_ulp>1) {
            std::printf("ulp: %d\nx: %a\nlx: %a\ncx: %a\n", err_ulp, x, lx, cx);
        }
    }
}

int main() {
    run_check(0.003, 2300, 1000000);
}
