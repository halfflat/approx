#include <cmath>
#include <cstring>
#include <stdexcept>
#include <iostream>

#ifdef USE_ARB
#include <arb.h>
#include <arb_hypgeom.h>
#endif

#include "ref_functions.h"

unary_fn<double> std_exp_double = {
    [](double x) { return std::exp(x); }, "std::exp"
};

unary_fn<double> std_expm1_double = {
    [](double x) { return std::expm1(x); }, "std::expm1"
};

unary_fn<double> std_log_double = {
    [](double x) { return std::log(x); }, "std::log"
};

unary_fn<double> std_log1p_double = {
    [](double x) { return std::log1p(x); }, "std::log1p"
};

#ifdef USE_ARB

// Wrapper around arb_t; manages memory, apply method
// evaluates up to double precision bits + 1 (i.e. 53).

struct arb {
    static constexpr slong prec = 53;

    arb() { init(); }
    arb(const arb_t& x) { init(); (*this)=x; }
    arb(const arb& x) { init(); (*this)=x; }
    arb(arb&& x) {
        std::memcpy(ptr(), x.ptr(), sizeof(arb_struct));
        x.allocated = false;
    }

    explicit arb(double x) { init(); (*this)=x; }

    arb& operator=(const arb& x) { return arb_set(a, x.a), *this; }
    arb& operator=(const arb_t& x) { return arb_set(a, x), *this; }
    arb& operator=(double x) { return arb_set_d(a, x), *this; }

    arb& operator=(arb&& x) {
        std::memcpy(ptr(), x.ptr(), sizeof(arb_struct));
        x.allocated = false;
        return *this;
    }

    arb_t& get() { return a; }
    const arb_t& get() const { return a; }

    explicit operator double() const {
        return arf_get_d(arb_midref(a), ARF_RND_NEAR);
    }

    double radius() const {
        return mag_get_d(arb_radref(a));
    }

    ~arb() { if (allocated) arb_clear(a); }

    arb apply(void (*arbfn)(arb_t, const arb_t, slong)) const {
        const slong maxprec = 20000;
        arb y;
        for (slong p = prec+2; p<maxprec; p += p/2) {
            (*arbfn)(y.a, a, p);
            if (arb_rel_accuracy_bits(y.a)>=prec) return y;
        }
        throw std::runtime_error("lost precision");
    }

private:
    arb_t a;
    bool allocated = false;

    void init() { arb_init(a); allocated = true; }
    arb_struct* ptr() { return &a[0]; }
    const arb_struct* ptr() const { return &a[0]; }
};

unary_fn<double> arb_exp_double = {
    [](double x) {
        return static_cast<double>(arb(x).apply(arb_exp));
    },
    "arb_exp"
};

unary_fn<double> arb_expm1_double = {
    [](double x) {
        return static_cast<double>(arb(x).apply(arb_expm1));
    },
    "arb_expm1"
};

unary_fn<double> arb_log_double = {
    [](double x) {
        return static_cast<double>(arb(x).apply(arb_log));
    },
    "arb_log"
};

unary_fn<double> arb_log1p_double = {
    [](double x) {
        return static_cast<double>(arb(x).apply(arb_log1p));
    },
    "arb_log1p"
};

void exprel(arb_t res, const arb_t x, slong prec) {
    arb_t one, two;
    arb_init(one);
    arb_init(two);

    arb_one(one);
    arb_set_si(two, 2);
    arb_hypgeom_1f1(res, one, two, x, false, prec);
}

unary_fn<double> arb_exprel_double = {
    [](double x) {
        return static_cast<double>(arb(x).apply(exprel));
    },
    "arb_exprel"
};

#endif // def USE_ARB
