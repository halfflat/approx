#pragma once

// Compare reference function and evaluation function across sequence of values;
// return histogram of ULP absolute error and worst case.

#include <cfloat>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <iomanip>
#include <map>
#include <utility>

template <typename V>
struct ulp_check_result {
    struct { V x, ref, eval, err, err_ulp; } worst;
    std::map<std::uintmax_t, std::uintmax_t> ulp_hist;
    std::uintmax_t count;
};

template <typename Generator, typename RefFn, typename EvalFn, typename V = decltype(std::declval<Generator>()())>
ulp_check_result<V> ulp_check(std::size_t n, RefFn ref, EvalFn eval, Generator gen) {
    ulp_check_result<V> h;
    h.worst.err_ulp = 0;
    h.worst.err = 0;
    h.count = 0;

    for (std::size_t i = 0; i<n; ++i) {
        V x = gen();
        V r = ref(x);
        V y = eval(x);

        V err = std::abs(y-r);
        V ulp = std::scalbn(V(1), std::logb(r)-52);
        V err_ulp = err/ulp;

        if (i==0 || err_ulp>h.worst.err_ulp) {
            h.worst = {x, r, y, err, err_ulp};
        }

        std::uintmax_t key = err_ulp>UINTMAX_MAX? UINTMAX_MAX: err_ulp;
        ++h.ulp_hist[key];
        ++h.count;
    }

    return h;
}

template <typename V>
std::ostream& pretty_print(std::ostream& out, const ulp_check_result<V>& h, int maxbars = 80, int barwidth = 80) {
    using uintmax = std::uintmax_t;

    if (h.ulp_hist.empty()) return out;

    double maxc = 0;
    std::size_t nbin = h.ulp_hist.size();

    for (auto& kv: h.ulp_hist) {
        if (kv.second>maxc) maxc = kv.second;
    }

    // search heuristic for a reasonable cut off:
    double threshold = 0;
    if (nbin > maxbars) {
        const double alpha = 0.7;
        double tmin = maxc/barwidth;
        for (double t = maxc*alpha; t>tmin; t *= alpha) {
            std::size_t count = 0;
            for (auto& kv: h.ulp_hist) {
                if (kv.second>t) ++count;
            }
            if (count>maxbars) break;
        }
        threshold = tmin/alpha;
    }

    int margin = std::ceil(std::log10(h.ulp_hist.end()->first));
    double xinc = maxc/barwidth;

    auto emit = [&](std::pair<uintmax, uintmax> kv) {
        out << std::setw(margin) << kv.first << '|';
        int percent = std::round(100.0*kv.second/h.count);
        out << std::setw(3) << percent << "%  ";
        int x = kv.second/xinc;
        if (x) out << std::string(x, '*');
        else if (kv.second) out << '.';
        out << '\n';
    };

    uintmax remaining = h.count;
    int lines = 1;
    bool trailing = false;

    uintmax i = h.ulp_hist.begin()->first;
    for (auto& kv: h.ulp_hist) {
        if (kv.second<threshold) {
            trailing = true;
            continue;
        }

        if (kv.first > i) {
            if (threshold==0 && kv.first < i+4)
                while (i<kv.first) emit({i++, uintmax(0)});
            else
                out << "//\n";
        }

        trailing = false;
        emit(kv);
        i = kv.first+1;
        ++lines;

        if (lines>maxbars) {
            trailing = true;
            break;
        }
    }

    if (trailing) out << "+++\n";

    out << "worst case:\n"
        << std::setprecision(17)
        << "x:       " << h.worst.x << "\n"
        << "ref(x):  " << h.worst.ref << "\n"
        << "eval(x): " << h.worst.eval << "\n"
        << "abserr:  " << h.worst.err << "\n"
        << "ulperr:  " << h.worst.err_ulp << "\n";

    return out;
}

