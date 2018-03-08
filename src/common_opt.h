#pragma once

#include <stdexcept>

struct common_opt {
    unsigned N = 10000;     // number of trials per interval
    unsigned seed = 12345;  // RNG seed
    bool raw = false;
    bool help = false;
};

common_opt common_options(int& argc, char** argv);
extern const char* common_option_summary;
