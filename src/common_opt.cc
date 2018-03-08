#include <stdexcept>

#include "common_opt.h"

const char* common_option_summary =
" -n COUNT       run N samples per interval\n"
" -s SEED        use integer SEED to seed each interval\n"
" -r             print raw values instead of summary\n"
" -h             display this help message\n";

common_opt common_options(int& argc, char** argv) {
    auto skip = [&]() { ++argv; };
    auto consume = [&]() {
        --argc;
        for (char** p = argv; *p; ++p) *p = p[1];
    };

    common_opt opt;
    enum parse_state { opt_none, opt_n, opt_s } state = opt_none;

    ++argv;
    while (*argv) {
        switch (state) {
        case opt_none:
            switch (**argv=='-'? (*argv)[1]: 0) {
            case 'n':
                consume(), state = opt_n;
                break;
            case 's':
                consume(), state = opt_s;
                break;
            case 'r':
                consume();
                opt.raw = true;
                break;
            case 'h':
                consume();
                opt.help = true;
                break;
            default:
                skip();
            }
            break;

        case opt_n:
            opt.N = std::strtoul(*argv, 0, 0);
            consume(), state = opt_none;
            break;

        case opt_s:
            opt.seed = std::strtoul(*argv, 0, 0);
            consume(), state = opt_none;
            break;
        }
    }

    return opt;
}
