#include "QComputations_SINGLE.hpp"
#include <complex>

bool is_hermit(const QComputations::Matrix<COMPLEX>& H) {
    bool res = true;

    for (int i = 0; res and i < H.n(); i++) {
        for (int j = i + 1; res and j < H.m(); j++) {
            if (std::abs(H[i][j] - std::conj(H[j][i])) >= QComputations::QConfig::instance().eps()) {
                res = false;
                std::cout << i << ' ' << j << std::endl;
                std::cout << H[i][j] << " " << H[j][i] << std::endl;
            }
        }
    }

    return res;
}

namespace QComputations {
    void prepare_state(TCH_State& state, double energy_diff = 1) {
        state.set_w(1 + energy_diff, 2);
        state.set_w(1, 1);
        state.set_g_level(COMPLEX(0.01 * std::sqrt(1 + energy_diff), 0), 2, 0);
        state.set_g_level(COMPLEX(0.01, 0), 2, 1);
        
        for (int i = 0; i < state.m(0); i++) {
            state.set_move(true, i, 0, 0, 2);
            state.set_move(true, i, 0, 1, 2);
            state.set_move(false, i, 0, 0, 1);
            state.set_atom(2, i, 0);
        }

        for (int j = 1; j < state.cavities_count(); j++) {
            for (int i = 0; i < state.m(j); i++) {
                state.set_move(true, i, j, 1, 2);
                state.set_move(false, i, j, 0, 1);
                state.set_move(false, i, j, 0, 2);
                state.set_atom(2, i, j);
            }
        }
        //state.set_n(1, 0, 0, 2);
    }
}

int main(void) {
    using namespace QComputations;
    constexpr int end = 10000;
    constexpr int steps_count = 50000;
    constexpr int n_start = 1;
    constexpr int n_end = 2;
    const std::string dir = "new_res/";
    bool is_common = true;
    bool is_all_ones = false;
    bool is_all_connected = false;
    int levels_count = 3;
    QConfig::instance().set_max_photons(10);

    auto time_vec = linspace(0, end, steps_count);
    /*
    {
        std::vector<size_t> grid = {1, 0};
        TCH_State state(grid);
        //prepare_state(state);
        state.set_n(1, 1);
        state.set_waveguide(0, 1, 0.01, 0);
        //std::cout << state.g(0, 0, 1, 2) << std::endl;
        H_TCH H(state);

        show_basis(H.get_basis());

        H.show();

        auto probs = quantum_master_equation(state, H, time_vec);

        make_probs_files(H, probs, time_vec, H.get_basis(), dir + "one_level");
    }
    {
        std::vector<size_t> grid = {1};
        TCH_State state(grid, levels_count);
        prepare_state(state);
        state.set_n(4, 0, 1, 2);
        //std::cout << state.g(0, 0, 1, 2) << std::endl;
        H_TCH H(state);

        show_basis(H.get_basis());

        H.show();

        auto probs = quantum_master_equation(state, H, time_vec);
        auto p = probs_to_qudits(probs, H.get_basis(), {state.ph_types_count()});

        make_probs_files(H, p.first, time_vec, p.second, dir + "one_cavity_reduced");
        make_probs_files(H, probs, time_vec, H.get_basis(), dir + "one_cavity");
    }
    */
    size_t n_atoms = 1;
    auto energy_diff = linspace(double(1), double(1), 1);
    auto waveguide_diff = linspace(double(0.001), double(0.1), 32);
    auto gamma_diff = linspace(double(0.001), double(0.1), 32);
    Matrix<double> w_time(C_STYLE, waveguide_diff.size(), gamma_diff.size());
    Matrix<double> w_prob_diff(C_STYLE, waveguide_diff.size(), gamma_diff.size());
    for (auto en: energy_diff) {
        std::cout << en << std::endl;
    }
    {
        std::vector<size_t> grid = {n_atoms};
        TCH_State state(grid, levels_count);
        prepare_state(state);
        //std::cout << state.g(0, 0, 1, 2) << std::endl;
        state.set_leak_for_cavity(0, 0.02);
        H_TCH H(state);

        show_basis(H.get_basis());

        //H.show();
        std::cout << H.size() << std::endl;

        //auto probs = schrodinger(state, H, time_vec);
        auto probs = quantum_master_equation(state, H, time_vec);
        //probs.show();
        auto p = probs_to_qudits(probs, H.get_basis(), {state.ph_types_count()});
        std::vector<std::string> st_names = {"|S>", "|A>", "|B>"};

        make_probs_files(H, p.first, time_vec, st_names, dir + "start_reduced");
        make_probs_files(H, probs, time_vec, H.get_basis(), dir + "start");
    }

    for (int n = n_start; is_common and n < n_end; n++) {
        #pragma omp parallel for
        for (int w = 0; w < waveguide_diff.size(); w++) {
            for (int g = 0; g < gamma_diff.size(); g++) {
                //std::vector<size_t> grid_config(n + 1, 0);
                std::vector<size_t> grid_config(2, n_atoms);
                //grid_config[0] = n_atoms;
                TCH_State state(grid_config, levels_count);
                prepare_state(state);
                //for (int i = 1; i < grid_config.size(); i++) {
                //    state.set_n(1, i, 1, 2);
                //    state.set_waveguide(0, i, 0.005, 0);
                //}

                state.set_leak_for_cavity(0, gamma_diff[g]);
                //state.set_n(n, 1, 1, 2);
                state.set_waveguide(0, 1, waveguide_diff[w], 0);

                H_TCH H(state);

                if (!is_hermit(H.get_matrix())) {
                    std::cerr << "NOT HERMIT!\n";
                }

                //std::cout << "H size = " << H.size() << std::endl;
                std::cout << w << " " << g << std::endl;

                /*
                if (w < 2) {
                    show_basis(H.get_basis());
                    H.show();
                }
                */

                //auto time_vec = linspace(0, int(double(M_PI)/(0.005)), int(double(M_PI)/(0.005)));
                auto probs = quantum_master_equation(state, H, time_vec);

                //make_probs_files(H, probs, time_vec, H.get_basis(), dir + "n=" + std::to_string(n) + "_common_all");

                for (int i = 0; i < state.m(0); i++) {
                    auto p = probs_to_qudits(probs, H.get_basis(), {state.ph_types_count() + i});
                    auto probs_in_cavity = p.first;
                    size_t target_index = 0;
                    for (int t = probs_in_cavity.m(); t >= 0; --t) {
                        if (probs_in_cavity[1][t] < probs_in_cavity[0][t] || probs_in_cavity[1][t] < probs_in_cavity[2][t]) {
                            target_index = t + 1;
                            break;
                        }
                    }
                    w_time[w][g] = target_index;
                    w_prob_diff[w][g] = probs_in_cavity[1][probs_in_cavity.m() - 1] - probs_in_cavity[2][probs_in_cavity.m() - 1];
                    auto basis = p.second;
                    std::vector<std::string> st_names = {"|S>", "|A>", "|B>"};

                    //make_probs_files(H, probs_in_cavity, time_vec, basis, dir + "n=" + std::to_string(n) + "_common atom = " + std::to_string(i));
                    //make_probs_files(H, probs_in_cavity, time_vec, st_names, dir + "n_atoms=" + std::to_string(state.m(0)) + "_wave_diff=" + std::to_string(waveguide_diff[w]) + " gamma = " + std::to_string(gamma_diff[g]) + " atom= " + std::to_string(i));
                    //make_probs_files(H, probs_in_cavity, time_vec, basis, dir + "qme_test");
                }

                /*
                probs = schrodinger(state, H, time_vec);

                for (int i = 0; i < state.m(0); i++) {
                    auto p = probs_to_qudits(probs, H.get_basis(), {state.ph_types_count() + i});
                    auto probs_in_cavity = p.first;
                    auto basis = p.second;

                    //make_probs_files(H, probs_in_cavity, time_vec, basis, dir + "n=" + std::to_string(n) + "_common atom = " + std::to_string(i));
                    //make_probs_files(H, probs_in_cavity, time_vec, basis, dir + "shr_test");
                }
                */
            }
        }
    }

    w_time.write_to_csv_file("wave_time.csv");
    w_prob_diff.write_to_csv_file("wave_prob_diff.csv");

    for (int n = n_start; is_all_connected and n < n_end; n++) {
        std::vector<size_t> grid_config(n + 1, 0);
        TCH_State state(grid_config);
        //state.set_n(1, 0);
        for (int i = 1; i < grid_config.size(); i++) {
            state.set_n(1, i);
            state.set_waveguide(0, i, 0.01, 0);
        }

        for (int i = 1; i < grid_config.size(); i++) {
            for (int j = i + 1; j < grid_config.size(); j++) {
                state.set_waveguide(i, j, 0.01, 0);
            }
        }

        state.set_leak_for_cavity(0, 0.01);

        H_TCH H(state);
        if (!is_hermit(H.get_matrix())) {
            std::cerr << "NOT HERMIT!\n";
        }

        std::cout << "H size = " << H.size() << std::endl;

        //show_basis(H.get_basis());
        //H.show();

        auto probs = quantum_master_equation(state, H, time_vec);

        make_probs_files(H, probs, time_vec, H.get_basis(), "res/n=" + std::to_string(n) + "_all_connected");
    }

    for (int n = n_start; is_all_ones and n < n_end; n++) {
        std::vector<size_t> grid_config(n + 1, 1);
        TCH_State state(grid_config);
        state.set_n(1, 0);
        for (int i = 1; i < grid_config.size(); i++) {
            state.set_n(1, i);
            state.set_waveguide(0, i, 0.01, 0);
        }

        state.set_leak_for_cavity(0, 0.01);

        H_TCH H(state);

        if (!is_hermit(H.get_matrix())) {
            std::cerr << "NOT HERMIT!\n";
        }

        std::cout << "H size = " << H.size() << std::endl;

        //show_basis(H.get_basis());
        //H.show();

        auto probs = quantum_master_equation(state, H, time_vec);

        make_probs_files(H, probs, time_vec, H.get_basis(), "res/n=" + std::to_string(n) + "_all_ones");
    }

    return 0;
}