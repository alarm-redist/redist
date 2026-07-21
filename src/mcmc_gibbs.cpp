#include "mcmc_gibbs.h"
#include "make_swaps_helper.h"
#include "map_calc.h"

/*
 * Helper function to iterate over constraints and apply them
 */
double add_constraint(const std::string &name, Rcpp::List constraints, std::vector<int> districts,
                      Rcpp::NumericVector &psi_vec, std::function<double(Rcpp::List, int)> fn_constr) {
    if (!constraints.containsElementNamed(name.c_str()))
        return 0;

    Rcpp::List constr = constraints[name];
    double val = 0;
    double psi = 0.0;
    for (int i = 0; i < constr.size(); i++) {
        Rcpp::List constr_inst = constr[i];
        double strength = constr_inst["strength"];
        if (strength != 0) {
            for (int dist : districts) {
                psi = fn_constr(constr_inst, dist);
                psi_vec[name] = psi + psi_vec[name];
                val += strength * (psi);
            }
        }
    }
    return val;
}

/*
 * Add specific constraint weights & return the cumulative weight vector
 */
double calc_gibbs_tgt(const arma::subview_col<arma::uword> &plan, int n_distr, int V,
                      std::vector<int> districts, Rcpp::NumericVector &psi_vec, const std::vector<unsigned int> &pop,
                      double parity, const Graph &g, Rcpp::List constraints) {
    if (constraints.size() == 0)
        return 0.0;
    double log_tgt = 0;
    double n_consider = (double)districts.size();

    log_tgt += add_constraint("pop_dev", constraints, districts, psi_vec,
                              [&](Rcpp::List l, int distr) -> double {
                                  return eval_pop_dev(plan, distr, distr, pop, parity);
                              });

    std::vector<int> distr_dummy = {1};
    log_tgt += add_constraint(
        "splits", constraints, distr_dummy, psi_vec, [&](Rcpp::List l, int distr) -> double {
            return eval_splits(plan, distr, Rcpp::as<arma::uvec>(l["admin"]), l["n"], false);
        });
    log_tgt += add_constraint(
        "multisplits", constraints, distr_dummy, psi_vec, [&](Rcpp::List l, int distr) -> double {
            return eval_multisplits(plan, distr, Rcpp::as<arma::uvec>(l["admin"]), l["n"], false);
        });
    log_tgt += add_constraint(
        "total_splits", constraints, distr_dummy, psi_vec, [&](Rcpp::List l, int distr) -> double {
            return eval_total_splits(plan, distr, Rcpp::as<arma::uvec>(l["admin"]), l["n"], false);
        });

    log_tgt += add_constraint(
        "segregation", constraints, districts, psi_vec, [&](Rcpp::List l, int distr) -> double {
            return eval_segregation(plan, distr, distr, V, Rcpp::as<arma::uvec>(l["group_pop"]),
                                    Rcpp::as<arma::uvec>(l["total_pop"]));
        });

    log_tgt += add_constraint(
        "grp_pow", constraints, districts, psi_vec, [&](Rcpp::List l, int distr) -> double {
            return eval_grp_pow(plan, V, distr, distr, Rcpp::as<arma::uvec>(l["group_pop"]),
                                Rcpp::as<arma::uvec>(l["total_pop"]), Rcpp::as<double>(l["tgt_group"]),
                                Rcpp::as<double>(l["tgt_other"]), Rcpp::as<double>(l["pow"]));
        });

    log_tgt += add_constraint(
        "grp_hinge", constraints, districts, psi_vec, [&](Rcpp::List l, int distr) -> double {
            return eval_grp_hinge(plan, V, distr, distr, Rcpp::as<arma::vec>(l["tgts_group"]),
                                  Rcpp::as<arma::uvec>(l["group_pop"]), Rcpp::as<arma::uvec>(l["total_pop"]));
        });

    log_tgt += add_constraint(
        "grp_inv_hinge", constraints, districts, psi_vec, [&](Rcpp::List l, int distr) -> double {
            return eval_grp_hinge(plan, V, distr, distr, Rcpp::as<arma::vec>(l["tgts_group"]),
                                  Rcpp::as<arma::uvec>(l["group_pop"]), Rcpp::as<arma::uvec>(l["total_pop"]));
        });

    log_tgt += add_constraint("compet", constraints, districts, psi_vec,
                              [&](Rcpp::List l, int distr) -> double {
                                  arma::uvec dvote = l["dvote"];
                                  arma::uvec total = dvote + Rcpp::as<arma::uvec>(l["rvote"]);
                                  return eval_grp_pow(plan, V, distr, distr, dvote, total, 0.5,
                                                      0.5, Rcpp::as<double>(l["pow"]));
                              });

    log_tgt += add_constraint(
        "status_quo", constraints, districts, psi_vec, [&](Rcpp::List l, int distr) -> double {
            return eval_sq_entropy(plan, Rcpp::as<std::vector<unsigned int>>(l["current"]), distr, distr, pop, n_distr,
                                   Rcpp::as<int>(l["n_current"]), V);
        });

    log_tgt += add_constraint(
        "incumbency", constraints, districts, psi_vec, [&](Rcpp::List l, int distr) -> double {
            return eval_inc(plan, distr, distr, Rcpp::as<arma::uvec>(l["incumbents"]));
        });

    log_tgt += add_constraint(
        "polsby", constraints, districts, psi_vec, [&](Rcpp::List l, int distr) -> double {
            return eval_polsby(plan, distr, distr, V, Rcpp::as<arma::ivec>(l["from"]), Rcpp::as<arma::ivec>(l["to"]),
                               Rcpp::as<arma::vec>(l["area"]), Rcpp::as<arma::vec>(l["perimeter"]));
        });

    log_tgt += add_constraint(
        "fry_hold", constraints, districts, psi_vec, [&](Rcpp::List l, int distr) -> double {
            return eval_fry_hold(plan, distr, Rcpp::as<arma::uvec>(l["total_pop"]), Rcpp::as<arma::mat>(l["ssdmat"]),
                                 Rcpp::as<double>(l["denominator"]));
        });

    log_tgt += add_constraint(
        "log_st", constraints, districts, psi_vec, [&](Rcpp::List l, int distr) -> double {
            return eval_log_st(plan, g, Rcpp::as<arma::uvec>(l["admin"]), n_distr) / n_consider;
        });

    log_tgt += add_constraint(
        "edges_removed", constraints, districts, psi_vec,
        [&](Rcpp::List l, int distr) -> double { return eval_er(plan, g, n_distr) / n_consider; });

    log_tgt += add_constraint(
        "qps", constraints, districts, psi_vec, [&](Rcpp::List l, int distr) -> double {
            return eval_qps(plan, distr, Rcpp::as<arma::uvec>(l["total_pop"]), Rcpp::as<arma::uvec>(l["cities"]),
                            Rcpp::as<int>(l["n_city"]), n_distr);
        });

    log_tgt += add_constraint("custom", constraints, districts, psi_vec,
                              [&](Rcpp::List l, int distr) -> double {
                                  Rcpp::Function fn = l["fn"];
                                  return Rcpp::as<Rcpp::NumericVector>(fn(plan, distr))[0];
                              });

    return log_tgt;
}