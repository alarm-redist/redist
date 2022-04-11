#include "mcmc_gibbs.h"

/*
 * Helper function to iterate over constraints and apply them
 */
double add_constraint(const std::string& name, List constraints,
                      std::vector<int> districts, NumericVector &psi_vec,
                      std::function<double(List, int)> fn_constr) {
    if (!constraints.containsElementNamed(name.c_str())) return 0;

    List constr = constraints[name];
    double val = 0;
    double psi = 0.0;
    for (int i = 0; i < constr.size(); i++) {
        List constr_inst = constr[i];
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
double calc_gibbs_tgt(const subview_col<uword> &plan, int n_distr, int V,
                      std::vector<int> districts, NumericVector &psi_vec, const uvec &pop,
                      double parity, const Graph &g, List constraints) {
    if (constraints.size() == 0) return 0.0;
    double log_tgt = 0;
    double n_consider = (double) districts.size();

    log_tgt += add_constraint("pop_dev", constraints, districts, psi_vec,
                              [&] (List l, int distr) -> double {
                                  return eval_pop_dev(plan, distr,
                                                      pop, parity);
                              });

    log_tgt += add_constraint("splits", constraints, districts, psi_vec,
                              [&] (List l, int distr) -> double {
                                  return eval_splits(plan, distr, as<uvec>(l["admin"]), l["n"], false);
                              });

    log_tgt += add_constraint("multisplits", constraints, districts, psi_vec,
                              [&] (List l, int distr) -> double {
                                  return eval_multisplits(plan, distr, as<uvec>(l["admin"]), l["n"], false);
                              });

    log_tgt += add_constraint("total_splits", constraints, districts, psi_vec,
                              [&] (List l, int distr) -> double {
                                  return eval_total_splits(plan, distr, as<uvec>(l["admin"]), l["n"]);
                              });

    log_tgt += add_constraint("segregation", constraints, districts, psi_vec,
                              [&] (List l, int distr) -> double {
                                  return eval_segregation(plan, distr, as<uvec>(l["group_pop"]), as<uvec>(l["total_pop"]));
                              });

    log_tgt += add_constraint("grp_pow", constraints, districts, psi_vec,
                              [&] (List l, int distr) -> double {
                                  return eval_grp_pow(plan, distr, as<uvec>(l["group_pop"]),
                                                      as<uvec>(l["total_pop"]), as<double>(l["tgt_group"]),
                                                      as<double>(l["tgt_other"]), as<double>(l["pow"]));
                              });

    log_tgt += add_constraint("grp_hinge", constraints, districts, psi_vec,
                              [&] (List l, int distr) -> double {
                                  return eval_grp_hinge(plan, distr, as<vec>(l["tgts_group"]),
                                                        as<uvec>(l["group_pop"]), as<uvec>(l["total_pop"]));
                              });

    log_tgt += add_constraint("grp_inv_hinge", constraints, districts, psi_vec,
                              [&] (List l, int distr) -> double {
                                  return eval_grp_hinge(plan, distr, as<vec>(l["tgts_group"]),
                                                            as<uvec>(l["group_pop"]), as<uvec>(l["total_pop"]));
                              });

    log_tgt += add_constraint("compet", constraints, districts, psi_vec,
                              [&] (List l, int distr) -> double {
                                  uvec dvote = l["dvote"];
                                  uvec total = dvote + as<uvec>(l["rvote"]);
                                  return eval_grp_pow(plan, distr, dvote, total, 0.5, 0.5,
                                                      as<double>(l["pow"]));
                              });

    log_tgt += add_constraint("status_quo", constraints, districts, psi_vec,
                              [&] (List l, int distr) -> double {
                                  return eval_sq_entropy(plan, as<uvec>(l["current"]), distr,
                                                         pop, n_distr, as<int>(l["n_current"]), V);
                              });

    log_tgt += add_constraint("incumbency", constraints, districts, psi_vec,
                              [&] (List l, int distr) -> double {
                                  return eval_inc(plan, distr, as<uvec>(l["incumbents"]));
                              });


    log_tgt += add_constraint("polsby", constraints, districts, psi_vec,
                              [&] (List l, int distr) -> double {
                                  return eval_polsby(plan, distr, as<ivec>(l["from"]),
                                                     as<ivec>(l["to"]), as<vec>(l["area"]),
                                                     as<vec>(l["perimeter"]));
                              });

    log_tgt += add_constraint("fry_hold", constraints, districts, psi_vec,
                              [&] (List l, int distr) -> double {
                                  return eval_fry_hold(plan, distr, as<uvec>(l["total_pop"]),
                                                       as<mat>(l["ssdmat"]),
                                                       as<double>(l["denominator"]));
                              });

    log_tgt += add_constraint("log_st", constraints, districts, psi_vec,
                              [&] (List l, int distr) -> double {
                                  return eval_log_st(plan, g,
                                                     as<uvec>(l["admin"]),
                                                     n_distr) / n_consider;
                              });

    log_tgt += add_constraint("edges_removed", constraints, districts, psi_vec,
                              [&] (List l, int distr) -> double {
                                  return eval_er(plan, g, n_distr) / n_consider;
                              });

    log_tgt += add_constraint("qps", constraints, districts, psi_vec,
                              [&] (List l, int distr) -> double {
                                  return eval_qps(plan, distr, as<uvec>(l["total_pop"]),
                                                  as<uvec>(l["cities"]), as<int>(l["n_city"]),
                                                  n_distr);
                              });

    log_tgt += add_constraint("custom", constraints, districts, psi_vec,
                              [&] (List l, int distr) -> double {
                                  Function fn = l["fn"];
                                  return as<NumericVector>(fn(plan, distr))[0];
                              });


    return log_tgt;
}
