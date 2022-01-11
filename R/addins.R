custom_constraint <- function() {
    if (!requireNamespace('rstudioapi', quietly = TRUE)) {
        cli::cli_abort('{.pkg rstudioapi} required to use add ins.')
    }

    rstudioapi::insertText('
    add_constr_custom(
      strength = 1,
      fn = function(plan, distr) {
        0
      }
    )')
}
