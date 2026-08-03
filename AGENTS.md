## This package

<!-- Insert package-specific content here. -->

## Analysis posture

- Correctness comes before cleverness or volume.
- Code should work on real inputs, not only toy cases or happy paths.
- Do not declare work complete while the script, report, model, figure, table, package, or workflow still fails in practice.
- Stay within the requested scope.
- When reimplementing an existing method, compare carefully against the reference and identify meaningful behavioral differences.
- Prefer readable, explicit steps over dense abstractions that hide what the work is doing.

## Package development

### Key commands

```r
# Execute code
devtools::load_all()
code

# Run tests
devtools::test()
devtools::test(filter = '^{name}')
devtools::test_active_file('R/{name}.R')
devtools::test_active_file('R/{name}.R', desc = 'blah')

# Check test coverage
devtools::test_coverage()
devtools::test_coverage_active_file('R/{name}.R')

# Update and check documentation
devtools::document()
pkgdown::check_pkgdown()

# Run R CMD check
devtools::check()
```

### Formatting and linting

```sh
air format .
jarl check .
```

- Always run `air format .` and `jarl check .` after generating or changing code.

### Running R

- Use `Rscript` to run standalone scripts.
- On Windows, do not use `Rscript -e '...'` for multiline or complex R code because it will segfault.
- For multiline or complex code, write a temporary `.R` file and run `Rscript path/to/file.R`.
- Simple one-liners such as `Rscript -e 'devtools::test()'` are acceptable.
- R 4.6.1 and Rtools 4.5 are correctly paired; there is no Rtools 4.6.

## Code style

- Follow the tidyverse style guide.
- Use idiomatic R and established, well-tested packages instead of unnecessary custom wrappers.
- Use `<-` for assignment.
- Use single quotes for strings.
- Use the native pipe, `|>`, not `%>%`.
- Use braces for every `if` statement, including single-line bodies.
- Do not assign from an `if` expression such as `x <- if (...) { ... }`.
- Use `\(x) x + 1` for single-line anonymous functions.
- Use `function(x) { ... }` for multiline anonymous functions.
- Use `snake_case` for variables and columns.
- Function names must not begin with a period.
- Use `seq_len(n)`, not `1:n`.
- Use the `L` suffix for logically integer literals, such as `1L` and `8L`.
- Use `glue::glue()` for string construction unless exact formatting calls for `sprintf()`.
- Use `glue::glue_sql()` for SQL strings.
- Use `fs` for file-system operations in packages or when path manipulation is nontrivial.
- Plain strings are sufficient for simple paths.
- Avoid deprecated tidyverse functions, including `recode()`, `transmute()`, and purrr's `_dfr` and `_dfc` variants.
- Prefer `imap()` followed by `list_rbind()` over `imap_dfr()`.
- Use `anyNA(x)`, not `any(is.na(x))`.
- Use `nlevels(x)`, not `length(levels(x))`.
- Do not use `.by` in `summarize()` or `mutate()`.
- Use explicit `group_by()`, then `.groups` in `summarize()` or `ungroup()` as appropriate.

## Functions and abstractions

- Use `match.arg()` with a character-vector default for arguments that select among named options.
- Use booleans only for flags that modify one algorithm, not to choose among algorithms.
- Follow base R naming conventions for parameters, such as `method` rather than `algorithm_mode`.
- Use a consistent verb-object naming pattern for helpers.
- Internal functions do not need pseudo-private prefixes.
- Name files for their primary function or responsibility.
- Put nontrivial shared helpers in clearly named utility files.
- Do not hide an entire analysis behind helper layers when the steps matter.
- Do not create wrappers whose main effect is to hide behavior or swallow errors.
- Remove unreachable branches and duplicate functions.

## Null, NA, and missing values

- Use base R `%||%` for `NULL` fallbacks.
- `%||%` does not coalesce `NA`; guard missing values separately.

## Messages and errors

- Use `cli::cli_warn()` and `cli::cli_inform()` instead of `warning()`, `message()`, or `cat()` for user-facing output.
- Reserve `cli::cli_warn()` for actual problems.
- Use `cli::cli_inform()` for status updates.
- Add `.frequency = 'once'` to `cli::cli_inform()` calls inside loops.
- Do not suppress errors needed to diagnose failures.

## Data, SQL, and performance

- Keep work in SQL when possible rather than materializing data in R.
- Treat `DBI::dbGetQuery()`, `collect()`, and tibble construction from SQL output as materialization points that require scrutiny.
- If R immediately groups, counts, filters, or joins rows returned by SQL, that work probably belongs in SQL.
- Avoid unnecessary materialization, copying, and conversion.
- Prefer vectorized or backend-aware operations for large data.
- Base performance changes on measurement, not speculation.
- Preserve the native data structure in use, such as tibbles, `sf` objects, database tables, or package-specific objects.

## ggplot2

- Do not put `aes()` in the initial `ggplot()` call.
- Reshape data with `pivot_*()` before plotting instead of creating awkward mappings in `aes()`.
- Use facets to compare related quantities rather than rescaling them awkwardly.
- Use `scale_*()` for axis and legend labels.
- Use `labs()` for titles and subtitles.
- Do not add subtitles or captions to figures intended for academic articles.
- Keep themes and scales minimal unless they serve a clear purpose.
- Exported package plotting functions must not set a default theme.
- Do not hardcode aesthetic values that override defaults without a specific reason.
- Do not select one-off colors without a reason.
- Extract repeated theme code into a named `theme_*()` helper.
- Do not resize individual geoms or theme elements away from defaults unless solving a specific problem.
- Do not map color, fill, or shape to a variable already represented by a facet or categorical axis.

## Testing

- Tests for `R/{name}.R` belong in `tests/testthat/test-{name}.R`.
- All new code should have an accompanying test.
- Put new tests next to similar existing tests.
- Keep tests minimal and comments sparse.
- Never put executable test code outside a `test_that()` block in a `test-{name}.R` file.
- Put shared setup code in `tests/testthat/helper.R` or a clearly named helper file.
- Avoid `expect_true()` and `expect_false()` when a more specific expectation would produce a better failure message.
- Consider newer expectations such as `expect_all_true()`, `expect_all_equal()`, and `expect_r6_class()`.
- Use `expect_error()` or `expect_warning()` only when the condition has a known class; otherwise prefer snapshots.
- Avoid the `.package` argument to `local_mocked_bindings()` because it modifies another package's namespace.
- Instead, create a mockable function in the current package.
- Tests should verify package behavior, not basic facts about R or fixture data.
- Do not put `expect_*()` calls inside loops.
- Consolidate near-duplicate tests.
- Do not specify default arguments explicitly in tests unless the default itself is under test.
- Add regression tests for real failure modes and edge cases.
- Do not add regression tests for code that has been removed.
- Use `skip_if()` when a test depends on a suggested package.
- Delete tests that add runtime without testing meaningful behavior.

## Documentation

- Every user-facing function should be exported and have roxygen2 documentation.
- Internal functions should not have roxygen2 documentation.
- Give each exported function its own roxygen block.
- Do not use `@describeIn`, shared `@name` blocks, or `@inheritParams`.
- Each exported function should have its own title, `@description`, `@param`, `@return`, and useful runnable `@examples`.
- Wrap roxygen2 comments at 80 characters.
- Do not use `\dontrun{}` or `\donttest{}` in examples.
- Use `tempfile()`, small real examples, and explicit `requireNamespace()` guards for optional packages.
- Add every new non-internal documentation topic to `_pkgdown.yml`.
- Run `devtools::document()` after changing roxygen2 comments.
- Run `pkgdown::check_pkgdown()` to verify that every topic appears in the reference index.
- Keep roxygen, pkgdown, NEWS, vignettes, and code synchronized.
- README and vignette examples should demonstrate real workflows rather than fake data or function listings.

## Package dependencies and compatibility

- Keep `DESCRIPTION` clean.
- Put packages used unconditionally in `Imports`.
- Put packages used only in tests, examples, or documentation in `Suggests`.
- Once a package is in `Imports`, remove `requireNamespace()` guards for it.
- Do not list base or recommended packages such as `stats` or `utils` without a specific reason.
- Use `devtools::check()` rather than raw commands that leave `Rcheck` directories behind.
- For packages that are not yet public, prefer correcting the API directly instead of adding compatibility shims.

## `NEWS.md`

- Add a bullet for every user-facing change.
- Do not add bullets for small documentation changes, internal refactors, or bugs introduced and fixed within the current development version.
- Briefly describe the change from the user's perspective.
- Mention a related issue in parentheses when applicable.
- Keep each bullet on one physical line.
- Put a function name early in the bullet when the change concerns that function.
- Include a GitHub username only when the pull request author is not a package author.
- Sort function-related bullets alphabetically by function name.
- Put bullets without function names first.
- In prerelease packages, avoid release-note framing such as 'legacy', 'new', 'now', or 'fixed'; describe current feature coverage instead.

## Comments

- Comments should explain why: a constraint, invariant, analytic decision, or workaround.
- Do not narrate what the code already states.
- Preserve comments that explain algorithmic reasoning or correspondence with a reference implementation.
- Do not use section headers in package code.

## Git

- Never touch git.

## Writing

- Use sentence case for headings.
- Use US English.

## Proofreading

When the user asks for proofreading:

- Act as an expert proofreader and editor focused on clear, engaging, well-structured writing.
- Work paragraph by paragraph.
- Begin with a TODO list containing one item for each top-level section.
- Fix spelling, grammar, and minor problems without asking first.
- Mark unclear, confusing, or ambiguous sentences with a `FIXME` comment.
- Report only what changed.
