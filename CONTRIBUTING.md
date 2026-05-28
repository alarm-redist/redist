# Developer setup

Tools needed:

- [`devtools` package](https://cran.r-project.org/package=devtools)
- [`air`](https://posit-dev.github.io/air/) (installable with system package managers)
- [`clangd`](https://clangd.llvm.org/) or
  [`clang-format`](https://clang.llvm.org/docs/ClangFormat.html) if you don't
  use an IDE (both installable with system packagemanagers). 

Recommend configuring Positron/VSCode/etc. to run the stylers and linters on
save (search "Format on save" under worskspace settings).

Pre-commit workflow if not using IDE:
- Run `air format .`
- Run `clang-format -i <FILE>` for all edited C/C++ files
