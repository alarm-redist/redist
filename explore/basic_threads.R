library(Rcpp)

Rcpp::cppFunction("
void func(int n) {
    RcppThread::ProgressBar bar(n, 1);
    RcppThread::parallelFor(0, n, [&] (int i) {
        std::this_thread::sleep_for(std::chrono::milliseconds(100));
        bar++;
    });
}
", depends="RcppThread", plugins="cpp11")

cat("start\n")
func(500)
cat("done\n")
