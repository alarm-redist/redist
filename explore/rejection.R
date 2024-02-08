# goal: unbiasedly estimate mean of geometric distribution
library(purrr)

N = 2e4

est_mle <- function(x, n) n / sum(x)
est_mle_db <- function(x, n) {
    phat = n / sum(x)
    phat = phat*(1 - (1-phat) / n)
    phat
}
rgeo <- function(n, p=0.25) 1 + rgeom(n, p)

rN <- function(n, p=0.25) 1 + rgeom(n, p)
excN <- function(q, p=0.25) pgeom(q - 2, p, lower.tail=FALSE)

ests_mle <- map_dbl(rep(2, N), \(n) est_mle(rgeo(n), n))
ests_mle_db <- map_dbl(rep(2, N), \(n) est_mle_db(rgeo(n), n))
# ests_mle <- map_dbl(rN(N), \(n) est_mle(rgeo(n), n))
summary(ests_mle)
summary(ests_mle_db)
t.test(ests_mle, mu=0.25)
t.test(ests_mle_db, mu=0.25)

ests_db <- map_dbl(rN(N), \(n) {
    x = rgeo(n)
    s = c(0, seq_len(n) / cumsum(x)) # starting at 0
    sum(diff(s) / excN(1:n))
})
summary(ests_db)
t.test(ests_db, mu=0.25)
summary(ests_db[abs(ests_db) <= 4])
t.test(ests_db[abs(ests_db) <= 50], mu=0.25)
# summary(ests_db[0 <= ests_db & ests_db <= 1])
hist(ests_db)


# get some stats on USTs

# data(iowa)
# ia = redist_map(iowa, existing_plan = cd_2010, pop_tol = 0.01)
mapa = alarmdata::alarm_50state_map("NC")
V = nrow(mapa)

ks = map_int(1:100, function(i) {
    ust = sample_ust(mapa$adj, mapa$pop, attr(mapa, "pop_bounds")[1],
                     attr(mapa, "pop_bounds")[2], rep(1, V), rep(FALSE, V))
    pop_below = map_dbl(seq_along(ust), ~ tree_pop(ust, . - 1, mapa$pop, rep(0, V), rep(0, V)))
    devs = abs(pop_below / get_target(mapa) - 1)
    sum(devs <= 0.005)
}, .progress=TRUE)
plot(table(ks))
mean(ks > 0)

# prob of any possible edge vs prob of selecting from top k, k=1..10
map_dbl(1:10, ~ mean(ks > 0) / mean(pmin(ks / ., 1))) |> plot(type='b')
plot(table(round(ks/7, 2)))

# rough prob of successful split
p = mean(ks / 7)
# hierarchical splitting process
splits = list(14, c(7, 7), rep(c(3, 4), 2), rep(c(2, 1, 2, 2), 2)) |>
    map_int(~ sum(. > 1))

# old way: n-1 splits, each taking 1/p USTs to succeed
ps = 0.5 * p^seq(1, by=0.1, length.out=13)
# 13 * (1/p)
sum(1/ps)
# new way (a): wait til a successful split for all sibling splits
sum(splits/ps[seq_along(splits)]^splits)
# new way (b): optimally determine the number of parallel tries per set of siblings
map2_dbl(splits, ps[seq_along(splits)], \(k, p) {
    res = optimize(\(n) k * n / (1 - (1-p)^n)^k, c(1, 100))
    n = round(res$minimum)
    print(n)
    k * n / (1 - (1-p)^n)^k
}) |>
    sum()
