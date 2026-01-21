# data
n <- 2520L
n <- 840
x <- tdaunif::sample_circle(n = n, sd = .1)
y <- tdaunif::sample_torus_flat(n = n, ar = 1.5, sd = .1)
z <- tdaunif::sample_projective_plane(n = n, sd = .1)
readLines("sandbox/o3_4096.txt") |>
  head(n = n) |> 
  sapply(strsplit, split = "\t") |>
  unname() |>
  lapply(as.double) |>
  do.call(what = rbind) |> 
  tdaunif::add_noise(sd = .1) ->
  w
save(x,y,z,w,file="sandbox/sample.rda")

# test on v 0.5.0
pak::pkg_install("tdaverse/ripserr", ask = FALSE)
.rs.restartR()

xt0 <- bench::bench_time(xd0 <- ripserr::vietoris_rips(x))
yt0 <- bench::bench_time(yd0 <- ripserr::vietoris_rips(y))
zt0 <- bench::bench_time(zd0 <- ripserr::vietoris_rips(z))

# patch
devtools::install()
.rs.restartR()

xt9 <- bench::bench_time(xd9 <- ripserr::vietoris_rips(x))
yt9 <- bench::bench_time(yd9 <- ripserr::vietoris_rips(y))
zt9 <- bench::bench_time(zd9 <- ripserr::vietoris_rips(z))

# # test on v 1.0.0
# pak::pkg_install("tdaverse/ripserr@ripser-upgrade", ask = FALSE)
# .rs.restartR()
# 
# at1 <- bench::bench_time(ad1 <- ripserr::vietoris_rips(a))
# bt1 <- bench::bench_time(bd1 <- ripserr::vietoris_rips(b))
# ct1 <- bench::bench_time(cd1 <- ripserr::vietoris_rips(c))
# xt1 <- bench::bench_time(xd1 <- ripserr::vietoris_rips(x))
# yt1 <- bench::bench_time(yd1 <- ripserr::vietoris_rips(y))
# zt1 <- bench::bench_time(zd1 <- ripserr::vietoris_rips(z))

# compare output
xd0 |> as.data.frame() |> subset(subset = dimension == 1) |> head()
xd9 |> as.data.frame() |> subset(subset = dimension == 1) |> head()
xd0 |> as.data.frame() |> subset(subset = dimension == 1) |> tail()
xd9 |> as.data.frame() |> subset(subset = dimension == 1) |> tail()
range(as.matrix(as.data.frame(xd0)) -
        as.matrix(subset(as.data.frame(xd9), subset = is.finite(death))))
yd0 |> as.data.frame() |> subset(subset = dimension == 1) |> head()
yd9 |> as.data.frame() |> subset(subset = dimension == 1) |> head()
range(as.matrix(as.data.frame(yd0)) -
        as.matrix(subset(as.data.frame(yd9), subset = is.finite(death))))

# compare runtimes
xt0
xt9
yt0
yt9
zt0
zt9

stop()

# v0.1.1
load("sandbox/sample.rda")
m0 <- bench::mark(
  ripserr:::ripser_cpp_dist(dist(w), dim = 2L, thresh = -1, p = 2L),
  ripserr:::ripser_cpp_dist(dist(x), dim = 2L, thresh = -1, p = 2L),
  ripserr:::ripser_cpp_dist(dist(y), dim = 2L, thresh = -1, p = 2L),
  ripserr:::ripser_cpp_dist(dist(z), dim = 2L, thresh = -1, p = 2L),
  check = FALSE, iterations = 3L
)
save(m0, file = "sandbox/sample-v0.rda")

# v1.0.0
load("sandbox/sample.rda")
m1 <- bench::mark(
  ripserr:::ripser_cpp_dist(dist(w), dim = 2L, thresh = Inf, ratio = 1, p = 2L),
  ripserr:::ripser_cpp_dist(dist(x), dim = 2L, thresh = Inf, ratio = 1, p = 2L),
  ripserr:::ripser_cpp_dist(dist(y), dim = 2L, thresh = Inf, ratio = 1, p = 2L),
  ripserr:::ripser_cpp_dist(dist(z), dim = 2L, thresh = Inf, ratio = 1, p = 2L),
  check = FALSE, iterations = 3L
)
save(m1, file = "sandbox/sample-v1.rda")
