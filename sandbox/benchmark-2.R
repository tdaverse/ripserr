library(tidyverse)

load("sandbox/sample.rda")

# plot samples
list(S1 = x, T2 = y, RP2 = z, O3 = w) |> 
  enframe(name = "manifold", value = "sample") |> 
  print() -> bench_samples
for (i in seq(nrow(bench_samples))) {
  len <- if (i == 4L) 1440 else 480
  jpeg(
    paste0("sandbox/benchmark-sample-", bench_samples$manifold[[i]], ".jpg"),
    width = len, height = len
  )
  pairs(bench_samples$sample[[i]], asp = 1, pch = 19, cex = .5)
  dev.off()
}

load("sandbox/sample-v0.rda")
load("sandbox/sample-v1.rda")

# plot runtimes
bind_rows(mutate(m0, version = "0.1.1"), mutate(m1, version = "1.0.0")) |> 
  bind_cols(manifold = rep(c("O3", "S1", "T2", "RP2"), times = 2)) |> 
  mutate(across(c(version, manifold), fct_inorder)) |> 
  mutate(manifold = fct_relevel(manifold, "O3", after = Inf)) |> 
  select(version, manifold, time) |> 
  unnest(time) |> 
  ggplot(aes(x = time, y = manifold, color = version)) +
  geom_vline(xintercept = 0) +
  geom_jitter(width = 0, height = .1) +
  scale_x_continuous(
    limits = c(0, NA),
    labels = scales::unit_format(unit = "s")
  ) +
  scale_y_discrete(limits = rev) +
  scale_color_brewer(type = "qual", palette = 2) +
  labs(x = "Runtime", y = "Manifold", color = "Version") ->
  ripserr_benchmark_plot
ggsave(
  "sandbox/benchmark-plot.jpg", ripserr_benchmark_plot,
  width = 8, height = 3
)
