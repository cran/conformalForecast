## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  warning = FALSE,
  message = FALSE,
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(conformalForecast)
library(forecast)
library(ggplot2)
library(dplyr)
library(tibble)
library(tsibble)

## ----data, fig.width = 7, fig.height = 3--------------------------------------
set.seed(0)
series <- arima.sim(n = 1000, list(ar = c(0.8, -0.5)), sd = sqrt(1))
autoplot(series) +
  labs(
    title = "Time series generated from an AR(2) model",
    ylab = ""
  ) +
  theme_bw()

## ----cv-----------------------------------------------------------------------
far2 <- function(x, h, level) {
  Arima(x, order = c(2, 0, 0)) |> forecast(h = h, level)
}
fc <- cvforecast(series, forecastfun = far2, h = 3, level = c(80, 95),
                 forward = TRUE, window = 100, initial = 1)
summary(fc)

## ----cvplot, fig.width = 7, fig.height = 3------------------------------------
fc |>
  autoplot() +
  labs(
    title = "Forecasts produced using an AR(2) model",
    ylab = ""
  ) +
  theme_bw()

## ----cvinfo-------------------------------------------------------------------
(fc_score <- accuracy(fc, byhorizon = TRUE))
(fc_cov <- coverage(fc, window = 100, level = 95))
(fc_wid <- width(fc, window = 100, level = 95, includemedian = TRUE))

## ----scp----------------------------------------------------------------------
scpfc <- scp(fc, symmetric = FALSE, ncal = 100, rolling = TRUE,
             weightfun = NULL, kess = FALSE, quantiletype = 1)

(scpfc_score <- accuracy(scpfc, byhorizon = TRUE))
(scpfc_cov <- coverage(scpfc, window = 100, level = 95))
(scpfc_wid <- width(scpfc, window = 100, level = 95, includemedian = TRUE))

## ----scp-exp------------------------------------------------------------------
expweight <- function(n) 0.99^{n+1-(1:n)}
scpfc_exp <- scp(fc, symmetric = FALSE, ncal = 100, rolling = TRUE,
                 weightfun = expweight, kess = FALSE, quantiletype = 1)

(scpfc_exp_score <- accuracy(scpfc_exp, byhorizon = TRUE))
(scpfc_exp_cov <- coverage(scpfc_exp, window = 100, level = 95))
(scpfc_exp_wid <- width(scpfc_exp, window = 100, level = 95, includemedian = TRUE))

## ----acp----------------------------------------------------------------------
acpfc <- acp(fc, symmetric = FALSE, gamma = 0.005, ncal = 100, rolling = TRUE)

(acpfc_score <- accuracy(acpfc, byhorizon = TRUE))
(acpfc_cov <- coverage(acpfc, window = 100, level = 95))
(acpfc_wid <- width(acpfc, window = 100, level = 95, includemedian = TRUE))

## ----pid-setup----------------------------------------------------------------
# PID setup
Tg <- 1000; delta <- 0.01
Csat <- 2 / pi * (ceiling(log(Tg) * delta) - 1 / log(Tg))
KI <- 2
lr <- 0.1

## ----pid-nsf------------------------------------------------------------------
# PID without scorecaster
pidfc_nsf <- pid(fc, symmetric = FALSE, ncal = 100, rolling = TRUE,
                 integrate = TRUE, scorecast = FALSE,
                 lr = lr, KI = KI, Csat = Csat)

(pidfc_nsf_score <- accuracy(pidfc_nsf, byhorizon = TRUE))
(pidfc_nsf_cov <- coverage(pidfc_nsf, window = 100, level = 95))
(pidfc_nsf_wid <- width(pidfc_nsf, window = 100, level = 95, includemedian = TRUE))

## ----pid----------------------------------------------------------------------
# PID with a Naive method as the scorecaster
naivefun <- function(x, h) {
  naive(x) |> forecast(h = h)
}
pidfc <- pid(fc, symmetric = FALSE, ncal = 100, rolling = TRUE,
             integrate = TRUE, scorecast = TRUE, scorecastfun = naivefun,
             lr = lr, KI = KI, Csat = Csat)

(pidfc_score <- accuracy(pidfc, byhorizon = TRUE))
(pidfc_cov <- coverage(pidfc, window = 100, level = 95))
(pidfc_wid <- width(pidfc, window = 100, level = 95, includemedian = TRUE))

## ----acmcp--------------------------------------------------------------------
acmcpfc <- acmcp(fc, ncal = 100, rolling = TRUE, integrate = TRUE, scorecast = TRUE,
             lr = lr, KI = KI, Csat = Csat)

(acmcpfc_score <- accuracy(acmcpfc, byhorizon = TRUE))
(acmcpfc_cov <- coverage(acmcpfc, window = 100, level = 95))
(acmcpfc_wid <- width(acmcpfc, window = 100, level = 95, includemedian = TRUE))

## ----covplot, fig.width = 7, fig.height = 5-----------------------------------
acmcpfc_cov$rollmean |>
  as_tsibble() |>
  mutate(horizon = key, coverage = value) |>
  update_tsibble(key = horizon) |>
  select(-c(key, value)) |>
  ggplot(aes(x = index, y = coverage, group = horizon)) +
  geom_line() +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "blue") +
  facet_grid(horizon~., scales = "free_y") +
  xlab("Time") +
  ylab("Rolling mean coverage for AcMCP") +
  theme_bw()

## ----widplot, fig.width = 7, fig.height = 5-----------------------------------
acmcpfc_wid$rollmean |>
  as_tsibble() |>
  mutate(horizon = key, width = value) |>
  update_tsibble(key = horizon) |>
  select(-c(key, value)) |>
  ggplot(aes(x = index, y = width, group = horizon)) +
  geom_line() +
  facet_grid(horizon~., scales = "free_y") +
  xlab("Time") +
  ylab("Rolling mean width for AcMCP") +
  theme_bw()

## ----bind-cov, fig.width = 7, fig.height = 5----------------------------------
candidates <- c("fc", "scpfc", "scpfc_exp", "acpfc", "pidfc_nsf", "pidfc", "acmcpfc")
methods <- c("AR", "SCP", "WCP", "ACP", "PI", "PID", "AcMCP")
for (i in 1:length(candidates)) {
  out <- get(paste0(candidates[i], "_cov"))
  out_pivot <- out$rollmean |>
    as_tsibble() |>
    mutate(horizon = key, coverage = value) |>
    update_tsibble(key = horizon) |>
    select(-c(key, value)) |>
    mutate(method = methods[i]) |>
    as_tibble()
  assign(paste0(methods[i], "_cov"), out_pivot)
}
cov <- bind_rows(mget(paste0(methods, "_cov")))

cols <- c(
  "AR" = "black",
  "SCP" = "yellow",
  "WCP" = "#fa9200",
  "ACP" = "green",
  "PI" = "blue",
  "PID" = "purple",
  "AcMCP" = "red"
)
cov |>
  as_tsibble(index = index, key = c(horizon, method)) |>
  mutate(method = factor(method, levels = methods)) |>
  ggplot(aes(x = index, y = coverage, group = method, colour = method)) +
  geom_line(size = 0.8, alpha = 0.8) +
  scale_colour_manual(values = cols) +
  geom_hline(yintercept = 0.95, linetype = "dashed", colour = "gray") +
  facet_grid(horizon~.) +
  xlab("Time") +
  ylab("Rolling mean coverage") +
  theme_bw()

## ----bind-covmean-------------------------------------------------------------
cov_mean <- lapply(1:length(candidates), function(i) {
  out_cov <- get(paste0(candidates[i], "_cov"))
  out_score <- get(paste0(candidates[i], "_score"))
  out_mean <- data.frame(
      method = methods[i],
      covmean = as.vector(out_cov$mean),
      winkler = as.vector(out_score[, "Winkler_95"]),
      msis = as.vector(out_score[,"MSIS_95"])
    ) |>
    as_tibble() |>
    rownames_to_column("horizon") |>
    mutate(horizon = paste0("h=", horizon))
  out_mean
})
cov_mean <- do.call(bind_rows, cov_mean) |>
  mutate(method = factor(method, levels = methods)) |>
  mutate(covdiff = covmean - 0.95) |>
  arrange(horizon, method)
print(cov_mean, n = nrow(cov_mean))

## ----bind-wid, fig.width = 7, fig.height = 5----------------------------------
for (i in 1:length(candidates)) {
  out <- get(paste0(candidates[i], "_wid"))
  out_pivot <- out$rollmean |>
    as_tsibble() |>
    mutate(horizon = key, width = value) |>
    update_tsibble(key = horizon) |>
    select(-c(key, value)) |>
    mutate(method = methods[i]) |>
    as_tibble()
  assign(paste0(methods[i], "_wid"), out_pivot)
}
wid <- bind_rows(mget(paste0(methods, "_wid")))

wid |>
  as_tsibble(index = index, key = c(horizon, method)) |>
  mutate(method = factor(method, levels = methods)) |>
  ggplot(aes(x = index, y = width, group = method, colour = method)) +
  geom_line(size = 0.8, alpha = 0.8) +
  scale_colour_manual(values = cols) +
  facet_grid(horizon~.) +
  xlab("Time") +
  ylab("Rolling mean width") +
  theme_bw()

