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

## ----sim-data, fig.width = 7, fig.height = 3----------------------------------
set.seed(101)
n0 <- 200
n_batch <- 20
n_batches <- 5
total_n <- n0 + n_batch * n_batches

series_full <- arima.sim(n = total_n, list(ar = c(0.8, -0.5)), sd = sqrt(1))
old_series <- ts(series_full[1:n0])
batches <- split(
  as.numeric(series_full[(n0 + 1):total_n]),
  rep(seq_len(n_batches), each = n_batch)
)

autoplot(series_full) +
  geom_vline(xintercept = n0 + n_batch * (0:n_batches), linetype = "dashed",
             colour = "grey60") +
  labs(
    title = "Simulated AR(2) series arriving in five batches",
    subtitle = "Dashed lines mark batch boundaries; the first segment is the initial history",
    x = "Time", y = NULL
  ) +
  theme_bw()

## ----forecastfun--------------------------------------------------------------
far2 <- function(x, h, level) {
  Arima(x, order = c(2, 0, 0)) |> forecast(h = h, level = level)
}
naivefun <- function(x, h) {
  naive(x) |> forecast(h = h)
}

## ----baseline-----------------------------------------------------------------
h <- 3
level <- 95
window_len <- 50
ncal <- 40
Tg <- 200
KI <- 2

fc0 <- cvforecast(old_series, forecastfun = far2, h = h, level = level,
                  window = window_len, initial = 1)
cps_list <- NULL
cps_list[["scp"]] <- conformal(fc0, method = "scp", ncal = ncal, rolling = TRUE)
cps_list[["acp"]] <- conformal(fc0, method = "acp", ncal = ncal, rolling = TRUE)
cps_list[["pid"]] <- conformal(fc0, method = "pid", ncal = ncal, rolling = TRUE, 
                               scorecastfun = naivefun, Tg = Tg, KI = KI)
cps_list[["acmcp"]] <- conformal(fc0, method = "acmcp", ncal = ncal, rolling = TRUE, 
                                 KI = KI, Tg = Tg)
scp_initial <- cps_list$scp

## ----update-loop--------------------------------------------------------------
timings_update <- matrix(
  0, nrow = n_batches, ncol = length(cps_list),
  dimnames = list(paste("batch", seq_len(n_batches)), names(cps_list))
)

for (b in seq_len(n_batches)) {
  for (m in names(cps_list)) {
    timings_update[b, m] <- system.time(
      cps_list[[m]] <- update(cps_list[[m]], new_data = batches[[b]], forecastfun = far2)
    )["elapsed"]
  }
}
round(timings_update, 3)

## ----refit-loop---------------------------------------------------------------
current_series <- old_series
timings_full <- matrix(
  0, nrow = n_batches, ncol = length(cps_list) + 1L,
  dimnames = list(paste("batch", seq_len(n_batches)),
                  c("cvforecast", names(cps_list)))
)

for (b in seq_len(n_batches)) {
  current_series <- ts(c(current_series, batches[[b]]))
  timings_full[b, "cvforecast"] <- system.time(
    fc_full <- cvforecast(current_series, forecastfun = far2, h = h, level = level,
                          window = window_len, initial = 1)
  )["elapsed"]
  timings_full[b, "scp"] <- system.time(
    conformal(fc_full, method = "scp", ncal = ncal, rolling = TRUE)
  )["elapsed"]
  timings_full[b, "acp"] <- system.time(
    conformal(fc_full, method = "acp", ncal = ncal, rolling = TRUE)
  )["elapsed"]
  timings_full[b, "pid"] <- system.time(
    conformal(fc_full, method = "pid", ncal = ncal, rolling = TRUE, 
                               scorecastfun = naivefun, Tg = Tg, KI = KI)
  )["elapsed"]
  timings_full[b, "acmcp"] <- system.time(
    conformal(fc_full, method = "acmcp", ncal = ncal, rolling = TRUE, 
                                 KI = KI, Tg = Tg)
  )["elapsed"]
}
round(timings_full, 3)

## ----xreg-setup---------------------------------------------------------------
set.seed(202)
n_x <- 60
xreg_full <- matrix(rnorm(n_x + 10), ncol = 1)
y_x <- ts(as.numeric(arima.sim(n = n_x, list(ar = 0.5))) + 2 * xreg_full[1:n_x, 1])

lm_fun <- function(x, h, level, xreg, newxreg) {
  Arima(x, order = c(1, 0, 0), xreg = xreg) |>
    forecast(h = h, level = level, xreg = newxreg)
}

fc_x <- cvforecast(y_x, forecastfun = lm_fun, h = 2, level = 95, window = 30,
                   initial = 1, xreg = xreg_full[1:(n_x + 2), , drop = FALSE])
cp_x <- scp(fc_x, ncal = 20, rolling = TRUE)

## ----xreg-update--------------------------------------------------------------
new_y <- as.numeric(arima.sim(n = 5, list(ar = 0.5))) + 2 * xreg_full[(n_x + 1):(n_x + 5), 1]
new_xreg <- xreg_full[(n_x + 3):(n_x + 7), , drop = FALSE]

cp_x_upd <- update(cp_x, new_data = new_y, forecastfun = lm_fun, new_xreg = new_xreg)

## ----forward-refresh----------------------------------------------------------
scp_initial$mean
scp_upd <- update(scp_initial, new_data = batches[[1]], forecastfun = far2)
scp_upd$mean

