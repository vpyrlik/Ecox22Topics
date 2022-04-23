#' Продивнутая эконометрика (апрель-июнь 2022)
#' Компьютерная практика - 1
#' 
#'                       ------------------------- 
#'                      |       26.04.2022        |
#'                      |                         |
#'                      |  ПЕРВОЕ ПРИМЕНЕНИЕ ОММ  |
#'                       ------------------------- 

# УСТАНОВКА ПАКЕТОВ ===========================================================
# instal.packages(
#   c(
#     "ggplot2",
#     "dplyr",
#     "tidyr",
#     "lubridate",
#     "stats"
#   )
# )

# ПАКЕТЫ ======================================================================
#' для графики
library(dplyr)
library(tidyr)
library(ggplot2)
library(lubridate)

# ОПИСАНИЕ МОДЕЛИ =============================================================
#' исходная модель:
#' x ~ N(mu, sigma2)
#' условия на моменты:
#' E(x-mu) = 0
#' E(x-mu)^2 = sigma2
#' E(x-mu)^3 =0
#' E(x-mu)^4 = 3sigma2^2

#' параметры:
#' q = (m, s); m - аналог mu, s - аналог sigma2

# ДАННЫЕ: часть 1 - симуляция =================================================
mu <- -1.2
sigma2 <- 2.7
n <- c(8, 20, 50, 250, 10^4)
X <- sapply(n, function(i) rnorm(i, mu, sqrt(sigma2)))

# ФУНКЦИИ =====================================================================
#' функция условий на моменты
f <- function(x,q) {
  m <- q[1]
  s <- q[2]
  c(
    x-m,
    (x-m)^2-s,
    (x-m)^3,
    (x-m)^4-3*s^2
  )
}

#' выборочное среднее условий на моменты
fbar <- function(dta,q) {
  rowMeans(sapply(dta, function(x) f(x,q)))
}

#' целевая функция ОММ
Sn <- function(dta,q,Wn) {
  fbar_vec <- fbar(dta,q)
  (matrix(fbar_vec, nrow=1) %*% Wn %*% matrix(fbar_vec, ncol=1))[1,1]
}

#' целевая функция первого шага эффективного ОММ
SnI <- function(dta,q){Sn(dta,q,diag(c(1,1,1,1)))}

# ПРИНЦИП ОММ - ГРАФИЧЕСКАЯ ИЛЛЮСТРАЦИЯ =======================================
#' сетка для значений параметров
m_grid <- seq(-2, 2, by = 0.01)
s_grid <- seq(1, 5, by = 0.01)
q_grid <- expand.grid(m_grid, s_grid)

#' анализ SnI для одной из выборок
i <- 3
SnI_on_grid <- 
  q_grid %>% 
  rename(
    m = Var1,
    s = Var2
  ) %>% 
  group_by(
    m,
    s
  ) %>% 
  mutate(
    S = SnI(X[[i]], c(m,s))
  )

SnI_on_grid %>% 
  ggplot() +
  geom_raster(
    aes(
      x = m,
      y = s,
      fill = S
    )
  ) +
  geom_contour(
    aes(
      x = m,
      y = s,
      z = S
    ),
    bins = 200,
    color = "yellow",
    alpha = 0.5
  ) +
  geom_point(
    aes(
      x = mu,
      y = sigma2
    ),
    size=2,
    color = "red"
  ) + 
  geom_point(
    aes(
      x = SnI_on_grid$m[which.min(SnI_on_grid$S)],
      y = SnI_on_grid$s[which.min(SnI_on_grid$S)]
    ),
    color = "yellow",
    size = 2
  ) + 
  annotate(
    "text",
    label = floor(1000*min(SnI_on_grid$S))/1000,
    x = SnI_on_grid$m[which.min(SnI_on_grid$S)],
    y = SnI_on_grid$s[which.min(SnI_on_grid$S)]-0.05,
    color = "yellow"
  ) +
  annotate(
    "text",
    label = floor(1000*SnI(X[[i]],c(mu,sigma2)))/1000,
    x = mu,
    y = sigma2+0.05,
    color = "red"
  )

# ЭФФЕКТИВНЫЙ (ДВУХШАГОВЫЙ) ОММ ===============================================
#' для одной из выборок
i <- 3
m_grid <- seq(-2, 0, by = 0.025)
s_grid <- seq(0.25, 3, by = 0.025)
q_grid <- expand.grid(m_grid, s_grid)
SnI_on_grid <- 
  q_grid %>% 
  rename(
    m = Var1,
    s = Var2
  ) %>% 
  group_by(
    m,
    s
  ) %>% 
  mutate(
    S = SnI(X[[i]], c(m,s))
  )

#' оценка первого шага
q1 <- c(m = SnI_on_grid$m[which.min(SnI_on_grid$S)],
        s = SnI_on_grid$s[which.min(SnI_on_grid$S)])

#' оценка оптимальной весовой матрицы
Wnopt <- solve(n[i]^-1 * crossprod(t(sapply(X[[i]],function(x) f(x,q1)))))

#' целевая функция второго шага
SnII <- function(dta,q) {Sn(dta,q,Wnopt)}

#' графическая "оптимизация"
SnII_on_grid <- 
  q_grid %>% 
  rename(
    m = Var1,
    s = Var2
  ) %>% 
  group_by(
    m,
    s
  ) %>% 
  mutate(
    S = SnII(X[[i]], c(m,s))
  )

SnII_on_grid %>% 
  ggplot() +
  geom_raster(
    aes(
      x = m,
      y = s,
      fill = S
    )
  ) +
  geom_contour(
    aes(
      x = m,
      y = s,
      z = S
    ),
    bins = 200,
    color = "cyan",
    alpha = 0.25
  ) +
  geom_point(
    aes(
      x = mu,
      y = sigma2
    ),
    size=2,
    color = "red"
  ) + 
  geom_point(
    aes(
      x = SnI_on_grid$m[which.min(SnI_on_grid$S)],
      y = SnI_on_grid$s[which.min(SnI_on_grid$S)]
    ),
    color = "yellow",
    size = 2
  ) + 
  geom_point(
    aes(
      x = SnII_on_grid$m[which.min(SnII_on_grid$S)],
      y = SnII_on_grid$s[which.min(SnII_on_grid$S)]
    ),
    color = "cyan",
    size = 2
  ) + 
  annotate(
    "text",
    label = floor(1000*min(SnI_on_grid$S))/1000,
    x = SnI_on_grid$m[which.min(SnI_on_grid$S)],
    y = SnI_on_grid$s[which.min(SnI_on_grid$S)]-0.05,
    color = "yellow"
  ) +
  annotate(
    "text",
    label = floor(1000*min(SnII_on_grid$S))/1000,
    x = SnII_on_grid$m[which.min(SnII_on_grid$S)],
    y = SnII_on_grid$s[which.min(SnII_on_grid$S)]-0.05,
    color = "cyan"
  ) +
  annotate(
    "text",
    label = floor(1000*SnI(X[[i]],c(mu,sigma2)))/1000,
    x = mu,
    y = sigma2+0.05,
    color = "red"
  )

# ДАННЫЕ: часть 2 - доходность финансового актива =============================
dta <- 
  read.csv(
    "SBERP_180101_201231.csv",
    sep = ";",
    header = TRUE
  ) %>% 
  rename(
    DATE = X.DATE.,
    TIME = X.TIME.,
    CLOSE = X.CLOSE.
  ) %>% 
  mutate(
    DATETIME = as_datetime(paste(DATE,TIME), format = "%d/%m/%y %H:%M"),
    RETURN =log(CLOSE/lag(CLOSE)),
  ) %>% 
  select(
    DATETIME,
    CLOSE,
    RETURN
  ) %>% 
  filter(
    complete.cases(.)
  )

dta <-
  dta %>% 
  mutate(
    NORMAL = sapply(RETURN, function(r) dnorm(r, mean(dta$RETURN), sd(dta$RETURN)))
  )

#' обзор
dta %>% 
  select(
    -NORMAL
  ) %>% 
  pivot_longer(
    -DATETIME,
    values_to = "VALUE",
    names_to = "SERIES"
  ) %>% 
  ggplot() +
  geom_line(
    aes(
      x = DATETIME,
      y = VALUE,
      group = SERIES
    )
  ) + 
  facet_wrap(
    .~SERIES,
    scales = "free_y"
  )

dta %>%
  ggplot() +
  geom_density(
    aes(
      x = RETURN
    ),
    size = 1.2
  ) +
  geom_line(
    aes(
      x = RETURN,
      y = NORMAL
    ),
    color = "red",
    size = 1.1,
    alpha = 0.5
  ) + 
  scale_x_continuous(
    limits = c(-0.015,0.015)
  ) +
  scale_y_continuous(
    name = "DENSITY"
  )

#dta <- dta %>% filter(DATETIME >= "2020-10-01")

## Эффективный ОММ ============================================================
# ЭФФЕКТИВНЫЙ (ДВУХШАГОВЫЙ) ОММ ===============================================
#' для одной из выборок
m_grid <- seq(-0.001, 0.001, by = 0.0001)
s_grid <- seq(10^(-6), 10^(-5), by = 10^(-7))
q_grid <- expand.grid(m_grid, s_grid)
SnI_on_grid <- 
  q_grid %>% 
  rename(
    m = Var1,
    s = Var2
  ) %>% 
  group_by(
    m,
    s
  ) %>% 
  mutate(
    S = SnI(dta$RETURN, c(m,s))
  )

#' оценка первого шага
q1 <- c(m = SnI_on_grid$m[which.min(SnI_on_grid$S)],
        s = SnI_on_grid$s[which.min(SnI_on_grid$S)])

#' оценка оптимальной весовой матрицы
Wnopt <- solve(nrow(dta)^-1 * crossprod(t(sapply(dta$RETURN,function(x) f(x,q1)))))

#' целевая функция второго шага
SnII <- function(dta,q) {Sn(dta,q,Wnopt)}

#' графическая "оптимизация"
SnII_on_grid <- 
  q_grid %>% 
  rename(
    m = Var1,
    s = Var2
  ) %>% 
  group_by(
    m,
    s
  ) %>% 
  mutate(
    S = SnII(dta$RETURN, c(m,s))
  )

SnII_on_grid %>% 
  ggplot() +
  geom_raster(
    aes(
      x = m,
      y = s,
      fill = S
    )
  ) +
  geom_contour(
    aes(
      x = m,
      y = s,
      z = S
    ),
    bins = 100,
    color = "cyan",
    alpha = 0.25
  ) +
  geom_point(
    aes(
      x = SnI_on_grid$m[which.min(SnI_on_grid$S)],
      y = SnI_on_grid$s[which.min(SnI_on_grid$S)]
    ),
    color = "yellow",
    size = 2
  ) + 
  geom_point(
    aes(
      x = SnII_on_grid$m[which.min(SnII_on_grid$S)],
      y = SnII_on_grid$s[which.min(SnII_on_grid$S)]
    ),
    color = "cyan",
    size = 2
  ) + 
  annotate(
    "text",
    label = floor(10^6*SnII_on_grid$S[which.min(SnI_on_grid$S)])/10^6,
    x = SnI_on_grid$m[which.min(SnI_on_grid$S)],
    y = SnI_on_grid$s[which.min(SnI_on_grid$S)],
    color = "yellow"
  ) +
  annotate(
    "text",
    label = floor(10^6*min(SnII_on_grid$S))/10^6,
    x = SnII_on_grid$m[which.min(SnII_on_grid$S)],
    y = SnII_on_grid$s[which.min(SnII_on_grid$S)],
    color = "cyan"
  )

## Проверка модели: J-тест ====================================================
J <- min(SnII_on_grid$S)*nrow(dta)
J_critical <- qchisq(0.95, 4-2)
P_value <- 1 - pchisq(J,4-2)



