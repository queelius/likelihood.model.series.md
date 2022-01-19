library(masked.data)
library(tidyverse)
library(ggplot2)

mle <- point(md_mle_exp_series_m0(exp_series_data_1))
kloglik <- md_kloglike_exp_series_m0(exp_series_data_1)
profile <- function(x,y) { loglik(c(x,y,mle[3])) }

# define the drawing grid
# create all possible combinations of x and y
inputs <- dplyr::tibble(x=seq(mle[1]-1,mle[1]+1,0.05),y=seq(mle[2]-1,mle[2]+1,0.05)) %>% purrr::cross_df()
n <- nrow(inputs)
surf <- inputs
surf$z <- numeric(n)
for (i in 1:n)
    surf$z[i] <- profile(surf$x[i],surf$y[i])
surf %>% ggplot(aes(x=x,y=y,z=z)) + geom_contour()


plot_ly(x=surf$x,y=surf$y,z=surf$z,type="scatter3d",mode="markers",color=surf$z)
