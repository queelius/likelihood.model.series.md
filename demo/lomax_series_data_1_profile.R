library(masked.data)
library(tidyverse)
library(ggplot2)
library(plotly)

#est <- md_mle_exp_series_m0(exp_series_data_2)
theta <- c(3,4,5,2,3,4)
kloglik <- md_kloglike_lomax_series_m0_ref(lomax_series_data_1)
profile2 <- function(x,y) { kloglik(c(x,y,theta[3],theta[4],theta[5],theta[6])) }

# define the drawing grid
# create all possible combinations of x and y
data <- dplyr::tibble(x=seq(theta[1]-2,theta[1]+2,0.1),
                      y=seq(theta[2]-2,theta[2]+2,0.1)) %>% purrr::cross_df()
n <- nrow(data)
data$z <- numeric(n)
for (i in 1:n)
    data$z[i] <- profile2(data$x[i],data$y[i])


data %>% ggplot(aes(x=x,y=y,z=z)) + geom_contour()

plot_ly(x=data$x,y=data$y,z=data$z,type="scatter3d",mode="markers") %>%
    add_trace(x=mle[1],y=mle[2],z=est$max_kloglike)

plot_ly(x=data$x,
        y=data$y,
        z=data$z,
        type="contour",
        contours=list(coloring='heatmap')) %>% colorbar(title="log-likelihood")

profile1 <- function(x) { kloglik(c(x,theta[2],theta[3])) }

data <- dplyr::tibble(x=seq(0,mle[1]+4,0.01))
n <- nrow(data)
data$y <- numeric(n)
for (i in 1:n)
    data$y[i] <- profile1(data$x[i])

data %>% ggplot(aes(x=x,y=y)) + geom_line() + geom_point(aes(x=mle[1],profile1(mle[1])),color="blue") +
    geom_text(label=paste("(",toString(round(mle,2)),")",sep=""),aes(x=mle[1],y=100+profile1(mle[1])))

mle_text <- paste("theta = (",toString(round(mle,2)),")'",sep="")
data %>% ggplot(aes(x=x,y=y)) + geom_line() + geom_point(aes(x=mle[1],profile1(mle[1])),color="blue") +
    annotate("text",label=mle_text,x=mle[1],y=100+profile1(mle[1]))


