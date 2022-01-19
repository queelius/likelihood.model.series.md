library(masked.data)
library(tidyverse)
library(ggplot2)
library(plotly)

lambda <- c(1,1.5,.75)
kappa <- c(2,1.5,2.5)
theta <- c(lambda,kappa)

kloglik <- md_kloglike_lomax_series_m0_ref(lomax_series_data_2)
profile2 <- function(x,y) { kloglik(c(x,y,lambda[3],kappa[1],kappa[2],kappa[3])) }

data <- dplyr::tibble(lam1=seq(lambda[1]-.1,lambda[1]+.1,.05),
                      lam2=seq(lambda[2]-.1,lambda[2]+.1,.05),
                      lam3=seq(lambda[3]-.1,lambda[3]+.1,.05),
                      kap1=seq(kappa[1]-.1,kappa[1]+.1,.05),
                      kap2=seq(kappa[2]-.1,kappa[2]+.1,.05),
                      kap3=seq(kappa[3]-.1,kappa[3]+.1,.05)) %>% purrr::cross_df()

n <- nrow(data)
data$z <- numeric(n)
for (i in 1:n)
    data$z[i] <- kloglik(data[i,])


# define the drawing grid
# create all possible combinations of x and y
data <- dplyr::tibble(x=seq(theta[1]-1,theta[1]+1,0.025),
                      y=seq(theta[2]-1,theta[2]+1,0.025)) %>% purrr::cross_df()
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

mle_text <- paste("𝜃 = (",toString(round(mle,2)),")'",sep="")
data %>% ggplot(aes(x=x,y=y)) + geom_line() + geom_point(aes(x=mle[1],profile1(mle[1])),color="blue") +
    annotate("text",label=mle_text,x=mle[1],y=100+profile1(mle[1]))


