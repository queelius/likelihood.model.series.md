#' Function to evaluate performance
#' @export
eval_perf <- function(estimates, true_params)
{
    params.est <- sapply(estimates, point)
    # Calculate bias
    bias <- rowMeans(params.est - true_params)
    # Calculate mean squared error (MSE)
    mse <- mean((params.est - true_params)^2)
    # Calculate coverage probability
    ci_list <- lapply(estimates, confint)
    count <- 0L
    for (ci in ci_list)
    {
        if (ci[, 1] <= true_params && ci[, 2] >= true_params)
            count <- count + 1L
    }
    list(bias = matrix(bias), mse = mse, coverage = count / length(estimates))
}

#' Function to aggregate results
#' @export
aggregate_results <- function(results) {
    results_summary <- results %>%
        group_by(n, p) %>%
        summarise(
            bias = mean(bias),
            mse = mean(mse),
            coverage = mean(coverage),
            .groups = "drop"
        )

    return(results_summary)
}


#' Function to plot results
#' @export
plot_results <- function(results_summary) {
    library(ggplot2)

    # Plot bias
    ggplot(results_summary, aes(x = n, y = bias, color = factor(p))) +
        geom_line() +
        labs(x = "Sample size", y = "Bias", color = "Probability p") +
        theme_minimal() +
        ggtitle("Bias vs. Sample Size")

    # Plot MSE
    ggplot(results_summary, aes(x = n, y = mse, color = factor(p))) +
        geom_line() +
        labs(x = "Sample size", y = "Mean Squared Error", color = "Probability p") +
        theme_minimal() +
        ggtitle("Mean Squared Error vs. Sample Size")

    # Plot coverage probability
    ggplot(results_summary, aes(x = n, y = coverage, color = factor(p))) +
        geom_line() +
        labs(x = "Sample size", y = "Coverage Probability", color = "Probability p") +
        theme_minimal() +
        ggtitle("Coverage Probability vs. Sample Size")
}
