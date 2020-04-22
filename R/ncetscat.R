#' Eliminating the Unmeasured Confounders and Estimating Causal Effect
#' for Categorical Outcome.
#'
#' @param data an optional data frame containing the variables in the model.
#' @param x1_name the name of pre-outcome exposure
#' @param x3_name the name of post-outcome exposure
#' @param y_name the name of outcome
#' @param boots_no the number of bootstrap
#'
#'
#' @return  \code{coefficient} the casual effect coefficients
#' and its 95%CI for categorical variables.

#' @examples
#' u1 <- rbinom(1000,1,0.5)
#' px1 <- exp(0.8*u1)/(1+exp(0.8*u1))
#' x1 <- apply(matrix(px1,nrow=1),2,rbinom,n=1,size=1)
#' px3 <- exp(0.8*u1)/(1+exp(0.8*u1))
#' x3 <- apply(matrix(px3,nrow=1),2,rbinom,n=1,size=1)
#' py<- 0.2*x1+0.4*u1
#' y<- apply(matrix(py,nrow=1),2,rbinom,n=1,size=1)
#'
#' data <- data.frame(x1,x3,y)
#' model <- ncets_cat ( data = data,x1_name = 'x1',
#' x3_name = 'x3', y_name ='y',boots_no = 1000)
#' model
#'



ncets_cat <- function(data = data, x1_name = "x1", x3_name = "x3",
    y_name = "y", boots_no = 1000) {

    data <- na.omit(data)
    new_data <- data[, c(y_name, x1_name, x3_name)]
    colnames(new_data) <- c("y", "x1", "x3")

    cat <- unique(new_data$x1)
    cat <- cat[order(cat)]
    n_length <- length(cat) - 1
    comdata <- min(cat)

    all_result_nce <- NULL
    for (i in 1:n_length) {

        nce_result <- model_nce(new_data, cat[i + 1],
            comdata)

        result_bb <- rep(NA, boots_no)
        for (j in 1:boots_no) {
            chou <- sample(1:dim(new_data)[1], dim(new_data)[1],
                replace = T)
            boost_sample <- new_data[chou, ]
            result_bb[j] <- model_nce(boost_sample, cat[i +
                1], comdata)
        }

        qua1 <- quantile(result_bb, probs = c(0.025,
            0.975))

        result_all <- c(nce_result, qua1)
        result_all <- t(data.frame(result_all))
        rownames(result_all) <- paste0(cat[i + 1], " v.s.",
            comdata)

        all_result_nce <- rbind(all_result_nce, result_all)

    }

    colnames(all_result_nce) <- c("Estimation", "2.5%",
        "97.5%")


    result_end <- list()
    result_end$coefficient <- all_result_nce

    return(result_end)
}


model_nce <- function(data, big, small) {
    p1 <- data[data$x1 == big & data$x3 == small, ]
    coef1pav <- sum(p1$y) / dim(p1)[1]

    p2 <- data[data$x1 == small & data$x3 == big, ]
    coef2pav <- sum(p2$y) / dim(p2)[1]

    effectpav <- coef1pav - coef2pav

    return(effectpav)
}
