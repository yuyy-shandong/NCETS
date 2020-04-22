#' Eliminating the Unmeasured Confounders and Estimating Causal Effect
#' for Continuous Outcome.
#'
#' @param formula an object of class 'formula' :
#' a symbolic description of the model to be fitted
#' @param data an optional data frame containing the variables in the model.
#' @param sdmethod Choosing a method for estimating variance
#' @param x1_name the name of pre-outcome exposure
#' @param x3_name the name of post-outcome exposure
#' @param boots_no the number of bootstrap
#'
#'
#' @return  \code{formula} the model
#' @return  \code{bootstrap} the casual effect coefficients
#' when the variance is calculated by bootstrap
#' @return  \code{normal} the casual effect coefficients
#' when the variance is calculated in theory
#' @return  \code{quant} the confidence interval
#' calculated by bootstrap
#' @examples
#'   u <- rnorm(1000,0,1)
#'   x1 <- 0.5*u +rnorm(1000,0,1)
#'   x3 <- 0.5*u +rnorm(1000,0,1)
#'   y <- 2*x1 + 1*u +rnorm(1000,0,1)
#'   data <- data.frame(u,x1,x3,y)
#'   model <- ncets_con (y~x1+x3, data = data, sdmethod ='normal',
#'   x1_name = 'x1',x3_name = 'x3',boots_no = NULL)
#'   model
#'


ncets_con <- function(formula, data = data, sdmethod = c("normal",
    "bootstrap", "all"), x1_name = "x1", x3_name = "x3",
    boots_no = NULL) {
    formula <- as.formula(formula)
    data <- na.omit(data)
    model_pav <- lm(formula, data = data)
    coef1_pav <- summary(model_pav)$coefficients[2, 1]
    coef2_pav <- summary(model_pav)$coefficients[3, 1]
    result1_pav <- coef1_pav - coef2_pav


    if (sdmethod == "bootstrap") {
        result_pav_b <- NULL
        for (k in 1:boots_no) {
            sample_boot <- data[sample(1:dim(data)[1],
                dim(data), replace = T), ]
            model_pav_b <- lm(formula = formula, data = sample_boot)
            pav_b_1 <- summary(model_pav_b)$coefficients[2,
                1]
            pav_b_2 <- summary(model_pav_b)$coefficients[3,
                1]

            result_pav_b[k] <- pav_b_1 - pav_b_2

        }

        sd_pav_p <- sd(result_pav_b)
        p1_b <- pnorm(result1_pav, mean = 0, sd = sd_pav_p,
            lower.tail = TRUE, log.p = FALSE)
        p2_b <- 1 - pnorm(result1_pav, mean = 0, sd = sd_pav_p,
            lower.tail = TRUE, log.p = FALSE)
        p_b <- 2 * min(p1_b, p2_b)
        conup_b <- qnorm(0.975, mean = result1_pav, sd = sd_pav_p,
            lower.tail = TRUE, log.p = FALSE)
        conlower_b <- qnorm(0.025, mean = result1_pav,
            sd = sd_pav_p, lower.tail = TRUE, log.p = FALSE)
        conup_b_q <- quantile(result_pav_b, c(0.025,
            0.975))[2]
        conlower_b_q <- quantile(result_pav_b, c(0.025,
            0.975))[1]
        result_boot <- list()
        result_boot$formula <- formula
        boos_coef <- cbind(result1_pav, sd_pav_p, p_b,
            conlower_b, conup_b)
        colnames(boos_coef) <- c("Estimation", "sd",
            "P_value", "2.5%", "97.5%")
        result_boot$bootstrap <- boos_coef
        quant <- cbind(conlower_b_q, conup_b_q)
        colnames(quant) <- c("2.5%", "97.5%")
        result_boot$quant <- quant



        return(result_boot)

    } else {
        if (sdmethod == "normal") {
            covqian <- cbind(data[, x1_name, with = F],
                data[, x3_name, with = F])
            cov <- cbind(1, covqian)
            cov <- as.matrix(cov)
            r <- qr(cov)$rank
            n <- nrow(cov)
            resivar <- sum(summary(model_pav)$residuals^2) / (n -
                r)
            cov2 <- solve(t(cov) %*% cov)
            c <- cbind(0, 1, -1)
            varco <- c %*% cov2 %*% t(c)
            sd_pav <- (resivar * varco) ^ (1 / 2)
            z_pav <- result1_pav / sd_pav

            p1_pav <- pnorm(z_pav, mean = 0, sd = 1,
                lower.tail = TRUE, log.p = FALSE)
            p2_pav <- 1 - pnorm(z_pav, mean = 0, sd = 1,
                lower.tail = TRUE, log.p = FALSE)

            p_pav <- 2 * min(p1_pav, p2_pav)

            conup_pav <- qnorm(0.975, result1_pav, sd_pav)
            conlower_pav <- qnorm(0.025, result1_pav,
                sd_pav)


            result_normal <- list()
            result_normal$formula <- formula
            coef <- cbind(result1_pav, sd_pav, z_pav,
                p_pav, conlower_pav, conup_pav)
            colnames(coef) <- c("Estimation", "sd", "Z",
                "P_value", "2.5%", "97.5%")
            result_normal$normal <- coef

            return(result_normal)

        } else {
            if (sdmethod == "all") {
                result_pav_b <- NULL
                for (k in 1:boots_no) {
                  sample_boot <- data[sample(1:dim(data)[1],
                    dim(data), replace = T), ]
                  model_pav_b <- lm(formula, data = sample_boot)
                  pav_b_1 <- summary(model_pav_b)$coefficients[2,
                    1]
                  pav_b_2 <- summary(model_pav_b)$coefficients[3,
                    1]

                  result_pav_b[k] <- pav_b_1 - pav_b_2

                }

                sd_pav_p <- sd(result_pav_b)
                p1_b <- pnorm(result1_pav, mean = 0,
                  sd = sd_pav_p, lower.tail = TRUE, log.p = FALSE)
                p2_b <- 1 - pnorm(result1_pav, mean = 0,
                  sd = sd_pav_p, lower.tail = TRUE, log.p = FALSE)
                p_b <- 2 * min(p1_b, p2_b)
                conup_b <- qnorm(0.975, mean = result1_pav,
                  sd = sd_pav_p, lower.tail = TRUE, log.p = FALSE)
                conlower_b <- qnorm(0.025, mean = result1_pav,
                  sd = sd_pav_p, lower.tail = TRUE, log.p = FALSE)
                conup_b_q <- quantile(result_pav_b, c(0.025,
                  0.975))[2]
                conlower_b_q <- quantile(result_pav_b,
                  c(0.025, 0.975))[1]
                result_boot <- list()
                result_boot$formula <- formula
                boos_coef <- cbind(result1_pav, sd_pav_p,
                  p_b, conlower_b, conup_b)
                colnames(boos_coef) <- c("Estimation",
                  "sd", "P_value", "2.5%", "97.5%")
                result_boot$coefficients <- boos_coef
                quant <- cbind(conlower_b_q, conup_b_q)
                colnames(quant) <- c("2.5%", "97.5%")
                result_boot$quant <- quant


                covqian <- cbind(data[, x1_name, with = F],
                  data[, x3_name, with = F])
                cov <- cbind(1, covqian)
                cov <- as.matrix(cov)
                r <- qr(cov)$rank
                n <- nrow(cov)
                resivar <- sum(summary(model_pav)$residuals^2) / (n -
                  r)
                cov2 <- solve(t(cov) %*% cov)
                c <- cbind(0, 1, -1)
                varco <- c %*% cov2 %*% t(c)
                sd_pav <- (resivar * varco) ^ (1 / 2)
                z_pav <- result1_pav / sd_pav

                p1_pav <- pnorm(z_pav, mean = 0, sd = 1,
                  lower.tail = TRUE, log.p = FALSE)
                p2_pav <- 1 - pnorm(z_pav, mean = 0,
                  sd = 1, lower.tail = TRUE, log.p = FALSE)

                p_pav <- 2 * min(p1_pav, p2_pav)

                conup_pav <- qnorm(0.975, result1_pav,
                  sd_pav)
                conlower_pav <- qnorm(0.025, result1_pav,
                  sd_pav)


                result_normal <- list()
                result_normal$formula <- formula
                coef <- cbind(result1_pav, sd_pav, z_pav,
                  p_pav, conlower_pav, conup_pav)
                colnames(coef) <- c("Estimation", "sd",
                  "Z", "P_value", "2.5%", "97.5%")
                result_normal$coefficients <- coef

                result_all <- list()
                result_all$formula <- formula
                result_all$bootstrap <- boos_coef
                result_all$normal <- coef
                result_all$quant <- quant


                return(result_all)
            }

        }

    }


}
