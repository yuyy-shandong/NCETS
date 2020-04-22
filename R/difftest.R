#' Diff_test for test causal effect among time-varying Exposure.
#'
#' @param data an optional data frame containing the variables in the model.
#' @param x1_name the name of pre-outcome exposure
#' @param x3_name the name of post-outcome exposure
#' @param x2_name the name of exposure between pre-outcome exposure
#' and post-outcome exposure
#' @param boots_no the number of bootstrap
#'
#'
#' @return  \code{diff_test} the result of Diff and its 95%CI .

#' @examples
#' a <- rnorm(1000,0,1)
#' b <- a+rnorm(1000,0,1)
#' c <- b+rnorm(1000,0,1)
#' hhh <- data.frame(a,b,c)

#' hhh$x1 <- ifelse(a<(-0.5),0,ifelse(a>0.5,2,1))

#' hhh$x2 <- ifelse(b<(-0.5),0,ifelse(b>0.5,2,1))

#' hhh$x3 <- ifelse(c<(-0.5),0,ifelse(c>0.5,2,1))

#' con_result <- diff_test (data = hhh, expo_name1 =  'x1', expo_name2 =  'x2',
#' expo_name3 =  'x3', continuous = F, boot_no =1000)

#'
#'

diff_test <- function(data = data, expo_name1 = "x1",
    expo_name2 = "x2", expo_name3 = "x3", continuous = T,
    boot_no = 1000) {


    if (continuous == T) {

        con_result <- diff_test_pre(data = data, expo_name1 = expo_name1,
            expo_name2 = expo_name2, expo_name3 = expo_name3,
            continuous = T)

        cor_boot <- rep(NA, boot_no)
        for (i in 1:boot_no) {
            sample <- data[sample(1:dim(data)[1], dim(data)[1],
                replace = T), ]

            cor_boot[i] <- diff_test_pre(data = sample,
                expo_name1 = expo_name1, expo_name2 = expo_name2,
                expo_name3 = expo_name3, continuous = T)

        }
        qu_con <- quantile(cor_boot, probs = c(0.025,
            0.975), na.rm = T)
        con_result_all <- c(con_result, qu_con)
        con_result_all <- t(data.frame(con_result_all))
        colnames(con_result_all) <- c("Diff", "2.5%",
            "97.5")
        rownames(con_result_all) <- NULL

        result_end <- list()
        result_end$diff_test <- con_result_all


        return(result_end)
    } else {
        cat_result <- diff_test_pre(data = data, expo_name1 = expo_name1,
            expo_name2 = expo_name2, expo_name3 = expo_name3,
            continuous = F)
        cat_result <- data.frame(cat_result)
        cha_boot <- NULL

        for (i in 1:boot_no) {
            sample <- data[sample(1:dim(data)[1], dim(data)[1],
                replace = T), ]

            cha_boot_fen <- diff_test_pre(data = sample,
                expo_name1 = expo_name1, expo_name2 = expo_name2,
                expo_name3 = expo_name3, continuous = F)

            cha_boot <- rbind(cha_boot, cha_boot_fen)


        }
        cha_boot <- data.frame(cha_boot)

        test_result <- apply(cha_boot, 2, qu_fun)
        test_result <- t(data.frame(test_result))
        cha_end_hh <- cbind(cat_result, test_result)
        cha_end_hh <- data.frame(cha_end_hh)
        colnames(cha_end_hh) <- c("Diff", "2.5%", "97.5%")
        result_end <- list()
        result_end$diff_test <- cha_end_hh
        return(result_end)
    }
}


diff_test_pre <- function(data = data, expo_name1 = "x1",
    expo_name2 = "x2", expo_name3 = "x3", continuous = T) {
    data <- na.omit(data)

    if (continuous == T) {
        cor12 <- cor(data[, expo_name1], data[, expo_name2])
        cor23 <- cor(data[, expo_name2], data[, expo_name3])
        z_cor <- cor12 - cor23

        return(z_cor)

    } else {
        new_data <- data[, c(expo_name1, expo_name2,
            expo_name3)]
        colnames(new_data) <- c("x1", "x2", "x3")

        cat <- unique(new_data$x1)
        cat <- cat[order(cat)]
        n_length <- length(cat) - 1
        comdata <- min(cat)

        z_cha <- rep(NA, n_length)
        for (i in 1:n_length) {
            p11_qian <- sum(new_data$x1 == cat[i + 1] &
                new_data$x2 == cat[i + 1]) / length(new_data$x1)
            p10_qian <- sum(new_data$x1 == comdata &
                new_data$x2 == cat[i + 1]) / length(new_data$x1)
            p01_qian <- sum(new_data$x1 == cat[i + 1] &
                new_data$x2 == comdata) / length(new_data$x1)
            p00_qian <- sum(new_data$x1 == comdata &
                new_data$x2 == comdata) / length(new_data$x1)

            p11_hou <- sum(new_data$x2 == cat[i + 1] &
                new_data$x3 == cat[i + 1]) / length(new_data$x3)
            p10_hou <- sum(new_data$x2 == comdata & new_data$x3 ==
                cat[i + 1]) / length(new_data$x3)
            p01_hou <- sum(new_data$x2 == cat[i + 1] &
                new_data$x3 == comdata) / length(new_data$x3)
            p00_hou <- sum(new_data$x2 == comdata & new_data$x3 ==
                comdata) / length(new_data$x3)


            cha_qian_11 <- p11_qian / (p11_qian + p10_qian)
            cha_qian_10 <- p01_qian / (p01_qian + p00_qian)

            cha1 <- cha_qian_11 - cha_qian_10

            cha_hou_11 <- p11_hou / (p11_hou + p10_hou)
            cha_hou_10 <- p01_hou / (p01_hou + p00_hou)

            cha2 <- cha_hou_11 - cha_hou_10

            z_cha[i] <- cha1 - cha2
            names(z_cha)[i] <- paste0(cat[i + 1], " v.s.",
                comdata)
        }

        return(z_cha)
    }
}


qu_fun <- function(data) {
    qua <- quantile(data, probs = c(0.025, 0.975), na.rm = T)
    return(qua)
}
