########################################################################
### MR script.                                                       ###
### This script performs four different types of MR analysis.        ###
### It will provide the Inverse variance weighted analysis, the one  ###
### most people do as default. It will also do the Eggers, median    ###
### and weighted median that Stephen Burgess has proposed.           ###
###                                                                  ###
### In addition this it will also draw funnel and forest plots for   ###
### the variants included in the model which are additional          ###
### sensitivity checks that are useful to perform on Mendelian       ###
### randomisation results.                                           ###
###                                                                  ###
### Finally, this script will also generate a dosage plot that       ###
### includes the lines of best fit based on each method to allow for ###
### a visual assessment of the coherence of the estimated effect.    ###
###                                                                  ###
### The script is set up to read a file which contains beta and SEs  ###
### for the associations with the SNP and the outcomes of interest   ###
### and the proposed middle factors. It is assumed that the input    ###
### file will have a first column called SNP, with subsequent cols   ###
### names Beta.XXX, SE.XXX where XXX indicates the code for the      ###
### variables in question. These Betas need to be aligned. Replace   ###
### the strings in the "inputs" and "outputs" vector with the codes. ###
###                                                                  ###
###     Written by Felix 16/09/2015                                  ###
########################################################################

MR <- function (test_data, inputs, outputs, MR_results=NULL, nb=1000)
{
## Set-up
    if (require("metafor")) {
        print("metafor is loaded correctly")
    }
    else {
        print("Trying to install metafor")
        install.packages("metafor")
        if (require(metafor)) {
            print("metafor installed and loaded")
        }
        else {
            stop("Could not install metafor")
        }
    }
    require(ggplot2)
    if (require("ggplot2")) {
        print("ggplot2 is loaded correctly")
    }
    else {
        print("Trying to install ggplot2")
        install.packages("ggplot2")
        if (require(ggplot2)) {
            print("ggplot2 installed and loaded")
        }
        else {
            stop("Could not install ggplot2")
        }
    }

## Function for weighted medians
    weighted.median <- function(betaIV.in, weights.in) {
        betaIV.order = betaIV.in[order(betaIV.in)]
        weights.order = weights.in[order(betaIV.in)]
        weights.sum = cumsum(weights.order) - 0.5 * weights.order
        weights.sum = weights.sum/sum(weights.order)
        below = max(which(weights.sum < 0.5))
        weighted.est = betaIV.order[below] + (betaIV.order[below + 
            1] - betaIV.order[below]) * (0.5 - weights.sum[below])/(weights.sum[below + 
            1] - weights.sum[below])
        return(weighted.est)
    }

## Bootstrapping the weighted median to get SEs
    weighted.median.boot = function(betaXG.in, betaYG.in, sebetaXG.in, 
        sebetaYG.in, weights.in) {
        med = NULL
        for (i in 1:nb) {
            betaXG.boot = rnorm(length(betaXG.in), mean = betaXG.in, 
                sd = sebetaXG.in)
            betaYG.boot = rnorm(length(betaYG.in), mean = betaYG.in, 
                sd = sebetaYG.in)
            betaIV.boot = betaYG.boot/betaXG.boot
            med[i] = weighted.median(betaIV.boot, weights.in)
        }
        return(sd(med))
    }

## Set up results list
    n1 <- length(inputs)
    n2 <- length(outputs)
    n <- n1 * n2
    l <- vector("list", n)
    k = 1

## Loops
    for (j in inputs) {

       ## Defining input variable names
        BetaIn <- eval(parse(text = paste("test_data$Beta.", 
            j, sep = "")))
        SEIn <- eval(parse(text = paste("test_data$SE.", j, sep = "")))
        for (i in outputs) {

            ## Defining outcome variable names
            BetaOut <- eval(parse(text = paste("test_data$Beta.", 
                i, sep = "")))
            SEOut <- eval(parse(text = paste("test_data$SE.", 
                i, sep = "")))
            ColOut <- (parse(text = paste(j, ".", i, sep = "")))
            l[k] <- paste(j, ".", i, sep = "")
            k = k + 1

            ## Defining graph output options
            graph1 <- paste(j, "_", i, "_funnel.png", sep = "")
            graph2 <- paste(j, "_", i, "_forest.png", sep = "")
            graph3 <- paste(j, "_", i, "_dosage.png", sep = "")
            graph_title <- paste("Effect of SNPs for", j, "on", 
                i, sep = " ")
            x_title <- paste("Effect on", j, sep = " ")
            y_title <- paste("Effect on", i, sep = " ")

            ## IVW analysis
            betaIVW = summary(lm(BetaOut ~ BetaIn - 1, weights = SEOut^-2))$coef[1, 
                1]
            sebetaIVW = summary(lm(BetaOut ~ BetaIn - 1, weights = SEOut^-2))$coef[1, 
                2]/min(summary(lm(BetaOut ~ BetaIn - 1, weights = SEOut^-2))$sigma, 
                1)

            ## Egger's analysis
            betaEGGER = summary(lm(BetaOut ~ BetaIn, weights = SEOut^-2))$coef[2, 
                1]
            sebetaEGGER = summary(lm(BetaOut ~ BetaIn, weights = SEOut^-2))$coef[2, 
                2]/min(summary(lm(BetaOut ~ BetaIn, weights = SEOut^-2))$sigma, 
                1)

            ## Recording the intercept - used as a test of directional plietropy.
            interEGGER = summary(lm(BetaOut ~ BetaIn, weights = SEOut^-2))$coef[1, 
                1]
            seinterEGGER = summary(lm(BetaOut ~ BetaIn, weights = SEOut^-2))$coef[1, 
                2]

            ## Weighted median and penalised weighted median analysis
            betaIV = BetaOut/BetaIn
            weights = (SEOut/BetaIn)^-2
            betaIVW = sum(BetaOut * BetaIn * SEOut^-2)/sum(BetaIn^2 * 
                SEOut^-2)
            penalty = pchisq(weights * (betaIV - betaIVW)^2, 
                df = 1, lower.tail = FALSE)
            pen.weights = weights * pmin(1, penalty * 20)
            betaWM = weighted.median(betaIV, weights)
            sebetaWM = weighted.median.boot(BetaIn, BetaOut, 
                SEIn, SEOut, weights)
            betaPWM = weighted.median(betaIV, pen.weights)
            sebetaPWM = weighted.median.boot(BetaIn, BetaOut, 
                SEIn, SEOut, pen.weights)

            ## Generating funnel and forest plot data
            res <- rma(betaIV, (1/sqrt(weights)), method = "FE")
            CochQ = res$QE
            CochQp = res$QEp

            ## Funnel plots        
#           png(file = graph1)
            funnel(res, main = paste(graph_title), xlab = y_title)
            abline(v = 0, lty = "dashed", col = "red")
#           dev.off()

            ## Forest plots
#           png(file = graph2)
            forest(res, main = paste(graph_title), addfit = FALSE, 
                slab = test_data$SNP, xlab = y_title)
#           dev.off()

            ## Dosage plots
            dosage_data <- data.frame(BetaIn, SEIn, BetaOut, 
                SEOut)
            dosage_data$inLCI <- dosage_data$BetaIn + (1.96 * 
                dosage_data$SEIn)
            dosage_data$inUCI <- dosage_data$BetaIn - (1.96 * 
                dosage_data$SEIn)
            dosage_data$outLCI <- BetaOut - (1.96 * SEOut)
            dosage_data$outUCI <- BetaOut + (1.96 * SEOut)
            dose_graph <- ggplot(data = dosage_data, aes(x = BetaIn, 
                y = BetaOut)) + geom_point() + theme_bw() + geom_errorbar(aes(ymin = outLCI, 
                ymax = outUCI)) + geom_errorbarh(aes(xmin = inLCI, 
                xmax = inUCI)) + geom_abline(intercept = 0, slope = 0, 
                size = 1) + geom_abline(intercept = 0, slope = betaIVW, 
                size = 1, colour = "red") + geom_abline(intercept = 0, 
                slope = betaEGGER, size = 1, colour = "blue") + 
                geom_abline(intercept = 0, slope = betaWM, size = 1, 
                  colour = "green") + geom_abline(intercept = 0, 
                slope = betaPWM, size = 1, colour = "orange") + 
                geom_vline(xintercept = 0, size = 1) + ggtitle(graph_title) + 
                xlab(x_title) + ylab(y_title)
#           ggsave(filename = graph3, plot = dose_graph)
            plot(dose_graph)
            rm(dosage_data)

            ##Setting up vector of with-in loop results
            assign(paste(ColOut), c(paste(ColOut), betaIVW, sebetaIVW, 
                CochQ, CochQp, betaEGGER, sebetaEGGER, interEGGER, 
                seinterEGGER, betaWM, sebetaWM, betaPWM, sebetaPWM))
        }
    }

#Formatting results table
    results <- c("IV", "betaIVW", "sebetaIVW", "CochQ", "CochQp", 
        "betaEGGER", "sebetaEGGER", "interEGGER", "seinterEGGER", 
        "betaWM", "sebetaWM", "betaPWM", "sebetaPWM")
    for (p in l) {
        q <- eval(parse(text = paste(p)))
        results <- rbind(results, q)
    }
    if (!is.null(MR_results)) write.table(results, file = MR_results,
         sep = ",", row.names = FALSE, col.names = FALSE)
    invisible(results)
}

load("test_data.rda")
all_outputs = c("AFB", "ALS", "NoC", "Hap", "Int", "NuC", "Men", 
            "P16", "RsT", "Tow", "Irr", "Alc", "Smo", "SFN", "BMI", 
            "Une", "Hgt", "NSP")
m <- MR (test_data, "afs", all_outputs, "MR_afs.csv")
