#' @title \code{getPriors}: Function to obtain priors from a fitted \code{\link[INLA]{inla}} object.
#' @description This function is used to obtain the prior distributions relevant to models fit with intSDM, and return them in a readable format.
#' @param model A \code{\link[INLA]{inla}} model object.
#'
#' @return The function returns a list with two items: fixedEffects describing the priors used for the fixed effects, and randomEffects describing the priors used for the random effects.
#' @export

getPriors <- function(model) {

  if (!inherits(model, 'inla')) stop ('model must be an inla object')

   #If not empty
  fixedEffects <- t(setNames(data.frame(sapply(model$all.hyper$linear,
                                               function(x) c(mean = x$prior.mean, prec = x$prior.prec))),
                             sapply(model$all.hyper$linear,
                                    function(x) x$label)))

  if (!is.null(model$summary.hyperpar)) {

   randomIn <- unlist(lapply(sapply(model$all.hyper$random, function(x) x[['hyperid']]),
                              function(y) if (is.null(y)) NA else y))

   names(model$all.hyper$random) <- randomIn

   randomList <- list()

   for (eff in randomIn) {
      ##Need to do for group model too
     if (!is.na(eff)) {

       if (inherits(model$bru_info$model$effects[[eff]]$main$model, 'inla.spde2')) {

         if (identical(model$bru_info$model$effects[[eff]]$main$model$model, 'pcmatern')) {

           #not fixed
           if (any(grepl(paste('Range for', eff), row.names(model$summary.hyperpar)))) {

             Range1 <- exp(model$bru_info$model$effects[[eff]]$main$model$f$hyper.default$theta1$initial -1)
             Range2 <- exp(-1 * model$bru_info$model$effects[[eff]]$main$model$f$hyper.default$theta1$param[1]/Range1)


           } else {
             #Fixed
             Range1 <- exp(model$bru_info$model$effects[[eff]]$main$model$f$hyper.default$theta1$initial)
             Range1 <- NA
           }

           #not fixed
           if (any(grepl(paste('Stdev for', eff), row.names(model$summary.hyperpar)))) {

             Sigma1 <- exp(model$bru_info$model$effects[[eff]]$main$model$f$hyper.default$theta2$initial + 1)
             Sigma2 <- exp(-1 * model$bru_info$model$effects[[eff]]$main$model$f$hyper.default$theta1$param[2]*Sigma1)


           } else {
             #Fixed
             Sigma1 <- exp(model$bru_info$model$effects[[eff]]$main$model$f$hyper.default$theta2$initial)
             Sigma2 <- NA
           }

           randomList[[eff]] <- list(prior = 'pcmatern', values = list(Range = c(range = Range1, prob = Range2),
                                                                       StDev = c(sigma = Sigma1, prob = Sigma2)))


         } else {

           #MVNORM prior for the thetas

           nTheta <- sum(grepl('Theta', row.names(model$summary.hyperpar)) & grepl(eff, row.names(model$summary.hyperpar)))

           meanVec <- setNames(model$bru_info$model$effects[[eff]]$main$model$f$hyper.default$theta1$param[1:nTheta], paste0('Theta', 1:nTheta))
           precMat <- matrix(model$bru_info$model$effects[[eff]]$main$model$f$hyper.default$theta1$param[(nTheta+1):length(model$bru_info$model$effects[[eff]]$main$model$f$hyper.default$theta1$param)],
                             nrow = nTheta)

           row.names(precMat) <- colnames(precMat) <- paste0('Theta', 1:nTheta)


           randomList[[eff]] <- list(prior = 'mvnorm', values = list(meanVector = meanVec,
                                                                     precMatrix = precMat))

         }


       }
       else {

           lengthPar <- length(model$all.hyper$random[[eff]]$hyper)
           parList <- vector(mode = 'list', length = lengthPar)

           for (k in 1:lengthPar) {

             parList[[k]] <- list(name =  model$all.hyper$random[[eff]]$hyper[[k]]$name,
                                  prior =  model$all.hyper$random[[eff]]$hyper[[k]]$prior,
                                  param = model$all.hyper$random[[eff]]$hyper[[k]]$param,
                                  fixed = model$all.hyper$random[[eff]]$hyper[[k]]$fixed)

           }

           randomList[[eff]] <- parList

       }

       if (!is.null(model$bru_info$model$effects[[eff]]$group)) {

         lengthGroup <- length(model$all.hyper$random[[eff]]$group.hyper)
         groupList <- vector(mode = 'list', length = lengthGroup)

         for (j in 1:lengthGroup) {

           groupList[[j]] <- list(name = model$all.hyper$random[[eff]]$group.hyper$theta$name,
                                  prior = model$all.hyper$random[[eff]]$group.hyper$theta$prior,
                                  param = model$all.hyper$random[[eff]]$group.hyper$theta$param,
                                  fixed = model$all.hyper$random[[eff]]$group.hyper$theta$fixed)

         }

         randomList[[eff]][['group']] <- groupList

       }


     }


   }




  }

  return(list(fixedEffects = fixedEffects, randomEffects = randomList))

}
