RefTM <- function(sc_data, ref_data, k1 = 5, k2 = NULL, workflow = "LDA", covariate = NULL){
##RefTM-LDA workflow
  if(workflow == "LDA"){
    model_ref <- RefTM.LDA(ref_data, k = k1, k0 = 0, bulk_beta = NULL)
    print("Finished modeling the reference data.")
    print(k2)
    if(!is.null(k2)){
      print(k2)
      sprintf("Number of topics is fixed to be %s.", k1+5)
      model_sc <- RefTM.LDA(t(sc_data), k = k1+k2, method = "VEM", k0 = k1, bulk_beta = model_ref@beta, control = list(estimate.alpha = F))
    }else{
      print("Select the best model.")
      k <- k1 + c(0,2,5,10,15,20)
      models <- list()
      perplex <- c()
      for (i in 1:6) {
        models[[i]] <- RefTM.LDA(t(sc_data), k = k[i], method = "VEM", k0 = k1, bulk_beta = model_ref@beta, control = list(estimate.alpha = F))
        perplex[i] <- perplexity(models[[i]])
      }
      log.lik <- data.frame(topics = k, LL = perplex)
      log.lik$first_derivative <- c(-Inf, (diff(log.lik$LL) / diff(log.lik$topics)))
      log.lik$second_derivative <- c(-Inf, -Inf, diff(log.lik$first_derivative)[-1]/diff(log.lik$topics[-1]))
      par(bty = 'n')
      plot(log.lik$topics, log.lik$LL, xlab="Number of topics", ylab="perplexity", type='o', pch=16, col='black', main='Model selection')
      points(log.lik$topics[which(log.lik$second_derivative == max(log.lik$second_derivative))], log.lik$LL[which(log.lik$second_derivative == max(log.lik$second_derivative))], pch=4, col='red', lwd = 7)
      model_sc <- models[[which(log.lik$second_derivative == max(log.lik$second_derivative))]]
    }
  }else

##RefTM-STM workflow
  if(workflow == "STM"){
    if(is.null(covariate)){
      stop("Please use LDA workflow.")
    }
    ref_data  = ref_data[,which(colSums(ref_data) != 0)]
    sc_data = sc_data[which(colSums(ref_data) != 0),]
    getCorpus = function(mat){
      corpus = list()
      for (i in 1:dim(mat)[1]) {
        term = which(mat[i,] != 0)
        corpus[[i]] = rbind(term,mat[i,term])
      }
      return(corpus)
    }
    model_ref <- RefTM.STM(getCorpus(ref_data), as.character(1:dim(ref_data)[2]), K = k1)
    print("Finished modeling the reference data.")
    if(!is.null(k2)){
      sprintf("Number of topics is fixed to be %s.", k1+5)
      model_sc <- RefTM.STM(getCorpus(t(sc_data)), as.character(1:dim(sc_data)[1]), K = k1+5, prevalence = ~covariate, refBeta = exp(model_ref$beta$logbeta[[1]]))
    }else{
      print("Select the best model.")
      perp = function(model){
        exp(-model$convergence$bound[model$convergence$its]/sum(model$settings$dim$wcounts$x))
      }
      k <- k1 + c(0,2,5,10,15,20)
      models <- list()
      perplex <- c()
      for (i in 1:6) {
        models[[i]] <- RefTM.STM(getCorpus(t(sc_data)), as.character(1:dim(sc_data)[1]), K = k[i], prevalence = ~covariate, refBeta = exp(model_ref$beta$logbeta[[1]]))
        perplex[i] <- perp(models[[i]])
      }
      log.lik <- data.frame(topics = k, LL = perplex)
      log.lik$first_derivative <- c(-Inf, (diff(log.lik$LL) / diff(log.lik$topics)))
      log.lik$second_derivative <- c(-Inf, -Inf, diff(log.lik$first_derivative)[-1]/diff(log.lik$topics[-1]))
      par(bty = 'n')
      plot(log.lik$topics, log.lik$LL, xlab="Number of topics", ylab="perplexity", type='o', pch=16, col='black', main='Model selection')
      points(log.lik$topics[which(log.lik$second_derivative == max(log.lik$second_derivative))], log.lik$LL[which(log.lik$second_derivative == max(log.lik$second_derivative))], pch=4, col='red', lwd = 7)
      model_sc <- models[[which(log.lik$second_derivative == max(log.lik$second_derivative))]]
    }
  }
  return(model_sc = model_sc)
}



