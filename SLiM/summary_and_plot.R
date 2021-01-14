get_gen_means <- function(res){
  # take mean at each generation over replicates
  mean_res = aggregate(res[, c('expected_Ne','realized_Ne', 'expected_F', 'realized_F' )], list(res[,'g']), mean)
  return(mean_res)
}

plot_tempF <- function(mean_res, title=NULL){
  # plot observed vs expected temporal F
  plot(mean_res$Group.1, mean_res$realized_F, type="p", lwd=3, col="red", xlab="Generation", ylab="Temporal F")
  lines(mean_res$Group.1, mean_res$expected_F,type="l",lty=3,lwd=3,col="blue")
  legend(0.15*max(mean_res$Group.1), 0.8*max(mean_res$realized_F),
         legend=c("Observed","Expected"),
         col=c("red", "blue"), lty=1:2, cex=1.5)
  title(title)
}

res_SLIM = data.matrix(read.table("./res.txt", sep = ",", header = FALSE))
colnames(res_SLIM) = c('rep', 'g', 'expected_Ne','realized_Ne', 'expected_F', 'realized_F')
plot_tempF(get_gen_means(res_SLIM), title = 'res_SLIM')
