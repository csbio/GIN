# author: Arshia Z. Hassan. hassa418@umn.edu

options(java.parameters = "- Xmx1024m")
packages <- c("stringi","stringr","rstudioapi","RColorBrewer"
              )
for (p in packages) {
  library(p, character.only = TRUE)
}

# function to generate data for precision-recall curve -- ( modified fFLEX package function GenerateDataForPerfCurve)
GenerateDataForPerfCurve_mod <- function(value.predicted, value.true, neg.to.pos = FALSE, x.axis = 'sensitivity', y.axis = 'precision'){
  
  # Sorting direction
  if (neg.to.pos == FALSE){
    indices <- order(value.predicted, decreasing = TRUE)
  } else{
    indices <- order(value.predicted)
  }
  value.true <- value.true[indices]
  value.predicted <- value.predicted[indices]
  
  # Calculate basic elements
  num.real.true <- sum(value.true)
  num.predicted.true <- 1 : length(value.true) # Predicted Positive <= Threshold
  
  TP <- cumsum(value.true)
  FP <- num.predicted.true - TP
  FN <- num.real.true - TP
  TN <- length(value.true) - (TP + FP + FN)
  
  # *** Get the indices of unique predicted values (last occurrence)
  # For cases when we have the same prediction for many elements
  # unique.index <- which(!duplicated(value.predicted)) # Gives index of first occurrence
  # unique.index <- c(unique.index - 1, length(TP)) # Now, it's the last occurrence.
  # unique.index <- unique.index[-1] # The first element is 0 anyway!
  # 
  # TP <- TP[unique.index]
  # FP <- FP[unique.index]
  # FN <- FN[unique.index]
  # TN <- TN[unique.index]
  
  #value.true <- value.true[unique.index]
  #value.predicted <- value.predicted[unique.index]
  
  precision <- TP / (TP + FP) # PPV
  sensitivity <- TP / (TP + FN)  # Recall / TPR
  FPR <- FP / (FP + TN) # FPR 
  
  ## *** TODO: What about the cases when we have very few unique precision and/or sensitivity values?
  
  # Find which to return for x and y
  switch(x.axis, TP = {x = TP}, FP = {x = FP}, TN = {x = TN}, FN = {x = FN},
         precision = {x = precision}, 
         FPR = {x = FPR}, 
         sensitivity = {x = sensitivity}, recall = {x = sensitivity}, TPR = {x = sensitivity})
  switch(y.axis, TP = {y = TP}, FP = {y = FP}, TN = {y = TN}, FN = {y = FN},
         precision = {y = precision}, 
         FPR = {y = FPR}, 
         sensitivity = {y = sensitivity}, recall = {y = sensitivity}, TPR = {y = sensitivity})
  
  # area under curve: according to trapizoidal approximation (make sense for roc or pr curve only: unit area)
  # https://www.r-bloggers.com/calculating-auc-the-area-under-a-roc-curve/
  # area of trapezoid = 0.5 * h * (b1 + b2) # diff in x = h (height) # b1, b2 are values on the y axis surrounding a trapezoid
  auc = abs(0.5 * sum((x[-1] - x[-length(x)]) * (y[-1] + y[-length(y)]))) # same formula as used in matlab (and gives the same value)
  
  # Or this works too
  # dx <- c(diff(x), 0)
  # dy <- c(diff(y), 0)
  # sum(y * dx) + sum(dx * dy)/2
  
  return(list(x = x, y = y, true = value.true, predicted = value.predicted, auc = auc ))
}

# wrapper function to generate and store precision-recall data and curve plot
generate_corum_similarity_plot<-function(corum_sim_list, output_folder, output_file_name_noext, curve_names, cols, title)
{
  
  profile_similarity_plot <- paste(output_file_name_noext,  sep = "")
  labs <- c('TP', 'Precision')
  plot_data = PlotPRSimilarity_mod(corum_sim_list, outfile.name = profile_similarity_plot, 
                   outfile.type = 'pdf', fig.labs = labs, fig.title = title,
                   subsample = F,
                   legend.names = curve_names, is.bgdline = T,
                   legend.color = cols, save.figure = TRUE)
  write.csv(plot_data,file.path(output_folder,paste0(output_file_name_noext,'.csv')),row.names=F)

  files <- c( 
    paste(profile_similarity_plot, '.pdf', sep = "")#, paste(profile_similarity_plot, '.png', sep = "")
  )
  for (f in files)
  {
    file.copy(f, file.path(output_folder, f), overwrite = TRUE)
    #file.remove(f)
  }
}

# function to generate precision-recall curve -- ( modified fFLEX package function PlotPRSimilarity)
PlotPRSimilarity_mod <- function(pred.ca, subsample = FALSE,
                             neg.to.pos = FALSE, 
                             type.plot = 'log', is.bgdline = FALSE,
                             provided.xlim = NULL, provided.ylim = NULL,
                             legend.names = NULL, legend.color = c('blue'),
                             legend.ltype = NULL, box.type = 'L',
                             fig.title = NULL, fig.labs = c('TP', 'Precision'),
                             save.figure = FALSE, 
                             outfile.name = 'test_PR_sim', outfile.type = 'pdf') {
  
  ## *** Check input data format
  if(class(pred.ca) != 'list'){
    stop ('A list of list is expected as input ... ')
  } else{ # check if we have the correct data inside the list
    
    corr_columns <- sum(grepl('true', names(pred.ca[[1]]))) + sum(grepl('predicted', names(pred.ca[[1]])))
    if(corr_columns != 2){
      stop (" All list element should contain a 'true' and a 'predicted' list ... ")
    }
  }
  
  if (length(legend.color) < length(pred.ca)){
    warning('Color not provided for each curve !! Making our own ')
    if (length(pred.ca) == 1){
      legend.color = c('blue')
    }
    else{
      legend.color = palette(rainbow(n = length(pred.ca)))
    }
  }
  
  
  ## *** Make the data to plot on the x (TP) and y (Precision) axis
  for (i in 1 : length(pred.ca)){
    
    test <- GenerateDataForPerfCurve_mod(value.predicted = pred.ca[[i]]$predicted, value.true = pred.ca[[i]]$true, neg.to.pos = neg.to.pos, x.axis = 'TP', y.axis = 'precision')
    
    ## ** If we want a faster plot (less points)
    if (subsample){
      
      # Take the last precision value for the same TP
      unique.index <- which(!duplicated(test$x))
      unique.index <- c(unique.index - 1, length(test$x))
      unique.index <- unique.index[-1]
      
      small_x <- test$x[unique.index]
      small_y <- test$y[unique.index]
      small_true <- test$true[unique.index]
      small_predicted <- test$predicted[unique.index]
      
      # Now subsample from the whole space (top 100, and then 100 from each 10^ range)
      highest.power <- ceiling(log10(length(small_x)))
      
      if (highest.power > 3){
        # All from 1 to 100
        keep.index <- c(1:100)
        
        # 100 from each interval (i.e 1001 to 10000)
        for (pow in 3 : (highest.power - 1)){
          keep.index <- c(keep.index, c(round(seq(10^(pow-1)+1, 10^pow, length = 100))))
        }
        max.available <- min(100, length(small_x) - (10^pow))
        keep.index <- c(keep.index, 
                        c(round(seq(10^pow+1, length(small_x), 
                                    length = max.available))))
        
        test$x <- small_x[keep.index]
        test$y <- small_y[keep.index]
        test$true <- small_true[keep.index]
        test$predicted <- small_predicted[keep.index]
        
      } else{ # If the number of points are < 1000, no need for subsampling
        test$x <- small_x
        test$y <- small_y
        test$true <- small_true
        test$predicted <- small_predicted
      }
    }
    
    if (i == 1){
      plot.data <- list(list(x = test$x, y = test$y, true= test$true, predicted = test$predicted))
    } else{
      plot.data <- append(plot.data, list(list(x = test$x, y = test$y, true= test$true, predicted = test$predicted)))
    }
  } # end for
  
  
  ## *** Calculate the max y-lim and x-lim for the plots (if not provided)
  #  In R, we have to this beforehand as we are going to plot multiple lines
  plot.xlim = -1
  plot.ylim = -1
  
  for (i in 1 : length(plot.data)) {
    tmp <- plot.data[[i]]
    
    ind.10 <- which(tmp$x == 9) # Remove data for TP < 10
    
    if (length(ind.10) > 0){ 
      # The first 10 points are excluded from plotting
      x <- max(tmp$x[-(1:ind.10[length(ind.10)])])
      y <- max(tmp$y[-(1:ind.10[length(ind.10)])])
    } else{ # What if it's not there?
      x <- max(tmp$x)
      y <- max(tmp$y)
    }
    
    if (x > plot.xlim){
      plot.xlim <- x
    }
    if (y > plot.ylim){
      plot.ylim <- y
    }
  }
  
  if(!is.null(provided.xlim)) {plot.xlim <- provided.xlim}
  if(!is.null(provided.ylim)) {plot.ylim <- provided.ylim}
  
  if (type.plot == 'log'){
    plot.xlim <- 10 ^ ceiling(log10(plot.xlim))  
  }
  
  plot.ylim <- round(plot.ylim + 0.01, 2)
  
  ## *** Save to an output file
  if (save.figure == TRUE){
    if(outfile.type == 'png'){
      png(paste(outfile.name, ".png", sep = ''), width = 7, height = 5, units="in", res = 300)
    }else{
      pdf(paste(outfile.name, ".pdf", sep = '') )   
    }
  }
  
  for (i in 1 : length(plot.data)) {
    # print(i)
    
    tmp <- plot.data[[i]]
    y.8 = tmp$y[max(which(tmp$predicted > .8))]
    y.5 = tmp$y[max(which(tmp$predicted > .5))]
    y.2 = tmp$y[max(which(tmp$predicted > .2))]
    y.0 = tmp$y[max(which(tmp$predicted > -1))]
    #auprc_baseline = sum(tmp$true)/length(tmp$true)
    print(y.8)
    print(y.5)
    print(y.2)
    print(y.0)
    #print(auprc_baseline)
    
    ind.10 <- which(tmp$x == 9) # Remove data for TP < 10
    
    if (length(ind.10) > 0){ # What if it's not there?
      # The first 10 points are excluded from plotting
      data.x.axis <- tmp$x[-(1:ind.10[length(ind.10)])]
      data.y.axis <- tmp$y[-(1:ind.10[length(ind.10)])]
    } else{
      data.x.axis <- tmp$x
      data.y.axis <- tmp$y
    }
    
    if (i == 1){ # For the first time
      if (type.plot == 'log'){
        
        plot(data.x.axis, data.y.axis, 
             xlim = c(10, plot.xlim), ylim = c(0, plot.ylim),
             log = 'x', type = 'l', main = fig.title,
             bty = box.type, # 'n' -> no box, nothing - all boxes
             xlab = fig.labs[1], ylab = fig.labs[2], 
             lwd = 2, col = legend.color[i], lty = legend.ltype[i],
             cex.lab = 1.4, cex.main = 1.4, cex.axis = 1.4,
             font.main = 1, xaxt="n") # Not plotting x axis tick labels
        
        # Adding customized xtick labels
        # x.axis.ticks <- seq(1, (floor(log10(plot.xlim))), 1)
        x.axis.ticks <- seq(1, log10(plot.xlim), 1)
        x.tick.labels <- as.expression(lapply(x.axis.ticks, function(E) bquote(10^.(E))))
        axis(1, at=10^(x.axis.ticks), 
             labels=as.expression(lapply(x.axis.ticks, function(E) bquote(10^.(E)))), 
             cex.axis = 1.4)
        abline(h=y.8, col=legend.color[i], lwd = 1)
        abline(h=y.5, col=legend.color[i], lwd = 1)
        abline(h=y.2, col=legend.color[i], lwd = 1)
        abline(h=y.0, col=legend.color[i], lwd = 1)
        #abline(h=auprc_baseline, col="black", lwd = 2, lty=2)
        
        #mtext('S>0.8',side = 4,at = y.8,adj = 0)
        text(plot.xlim, y.8+.005, 'S>0.8', xpd=NA)
        text(plot.xlim, y.5+.005, 'S>0.5', xpd=NA)
        text(plot.xlim, y.2-.005, 'S>0.2', xpd=NA)
        text(plot.xlim, y.0+.005, 'pos_GI', xpd=NA)
      }
      else{
        # Normal plot
        plot(data.x.axis, data.y.axis, xlim = c(10, plot.xlim), ylim = c(0, plot.ylim), 
             type = 'l', main = fig.title,
             bty = box.type, # 'n' -> no box, nothing - all boxes
             xlab = fig.labs[1], ylab = fig.labs[2], font.main = 1,
             lwd = 2, col = legend.color[i], lty = legend.ltype[i],
             cex.lab = 1.4, cex.main = 1.2, cex.axis = 1.4)
        abline(h=y.8, col=legend.color[i], lwd = 1)
        abline(h=y.5, col=legend.color[i], lwd = 1)
        abline(h=y.2, col=legend.color[i], lwd = 1)
        abline(h=y.0, col=legend.color[i], lwd = 1)
        #abline(h=auprc_baseline, col="black", lwd = 2, lty=2)
        
        #mtext('S>0.8',side = 4,at = y.8,adj = 0)
        text(plot.xlim, y.8+.005, 'S>0.8', xpd=NA)
        text(plot.xlim, y.5+.005, 'S>0.5', xpd=NA)
        text(plot.xlim, y.2-.005, 'S>0.2', xpd=NA)
        text(plot.xlim, y.0+.005, 'pos_GI', xpd=NA)
      }
    }
    else{ # For every other time
      lines(data.x.axis, data.y.axis, col=legend.color[i], lwd = 2, lty = legend.ltype[i]) 
      abline(h=y.8, col=legend.color[i], lwd = 1)
      abline(h=y.5, col=legend.color[i], lwd = 1)
      abline(h=y.2, col=legend.color[i], lwd = 1)
      abline(h=y.0, col=legend.color[i], lwd = 1)
      #abline(h=auprc_baseline, col="black", lwd = 2, lty=2)
      
      #mtext('S>0.8',side = 4,at = y.8,adj = 0)
      text(plot.xlim, y.8+.005, 'S>0.8', xpd=NA)
      text(plot.xlim, y.5+.005, 'S>0.5', xpd=NA)
      text(plot.xlim, y.2-.005, 'S>0.2', xpd=NA)
      text(plot.xlim, y.0+.005, 'pos_GI', xpd=NA)
    }
  }
  
  if (!is.null(legend.names)){
    legend("topright", legend = legend.names, 
           fill = legend.color, lwd = 2, lty = legend.ltype[i],
           bty = box.type, # to remove the box
           cex = 1.2, text.col = "black", horiz = F)
    
    # Printing some warning, just in case!
    if (length(pred.ca) > length(legend.names)){
      warning('Legend not provided for all curves !!')
    }
    if (length(legend.color) != length(legend.names)){
      warning('Number of legends and colors are not equal !!!')
    }
  }
  
  if (is.bgdline){
    abline(h = min(data.y.axis), col = 'black', lty=2)
  }
  
  if (save.figure == TRUE){
    dev.off()
  }
  return(test)
}


# load precision-recall data generated by 
# generate_PR_GO_data_flex_suppression_score.r (wrt GO-BP standard) 
# and generate_PR_PW_data_flex_suppression_score.r (wrt pathway standard) scripts.
# Generate PR-curve plots and data

input_folder<-'./'
output_folder<-'./'

curve_labels <- c('')
col_set <- c()
cols0<-brewer.pal(4, "PiYG")
col_set <- c(col_set,cols0[1])

load(file.path(input_folder,'suppressors_all_pairs_GO_metric_2_.9.Rdata'))
corum_sim <- list(list(true = go$true, predicted = go$predicted))
file_ext <- "all_pairs_GO_PR"
generate_corum_similarity_plot(corum_sim, output_folder, file_ext, curve_labels, col_set, "All gene pairs" )

load(file.path(input_folder,'suppressors_no_mito_GO_metric_2_.9.Rdata'))
corum_sim <- list(list(true = go_no_mito$true, predicted = go_no_mito$predicted))
file_ext <- "no_mito_GO_PR"
generate_corum_similarity_plot(corum_sim, output_folder, file_ext, curve_labels, col_set, "No mito. gene pairs" )

load(file.path(input_folder,'suppressors_all_pairs_PW_metric_2_.9.Rdata'))
corum_sim <- list(list(true = pathway$true, predicted = pathway$predicted))
file_ext <- "all_pairs_PW_PR"
generate_corum_similarity_plot(corum_sim, output_folder, file_ext, curve_labels, col_set, "All gene pairs" )

load(file.path(input_folder,'suppressors_no_mito_PW_metric_2_.9.Rdata'))
corum_sim <- list(list(true = pathway_no_mito$true, predicted = pathway_no_mito$predicted))
file_ext <- "no_mito_PW_PR"
generate_corum_similarity_plot(corum_sim, output_folder, file_ext, curve_labels, col_set, "No mito. gene pairs" )
