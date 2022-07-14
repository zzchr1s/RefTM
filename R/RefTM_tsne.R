RefTM_tsne = function(mat,label,donor = NULL,title = NULL,seed = 1,perplexity = 30,size = 2,legend = FALSE,palette = NULL){
  set.seed(seed)
  tsne.info = Rtsne::Rtsne(mat,perplexity = perplexity)
  colnames(tsne.info$Y) = c("tSNE_1","tSNE_2")
  if(is.null(donor) == TRUE){
    tsne.data = data.frame(label = label,tsne.info$Y)
    p = ggpubr::ggscatter(tsne.data,title = title,x = "tSNE_1", y = "tSNE_2",color = "label",size = size,xlab = NULL, ylab = NULL,palette = palette)+
      ggthemes::theme_base() + ggplot2::theme(plot.title = element_text(hjust = 0.5))+ ggplot2::theme(legend.position="left")

    if(legend == FALSE) p = p + ggplot2::theme(legend.position="none")
  }
  if(is.null(donor) == FALSE){
    tsne.data = data.frame(label = label,donor = donor, tsne.info$Y)
    p = ggpubr::ggscatter(tsne.data,title = title,x = "tSNE_1", y = "tSNE_2",color = "label",shape = "donor",size = size,xlab = NULL
                  , ylab = NULL,palette = palette)+ ggplot2::scale_shape_manual(values=c(16,8)) + ggthemes::theme_base() + ggplot2::theme(plot.title = element_text(hjust = 0.5))+
      ggplot2::theme(legend.position="left")
    if(legend == FALSE) p = p + ggplot2::theme(legend.position="none")
  }
  p
}
