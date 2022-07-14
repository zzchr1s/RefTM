RefTM_motif = function (peaks, label_true, label_est, sc_data, cluster_peak_num = 1000,
                       motif_num = 50)
{
  result <- list()
  metaData <- data.frame(true_label = (label_true), pred_label = label_est)
  metaDataInd = with(metaData, order(pred_label, true_label))
  resultPeakSelection <- RA3:::getClusterSpecificPvalue(data = sc_data,
                                                        cluster = label_est, offset = colSums(sc_data))
  pvaluematrix = resultPeakSelection$pvalue
  rownames(pvaluematrix) = peaks
  colData <- data.frame(true_label = (label_true))
  colData <- colData[metaDataInd, ]
  peakselected = c()
  for (i in 1:ncol(pvaluematrix)) {
    peakselected = c(peakselected, head(order(pvaluematrix[,
                                                           i], decreasing = FALSE), cluster_peak_num))
  }
  peakselected = unique(sort(peakselected))
  chr = c()
  p1 = c()
  p2 = c()
  for (i in 1:length(peaks[peakselected])) {
    strings = unlist(strsplit(peaks[peakselected][i], "_"))
    chr = c(chr, strings[1])
    p1 = c(p1, strings[2])
    p2 = c(p2, strings[3])
  }
  peakrange = data.frame(chr, p1, p2)
  peakrange$p1 = as.numeric(as.character(peakrange$p1))
  peakrange$p2 = as.numeric(as.character(peakrange$p2))
  peaks.gr = GenomicRanges::GRanges(peakrange[, 1], IRanges::IRanges(peakrange[,
                                                                               2], peakrange[, 3]))
  matrix_use = sc_data
  matrix_select = matrix_use[peakselected, ]
  matrix_select = matrix_select[, metaDataInd]
  colnames(matrix_select) = 1:ncol(matrix_select)
  frag_counts = SummarizedExperiment::SummarizedExperiment(assays = S4Vectors::SimpleList(counts = matrix_select),
                                                           rowRanges = peaks.gr, colData = colData)
  frag_counts = chromVAR::addGCBias(frag_counts, genome = BSgenome.Mmusculus.UCSC.mm9::BSgenome.Mmusculus.UCSC.mm9)
  motifs <- chromVAR::getJasparMotifs()
  motifs.matched = motifmatchr::matchMotifs(motifs, frag_counts,
                                            genome = BSgenome.Mmusculus.UCSC.mm9::BSgenome.Mmusculus.UCSC.mm9)
  dev = chromVAR::computeDeviations(object = frag_counts, annotations = motifs.matched)
  dev.scores = chromVAR::deviationScores(dev)
  variability = chromVAR::computeVariability(dev)
  top_motifs = variability$name[head(order(variability$variability,
                                           decreasing = TRUE), motif_num)]
  names(top_motifs) = rownames(variability[head(order(variability$variability,
                                                      decreasing = TRUE), motif_num), ])
  top_devs = dev.scores[which(rownames(dev.scores) %in% names(top_motifs)),
                        ]
  rownames(top_devs) = top_motifs[match(rownames(top_devs),
                                        names(top_motifs))]
  result$metaDataInd <- metaDataInd
  result$top_devs <- top_devs
  return(result)
}
