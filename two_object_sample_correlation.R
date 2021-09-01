
#' two_object_sample_correlation
#'
#' @description Prints a correlation heatmap with Pearson R values \
#' for two objects' samples split by each given COI
#' @note Using spearman correlation uses raw cell numbers instead of percentages of sample
#' @note Version 3: returns permutation test results included in pheatmap
#' @importMethodsFrom jmuOutlier
#' @param SOs A list of seurat objects to sample so1 and so2 from
#' @param COIs A condition of interest for each object to split heatmap by EX: COIs <- list(so1_coi, so2_coi)
#' @param statmethod The correlation statistic, either pearson or spearman
#' @param cellpcts Boolean, whether we want percentage of orig.ident or raw cell numbers
#' @param num_sim The number of simulations to run for the correlation permutation test for p-values
#' @param noplot Whether to output a figure or not
#' @param with_permtest Whether to add the p-values to the figure if noplot = F, or to the returned lists if noplot=T
#' @param norm_vector An external vector to normalize sample counts by, for instance, num CD45+ cells from larger object \
#' (only works with cellpcts=F)
#' @return NULL (prints the figure instead)
two_object_sample_correlation <- function(SOs, COIs, statmethod="pearson",
                                          cellpcts=T, num_sim = 10000,
                                          noplot = F, with_permtest = T,
                                          norm_vector = NULL) {
  require(jmuOutlier)

  so1 <- SOs[[1]]
  so2 <- SOs[[2]]
  coi1 <- COIs[[1]]
  coi2 <- COIs[[2]]

  if (statmethod != "pearson" && statmethod != "spearman") {
    print("statmethod must either be pearson or spearman, breaking...")
    return(NULL)
  }
  if (!(all(sort(unique(so1@meta.data$orig.ident)) == sort(unique(so2@meta.data$orig.ident))))) {
    print("object orig.ident's do not match...")
    return(NULL)
  } else {
    name.ordering = sort(unique(so1@meta.data$orig.ident))
  }
  if (so1@project.name == "SeuratProject") {
    name1 = "object1"
  } else {
    name1 = so1@project.name
  }
  if (so2@project.name == "SeuratProject") {
    name2 = "object2"
  } else {
    name2 = so2@project.name
  }
  cor.df <- data.frame(matrix(nrow = length(unique(so1@meta.data[,coi1])),
                              ncol = length(unique(so2@meta.data[,coi2]))))
  clusters1 <- sort(unique(so1@meta.data[,coi1]))
  clusters2 <- sort(unique(so2@meta.data[,coi2]))
  rownames(cor.df) <- paste0(name1, "_", clusters1)
  colnames(cor.df) <- paste0(name2, "_", clusters2)
  if (with_permtest) {
    cor.df.p <- cor.df #make a copy to use for p-values from permutation test
  }
  for (i in clusters1) {
    for (j in clusters2) {
      if (cellpcts) {
        cur.so1 <- table(factor(so1@meta.data$orig.ident[which(so1@meta.data[,coi1]==i)], levels=name.ordering)) /
          table(so1@meta.data$orig.ident)
        cur.so2 <- table(factor(so2@meta.data$orig.ident[which(so2@meta.data[,coi2]==j)], levels=name.ordering)) /
          table(so2@meta.data$orig.ident)
      } else {
        if (is.null(norm_vector)) {
          cur.so1 <- table(factor(so1@meta.data$orig.ident[which(so1@meta.data[,coi1]==i)], levels=name.ordering))
          cur.so2 <- table(factor(so2@meta.data$orig.ident[which(so2@meta.data[,coi2]==j)], levels=name.ordering))
        } else {
          if (length(intersect(names(norm_vector), name.ordering)) == length(name.ordering)) { #check if the same elements are in both
            norm_vector <- norm_vector[match(name.ordering, names(norm_vector))]
          } else {
            stop("norm_vector names don't match names provided by objects")
          }
          cur.so1 <- table(factor(so1@meta.data$orig.ident[which(so1@meta.data[,coi1]==i)], levels=name.ordering)) /
            norm_vector
          cur.so2 <- table(factor(so2@meta.data$orig.ident[which(so2@meta.data[,coi2]==j)], levels=name.ordering)) /
              norm_vector
        }
      }
      cur.so1 <- cur.so1[match(name.ordering, names(cur.so1))]
      cur.so2 <- cur.so2[match(name.ordering, names(cur.so2))]
      if (with_permtest) {
        set.seed(415)
        cur.permtest <- perm.cor.test(as.numeric(cur.so1), as.numeric(cur.so2), alternative = "two.sided", method = statmethod, num.sim = num_sim)
        cor.df.p[paste0(name1, "_", i),paste0(name2, "_", j)] <- cur.permtest$p
      }
      cor.df[paste0(name1, "_", i),paste0(name2, "_", j)] <- cor(cur.so1, cur.so2, method = statmethod)
    }
  }
  if (noplot) {
    if (with_permtest) {
      return(list("correlations"=cor.df, "permutation_test_pvalues"=cor.df.p))
    } else {
      return(list("correlations"=cor.df, "permutation_test_pvalues"=c()))
    }
  } else {
    if (with_permtest) {
      return(pheatmap(cor.df, display_numbers=cor.df.p, cluster_cols = F, cluster_rows = F, angle_col = 315))
    } else {
      return(pheatmap(cor.df, cluster_cols = F, cluster_rows = F, angle_col = 315))
    }
  }
}

