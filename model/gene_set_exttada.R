gene_set_exttada <- function(geneset, cases, samplenumber, reference, rankPercentileData = NA, cal_pval = FALSE, outputpath=NA) {
  # form a proper input for TADA
  mytada.data1 = data.frame(gene.id=geneset,
                            mut.cls0=numeric(length(geneset)),
                            mut.cls1=numeric(length(geneset)),
                            mut.cls2=numeric(length(geneset)),
                            dn.cls0=numeric(length(geneset)),
                            dn.cls1=numeric(length(geneset)),
                            dn.cls2=numeric(length(geneset)))
  cases1 = cases
  samplenumber1 = samplenumber
  for (i in 1:length(geneset)) {
    index = match(as.character(geneset[i]), reference$GeneID)
    # check whether we found the gene
    if (!is.na(index)) {
      mytada.data1$mut.cls0[i] = reference$Mu_Silent[index]
      mytada.data1$mut.cls1[i] = reference$Mu_LoF[index]
      mytada.data1$mut.cls2[i] = reference$Mu_Dmis_REVEL0.5[index]
      
    } else {
      mytada.data1$mut.cls0[i] = .Machine$double.eps
      mytada.data1$mut.cls1[i] = .Machine$double.eps
      mytada.data1$mut.cls2[i] = .Machine$double.eps
      
    }
    if (mytada.data1$mut.cls0[i] <= 0) {
      mytada.data1$mut.cls0[i] = .Machine$double.eps
    }
    if (mytada.data1$mut.cls1[i] <= 0) {
      mytada.data1$mut.cls1[i] = .Machine$double.eps
    }
    if (mytada.data1$mut.cls2[i] <= 0) {
      mytada.data1$mut.cls2[i] = .Machine$double.eps
    }
    gene_cases1 = cases1[cases1$GeneID==as.character(geneset[i]),]
    
    mytada.data1$dn.cls0[i] = sum(gene_cases1$vclass=="syn")
    mytada.data1$dn.cls1[i] = sum(gene_cases1$vclass=="LGD")
    mytada.data1$dn.cls2[i] = sum(gene_cases1$vclass=="Dmis")
    
  }
  # first do sanity check using syn mutation
  message("sanity check before correction")
  syn_oe1 = sum(mytada.data1$dn.cls0) / (sum(mytada.data1$mut.cls0) * samplenumber1 * 2)
  message(paste0("syn_oe1 = ", syn_oe1))
  message("adjusting mutation rate according to syn_oe")
  constant1 = sum(mytada.data1$dn.cls0) / (sum(mytada.data1$mut.cls0) * samplenumber1 * 2)
  mytada.data1[,2:4]=mytada.data1[,2:4]*constant1
  
  syn_oe1 = sum(mytada.data1$dn.cls0) / (sum(mytada.data1$mut.cls0) * samplenumber1 * 2)
  LGD_oe1 = sum(mytada.data1$dn.cls1) / (sum(mytada.data1$mut.cls1) * samplenumber1 * 2)
  Dmis_oe1 = sum(mytada.data1$dn.cls2) / (sum(mytada.data1$mut.cls2) * samplenumber1 * 2)
  
  message("sanity check after correction")
  message(paste0("syn_oe1 = ", syn_oe1))
  message(paste0("LGD_oe1 = ", LGD_oe1))
  message(paste0("Dmis_oe1 = ", Dmis_oe1))
  message("output geneset datafile")
  write.table(mytada.data1, file = "mytada.data1.csv", quote = FALSE, sep = ",", row.names = FALSE)
  
  # use extTADA to sample values of parameters
  dataDN1 = mytada.data1[,c('dn.cls1','dn.cls2')]
  colnames(dataDN1) = c('dn_LGD','dn_Dmis')
  mutRate1 = mytada.data1[,c('mut.cls1','mut.cls2')]
  colnames(mutRate1) = c('mut_LGD','mut_Dmis')
  
  rankPercentile = rankPercentileData$value[match(geneset, rankPercentileData$gene)]
  rankPercentile[is.na(rankPercentile)] = 0.5
  rankPercentile = dplyr::percent_rank(rankPercentile)
  options(mc.cores = parallel::detectCores())
  mcmcDD <- extTADAmcmc(modelName = DNextTADA,
                        dataDN = dataDN1,
                        mutRate = mutRate1,
                        rankPercentile = rankPercentile,
                        Ndn = rep(samplenumber, 2),
                        nIteration = 2000,
                        nChain = 4,
                        nCore = 4)
  # saveRDS(mcmcDD, file = "extTADA_modified_7.mcmcDD.RDS")
  # mcmcDD <- readRDS("extTADA_modified_7.mcmcDD.RDS")
  # options(repr.plot.width = 4, repr.plot.height = 3)
  # par(mfrow = c(1,2))
  #plotParHeatmap(c("pi0[1]", "hyperGammaMeanDN[1]"), mcmcResult = mcmcDD)
  #plotParHeatmap(c("pi0[2]", "hyperGammaMeanDN[2]"), mcmcResult = mcmcDD)
  
  mcmcDataFrame <- summary(mcmcDD)$summary
  # use mean as estimation, not mode
  pars0 = mcmcDataFrame[grep("hyper|pi|A|B|C", row.names(mcmcDataFrame)), ]
  
  # pars0 = estimatePars(pars = mcmcDD@sim$fnames_oi,
  #                      mcmcResult = mcmcDD)
  parsFDR <- list(A = pars0[,1][2],
                  B = pars0[,1][3],
                  C = pars0[,1][4],
                  pi0 = pars0[,1][1],
                  gammaMeanDN = pars0[,1][5:6],
                  betaDN = pars0[,1][7:8],
                  nfamily = rep(samplenumber, 2),
                  rankPercentile = rankPercentile)

  dataFDR <- calculateFDR(pars = parsFDR,
                          dnData = dataDN1,
                          mutData = mutRate1,
                          geneName = geneset)
  result = list("mcmcDD"=mcmcDD, "pars0"=pars0, "dataFDR"=dataFDR)
  return(result)
}
