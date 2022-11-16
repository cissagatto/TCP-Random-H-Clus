###############################################################################
# TEST COMMUNITIES PARTITIONS
# Copyright (C) 2022
#
# This code is free software: you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option)
# any later version. This code is distributed in the hope that it will be
# useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
# Public License for more details.
#
# Elaine Cecilia Gatto | Prof. Dr. Ricardo Cerri | Prof. Dr. Mauri Ferrandin
# Federal University of Sao Carlos (UFSCar: https://www2.ufscar.br/) Campus
# Sao Carlos Computer Department (DC: https://site.dc.ufscar.br/)
# Program of Post Graduation in Computer Science
# (PPG-CC: http://ppgcc.dc.ufscar.br/)
# Bioinformatics and Machine Learning Group
# (BIOMAL: http://www.biomal.ufscar.br/)
#
###############################################################################


###############################################################################
# SET WORKSAPCE                                                               #
###############################################################################
FolderRoot = "~/TCP-Random-H-Clus"
FolderScripts = paste(FolderRoot, "/R", sep="")


##############################################################################
# FUNCTION MOUNT HYBRID PARTITION
##############################################################################
silho.build.test <- function(parameters){

  f = 1
  buildParalel <- foreach(f = 1:parameters$Number.Folds) %dopar%{
    #while(f<=parameters$Number.Folds){

    cat("\n#=========================================================")
    cat("\n# Fold: ", f)
    cat("\n#=========================================================")

    ########################################################################
    FolderRoot = "~/TCP-Random-H-Clus"
    FolderScripts = paste(FolderRoot, "/R", sep="")

    setwd(FolderScripts)
    source("utils.R")

    setwd(FolderScripts)
    source("libraries.R")


    #########################################################################
    converteArff <- function(arg1, arg2, arg3){
      str = paste("java -jar ", parameters$FoldersfolderUtils,
                  "/R_csv_2_arff.jar ", arg1, " ", arg2, " ", arg3, sep="")
      print(system(str))
      cat("\n")
    }

    ##########################################################################
    constroiParticoes <- function(TotalParticoes){

      data <- data.frame(matrix(NA,    # Create empty data frame
                                nrow = TotalParticoes,
                                ncol = 2))

      names(data) = c("numberPartition", "numberGroup")

      i = 1
      a = 1
      while(i<=nrow(data)){
        data[i,1] = a + 1
        data[i,2] = a + 1
        i = i + 1
        a = a + 1
        gc()
      }

      return(data)

    }


    ###########################################################################
    FolderPSplit = paste(parameters$Folders$folderPartitions,
                         "/Split-", f, sep="")

    FolderSplitComm = paste(parameters$Folders$folderCommunities,
                            "/Split-", f, sep="")

    FolderTSplit = paste(parameters$Folders$folderTestSilho,
                         "/Split-", f, sep="")
    if(dir.exists(FolderTSplit)==FALSE){dir.create(FolderTSplit)}

    ########################################################################
    #cat("\nOpen Train file ", f, "\n")
    setwd(parameters$Folders$folderCVTR)
    nome_arq_tr = paste(parameters$Dataset.Name, "-Split-Tr-", f,
                        ".csv", sep="")
    arquivo_tr = data.frame(read.csv(nome_arq_tr))

    #####################################################################
    #cat("\nOpen Validation file ", f, "\n")
    setwd(parameters$Folders$folderCVVL)
    nome_arq_vl = paste(parameters$Dataset.Name, "-Split-Vl-", f,
                        ".csv", sep="")
    arquivo_vl = data.frame(read.csv(nome_arq_vl))

    ########################################################################
    #cat("\nOpen Test file ", f, "\n")
    setwd(parameters$Folders$folderCVTS)
    nome_arq_ts = paste(parameters$Dataset.Name, "-Split-Ts-", f,
                        ".csv", sep="")
    cat("\n\t\t", nome_arq_ts)
    arquivo_ts = data.frame(read.csv(nome_arq_ts))

    # juntando treino com validação
    arquivo_tr2 = rbind(arquivo_tr, arquivo_vl)

    #######################################################3
    all = data.frame(parameters$Choosed$all.methods.choosed)
    all.2 = filter(all, all$split == f)
    total = nrow(all.2)


    all = data.frame(parameters$Choosed$all.methods.choosed)
    all.2 = data.frame(filter(all, all$split == f))

    partitions = data.frame(parameters$Choosed$all.partitions.choosed)
    partitions = data.frame(filter(partitions, fold == f))

    best = data.frame(parameters$best.silhouette$all.silhouette)
    best = data.frame(filter(best, best$fold==f))

    num.part = best$part
    num.groups = best$part

    partitions2 = data.frame(filter(partitions, partitions$method == all.2$method))
    partitions2 = partitions2[,c(-1,-2)]
    labels = partitions2$labels
    groups = partitions2[,num.part]
    particao = data.frame(labels, groups)

    g = 1
    while(g<=num.groups){

      cat("\n#=========================================================")
      cat("\n# Group = ", g)
      cat("\n#=========================================================")

      ###############################################################
      FolderTestGroup = paste(FolderTSplit, "/Group-", g, sep="")
      if(dir.exists(FolderTestGroup)== FALSE){dir.create(FolderTestGroup) }


      ###############################################################
      #cat("\nSpecific Group: ", g, "\n")
      grupoEspecifico = filter(particao, groups == g)


      ################################################################
      cat("\nTRAIN: Mount Group ", g, "\n")
      atributos_tr = arquivo_tr2[parameters$Dataset.Info$AttStart:parameters$Dataset.Info$AttEnd]
      n_a = ncol(atributos_tr)
      classes_tr = select(arquivo_tr2, grupoEspecifico$labels)
      n_c = ncol(classes_tr)
      grupo_tr = cbind(atributos_tr, classes_tr)
      fim_tr = ncol(grupo_tr)


      #####################################################################
      #cat("\n\tTRAIN: Save Group", g, "\n")
      setwd(FolderTestGroup)
      nome_tr = paste(parameters$Dataset.Name, "-split-tr-", f, "-group-",
                      g, ".csv", sep="")
      write.csv(grupo_tr, nome_tr, row.names = FALSE)


      #####################################################################
      #cat("\n\tINICIO FIM TARGETS: ", g, "\n")
      inicio = parameters$Dataset.Info$LabelStart
      fim = fim_tr
      ifr = data.frame(inicio, fim)
      write.csv(ifr, "inicioFimRotulos.csv", row.names = FALSE)


      ####################################################################
      #cat("\n\tTRAIN: Convert Train CSV to ARFF ", g , "\n")
      nome_arquivo_2 = paste(parameters$Dataset.Name, "-split-tr-", f,
                             "-group-", g, ".arff", sep="")
      arg1Tr = nome_tr
      arg2Tr = nome_arquivo_2
      arg3Tr = paste(inicio, "-", fim, sep="")
      str = paste("java -jar ", parameters$Folders$folderUtils,
                  "/R_csv_2_arff.jar ", arg1Tr, " ", arg2Tr, " ",
                  arg3Tr, sep="")
      print(system(str))


      ##################################################################
      #cat("\n\tTRAIN: Verify and correct {0} and {1} ", g , "\n")
      arquivo = paste(FolderTestGroup, "/", arg2Tr, sep="")
      str0 = paste("sed -i 's/{0}/{0,1}/g;s/{1}/{0,1}/g' ", arquivo, sep="")
      print(system(str0))


      ######################################################################
      #cat("\n\tTEST: Mount Group: ", g, "\n")
      atributos_ts = arquivo_ts[parameters$Dataset.Info$AttStart:parameters$Dataset.Info$AttEnd]
      classes_ts = select(arquivo_ts, grupoEspecifico$labels)
      grupo_ts = cbind(atributos_ts, classes_ts)
      fim_ts = ncol(grupo_ts)


      ######################################################################
      #cat("\n\tTEST: Save Group ", g, "\n")
      setwd(FolderTestGroup)
      nome_ts = paste(parameters$Dataset.Name, "-split-ts-", f, "-group-",
                      g, ".csv", sep="")
      write.csv(grupo_ts, nome_ts, row.names = FALSE)


      ####################################################################
      #cat("\n\tTEST: Convert CSV to ARFF ", g , "\n")
      nome_arquivo_3 = paste(parameters$Dataset.Name, "-split-ts-", f,
                             "-group-", g , ".arff", sep="")
      arg1Ts = nome_ts
      arg2Ts = nome_arquivo_3
      arg3Ts = paste(inicio, "-", fim, sep="")
      str = paste("java -jar ", parameters$Folders$folderUtils,
                  "/R_csv_2_arff.jar ", arg1Ts, " ", arg2Ts, " ",
                  arg3Ts, sep="")
      system(str)


      #####################################################################
      #cat("\n\tTEST: Verify and correct {0} and {1} ", g , "\n")
      arquivo = paste(FolderTestGroup, "/", arg2Ts, sep="")
      str0 = paste("sed -i 's/{0}/{0,1}/g;s/{1}/{0,1}/g' ", arquivo, sep="")
      cat("\n")
      print(system(str0))
      cat("\n")


      #####################################################################
      #cat("\nCreating .s file for clus")
      if(inicio == fim){

        nome_config = paste(parameters$Dataset.Name, "-split-", f, "-group-",
                            g, ".s", sep="")
        sink(nome_config, type = "output")

        cat("[General]")
        cat("\nCompatibility = MLJ08")

        cat("\n\n[Data]")
        cat(paste("\nFile = ", nome_arquivo_2, sep=""))
        cat(paste("\nTestSet = ", nome_arquivo_3, sep=""))

        cat("\n\n[Attributes]")
        cat("\nReduceMemoryNominalAttrs = yes")

        cat("\n\n[Attributes]")
        cat(paste("\nTarget = ", fim, sep=""))
        cat("\nWeights = 1")

        cat("\n")
        cat("\n[Tree]")
        cat("\nHeuristic = VarianceReduction")
        cat("\nFTest = [0.001,0.005,0.01,0.05,0.1,0.125]")

        cat("\n\n[Model]")
        cat("\nMinimalWeight = 5.0")

        cat("\n\n[Output]")
        cat("\nWritePredictions = {Test}")
        cat("\n")
        sink()

        ###################################################################
        cat("\nExecute CLUS: ", g , "\n")
        nome_config2 = paste(FolderTestGroup, "/", nome_config, sep="")
        str = paste("java -jar ", parameters$Folders$folderUtils,
                    "/Clus.jar ", nome_config2, sep="")
        print(system(str))
        cat("\n")

      } else {

        nome_config = paste(parameters$Dataset.Name, "-split-", f,
                            "-group-", g, ".s", sep="")
        sink(nome_config, type = "output")

        cat("[General]")
        cat("\nCompatibility = MLJ08")

        cat("\n\n[Data]")
        cat(paste("\nFile = ", nome_arquivo_2, sep=""))
        cat(paste("\nTestSet = ", nome_arquivo_3, sep=""))

        cat("\n\n[Attributes]")
        cat("\nReduceMemoryNominalAttrs = yes")

        cat("\n\n[Attributes]")
        cat(paste("\nTarget = ", inicio, "-", fim, sep=""))
        cat("\nWeights = 1")

        cat("\n")
        cat("\n[Tree]")
        cat("\nHeuristic = VarianceReduction")
        cat("\nFTest = [0.001,0.005,0.01,0.05,0.1,0.125]")

        cat("\n\n[Model]")
        cat("\nMinimalWeight = 5.0")

        cat("\n\n[Output]")
        cat("\nWritePredictions = {Test}")
        cat("\n")
        sink()

        cat("\nExecute CLUS: ", g , "\n")
        nome_config2 = paste(FolderTestGroup, "/", nome_config, sep="")
        str = paste("java -jar ", parameters$Folders$folderUtils,
                    "/Clus.jar ", nome_config2, sep="")
        print(system(str))
        cat("\n")

      }

      ##################################################################
      #cat("\n\nOpen predictions")
      nomeDoArquivo = paste(FolderTestGroup, "/", parameters$Dataset.Name,
                            "-split-", f,"-group-", g,
                            ".test.pred.arff", sep="")
      predicoes = data.frame(foreign::read.arff(nomeDoArquivo))


      #####################################################################
      #cat("\nS\nPLIT PREDICTIS")
      if(inicio == fim){
        #cat("\n\nOnly one label in this group")

        ###################################################################
        #cat("\n\nSave Y_true")
        setwd(FolderTestGroup)
        classes = data.frame(predicoes[,1])
        names(classes) = colnames(predicoes)[1]
        write.csv(classes, "y_true.csv", row.names = FALSE)

        #################################################################
        #cat("\n\nSave Y_true")
        rot = paste("Pruned.p.", colnames(predicoes)[1], sep="")
        pred = data.frame(predicoes[,rot])
        names(pred) = colnames(predicoes)[1]
        setwd(FolderTestGroup)
        write.csv(pred, "y_predict.csv", row.names = FALSE)

        ####################################################################
        rotulos = c(colnames(classes))
        n_r = length(rotulos)
        gc()

      } else {

        ##############################################################
        #cat("\n\nMore than one label in this group")
        comeco = 1+(fim - inicio)


        ####################################################################
        cat("\n\nSave Y_true")
        classes = data.frame(predicoes[,1:comeco])
        setwd(FolderTestGroup)
        write.csv(classes, "y_true.csv", row.names = FALSE)


        ##################################################################
        cat("\n\nSave Y_true")
        rotulos = c(colnames(classes))
        n_r = length(rotulos)
        nomeColuna = c()
        t = 1
        while(t <= n_r){
          nomeColuna[t] = paste("Pruned.p.", rotulos[t], sep="")
          t = t + 1
          gc()
        }
        pred = data.frame(predicoes[nomeColuna])
        names(pred) = rotulos
        setwd(FolderTestGroup)
        write.csv(pred, "y_predict.csv", row.names = FALSE)
        gc()
      } # FIM DO ELSE

      # deleting files
      um = paste(parameters$Dataset.Name, "-split-", f, "-group-", g, ".model", sep="")
      dois = paste(parameters$Dataset.Name, "-split-", f, "-group-", g, ".s", sep="")
      tres = paste(parameters$Dataset.Name, "-split-tr-", f, "-group-", g, ".arff", sep="")
      quatro = paste(parameters$Dataset.Name, "-split-ts-", f, "-group-", g, ".arff", sep="")
      cinco = paste(parameters$Dataset.Name, "-split-tr-", f, "-group-", g, ".csv", sep="")
      seis = paste(parameters$Dataset.Name, "-split-ts-", f, "-group-", g, ".csv", sep="")
      sete = paste(parameters$Dataset.Name, "-split-", f, "-group-", g, ".out", sep="")
      oito = paste("Variance_RHE_1.csv")

      setwd(FolderTestGroup)
      unlink(um, recursive = TRUE)
      unlink(dois, recursive = TRUE)
      unlink(tres, recursive = TRUE)
      unlink(quatro, recursive = TRUE)
      unlink(cinco, recursive = TRUE)
      unlink(seis, recursive = TRUE)
      unlink(sete, recursive = TRUE)
      unlink(oito, recursive = TRUE)

      g = g + 1
      gc()
    } # fim do grupo

    #f = f + 1
    gc()
  } # fim do for each

  gc()
  cat("\n###########################################################")
  cat("\n# TEST SILHO: Build and Test Hybrid Partitions End        #")
  cat("\n###########################################################")
  cat("\n\n\n\n")
}


##############################################################################
# FUNCTION SPLITS PREDCTIONS HYBRIDS
#############################################################################
silho.gather.predicts <- function(parameters){

  f = 1
  gatherR <- foreach(f = 1:parameters$Number.Folds) %dopar%{
    #while(f<=parameters$Number.Folds){

    cat("\n#========================================================#")
    cat("\n# Fold: ", f, "                                          #")
    cat("\n#========================================================#")

    #########################################################################
    FolderRoot = "~/TCP-Random-H-Clus"
    FolderScripts = paste(FolderRoot, "/R", sep="")

    setwd(FolderScripts)
    source("utils.R")

    setwd(FolderScripts)
    source("libraries.R")

    ##########################################################################
    constroiParticoes <- function(TotalParticoes){

      data <- data.frame(matrix(NA,    # Create empty data frame
                                nrow = TotalParticoes,
                                ncol = 2))

      names(data) = c("numberPartition", "numberGroup")

      i = 1
      a = 1
      while(i<=nrow(data)){
        data[i,1] = a + 1
        data[i,2] = a + 1
        i = i + 1
        a = a + 1
        gc()
      }

      return(data)

    }

    ###################################################################

    # "/dev/shm/j-GpositiveGO/Partitions/Split-1"
    FolderPSplit = paste(parameters$Folders$folderPartitions,
                         "/Split-", f, sep="")

    # "/dev/shm/j-GpositiveGO/Test-Silho/Split-1"
    FolderTSplit = paste(parameters$Folders$folderTestSilho,
                         "/Split-", f, sep="")

    # "/dev/shm/j-GpositiveGO/Communities/Split-1"
    FolderCommSplit = paste(parameters$Folders$folderCommunities,
                            "/Split-", f, sep="")

    #######################################################3
    all = data.frame(parameters$Choosed$all.methods.choosed)
    all.2 = filter(all, all$split == f)
    total = nrow(all.2)

    all = data.frame(parameters$Choosed$all.methods.choosed)
    all.2 = data.frame(filter(all, all$split == f))

    partitions = data.frame(parameters$Choosed$all.partitions.choosed)
    partitions = data.frame(filter(partitions, fold == f))

    best = data.frame(parameters$best.silhouette$all.silhouette)
    best = data.frame(filter(best, best$fold==f))

    num.part = best$part
    num.groups = best$part


    ################################################################
    apagar = c(0)
    y_true = data.frame(apagar)
    y_pred = data.frame(apagar)

    # GROUP
    g = 1
    while(g<=num.groups){

      cat("\n#=========================================================")
      cat("\n# Group = ", g)
      cat("\n#=========================================================")

      FolderTestGroup = paste(FolderTSplit, "/Group-", g, sep="")

      #cat("\n\nGather y_true ", g)
      setwd(FolderTestGroup)
      #setwd(FolderTG)
      y_true_gr = data.frame(read.csv("y_true.csv"))
      y_true = cbind(y_true, y_true_gr)

      setwd(FolderTestGroup)
      #setwd(FolderTG)
      #cat("\n\nGather y_predict ", g)
      y_pred_gr = data.frame(read.csv("y_predict.csv"))
      y_pred = cbind(y_pred, y_pred_gr)

      #cat("\n\nDeleting files")
      unlink("y_true.csv", recursive = TRUE)
      unlink("y_predict.csv", recursive = TRUE)
      unlink("inicioFimRotulos.csv", recursive = TRUE)

      g = g + 1
      gc()
    } # FIM DO GRUPO

    #cat("\n\nSave files ", g, "\n")
    setwd(FolderTSplit)
    y_pred = y_pred[,-1]
    y_true = y_true[,-1]
    write.csv(y_pred, "y_predict.csv", row.names = FALSE)
    write.csv(y_true, "y_true.csv", row.names = FALSE)

    #f = f + 1
    gc()
  } # fim do foreach

  gc()
  cat("\n###############################################################")
  cat("\n# Gather Predicts: END                                        #")
  cat("\n###############################################################")
  cat("\n\n\n\n")

} # fim da função


#############################################################################
# FUNCTION EVALUATION HYBRID PARTITIONS
############################################################################
silho.evaluate.test <- function(parameters){

  f = 1
  avalParal <- foreach(f = 1:parameters$Number.Folds) %dopar%{
    #while(f<=parameters$Number.Folds){

    cat("\n#=========================================================")
    cat("\n# Fold: ", f)
    cat("\n#=========================================================")

    folders = list()

    ##############################################################
    FolderRoot = "~/TCP-Random-H-Clus"
    FolderScripts = paste(FolderRoot, "/R/", sep="")

    setwd(FolderScripts)
    source("utils.R")

    setwd(FolderScripts)
    source("libraries.R")

    #################################################################
    # "/dev/shm/j-GpositiveGO/Partitions/Split-1"
    FolderPSplit = paste(parameters$Folders$folderPartitions,
                         "/Split-", f, sep="")

    # "/dev/shm/j-GpositiveGO/Test-Silho/Split-1"
    FolderTSplit = paste(parameters$Folders$folderTestSilho,
                         "/Split-", f, sep="")

    # "/dev/shm/j-GpositiveGO/Communities/Split-1"
    FolderCommSplit = paste(parameters$Folders$folderCommunities,
                            "/Split-", f, sep="")

    #######################################################3
    all = data.frame(parameters$Choosed$all.methods.choosed)
    all.2 = filter(all, all$split == f)
    total = nrow(all.2)

    all = data.frame(parameters$Choosed$all.methods.choosed)
    all.2 = data.frame(filter(all, all$split == f))

    partitions = data.frame(parameters$Choosed$all.partitions.choosed)
    partitions = data.frame(filter(partitions, fold == f))

    best = data.frame(parameters$best.silhouette$all.silhouette)
    best = data.frame(filter(best, best$fold==f))

    #cat("\nData frame")
    apagar = c(0)
    confMatPartitions = data.frame(apagar)
    partitions = c()

    #cat("\nGet the true and predict lables")
    setwd(FolderTSplit)
    #setwd(FolderKnn)
    y_true = data.frame(read.csv("y_true.csv"))
    y_pred = data.frame(read.csv("y_predict.csv"))

    #cat("\nCompute measures multilabel")
    y_true2 = data.frame(sapply(y_true, function(x) as.numeric(as.character(x))))
    y_true3 = mldr_from_dataframe(y_true2 , labelIndices = seq(1,ncol(y_true2 )), name = "y_true2")
    y_pred2 = sapply(y_pred, function(x) as.numeric(as.character(x)))

    #cat("\nSave Confusion Matrix")
    setwd(FolderTSplit)
    #setwd(FolderKnn)
    salva3 = paste("Conf-Mat-Fold-", f, ".txt", sep="")
    sink(file=salva3, type="output")
    confmat = multilabel_confusion_matrix(y_true3, y_pred2)
    print(confmat)
    sink()

    #cat("\nCreating a data frame")
    confMatPart = multilabel_evaluate(confmat)
    confMatPart = data.frame(confMatPart)
    names(confMatPart) = paste("Fold-", f, sep="")
    namae = paste("Split-", f, "-Evaluated.csv", sep="")
    setwd(FolderTSplit)
    write.csv(confMatPart, namae)

    #cat("\nDelete files")
    setwd(FolderTSplit)
    unlink("y_true.csv", recursive = TRUE)
    unlink("y_predict.csv", recursive = TRUE)

    #f = f + 1
    gc()
  } # fim do for each

  gc()
  cat("\n############################################################")
  cat("\n# TEST SILHO: Evaluation Folds END                         #")
  cat("\n############################################################")
  cat("\n\n\n\n")
}


#################################################################
#
#################################################################
silho.build.data.frame <- function(){

  data <- data.frame(matrix(NA,    # Create empty data frame
                            nrow = 22,
                            ncol = 11))

  measures = c("accuracy", "average-precision", "clp", "coverage", "F1",
               "hamming-loss", "macro-AUC", "macro-F1", "macro-precision",
               "macro-recall", "margin-loss", "micro-AUC", "micro-F1",
               "micro-precision", "micro-recall", "mlp", "one-error",
               "precision", "ranking-loss", "recall", "subset-accuracy", "wlp")

  data$X1 = measures

  return(data)

}


##########################################################################
# FUNCTION GATHER EVALUATIONS
#######################################################################
silho.gather.evaluated <- function(parameters){

  # vector with names
  measures = c("accuracy","average-precision","clp","coverage","F1",
               "hamming-loss","macro-AUC", "macro-F1","macro-precision",
               "macro-recall","margin-loss","micro-AUC","micro-F1",
               "micro-precision","micro-recall","mlp","one-error",
               "precision","ranking-loss", "recall","subset-accuracy","wlp")

  # from fold = 1 to number_folders
  f = 1
  while(f<=parameters$Number.Folds){

    # data frame
    apagar = c(0)
    avaliadoFinal = data.frame(apagar, measures)
    avaliadoKnn = data.frame(apagar, measures)
    folds = c(0)
    threshold = c(0)
    nomesThreshold = c(0)
    nomesFolds = c(0)

    ######################################################################

    # "/dev/shm/j-GpositiveGO/Partitions/Split-1"
    FolderPSplit = paste(parameters$Folders$folderPartitions,
                         "/Split-", f, sep="")

    # "/dev/shm/j-GpositiveGO/Test-Silho/Split-1"
    FolderTSplit = paste(parameters$Folders$folderTestSilho,
                         "/Split-", f, sep="")

    #######################################################3
    all = data.frame(parameters$Choosed$all.methods.choosed)
    all.2 = filter(all, all$split == f)
    total = nrow(all.2)


    cat("\n#==========================")
    cat("\n# Fold \t", f)

    ####################################################################
    setwd(FolderTSplit)
    str = paste("Split-", f, "-Evaluated.csv", sep="")
    avaliado = data.frame(read.csv(str))
    names(avaliado)[1] = "medidas"
    avaliadoKnn = cbind(avaliadoKnn, avaliado[,2])
    #nomesThreshold[f] = paste("Fold-", f, sep="")
    #names(avaliadoKnn)[f+2] = nomesThreshold[f]
    #unlink(str, recursive = TRUE)

    avaliadoKnn = avaliadoKnn[,-1]
    setwd(FolderTSplit)
    write.csv(avaliadoKnn, paste("Evaluated-Fold-", f, ".csv", sep=""),
              row.names = FALSE)


    f = f + 1
    gc()

  } # end folds


  gc()
  cat("\n###############################################################")
  cat("\n# TEST SILHO: Gather Evaluations End                          #")
  cat("\n###############################################################")
  cat("\n\n\n\n")
}


#####################################################################
#
#####################################################################
silho.organize.evaluation <- function(parameters){

  measures = c("accuracy","average-precision","clp","coverage","F1",
               "hamming-loss","macro-AUC", "macro-F1","macro-precision",
               "macro-recall","margin-loss","micro-AUC","micro-F1",
               "micro-precision","micro-recall","mlp","one-error",
               "precision","ranking-loss", "recall","subset-accuracy","wlp")


  apagar = c(0)
  resultado = data.frame(apagar)
  nomes = c("")

  f = 1
  while(f<=10){
    cat("\n\tFold ", f)

    ######################################################################

    # "/dev/shm/j-GpositiveGO/Partitions/Split-1"
    FolderPSplit = paste(parameters$Folders$folderPartitions,
                         "/Split-", f, sep="")

    # "/dev/shm/j-GpositiveGO/Test-Silho/Split-1"
    FolderTSplit = paste(parameters$Folders$folderTestSilho,
                         "/Split-", f, sep="")

    setwd(FolderTSplit)
    nome = paste("Split-", f, "-Evaluated.csv",sep="")
    arquivo = data.frame(read.csv(nome))
    resultado = cbind(resultado, arquivo[,2])
    nomes[f] = paste("Fold-", f,sep="")

    f = f + 1
    gc()
  }

  resultado = resultado[,-1]
  names(resultado) = nomes
  resultado[is.na(resultado)] <- 0

  media = data.frame(apply(resultado, 1, mean))
  media = cbind(measures, media)
  names(media) = c("Measures", "Mean")

  resultado = cbind(measures, resultado)

  setwd(parameters$Folders$folderReports)
  write.csv(media, "Mean-10Folds.csv", row.names = FALSE)
  write.csv(resultado, "All-10folds.csv", row.names = FALSE)

  gc()
  cat("\n################################################################")
  cat("\n# TEST SILHO: Organize End                                     #")
  cat("\n################################################################")
  cat("\n\n\n\n")
}



##################################################################
#
##################################################################
testSilhouette <- function(parameters){


  cat("\n\n#############################################################")
  cat("\n# RUN TEST SILHO: build and test partitions                 #")
  cat("\n#############################################################\n\n")
  timeBuild = system.time(resBT <- silho.build.test(parameters))


  cat("\n\n############################################################")
  cat("\n# RUN TEST SILHO: Matrix Confusion                         #")
  cat("\n############################################################\n\n")
  timeSplit = system.time(resGather <- silho.gather.predicts(parameters))


  cat("\n\n#############################################################")
  cat("\n# RUN TEST SILHO: Evaluation Fold                           #")
  cat("\n#############################################################\n\n")
  timeAvalia = system.time(resEval <- silho.evaluate.test(parameters))


  cat("\n\n#############################################################")
  cat("\n# RUN TEST SILHO: Gather Evaluation                         #")
  cat("\n#############################################################\n\n")
  timeGather = system.time(resGE <- silho.gather.evaluated(parameters))


  cat("\n\n###########################################################")
  cat("\n# RUN TEST SILHO: Organize Evaluation                     #")
  cat("\n###########################################################\n\n")
  timeOrg = system.time(resOA <- silho.organize.evaluation(parameters))


  cat("\n\n############################################################")
  cat("\n# RUN TEST SILHO: Save Runtime                             #")
  cat("\n############################################################\n\n")
  Runtime = rbind(timeBuild,
                  timeSplit,
                  timeAvalia,
                  timeGather,
                  timeOrg)
  setwd(parameters$Folders$folderTestSilho)
  write.csv(Runtime, paste(parameters$Dataset.Name,
                           "-test-silho-runtime-1.csv", sep=""),
            row.names = FALSE)
}


#########################################################################
# Please, any errors, contact us: elainececiliagatto@gmail.com
# Thank you very much!
#########################################################################
