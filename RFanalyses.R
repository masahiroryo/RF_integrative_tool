#####################################################################
#
#  Random Forest integrative analysis tool
#
#                                    2019.01.30 Masahiro Ryo
#####################################################################


#####################################################################
#####                    Data info                              #####
#####################################################################

    #_____reading______
    df = data.frame(read.csv("data_example.csv", na='NA', skip = 1))

    #_____data handling______   
    list_responses  = colnames(df)[c(1,2)]
    list_predictors = colnames(df)[c(3:ncol(df))] #[c(3:ncol(df))]

    #_____transforming_______
    list_factortransform = c()
    for(i in list_factortransform)  df[,i] = as.factor(df[,i])
    
    #_____parameters______
    # RF
    RF_nperm = 20
    RF_ntree = 30
    RF_ncore = 6      # parallel computing
    RF_varlim = -1    # estimation of p-value only for top X variables if specified (default: -1)
    
    # stats
    significance = 0.01   # significance level
    Bonferroni = F        # Bonferroni correction
    
    # plot
    plot_R2_min = 0.01 # minimum threshold of R2 to visualize plots
    plot_to_pdf = F      # to visualize figures on the console or pdf
                         #  'partial dependence plots.pdf'

    
    
    
#================================================================================
# Hereafter do not change anything unless you understand the structure.
#================================================================================
    
    
    
    
        
    
    
#####################################################################
#####                    Setting                                #####
#####################################################################
source("Functions.R")
    
#####################################################################
#####                    Modelling                              #####
#####################################################################

      
#-------------------------------------------
#  Setting
#-------------------------------------------
    if(plot_to_pdf) pdf("partial dependence plots.pdf", paper="a4r", width=9.5, height=7, pointsize=10)  # start printing to pdf
    result_summary = numeric()
    result_summary_pval = numeric()
    iloop = 0
#-------------------------------------------
#  Modeling each response variable w/ loop
#-------------------------------------------
    for (response in list_responses){
          #-------------------------------------------
          #  Setting up before modeling
          #-------------------------------------------
          #_____setting______
          iloop = iloop + 1
          
          #_____NA removal______
          df.o = df
          if(any(is.na(df[response]))) df.o = df[-which(is.na(df[response])),]

          #_____predictor selection______
          formula = as.formula(paste(response, paste(list_predictors, collapse=" + "), sep=" ~ "))
          
          #-------------------------------------------
          #  Modeling and assessment
          #-------------------------------------------
          #_____modeling______ 
          #             NOTE) change nperm and ntree accordingly
          model.Hapfelmeier.RF = tryCatch(
            {model.Hapfelmeier.RF = RF_permutation(formula, data=df.o, nperm=RF_nperm, ntree=RF_ntree, 
                                                   ncore=RF_ncore, varlim=RF_varlim, alpha=significance)}
            , error = function(e){return(NA)})
          
          #_____evaluating______ 
          # vim: variable importance measure
          vim = numeric(length(list_predictors))
          names(vim) = c(list_predictors)
          
          # p-value (Ho: predictor X has no contribution to improve model performance)
          pvals = rep(NA, c(length(list_predictors)))
          names(pvals) = c(list_predictors)
          
          # model performance: fitting and OOB validation scores
          fitting_score = NA; validation_score = NA
          
          if(any(is.na(model.Hapfelmeier.RF$p.values)==FALSE)){
            if(Bonferroni){
              for(i in names(model.Hapfelmeier.RF$varimp.bonf)) {vim[i] = round(model.Hapfelmeier.RF$varimp.R2.bonf[i],3)}
              for(i in names(model.Hapfelmeier.RF$varimp.bonf)) {pvals[i] = model.Hapfelmeier.RF$p.values.bonf[i]}
              fitting_score = model.Hapfelmeier.RF$fitting.bonf[2]
              validation_score = model.Hapfelmeier.RF$validation.bonf[2]
            }else{
              for(i in names(model.Hapfelmeier.RF$varimp)) {vim[i] = round(model.Hapfelmeier.RF$varimp.R2[i],3)}
              for(i in names(model.Hapfelmeier.RF$varimp)) {pvals[i] = model.Hapfelmeier.RF$p.values[i]}
              fitting_score = model.Hapfelmeier.RF$fitting[2]
              validation_score = model.Hapfelmeier.RF$validation[2]
            }
          }
          
          #_____saving the result______
          result_summary = rbind(result_summary, c(fitting_score, validation_score, vim))
          result_summary_pval = rbind(result_summary_pval, c(fitting_score, validation_score,pvals))

          
          #_____data_summary_______
          los = rep("ns",length(pvals))
          for (s in c(0.05,0.01,0.001)) los[which(pvals<s)] = c("*", "**", "***")[which(s==c(0.05,0.01,0.001))]

          gplt_data = data.frame(varname = names(vim), vim = vim, p = pvals, sign = as.factor(los), type = sapply(df[list_predictors], class))
          gplt_data = gplt_data[order(-gplt_data$vim),]
          
          
                    
          #-------------------------------------------
          #  Partial dependence plot
          #-------------------------------------------
          
          #_____Modeling______
          #             NOTE)  only with selected predictors 
          # 
          if(is.na(fitting_score)==FALSE){
            data.tmp = data.frame(cbind(df.o[,response], df.o[,names(which(vim>0))]))
            colnames(data.tmp) = c(response, names(which(vim>0)))
            task = makeRegrTask(data = data.tmp, target = response)
            learner = makeLearner(cl="regr.cforest", predict.type = "response", ntree = RF_ntree, mtry = ceiling(sqrt(length(colnames(data.tmp)))))
            fit = mlr::train(learner, task)
            
            
            #_____saving as a plot panel______
            #
            jloop = 0
            trend = rep("positive",length(pvals))
            for(i in gplt_data$varname[which(gplt_data$sign!="ns" & gplt_data$vim >plot_R2_min)]){
              jloop = jloop + 1
              pdp = generatePartialDependenceData(fit, task,  i, individual=F)
              pdp_data = data.frame(cbind((pdp$data[,2]), unlist(pdp$data[,1])))
              colnames(pdp_data) = c("x","y") 
              if(gplt_data$type[which(gplt_data$varname==i)]!="factor"){
                 if(pdp_data$y[which(pdp_data$x==max(pdp_data$x))] < pdp_data$y[which(pdp_data$x==min(pdp_data$x))]){
                    gplt_data[i,"vim"] = -gplt_data[i,"vim"]
                 }
              }
              g =             
                  ggplot(data=pdp_data, aes(x = x, y=y)) +
                  xlab(i) + ylab(response) + 
                  geom_line() +
                  geom_point() 
              eval(parse(text=(paste("g_", jloop, "= g", sep=""))))  
            }
            color_list  = c("ns" = "#00000040", "*" = "#F5CD6DA0", "**" = "#EE825AA0", "***" = "#B84543A0")
            
#-------------------------------------------
#  Visualization
#-------------------------------------------
            #______________________________________________________
            # 3.2. visualizing the relative importance 
            
            g_0 = 
              ggplot(data=gplt_data[which(gplt_data$sign!="ns" & abs(gplt_data$vim) >plot_R2_min),], aes(x = reorder(varname, abs(vim)), y=vim, fill=sign)) +
              xlab("Predictors") + ylab("R2 [%]") + 
              geom_bar(stat="identity") +
              geom_text(aes(label = sign))+
              scale_fill_manual(values = color_list) + 
              theme(legend.position = 'none') +
              coord_flip() 
            
            if(jloop==1) print(g_0 | g_1)
            if(jloop==2) print(g_0 | (g_1 / g_2))
            if(jloop==3) print(g_0 | (g_1 / g_2 / g_3))
            if(jloop==4) print(g_0 | (g_1 / g_3)| (g_2 / g_4))
            if(jloop==5) print(g_0 | (g_1 / g_4)| (g_2 / g_5) | (g_3 / plot_spacer()))
            if(jloop==6) print(g_0 | (g_1 / g_4)| (g_2 / g_5) | (g_3 / g_6))
            if(jloop==7) print(g_0 | (g_1 / g_4 / g_7) | (g_2 / g_5 / plot_spacer()) | (g_3 / g_6 / plot_spacer()))
            if(jloop==8) print(g_0 | (g_1 / g_4 / g_7) | (g_2 / g_5 / g_8) | (g_3 / g_6 / plot_spacer()))
            if(jloop>=9) print(g_0 | (g_1 / g_4 / g_7) | (g_2 / g_5 / g_8) | (g_3 / g_6 / g_9))
            #if(jloop>=10) print(g_0 | (g_1 / g_6)| (g_2 / g_7)| (g_3 / g_8)| (g_4 / g_9)| (g_5 / g_10))
            Sys.sleep(0)
            
            suppressMessages(for(j in c(0:jloop)) try(eval(parse(text=(paste("remove(g_", j, ")", sep=""))))))
          } # end: if(is.na(fitting_score)==FALSE)
          
    } # end: for (response in list_responses)

    
#####################################################################
#####                    Output                                 #####
#####################################################################
    
    #_____Excel output & device off_____
    colnames(result_summary) = c("R2-fitting","R2-validation", list_predictors)
    colnames(result_summary_pval) = c("R2-fitting","R2-validation", list_predictors)
    
    write.csv(round(result_summary, digits = 3), "variable_importance.csv", row.names = TRUE)
    write.csv(round(result_summary_pval, digits = 5), "p-value.csv", row.names = TRUE)
    if(plot_to_pdf) dev.off() # finish printing to pdf


    