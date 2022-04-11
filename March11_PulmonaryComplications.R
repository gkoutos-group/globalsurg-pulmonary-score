#!/usr/bin/env Rscript

#not getting same results in slurm and here

#http://www.svds.com/learning-imbalanced-classes/
#https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0118432
#https://towardsdatascience.com/gaining-an-intuitive-understanding-of-precision-and-recall-3b9df37804a7
#https://github.com/stephenbates19/nestedcv/blob/master/R/core.R

################

#Description: 3 scripts
# - FinalPreprocessing.Rmd (processing of datasets, missing values and table 1)
# - March11_PulmonaryComplications.R (main manuscript)
# - FinalExtraApp.R (processing of results from previous script to fit what is needed in the app, including 
#bar plot and tables)

##########################################################
###################### Directories ######################
##########################################################

mainDir <- "/rds/projects/g/gkoutosg-variant-prediction/VictorSaisa/FinalCode"
NameRun <- 'ResultsApril2022'
subDir <-
  paste0(sub('\\..*', '', NameRun), format(Sys.time(), '_%Y%m%d_%H%M'))
dir.create(file.path(mainDir, subDir))
IncludeFigHere <- file.path(mainDir, subDir)

##########################################################
###################### Libraries ########################
##########################################################

.libPaths(c(.libPaths(), "/rds/homes/l/lxb732/R/library/4.0.0/haswell"))
#.libPaths(c(.libPaths(), "/rds/homes/l/lxb732/R/x86_64-pc-linux-gnu-library/4.0" ))
#.libPaths(c(.libPaths(), "/rds/bear-apps/2020a/EL8-cas/software/R/4.0.0-foss-2020a/lib64/R/library" ))

Sys.setenv(RSTUDIO_PANDOC = "/rds/bear-apps/2020a/EL8-cas/software/rstudio/1.2.5042-foss-2020a-Java-11-R-4.0.0/bin/pandoc")

library(workflowsets)
library(tidymodels)
library(data.table)
library(ggplot2)
library(dplyr)
library(finetune)
library(parsnip)
library(baguette)
library(readr)
library(DALEX)
library(fastshap) #### delete?
library(vip)
library(stringr)
library(rms)
library(gridExtra)
library(grid)
library(patchwork)
library(progress)
library(modEvA)
library(boot)
library(doParallel)
library(plotROC)
library(pROC)
library(patchwork)
library(gridExtra)
library(grid)
library(ggpubr)
library(finalfit)
library(pmsampsize)

#setwd("/rds/projects/g/gkoutosg-variant-prediction/VictorSaisa/FinalCode") #problems rendering Rmd if not wd?
rmarkdown::render("/rds/projects/g/gkoutosg-variant-prediction/VictorSaisa/FinalCode/FinalPreprocessing2.Rmd") #2 is the newdataset

.libPaths()
##########################################################
###################### Data Loading ######################
##########################################################

# Preprocessing
dataFull <- read.csv(paste0(IncludeFigHere,"/CombinedThreeFinalDatasetReconInElectiveLabelChangeNewVars.csv"))
#dataFull <- read_csv("/rds/projects/g/gkoutosg-variant-prediction/Laura/Magda/JamesCodeFuns/CombinedThreeFinalDatasetReconInElectiveLabelChangeNewVars.csv")  #withoutNewVars

table(dataFull$Type)

dim(dataFull)

dataFull <- dataFull %>%
  mutate(
    Anaesthetic_type = recode(
      Anaesthetic_type,
      "General anaesthetic" = "General",
      "Regional/local anaesthetic" = "Regional/local"
    )
  ) %>%
  filter(Type != "recon") %>%
  select(-c("X")) %>% #X1
  mutate_if(is.character, as.factor) %>%
  #filter(!(Death == "Yes" & Label == "Negative")) %>% #509 deleted
  #filter(!(is.na(Death) & Label == "Negative")) %>% #27 deleted
  select(-Death)

table(dataFull$Type)

dataFull <- dataFull %>%
  select(-c(Country_income, Community_SARSCoV2_risk)) #patienttype || IncomeComRisk deleted before
dataFull$Age <- relevel(dataFull$Age, ref = "30-39 years") #30-39 years #18-29 years
dataFull$Current_smoker <- relevel(dataFull$Current_smoker, ref = "No")
dataFull$RCRI <- relevel(dataFull$RCRI, ref = "0")
dataFull$Preoperative_SARSCoV2_test <- relevel(dataFull$Preoperative_SARSCoV2_test, ref = "Negative")
dataFull$Preexisting_respiratory_disease <- relevel(dataFull$Preexisting_respiratory_disease, ref = "No")
dataFull$Sex <- relevel(dataFull$Sex, ref = "Female")
dataFull$BUPA_operation_grade <- relevel(dataFull$BUPA_operation_grade, ref = "Minor")
dataFull$Anaesthetic_type <- relevel(dataFull$Anaesthetic_type, ref = "Regional/local")
dataFull$Country_income_SARSCoV2_risk <- relevel(dataFull$Country_income_SARSCoV2_risk, ref = "High income_Low")

table(dataFull$Type)

# Data separation

dataFull2 <- dataFull %>%
  mutate(.row = row_number())

dataFull3 <- dataFull2 %>%
  select(-Type)

indices <-
  list(analysis   = dataFull2$.row[dataFull2$Type == "Validation"] ,
       assessment = dataFull2$.row[dataFull2$Type == "Development"])

split <- make_splits(indices, dataFull3 %>% select(-.row))

WholeDevelopment <- training(split)

split2 <- initial_split(WholeDevelopment, prop = 0.75, strata = Label)

Test <- split2 %>% testing() #Testing split derivation

data_trn  <- split2 %>% training() #Training split derivation 

data_tst <- testing(split) #Validation 1 (SurgWeek)

##########################################################
###################### Table 2 selected vars ###########################
##########################################################

#data_trnTable <- data_trn
#
## edited version
##data_trnTable <- data_trnTable %>% 
##  mutate(BUPA_operation_grade = factor(BUPA_operation_grade, levels = c("Minor", "Major"))) should not be needed as ref before
#
#explanatory = names(data_trnTable) %>% setdiff(c(
#  "Current_smoker",
#  "Preoperative_SARSCoV2_test",
#  "Hospital_type",
#  "Label"
#))
#dependent = "Label"
##data_trn %>%
##  summary_factorlist(dependent, explanatory) -> t
##
##write.csv(t, (paste0(IncludeFigHere, "/Table2.csv")))
#
#data_trnTable %>%
#  finalfit(dependent, explanatory, column = TRUE)  %>% #metrics = TRUE
#  ff_remove_p() -> t
#
#write.csv(t, (paste0(IncludeFigHere, "/Table2Univariate.csv")))

##########################################################
###################### Table 2 all vars ###########################
##########################################################

dim(data_trn)
print(head(data_trn))

data_trnTable <- data_trn

s <- filter(data_trn, Age == "30-39 years")
table(s$Label)

# edited version
#data_trnTable <- data_trnTable %>% 
#  mutate(BUPA_operation_grade = factor(BUPA_operation_grade, levels = c("Minor", "Major")))

explanatory = names(data_trnTable) %>% setdiff(c(
  #"Current_smoker",
  #"Preoperative_SARSCoV2_test",
  #"Hospital_type",
  "Label"
))
dependent = "Label"
#data_trn %>%
#  summary_factorlist(dependent, explanatory) -> t
#
#write.csv(t, (paste0(IncludeFigHere, "/Table2.csv")))

data_trnTable %>%
  finalfit(dependent, explanatory, column = TRUE)  %>% #metrics = TRUE
  ff_remove_p() -> t

print(head(t))

write.csv(t, (paste0(IncludeFigHere, "/Table2UnivariateAll.csv")))

##########################################################
###################### Model #############################
##########################################################

lasso_reg_spec <-
  logistic_reg(penalty = tune(), mixture = 1) %>%
  set_engine("glmnet")

xgb_spec <-
  boost_tree(
    tree_depth = tune(),
    learn_rate = tune(),
    loss_reduction = tune(),
    min_n = tune(),
    sample_size = tune(),
    trees = 1000
  ) %>% #tree_depth = tune(),#learn_rate = tune(), loss_reduction = tune(),
  set_engine("xgboost") %>%
  set_mode("classification")

##########################################################
###################### Recepies ##########################
##########################################################

#tried smote too in previous iterations

AllRec <- recipe(Label ~ ., data = data_trn)
PreparedPreProc <- AllRec %>% prep()
AllRecJ <- juice(PreparedPreProc)

Standard <- AllRec %>%
  step_dummy(all_nominal(),-all_outcomes(), one_hot = FALSE) #%>%

PreparedPreProc <- Standard %>% prep()
StandardMatrix <- juice(PreparedPreProc)

Themis <- AllRec %>%
  step_dummy(all_nominal(),-all_outcomes(), one_hot = FALSE) %>%
  themis::step_smote(Label, over_ratio = 0.5) #before no over_ratio

PreparedPreProc <- Themis %>% prep()
ThemisMatrix <- juice(PreparedPreProc)

XGBoostRec <- AllRec %>%
  step_dummy(all_nominal(),-all_outcomes(), one_hot = FALSE)

PreparedPreProc <- XGBoostRec %>% prep()
XGBoostMat <- juice(PreparedPreProc)


recipe_list <-
  list(StandardMatrix = Standard,  ThemisMatrix = Themis)


model_list <-
  list(lasso = lasso_reg_spec , xgboost = xgb_spec)

model_set <-
  workflow_set(preproc = recipe_list,
               models = model_list,
               cross = T)

model_set <- model_set %>%
  anti_join(tibble(
    wflow_id = c("ThemisMatrix_lasso", "ThemisMatrix_xgboost"),
    #take smote analysis out
    by = "wflow_id"
  ))

#"StandardMatrix_lasso"   "StandardMatrix_xgboost" "ThemisMatrix_lasso"     "ThemisMatrix_xgboost"

##########################################################
###################### Run Models Racing #################
##########################################################

#print(Sys.getenv("SLURM_NTASKS"))

registerDoParallel(cores = (Sys.getenv("SLURM_NTASKS"))) #slurm
#registerDoParallel(cores=20)
#cl <- makePSOCKcluster(Sys.getenv("SLURM_NTASKS")) # is this really working?
#registerDoParallel(20)
#stopCluster(cl)

#xgboost_grid <- xgb_spec %>%
#  parameters() %>%
#  grid_max_entropy(size = 5)
#
#lasso_grid <- lasso_reg_spec  %>%
#  parameters() %>%
#  grid_max_entropy(size = 5)

#Grid <- cbind(xgboost_grid,lasso_grid )

CV <- 10

cv_splits <-
  rsample::vfold_cv(data_trn, strata = Label, v = CV) #go up to 20

multi_metric2 <-
  metric_set(roc_auc) #https://www.tmwr.org/performance.html #f_meas

get_model <- function(x) {
  extract_fit_parsnip(x) %>%
    vi()
}

race_ctrl <-
  control_race(
    save_pred = TRUE,
    parallel_over = "everything",
    save_workflow = TRUE,
    extract = get_model
  )

race_results <-
  model_set %>%
  workflow_map(
    "tune_race_anova",
    seed = 1503,
    resamples = cv_splits,
    grid = 10,
    #10
    control = race_ctrl,
    verbose = TRUE,
    metrics = multi_metric2
  )

#,
#c = multi_metric2

save(
  race_results,
  model_set,
  cv_splits,
  IncludeFigHere,
  multi_metric2,
  file = paste0(IncludeFigHere, "/Modeling2data3modsThree.RData")
)

#library(pROC)
#
#LASSO <- Pred %>%
#    filter(.config == "Preprocessor1_Model01") %>%
#    group_by(id) %>%
#    nest() %>%
#    mutate(roc = map(data, ~roc(.,Label, .pred_Positive))) %>%
#    mutate(ci = map(roc, ~pROC::ci(.,method = "bootstrap")))
#


##########################################################
######################### Model Performamce CV ###########
##########################################################

#https://www.tmwr.org/resampling.html

LASSO_Best <- select_best(race_results$result[[1]], metric = "roc_auc")$.config
XGBoost_Best <- select_best(race_results$result[[2]], metric = "roc_auc")$.config


library(cvAUC)

LASSO <- collect_predictions(race_results$result[[1]]) %>%
  filter(.config == LASSO_Best)

XGBoost <- collect_predictions(race_results$result[[2]]) %>%
  filter(.config == XGBoost_Best)

CICVROC <- function(LASSO) {
  predictionsL <- LASSO$.pred_Positive #attributes dont seem necessary
  labelsL <- as.numeric((LASSO$Label)) - 1
  folds2 <-
    lapply((unique(LASSO$id)), function(x) {
      LASSO[LASSO$id == x, ]$.row
    })
  
  res <- ci.cvAUC(
    predictions = predictionsL,
    labels = labelsL,
    folds = folds2,
    confidence = 0.95
  )
  
  LASSOPlot <- LASSO %>%
    group_by(id) %>%
    nest() %>%
    mutate(roc = map(data, ~ roc(., Label, .pred_Positive)))
  
  Plot <-
    lapply(1:length(LASSOPlot$id), function(x) {
      data.frame(
        Resample = x,
        sensitivities = LASSOPlot$roc[[x]]$sensitivities,
        specificities = LASSOPlot$roc[[x]]$specificities
      )
    }) %>%
    bind_rows()
  
  return(list(Plot, res))
  
}

PlotLASSO <- CICVROC(LASSO)

PlotXGBoost <- CICVROC(XGBoost)

Plot <- rbind(
  PlotLASSO[[1]] %>%
    add_column(Type = "LASSO") ,
  PlotXGBoost[[1]] %>%
    add_column(Type = "XGBOOST")
)


pdf(paste0(IncludeFigHere, "/CrossValCI.pdf"), 5, 7)

print(
  ggplot(Plot, aes(x = specificities, y = sensitivities)) +
    geom_path(aes(group = Resample, colour = Type), alpha = 0.1) +
    geom_smooth(aes(colour = Type), size = 1.3) +
    ylim(0, 1) +
    geom_abline(aes(slope = 1, intercept = 1)) +
    scale_x_reverse(limit = c(1, 0)) +
    theme_classic() +
    theme(legend.position = "bottom") +
    annotate(
      "text",
      x = 0.25,
      y = 0.25,
      label = paste0("LASSO: ",
        round(PlotLASSO[[2]]$cvAUC, 3),
        " [",
        round(PlotLASSO[[2]]$ci[1], 3),
        ",",
        round(PlotLASSO[[2]]$ci[2], 3),
        "]"
      )
    ) +
    annotate(
      "text",
      x = 0.25,
      y = 0.20,
      label = paste0(
        "XGBOOST: ",
        round(PlotXGBoost[[2]]$cvAUC, 3),
        " [",
        round(PlotXGBoost[[2]]$ci[1], 3),
        ",",
        round(PlotXGBoost[[2]]$ci[2], 3),
        "]"
      )
    ) +  labs(x = "Specificity", y = "Sensitivity")
)

dev.off()

##########################################################
######################### Variable Importance ############
##########################################################

#https://cran.r-project.org/web/packages/vip/vignettes/vip-introduction.pdf

a <- race_results$result

# extract from CV - boost

XgboostImp <-
  lapply(1:length(a[[2]][[".extracts"]]), function(x) {
    a[[2]][[".extracts"]][[x]][[".extracts"]][[1]]
  }) %>%
  bind_rows() %>% #last 1 should be a 3? -- what to do with that
  group_by(Variable) %>%
  dplyr::summarize(Mean = mean(Importance),
                   N = n(),
                   Sd = sd(Importance)) %>%
  filter(N == CV) %>%
  arrange(-Mean) %>%
  head(5)


XgboostImp$Variable <- XgboostImp$Variable  %>%
  str_replace_all(
    c(
      "Grade." = "Grade",
      "3." = " 3",
      "Grade4." = "Grade 4+",
      "ASA_grade" = "ASA",
      "BUPA_operation_grade" = "BUPA Grade ",
      "_" = ": "
    )
  )


XgboostImpPlot <-  ggplot(XgboostImp, aes(x = Variable, y = Mean)) +
  geom_line() +
  geom_pointrange(aes(ymin = Mean - Sd, ymax = Mean + Sd)) + coord_flip() + theme_bw() + labs(y = "Mean Importance", title = "XGBoost variable importance")

LassoImp <-
  lapply(1:length(a[[1]][[".extracts"]]), function(x) {
    a[[1]][[".extracts"]][[x]][[".extracts"]][[1]]
  }) %>%
  bind_rows() %>% #last 1 should be a 3? -- what to do with that
  filter(Importance != 0) %>%
  group_by(Variable, Sign) %>% #no sign
  dplyr::summarize(Mean = mean(Importance),
                   N = n(),
                   Sd = sd(Importance)) %>%
  filter(N == CV) %>%
  arrange(-Mean) %>%
  head(5)


LassoImp$Variable <- LassoImp$Variable  %>%
  str_replace_all(
    c(
      "Country_income_SARSCoV2_risk_Low.middle.income_High" = "Low/middle income and high SARSCoV2 risk",
      "Country_income_SARSCoV2_risk_Upper.middle.income_Low" = "Upper/middle income and low SARSCoV2 risk",
      "Age_X90..years" = " Age: 90+years",
      "ASA_grade_Grade.4." = "ASA: Grade 4+",
      "Compartment_Thoracic" = "Compartment: Thoracic"
    )
  ) # will have to adapt to what appers in final mansucript


LassoImpPlot <-
  ggplot(LassoImp, aes(x = Variable, y = Mean)) + #no color
  geom_line() +
  geom_pointrange(aes(ymin = Mean - Sd, ymax = Mean + Sd)) + coord_flip() + theme_bw() + labs(y = "Mean Importance", title = "LASSO variable importance")


pdf(paste0(IncludeFigHere, "/VarImportance.pdf"), 11, 6)
XgboostImpPlot + LassoImpPlot
dev.off()


##################################################################
######################### Model Choice (LASSO) and Development ###
##################################################################

best_results <-
  race_results %>%
  extract_workflow_set_result("StandardMatrix_lasso") %>%
  select_best(metric = "roc_auc")
best_results

get_model <- function(x) {
  pull_workflow_fit(x) %>% tidy(conf.int = TRUE) #tidy conf int not working?
}
ctrl <- control_resamples(extract = get_model, save_pred = TRUE)

bootstraps <- 2000 #######

set.seed(132)

Boots <-  bootstraps(data_trn, times = bootstraps)

# selected model
boosting_test_results <- race_results %>%
  extract_workflow("StandardMatrix_lasso") %>% 
  finalize_workflow(best_results) %>%
  fit_resamples(resamples = Boots,
                control = ctrl,
                metrics = multi_metric2)

collect_metrics(boosting_test_results)

all_coef <-
  lapply(c(1:bootstraps), function(i) {
    boosting_test_results$.extracts[[i]][[1]]
  }) %>%
  bind_rows(.id = "Boots") %>%
  group_by(term) %>%
  # summarise_at(vars(estimate), list(Min = min, Mean = mean, Max = max, Sd = sd, n = n())) %>%
  summarise(
    Min = min(estimate),
    Mean = mean(estimate),
    Max = max(estimate),
    sd = sd(estimate),
    n = n()
  ) %>%
  mutate(
    se.mpg = sd / sqrt(n),
    conf.low = Mean - qt(1 - (0.05 / 2), n - 1) * se.mpg,
    conf.high = Mean + qt(1 - (0.05 / 2), n - 1) * se.mpg
  ) %>%
  mutate(Category = str_extract(term, "[a-z]{1,4}"))

all_coef$term <- all_coef$term  %>%
  str_replace_all(
    c(
      "Head...Neck" = "Head_Neck",
      "Grade." = "Grade",
      "_X" = "_",
      "\\.y" = "y",
      "3." = "3",
      "\\." = "-",
      "90-y" = "90+y",
      "Grade4-" = "Grade4+"
    )
  )

save(all_coef, file = paste0(IncludeFigHere, "/all_coef.RData"))


##################################################################
######################### Check results ##########################
##################################################################

#IncludeFigHere <- "ExtractInfoNoDeath_20220405_1047"
#load(paste0("/rds/projects/g/gkoutosg-variant-prediction/VictorSaisa/FinalCode/", IncludeFigHere, "/all_coef.RData"))
#load(paste0("/rds/projects/g/gkoutosg-variant-prediction/VictorSaisa/FinalCode/", IncludeFigHere, "/Modeling2data3modsThree.RData"))

##################################################################
##################################################################
##################################################################


all_coef2 <- all_coef %>%
  select(term, Mean, conf.low , conf.high)

all_coef2$term[1] <- "Intercept"

all_coef2 <- all_coef2 %>%
  arrange(term, desc(term))


###############################################################
############################## Calculator  ###################
############################################################### 

print("Calculator")

Calculator <- all_coef2 %>%
  mutate(Other = round(Mean * 3, 0)) %>%
  mutate(Odds = exp(Mean)) %>%
  mutate(Interaction = ifelse(str_detect(term, "_x_"), "Yes", "No")) %>%
  mutate(FirstLets = case_when(Interaction == "No" ~ substr(term, 1, 4), TRUE ~ "Interaction")) %>%
  mutate_if(is.numeric, round, 3) 


CalculatorTable1 <- Calculator %>% 
  mutate_if(is.numeric,as.character) %>%
  mutate(Mean = paste0(Mean," (", conf.low, "-", conf.high, ")")) %>%
  select(-c(conf.low, conf.high))

CalculatorTable1$term <- CalculatorTable1$term %>%
  str_replace_all(
    c(
      "years" = " years", 
      "Grade" = "Grade ", 
      "RCRI_3" = "RCRI_3+", 
      "Not-performed" = "Not performed", 
      "Upper-middle-" = "Upper middle ", 
      "High-" = "High ", 
      "Low-middle-" = "Low middle ", 
      "Hospital-with-no-defined-pathway" ="Hospital with no defined pathway", 
      "Head_Neck" = "Head & Neck"
    ))

names(t)[c(1,2)] <- c("Vars", "Levels")

t <- t %>%
  replace_with_na(replace = list(Vars = ""))

tMerge <- t %>%
  fill(Vars) %>%
  mutate(term = paste0(Vars, "_", Levels)) %>%
  full_join(CalculatorTable1) %>%
  rename(Score = Other) %>%
  select(-c(term, Interaction:FirstLets)) %>%
   replace(is.na(.), "")

write.table(tMerge,
          paste0(IncludeFigHere, "/CalculatorTableManuscript.tsv"),
          row.names = FALSE)

CalculatorDelete <- Calculator %>%
  filter(Other != 0) #take out smoking, test and pathway

CalculatorDelete[CalculatorDelete$term == "Intercept", ]$Other <- 0

CalculatorCoeffs <- CalculatorDelete %>%
  select(term, Other, Mean, conf.low , conf.high) %>%
  arrange(term, desc(term))


###############################################################
############################## Functions Model Fitting  ##########
###############################################################

invlogit <- function(x) {
  1 / (1 + exp(-x))
}

Funs <- function(all_coef2, dataBaked, data, Calculator) {  #KEY!!!
  
  TrainingInteractions <- dataBaked %>%
    select(-Label)
  
  TrainingInteractionsNoLabel <-
    data.frame(Intercept = 1, TrainingInteractions) %>%
    select(sort(current_vars()))
  
  ff <- as.matrix(data.frame(TrainingInteractionsNoLabel))
  
  linpred <- ff  %*% all_coef2$Mean #intercept has a value!!
  
  .pred_Positive <- invlogit(linpred)
  .pred_Negative <- 1 - invlogit(linpred)
  
  train_probsDerv <-
    data.frame(data, linpred, .pred_Negative, .pred_Positive)
  
  
  if (Calculator == 1) {
    TrainingInteractions <- dataBaked %>%
      select(-Label)
    
    TrainingInteractionsNoLabel <-
      data.frame(Intercept = 1, TrainingInteractions) %>%
      select(sort(current_vars()))
    
    ff <- as.matrix(data.frame(TrainingInteractionsNoLabel))
    
    linpred2 <- ff  %*% all_coef2$Other #intercept 0!!
    
    train_probsDerv <- data.frame(train_probsDerv, linpred2)
    
    
  }
  
  return(train_probsDerv)
  
  
}


###############################################################
############################## Model Fitting Datasets  #########
###############################################################

#Derivation whole set

WholeDev <- bake(PreparedPreProc, new_data = WholeDevelopment)

names(WholeDev) <- names(WholeDev)  %>%
  str_replace_all(
    c(
      "Head...Neck" = "Head_Neck",
      "Grade." = "Grade",
      "_X" = "_",
      "\\.y" = "y",
      "3." = "3",
      "\\." = "-",
      "90-y" = "90+y",
      "Grade4-" = "Grade4+"
    )
  )

WholeDevf <- WholeDev %>%
  select(setdiff(CalculatorDelete$term, "Intercept"), "Label")

train_probsDervWholeCalc <-
  Funs(CalculatorCoeffs, WholeDevf, WholeDev, 1)


#Derivation training set

names(StandardMatrix)  <- names(StandardMatrix)   %>%
  str_replace_all(
    c(
      "Head...Neck" = "Head_Neck",
      "Grade." = "Grade",
      "_X" = "_",
      "\\.y" = "y",
      "3." = "3",
      "\\." = "-",
      "90-y" = "90+y",
      "Grade4-" = "Grade4+"
    )
  )

trainf <- StandardMatrix %>%
  select(setdiff(CalculatorDelete$term, "Intercept"), "Label")

train_probsDervCalc <- Funs(CalculatorCoeffs, trainf, data_trn, 1)

#-------#Validation Surg Week

uni_testBake <- bake(PreparedPreProc, new_data = data_tst)

names(uni_testBake) <- names(uni_testBake)  %>%
  str_replace_all(
    c(
      "Head...Neck" = "Head_Neck",
      "Grade." = "Grade",
      "_X" = "_",
      "\\.y" = "y",
      "3." = "3",
      "\\." = "-",
      "90-y" = "90+y",
      "Grade4-" = "Grade4+"
    )
  )

uni_testBakef <- uni_testBake %>%
  select(setdiff(CalculatorDelete$term, "Intercept"), "Label")

train_probsValCalc <-
  Funs(CalculatorCoeffs, uni_testBakef, data_tst, 1)

#----------#Recon

#recon <- read_csv( "/rds/projects/g/gkoutosg-variant-prediction/Laura/Magda/JamesCodeFuns/CombinedThreeFinalDatasetReconInElectiveLabelChangeNewVars.csv") %>% #withoutNewVars

recon <-  read.csv(paste0(IncludeFigHere,"/CombinedThreeFinalDatasetReconInElectiveLabelChangeNewVars.csv")) %>%
filter(Type == "recon") %>%
  mutate_if(is.character, as.factor) %>%
  mutate(
    Anaesthetic_type = recode(
      Anaesthetic_type,
      "General anaesthetic" = "General",
      "Regional/local anaesthetic" = "Regional/local"
    )
  ) %>%
  select(-X)  %>% #X1 before
  #filter(!(Death == "Yes" & Label == "Negative")) %>% #15 deleted
  #filter(!(is.na(Death) & Label == "Negative")) %>% #2 deleted
  select(-Death)


reconBake <- bake(PreparedPreProc, new_data = recon)

names(reconBake) <- names(reconBake)  %>%
  str_replace_all(
    c(
      "Head...Neck" = "Head_Neck",
      "Grade." = "Grade",
      "_X" = "_",
      "\\.y" = "y",
      "3." = "3",
      "\\." = "-",
      "90-y" = "90+y",
      "Grade4-" = "Grade4+"
    )
  )

reconBakef <- reconBake %>%
  select(setdiff(CalculatorDelete$term, "Intercept"), "Label")

train_probsCalc <- Funs(CalculatorCoeffs, reconBakef, recon, 1)

#-----------# Derivation just testing part

uni_TestBake <- bake(PreparedPreProc, new_data = Test)

names(uni_TestBake) <- names(uni_TestBake)  %>%
  str_replace_all(
    c(
      "Head...Neck" = "Head_Neck",
      "Grade." = "Grade",
      "_X" = "_",
      "\\.y" = "y",
      "3." = "3",
      "\\." = "-",
      "90-y" = "90+y",
      "Grade4-" = "Grade4+"
    )
  )

uni_TestBakef <- uni_TestBake %>%
  select(setdiff(CalculatorDelete$term, "Intercept"), "Label")

train_probsTestCalc <-
  Funs(CalculatorCoeffs, uni_TestBakef, Test, 1)


######## Legend from above. 

#data_trn ->train_probsDervCalc:  Derivation training 
#data_tst ->train_probsValCalc : Validation (SurgWeek)
#recon -> train_probsCalc: Recon
#Test -> train_probsTestCalc: Derivation Test
#WholeDevelopment -> train_probsDervWholeCalc: WholeDerivation


###############################################################
############################## Discrimination: AUC  ###########
###############################################################

#True calculator AUC - same as previous, just multiplied and added (ranking stays the same)

ROCWithCI <- function(train_probsValCalc) {
  obj <-
    roc(
      train_probsValCalc$Label,
      train_probsValCalc$linpred2,
      ci = TRUE,
      plot = FALSE
    )
  ciobj <-
    ci.se(obj, specificities = seq(0, 1, l = 25), method = "bootstrap")
  dat.ci <- data.frame(x = as.numeric(rownames(ciobj)),
                       lower = ciobj[, 1],
                       upper = ciobj[, 3])
  
  basicplot <-
    ggplot(train_probsValCalc, aes(d = Label, m = linpred2)) +
    geom_roc(n.cuts = 50, labels = FALSE) +
    style_roc(
      ylab = "Sensitivity",
      xlab = "1 - Specificity",
      minor.breaks = c(seq(0, 0.1, by = 0.01), seq(0.9, 1, by = 0.01)),
      guide = TRUE
    )
  
  dat.ci2 <- dat.ci %>%
    mutate(x = 1 - x)
  
  
  Vald <- basicplot +
    #ggtitle(paste0(title,"-", ChosenRun)) +
    annotate(
      "text",
      x = .7,
      y = .25,
      label = paste("AUC =", round(calc_auc(basicplot)$AUC, 4)),
      fontface = "bold"
    ) +
    annotate(
      "text",
      x = .7,
      y = .15,
      label = paste(gsub(
        "\\s*\\([^\\)]+\\)", "", capture.output(obj$ci)
      )),
      fontface = "bold"
    ) +
    geom_abline(
      slope = 1,
      intercept = 0,
      color = "grey",
      linetype = "dashed"
    ) +
    geom_ribbon(
      data = dat.ci2,
      aes(x = x, ymin = lower, ymax = upper),
      inherit.aes = FALSE,
      fill = "steelblue",
      alpha = 0.2
    )
  
  return(Vald)
}


Validation <- ROCWithCI(train_probsValCalc)
Derivation <- ROCWithCI(train_probsDervCalc)
Rec <- ROCWithCI(train_probsCalc)
DervTest <- ROCWithCI(train_probsTestCalc)
WholeDerivation <- ROCWithCI(train_probsDervWholeCalc)

pdf(paste0(IncludeFigHere, "/AUCTogether.pdf"), 15, 5)
print(Validation + Derivation + Rec + DervTest + WholeDerivation )
dev.off()


###############################################################
############################## Calibration  ##################
###############################################################


#Calibration Plot
Plot1Cal1 <- function(DerivationWPred, Title) {
  pred_var <- "m1_pred"
  model_name <- Title
  
  # The calibration plot
  g1 <-
    mutate(DerivationWPred, bin = ntile(get(pred_var), 10)) %>% #50
    # Bin prediction into 10ths
    group_by(bin) %>%
    mutate(
      n = n(),
      # Get ests and CIs
      bin_pred = mean(get(pred_var)),
      bin_prob = mean(as.numeric(as.factor(Label)) - 1),
      se = sqrt((bin_prob * (1 - bin_prob)) / n),
      ul = bin_prob + 1.96 * se,
      ll = bin_prob - 1.96 * se
    ) %>%
    ungroup() %>%
    select(bin, bin_pred, bin_prob, ll, ul) %>%
    unique()
  
  return(g1)
}


set.seed(132)

# Calibration ------------- transform calculator scores to probabilities

CalibrationEwen <- function(DervCal, ReconCal) {
  
  m1 <- glm(Label ~ linpred2, data = DervCal, family = binomial)
  
  # Predicted probabilities in validation data
  ReconCal$m1_pred <-
    predict(m1, newdata = ReconCal, type = "response")
  
  ## Calibration in the large ---------------------------------------------------------
  linear_pred_m1 = predict(m1, newdata = ReconCal)
  fit_linear_pred_m1 = glm(Label ~ 1,
                           offset = linear_pred_m1,
                           data = ReconCal,
                           family = "binomial")
  
  CITL <- fit_linear_pred_m1[["coefficients"]]
  
  print(CITL)
  
  ## Calibration slope -----------------------------------------------------------------
  glm(Label ~ linear_pred_m1, data = ReconCal, family = "binomial")
  
  ## Brier -----------------------------------------------------------------------------
  Brier <-
    rms::val.prob(ReconCal$m1_pred, as.numeric(ReconCal$Label) - 1,
                  pl = FALSE) %>% round(3)
  
  print(Brier[11:13])
  
  return(ReconCal)
  
}


ReconCal <- CalibrationEwen(train_probsDervCalc, train_probsCalc)
ValidationCal <- CalibrationEwen(train_probsDervCalc, train_probsValCalc)
TestCal <- CalibrationEwen(train_probsDervCalc, train_probsTestCalc)
WholeDervCal <- CalibrationEwen(train_probsDervCalc, train_probsDervWholeCalc)

g1a <- list(
    Recon = Plot1Cal1(ReconCal, "Calibration in RECON Set"),
    Val = Plot1Cal1(ValidationCal, "Calibration in Validation Set") ,
    Test =  Plot1Cal1(TestCal, "Calibration in Test Set")
  ) %>%
  bind_rows(.id = "Datasets") #Derv =  Plot1Cal1( DervCal, "Calibration in Derivation Set")

g1a$Datasets <- g1a$Datasets %>%
  str_replace_all(
    c(
      #  "Derv" = "DerivationTrain:SurgWeek" ,
      "Val" = "Validation 1:CovidSurg-Cancer",
      "Recon" = "Validation 2:RECON",
      "Test" = "DerivationTest:SurgWeek"
    )
  )


#Calibration Plot
g1 <- ggplot(g1a,
         aes(
           x = bin_pred,
           y = bin_prob,
           ymin = ll,
           ymax = ul,
           colour = Datasets
         )) +
  geom_pointrange(size = 0.3) +
  scale_y_continuous(limits = c(0, 0.18), breaks = seq(0, 1, by = 0.1)) +
  coord_cartesian(expand = FALSE, ylim = c(0, 0.18)) +
  scale_x_continuous(limits = c(0, 0.18), breaks = seq(0, 1, by = 0.1)) +
  geom_abline(lty = 2) + # 45 degree line indicating perfect calibration
  geom_line() +
  xlab("Predicted probability") +
  ylab("Observed probability") +
  theme_classic() + theme(legend.position = "bottom")


#Histogram plot

# The distribution plot
gVal <- ggplot(ValidationCal, aes(x = m1_pred, alpha = Label)) +
  geom_histogram(fill = "dark green" , bins = 10) +
  scale_alpha_manual(values = c("Negative" = 0.7, "Positive" = 1)) +
  #scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1)) +
  xlab("Predicted Probability") +
  ylab("") +
  scale_alpha_discrete(name = "PPC") +
  theme_minimal() +
  theme(panel.grid.minor = element_blank())

gRecon <- ggplot(ReconCal, aes(x = m1_pred, alpha = Label)) +
  geom_histogram(fill = "blue" , bins = 10) +
  scale_alpha_manual(values = c("Negative" = 0.7, "Positive" = 1)) +
  #scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1)) +
  xlab("Predicted Probability") +
  ylab("") +
  scale_alpha_discrete(name = "PPC") +
  theme_minimal() +
  #scale_y_continuous(breaks = c(0, 100)) +
  theme(panel.grid.minor = element_blank())

gTest <- ggplot(TestCal, aes(x = m1_pred, alpha = Label)) +
  geom_histogram(fill = "dark orange", bins = 10) +
  #scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1)) +
  scale_alpha_manual(values = c("Negative" = 0.2, "Positive" = 1)) +
  xlab("Predicted Probability") +
  ylab("") +
  scale_alpha_discrete(name = "PPC") +
  theme_minimal() +
  #scale_y_continuous(breaks = c(0, 100)) +
  theme(panel.grid.minor = element_blank())

# Combine them
g <- arrangeGrob(
    g1,
    gTest,
    gVal,
    gRecon,
    respect = TRUE,
    heights = c(1, 0.25, 0.25, 0.25),
    ncol = 1
  )
grid.newpage()
grid.draw(g)

pdf(paste0(IncludeFigHere, "/TogetherCalibrationFINAL4.pdf"), 6.5, 7)
grid.draw(g)
dev.off()

save(
  train_probsDervCalc,
  train_probsCalc,
  train_probsValCalc,
  train_probsDervWholeCalc,
  data_trn,
  data_tst,
  Test,
  recon,
  train_probsTestCalc,
  WholeDevelopment,
  dataFull,
  file = paste0(IncludeFigHere, "/dataCal.RData")
)

########################### Transportability Studies #############################

#https://www.sciencedirect.com/science/article/pii/S0895435614002753

########################### Thresholds #############################

###### extract table with thresholds

MortalityTable <- function(Probs, train_probsCalc) {
  #cut(train_probsCalc$linpred2, round(c(1,5,10,15), 2), include.lowest = TRUE)
  
  CaseFatalValMort <- train_probsCalc %>%
    mutate(z = cut(linpred2, round(Probs, 2), include.lowest = TRUE)) %>%
    group_by(z) %>%
    mutate(CaseFatality = length(Label[Label == "Positive"]) / (length(Label[Label == "Negative"]) + length(Label[Label == "Positive"]))) %>%
    select(z, CaseFatality, linpred2, Label) %>%
    mutate(Risk = paste0(z, " & ", round(CaseFatality, 3)))
  
  Probs <- Probs[-c(1, length(Probs))]
  
  MortalityTab <- list()
  
  Fin <- CaseFatalValMort
  
  for (i in seq_along(Probs)) {
   
     print(i)
    
    Cutoff <- Probs[i]
    
    Fin <- Fin %>%
      mutate(Var2 = case_when(linpred2 <= Cutoff ~ "NegativePred", TRUE ~ "PositivePred"))
    
    MortalityTab[[i]] <-
      data.frame(table(Fin$z, Fin$Label, Fin$Var2)) %>%
      group_by(Var2, Var3) %>%
      mutate(Sum = sum(Freq)) %>%
      select(-c(Var1, Freq)) %>%
      unique() %>%
      mutate(
        Variables = case_when(
          Var2 == "Negative" & Var3 == "NegativePred" ~ "TN",
          Var2 == "Negative" & Var3 == "PositivePred" ~ "FP",
          Var2 == "Positive" & Var3 == "NegativePred" ~ "FN",
          TRUE ~ "TP"
        )
      ) %>%
      ungroup() %>%
      select(c(Variables, Sum)) %>%
      pivot_wider(names_from = Variables, values_from = Sum) %>%
      mutate(
        Cutoff = paste0("<=", Cutoff),
        NumPats = data.frame(table(Fin$Var2)) %>% filter(Var1 == "NegativePred") %>% select(Freq) %>% as.numeric(),
        Sensitivity = TP / (TP + FN) * 100,
        Specificity = TN / (TN + FP) * 100,
        PPV = TP / (TP + FP) * 100,
        NPV = TN / (TN + FN) * 100,
        Mortality = FN / NumPats * 100
      )
  }
  
  MortalityTab5 <- bind_rows(MortalityTab)
  
  MortalityTab <- list()
  
  for (i in seq_along(Probs)) {
    print(i)
    
    Cutoff <- Probs[i]
    
    Fin <- Fin %>%
      mutate(Var2 = case_when(linpred2 >= Cutoff ~ "PositivePred", TRUE ~ "NegativePred"))
    
    MortalityTab[[i]] <-
      data.frame(table(Fin$z, Fin$Label, Fin$Var2)) %>%
      group_by(Var2, Var3) %>%
      mutate(Sum = sum(Freq)) %>%
      select(-c(Var1, Freq)) %>%
      unique() %>%
      mutate(
        Variables = case_when(
          Var2 == "Negative" & Var3 == "NegativePred" ~ "TN",
          Var2 == "Negative" & Var3 == "PositivePred" ~ "FP",
          Var2 == "Positive" & Var3 == "NegativePred" ~ "FN",
          TRUE ~ "TP"
        )
      ) %>%
      ungroup() %>%
      select(c(Variables, Sum)) %>%
      pivot_wider(names_from = Variables, values_from = Sum) %>%
      mutate(
        Cutoff = paste0(">=", Cutoff),
        NumPats = data.frame(table(Fin$Var2)) %>% filter(Var1 == "PositivePred") %>% select(Freq) %>% as.numeric(),
        Sensitivity = TP / (TP + FN) * 100,
        Specificity = TN / (TN + FP) * 100,
        PPV = TP / (TP + FP) * 100,
        NPV = TN / (TN + FN) * 100,
        Mortality = TP / NumPats * 100
      )
  }
  
  MortalityTab6 <- bind_rows(MortalityTab)
  
  qq <- rbind(MortalityTab5, MortalityTab6) %>%
    #filter(Cutoff %in% c("<=0.2", "<=0.4", "<=0.6", ">=0.6", ">=0.8")) %>%
    mutate(Percentage = round(NumPats / dim(train_probsCalc)[1] * 100, 2)) %>% #before dataA1?
    mutate(NumPer = paste0(NumPats, " (", Percentage, ")")) %>%
    select(-c(Percentage, NumPats)) %>%
    relocate(Cutoff,
             NumPer,
             TP,
             TN,
             FP,
             FN,
             Sensitivity,
             Specificity,
             PPV,
             NPV,
             Mortality) %>%
    dplyr::rename(
      'Cutoff Risk' = Cutoff,
      'No of patients (%)' = NumPer,
      'Sensitivity (%)' = Sensitivity,
      'Specificity (%)' = Specificity,
      'PPV (%)' = PPV,
      'NPV (%)' = NPV,
      'Pulmonary Complications (%)' = Mortality
    ) %>%
    mutate_if(is.numeric, round, 2)
  
  
  
}


ProbsRecon <- c(1, 5, 7, 9, 11, 13, 15) #thresholds
qqRecon <- MortalityTable(ProbsRecon, train_probsCalc)
ProbsVal <- c(1, 5, 7, 9, 11, 13, 15, 17, 18)
qqVal <- MortalityTable(ProbsVal, train_probsValCalc)
ProbsDerv <- c(1, 5, 7, 9, 11, 13, 15, 17, 18, 19, 20)
qqDerv <- MortalityTable(ProbsDerv, train_probsDervCalc) #This is just the trainig part not whole derivation!
ProbsDervTest <- c(1, 5, 7, 9, 11, 13, 15, 17, 19)
qqDervTest <- MortalityTable(ProbsDervTest, train_probsTestCalc) #This is just the trainig part not whole derivation!
ProbsWholeDervTest <- c(1, 5, 7, 9, 11, 13, 15,17,18,19,20)
qqDervWholeDervTest <- MortalityTable(ProbsWholeDervTest, train_probsDervWholeCalc) #This is just the trainig part not whole derivation!


myt <- ttheme_default(
  core = list(bg_params = list(fill = c(
    "transparent", "light blue"
  ))),
  colhead = list(
    fg_params = list(col = "white"),
    bg_params = list(fill = "black")
  )
)



pdf(paste0(IncludeFigHere, "/RiskTable2.pdf"), 12, 7)

print(grid.table(qqDerv, theme = myt, rows = NULL))
grid.newpage()
print(grid.table(qqVal, theme = myt, rows = NULL))
grid.newpage()
print(grid.table(qqRecon, theme = myt, rows = NULL))
grid.newpage()
print(grid.table(qqDervTest, theme = myt, rows = NULL))
grid.newpage()
print(grid.table(qqDervWholeDervTest, theme = myt, rows = NULL))

dev.off()

write.csv(qqDerv,
          paste0(IncludeFigHere, "/RiskTableDerv2.csv"),
          row.names = FALSE)
write.csv(qqVal,
          paste0(IncludeFigHere, "/RiskTableVal2.csv"),
          row.names = FALSE)
write.csv(qqRecon,
          paste0(IncludeFigHere, "/RiskTableRecon2.csv"),
          row.names = FALSE)
write.csv(qqDervTest,
          paste0(IncludeFigHere, "/RiskTableDervTest2.csv"),
          row.names = FALSE)
write.csv(qqDervWholeDervTest,
          paste0(IncludeFigHere, "/RiskTableDervWholeTest2.csv"),
          row.names = FALSE)



###############################################################
############################## Subgroups  ###################
###############################################################

#  minor v major, elective v emergency, high risk (ASA 3-5) v low risk (ASA 1-2), income group (high versus low-middle income (will have to combine UMIC, LMIC, LIC))


Subgroup <- list(
    train_probsTestCalc %>% add_column(Type = "Test"),
    train_probsCalc,
    train_probsValCalc %>% add_column(Type = "Validation")
  ) %>%
  bind_rows() %>%
  mutate(Risk = as.factor(ifelse(
    ASA_grade %in% c("Grade 3", "Grade 4+"), #was coded wrongly before with 3,4,5?
    "High Risk",
    "Low Risk"
  )),
  Country_income = as.factor(
    ifelse(
      Country_income_SARSCoV2_risk %in% c(
        "High income_Low",
        "High income_High",
        "Upper middle income_High",
        "Upper middle income_Low"
      ),
      "High_income",
      "Low_income"
    )
  ))

Subgroup$Type <- Subgroup$Type %>%
  str_replace_all(
    c(
      "Validation" = "Validation 1:CovidSurg-Cancer",
      "recon" = "Validation 2:RECON",
      "Test" = "DerivationTest:SurgWeek"
    )
  )


SubgroupsNames <- c("BUPA_operation_grade", "Risk", "Country_income")

annotate <- list()
Histogram <- list()

for (j in SubgroupsNames) {
  print(j)
  
  Histogram[[j]] <-
    ggplot(Subgroup, aes(x = linpred2, fill = !!sym(j))) +
    geom_bar(alpha = 6 / 10) + facet_wrap(. ~ Type) + theme_bw() + labs(x = "Score") +
    theme(legend.position = "bottom") + labs(y = "Number of patients")
  
  
  for (i in levels(Subgroup[[j]])) {
    print(paste0(i, "----", j))
    
    basicplot <-
      ggplot(Subgroup[Subgroup[j] == i, ], aes(d = Label, m = linpred2, color = Type)) +
      geom_roc(n.cuts = 50, labels = FALSE) +
      style_roc(
        ylab = "Sensitivity",
        xlab = "1 - Specificity",
        minor.breaks = c(seq(0, 0.1, by = 0.01), seq(0.9, 1, by = 0.01)),
        guide = TRUE
      )  + scale_colour_manual(
        values = c(
          "DerivationTest:SurgWeek" = "orange",
          "Validation 2:RECON" = "blue",
          "Validation 1:CovidSurg-Cancer" = "darkgreen"
        )
      )
    
    colors <- c("orange", "blue", "darkgreen")
    names(colors) <-
      c(
        "DerivationTest:SurgWeek",
        "Validation 2:RECON",
        "Validation 1:CovidSurg-Cancer"
      )
    
    Numbers <- c(3, 1, 2)
    names(Numbers) <-
      c(
        "DerivationTest:SurgWeek",
        "Validation 2:RECON",
        "Validation 1:CovidSurg-Cancer"
      )
    
    for (x in unique(Subgroup[Subgroup[j] == i, ]$Type)) {
     
       Data <- Subgroup[Subgroup[j] == i, ] %>%
        filter(Type == x)
      
      obj <- roc(Data$Label,
                 Data$linpred2,
                 ci = TRUE,
                 plot = FALSE)
      ciobj <-
        ci.se(obj,
              specificities = seq(0, 1, l = 25),
              method = "bootstrap")
      dat.ci <- data.frame(x = as.numeric(rownames(ciobj)),
                           lower = ciobj[, 1],
                           upper = ciobj[, 3])
      
      dat.ci2 <- dat.ci %>%
        mutate(x = 1 - x)
      
      
      basicplot <- basicplot +
        annotate(
          "text",
          x = .7 ,
          y = .10 + Numbers[[x]] / 6,
          label = paste("AUC =", round(obj[["ci"]][2], 4)),
          fontface = "bold",
          color = colors[[x]]
        ) +
        annotate(
          "text",
          x = .7,
          y = .02 + Numbers[[x]] / 6,
          label = paste(gsub(
            "\\s*\\([^\\)]+\\)", "", capture.output(obj$ci)
          )),
          fontface = "bold",
          color = colors[[x]]
        ) +
        geom_abline(
          slope = 1,
          intercept = 0,
          color = "grey",
          linetype = "dashed"
        ) +
        geom_ribbon(
          data = dat.ci2,
          aes(
            x = x,
            ymin = lower,
            ymax = upper
          ),
          inherit.aes = FALSE,
          fill = colors[[x]],
          alpha = 0.2
        )
      ######################
      
      
    }
    
    annotate[[i]] <- basicplot
    
  }
  
  
}


pdf(paste0(IncludeFigHere, "/Subgroups_BUPA_operationFINAL.pdf"),
    7,
    10)
print((annotate[[1]] / annotate[[2]]) +  Histogram[[1]] +  plot_layout(guides = 'collect')  &
        theme(legend.position = 'bottom')
)
dev.off()

pdf(paste0(IncludeFigHere, "/Subgroups_RiskFINAL.pdf"), 7, 10)
print((annotate[[3]] / annotate[[4]]) +  Histogram[[2]]  +  plot_layout(guides = 'collect')  &
        theme(legend.position = 'bottom')
)
dev.off()

pdf(paste0(IncludeFigHere, "/Subgroups_Country_incomeFINAL.pdf"), 7,10)

print((annotate[[5]] / annotate[[6]]) +  Histogram[[3]]  +  plot_layout(guides = 'collect') &
        theme(legend.position = 'bottom')
)

dev.off()

########################################################################
####################################### Shiny App ##################
########################################################################

source("FinalExtraApps.R") #Final table with risk groups and Extra results for shiny app


#Perfectt to explore everything!
#obj <- roc(train_probsDerv$Label, train_probsDerv$.pred_Positive, ci = TRUE) %>%
#coords(ret = "all")

########################################################################
####################################### Sample Size ##################
########################################################################

#first sample size
pmsampsize(
  type = "b",
  #rsquared = NA,
  parameters = 37 ,
  # shrinkage = 0.9,
  prevalence = 0.03,
  cstatistic = 0.7,
  seed = 123456
)

#Based on developed model

pmsampsize(
  type = "b",
  #rsquared = NA,
  parameters = 22 ,
  # shrinkage = 0.9,
  prevalence = 0.03,
  cstatistic = 0.77,
  seed = 123456
)

##############

print(sessionInfo())
writeLines(capture.output(sessionInfo()), "/rds/projects/g/gkoutosg-variant-prediction/VictorSaisa/FinalCode/SessionInfoForVictor.txt")

#############


#IncludeFigHere <- "FinalNoDeathCorrectRiskCI_20220405_2055"

ValidationTable <- read.csv(paste0(IncludeFigHere,"/RiskSmallTableValidation.csv"))
DerivationTable <- read.csv(paste0(IncludeFigHere,"/RiskSmallTableDerivationWhole.csv"))
ReconTable <- read.csv(paste0(IncludeFigHere,"/RiskSmallTableRecon.csv"))

Tables <- list(DerivationTable,ValidationTable, ReconTable) 
names(Tables) <- c("DerivationTable","ValidationTable", "ReconTable")

Tables2 <- Tables %>%
  bind_rows(.id = "Names") %>%
  select(-X) 

write.csv(Tables2,paste0(IncludeFigHere,"/CondensedSmallTable.csv"))

#- run calculator code and recepies

