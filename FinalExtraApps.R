

############################################################
#################### .RDS App Files #######################
############################################################


myt <- ttheme_default(
  core = list(
    bg_params=list(fill=c("transparent", "light blue"))),
  colhead = list(fg_params=list(col="white"),
                 bg_params=list(fill="black"))
)

Calculator2 <- CalculatorDelete %>%
  select(term,Mean,Odds,Other, conf.low ,conf.high) %>%
  arrange(-Other) %>%
  rename(   "Calculator Number" = "Other", 
            "Variable" = "term" )

saveRDS(Calculator2, paste0(IncludeFigHere,"/Calculator2.rds"))


myt <- ttheme_default(
  # Use hjust and x to left justify the text
  # Alternate the row fill colours
  core = list(
    bg_params=list(fill=c("transparent", "light blue"))),
  
  # Change column header to white text and red background
  colhead = list(fg_params=list(col="white"),
                 bg_params=list(fill="black"))
)

 
grid.newpage()

Calculator2Round <- Calculator2 %>%
  mutate_if(is.numeric, round, 2) %>%
  select(Variable, Mean, conf.low, conf.high, Odds, "Calculator Number") %>%
  arrange(-Mean)

names(Calculator2Round) <- c("variable name","mean betas",  "conf.low", "conf.high", "odds", "score")


pdf(paste0(IncludeFigHere, "/ImprovedVariablesDecimals.pdf"), 12,8)
grid.table(Calculator2Round,theme=myt, rows=NULL)
dev.off()


write.csv(Calculator2Round,paste0(IncludeFigHere,"/ImprovedVariablesDecimals.csv"), row.names = FALSE)

dummyNames <- list()
for (i in 1:dim(dataFull)[2]) {
  dummyNames[[i]] <- levels(dataFull[[i]])
}

names(dummyNames) <- names(dataFull)
saveRDS(dummyNames, paste0(IncludeFigHere, "/Vars.rds"))


CalcValues <- recipe(Label ~ ., data = data_trn) %>%
  step_dummy(all_nominal(), -all_outcomes(), one_hot = TRUE)

CalValuesPrep <- CalcValues %>% prep()

CalValuesNames<- names(juice(CalValuesPrep)) %>% 
  str_replace_all(
    c("Head...Neck" = "Head_Neck",
      "Grade." = "Grade",
      "_X" = "_",
      "\\.y" = "y",
      "\\." = "-",
      "90-y" = "90+y",
      "Grade4-" = "Grade4+",
      "3-" = "3"
    ))


#left_join0 <- function(x, y, fill = 0L, ...){
#  z <- full_join(x, y, ...)
#  new_cols <- setdiff(names(z), names(x))
#  z <- replace_na(z, setNames(as.list(rep(fill, length(new_cols))), new_cols))
#  z
#}
#
#aa <- left_join0(data.frame(term = CalValuesNames), Calculator) %>%
#  filter(term != "Label") 
#
#aa$term <- aa$term %>%
#  str_replace_all(
#    c("Grade" = "Grade ", 
#      "years" = " years", 
#      "RCRI_3" = "RCRI_3+", 
#      "Low-middle-income" = "Low middle income",
#      "Upper-middle-income" = "Upper middle income")
#  )
#
#
#saveRDS(aa,paste0(IncludeFigHere, "/CalculatorHotEncoding.rds"))

############################################################
#################### Shiny App Plot #######################
############################################################

ShinyAppPlot <- function(train_probsValCalc, Title){
  
  CountPerVal <- train_probsValCalc %>%
  group_by(linpred2) %>%
  add_tally() %>%
  ungroup() %>%
  group_by(linpred2, Label) %>%
  add_tally() %>%
  select(n,nn, linpred2, Label) %>%
  unique() %>%
  mutate(Div = nn/n*100) %>%
  filter(Label == "Positive") 

p1 <- ggbarplot(CountPerVal, x = "linpred2", y = "n",
                fill = "#0073C299",               # change fill color by cyl
                color = "white",            # Set bar border colors to white
                palette = "jco",            # jco journal color palett. see ?ggpar
                #sort.val = "asc",           # Sort the value in ascending order
                sort.by.groups = TRUE,      # Sort inside each group
                x.text.angle = 90           # Rotate vertically x axis texts
) + ylab("Number of patients") + xlab(" ")
#+ geom_vline(xintercept = 11, linetype="dotted", color = "#86868699", size=1) 


CountPerVal$linpred2 <- as.numeric(CountPerVal$linpred2)

CountPerVal <- CountPerVal %>%
  arrange(linpred2) %>%
  group_by(Label) %>%
  mutate(CumsumPositive = cumsum(nn)) #%>%
  #mutate(CumsumTotal = cumsum(n)) %>%


CountPerVal <- CountPerVal %>%
    mutate(TotalPositive  = CumsumPositive/CountPerVal[[dim(CountPerVal)[1],6]]*100) %>%
  mutate(across(is.numeric, ~ round(., 2))) %>%
  ungroup() %>%  select(-Label) 

#saveRDS(CountPerVal,paste0(IncludeFigHere, "/CountPerVal", Title, ".rds"))

p2 <- ggbarplot(CountPerVal, x = "linpred2", y = "Div",
                fill = "#EFC00099",              
                color = "black",            
                palette = "jco", 
                size = 0.001,
                #sort.val = "asc",
                label = TRUE, label.pos = "out",# Sort the value in ascending order
                lab.col = "#86868699",
                lab.size = 3,
                sort.by.groups = TRUE,
                lab.nb.digits = 1,# Sort inside each group
                x.text.angle = 90           # Rotate vertically x axis texts
)  + xlab("Predicted score") + 
  ylab("Pulmonary complication (%)") #+
#geom_line(data= CountPerVal2, aes(x = (linpred2), y = TotalPositive)) + geom_vline(xintercept = as.character(10), linetype="dotted", 
 #            color = "#86868699", size=1)

figure <- ggarrange(p1, p2 + font("x.text", size = 10),
                    ncol = 1, nrow = 2)

ExtraTab <- train_probsValCalc %>%
        mutate(
          ExtraTab = case_when(
            linpred2 <= 5 ~ "0-5",
            linpred2 <= 10 ~ "6-10",
            linpred2 <= 15 ~ "11-15",
            TRUE ~ "16+") %>% 
            factor(levels = c("0-5", "6-10", "11-15", "16+"))
        ) %>% 
        group_by(ExtraTab) %>% 
        count(Label) %>% 
        mutate(nn = sum(n)) %>% 
        filter(Label == "Positive") %>% 
        drop_na() %>% 
        ungroup() %>% 
        mutate(
          total_n = sum(n),
          total_nn = sum(nn)) %>% 
  bind_rows() %>% 
  group_by(ExtraTab) %>% 
  #summarise_all(mean) %>% 
  group_by(ExtraTab) %>% 
  #summarise_all(ceiling) %>% 
  mutate(
    n_out = paste0(nn, " (", (100*nn/total_nn) %>% round(2), ")"),
    m_out = paste0(n, " (", (100*n/nn) %>% round(2), ")")
  ) %>%
  select(ExtraTab, n_out,m_out) %>%
  rename("Risk group" = ExtraTab) %>%
  rename("Num Patients" = n_out) %>%
  rename("Num Positives" = m_out)


write.csv(ExtraTab, paste0(IncludeFigHere, "/RiskSmallTable", Title, ".csv"))

names(CountPerVal) <- c("Total patients", "Positives", "Score", "Percentage positive", "Cumsum Positive", "Total Positives")

write.csv(CountPerVal, paste0(IncludeFigHere, "/DescriptionShinyApp", Title, ".csv"))


Tables <- grid.arrange(tableGrob(CountPerVal), tableGrob(ExtraTab) ,
          ncol = 1, nrow = 2, 
          heights = c(0.6, 0.3))


pdf(paste0(IncludeFigHere, "/PlotShiny", Title, ".pdf"), 13,10)


print(grid.arrange(
  Tables,
  figure,
  ncol = 2,
  widths = c(1, 0.7),
  clip = FALSE
))


dev.off()



}


ShinyAppPlot( train_probsCalc, "Recon") 
ShinyAppPlot( train_probsValCalc, "Validation") 
ShinyAppPlot(train_probsDervCalc, "Derivation")
ShinyAppPlot(train_probsTestCalc, "DerivationTest")
ShinyAppPlot(train_probsDervWholeCalc, "DerivationWhole")
