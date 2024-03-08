library(readxl)
library(tidyr)
metabData<-read_xlsx("data_examples/metabolica_caseStudy/JCI71180sd2.xlsx", sheet = "ScaledImpData", skip = 1)
metabData2 <- as_tibble(metabData) %>% drop_na(6)  %>% 
  .[-1,colnames(.) =="...6" | !is.na(as.numeric(colnames(.)))] %>% 
  dplyr::rename("keggID"="...6") %>% distinct(keggID, .keep_all = T) %>%
  column_to_rownames("keggID")

