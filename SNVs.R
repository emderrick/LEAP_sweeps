library(tidyverse)
library(dplyr)

mag_list <- list("I4_MAG_00006", "I4_MAG_00065", "L2_MAG_00052", "L3_MAG_00058", "L4_MAG_00099",
                 "L7_MAG_00020", "L7_MAG_00028", "L7_MAG_00043", "L8_MAG_00011", "L8_MAG_00019", "L8_MAG_00042")

ponds_list <- list(I4_MAG_00006 = list("Control B at T2", "Control E at T2", "GBH A at T2", "GBH D at T2"), 
                   I4_MAG_00065 = list("Control A at T2", "Control B at T2", "Control C at T2", "Control D at T2", "Control E at T2", "GBH B at T2"), 
                   L2_MAG_00052 = list("Control A at T2", "Control B at T2", "Control D at T2", "Control E at T2", "GBH A at T2"), 
                   L3_MAG_00058 = list("Control C at T2", "Control D at T2", "GBH C at T2", "GBH D at T2"), 
                   L4_MAG_00099 = list("Control D at T2", "GBH A at T2", "GBH B at T2", "GBH C at T2", "GBH D at T2"), 
                   L7_MAG_00020 = list("Control A at T1", "Control A at T2", "Control C at T1", "Control C at T2", "Control D at T1", "Control D at T2", "GBH A at T1", "GBH A at T2", "GBH C at T1"), 
                   L7_MAG_00028 = list("Control E at T2", "GBH A at T2", "GBH B at T2", "GBH C at T2"), 
                   L7_MAG_00043 = list("Control D at T2", "GBH B at T2", "GBH C at T2"), 
                   L8_MAG_00011 = list("Control E at T2", "GBH A at T2", "GBH D at T2"), 
                   L8_MAG_00019 = list("Control E at T2", "GBH A at T2", "GBH B at T2", "GBH D at T2"), 
                   L8_MAG_00042 = list("Control A at T2", "Control C at T2", "Control D at T2", "Control E at T2", "GBH D at T2"))

all_MAG_SNVs <- read_csv("all_MAG_SNVs_subsamp.csv")

small_MAG_SNVs <- all_MAG_SNVs[, c('scaffold', 'length', 'groups', 'gene', 'mag', 'mag_length', 'new_name', 'final_ref_freq')]
small_MAG_SNVs_hor <- pivot_wider(small_MAG_SNVs, names_from = "new_name", values_from = "final_ref_freq")
small_MAG_SNVs_hor$all_mean <- rowMeans(small_MAG_SNVs_hor[,c('Control A at T2', 'Control B at T2', 'Control C at T2', 'Control D at T2', 'Control E at T2', 'GBH A at T2', 'GBH B at T2', 'GBH C at T2', 'GBH D at T2')], na.rm = T)
small_MAG_SNVs_hor$control_mean <- rowMeans(small_MAG_SNVs_hor[,c('Control A at T2', 'Control B at T2', 'Control C at T2', 'Control D at T2', 'Control E at T2')], na.rm = T)
small_MAG_SNVs_hor$GBH_mean <- rowMeans(small_MAG_SNVs_hor[,c('GBH A at T2', 'GBH B at T2', 'GBH C at T2', 'GBH D at T2')], na.rm = T)
small_MAG_SNVs_hor$abs_val <- abs(small_MAG_SNVs_hor$control_mean - small_MAG_SNVs_hor$GBH_mean)
small_MAG_SNVs_hor[sapply(small_MAG_SNVs_hor, is.nan)] <- NA
write.csv(small_MAG_SNVs_hor, "small_MAG_SNVs_hor_means_subsamp.csv", row.names = F)

threshold_snvs <- subset(small_MAG_SNVs_hor, abs_val >= 0.5)
threshold_snvs$control <- with(threshold_snvs, ifelse(control_mean < 0.5, 'low', 'high'))
threshold_snvs$control <- with(threshold_snvs, ifelse(control_mean == 0.5, 0.5, control))

threshold_snvs$control_ref <- with(threshold_snvs, ifelse(control == "high", 
                                                         apply(threshold_snvs[,c('Control A at T2', 'Control B at T2', 'Control C at T2', 'Control D at T2', 'Control E at T2')], 1, min, na.rm = T), 
                                                         apply(threshold_snvs[,c('Control A at T2', 'Control B at T2', 'Control C at T2', 'Control D at T2', 'Control E at T2')], 1, max, na.rm = T)))
threshold_snvs$control_ref <- with(threshold_snvs, ifelse(control == 0.5, ifelse(GBH_mean >= 0.5, 
                                                                              apply(threshold_snvs[,c('Control A at T2', 'Control B at T2', 'Control C at T2', 'Control D at T2', 'Control E at T2')], 1, max, na.rm = T), 
                                                                              apply(threshold_snvs[,c('Control A at T2', 'Control B at T2', 'Control C at T2', 'Control D at T2', 'Control E at T2')], 1, min, na.rm = T)), control_ref))
threshold_snvs$GBH_ref <- with(threshold_snvs, ifelse(control == "high", 
                                                      apply(threshold_snvs[,c('GBH A at T2', 'GBH B at T2', 'GBH C at T2', 'GBH D at T2')], 1, max, na.rm = T), 
                                                      apply(threshold_snvs[,c('GBH A at T2', 'GBH B at T2', 'GBH C at T2', 'GBH D at T2')], 1, min, na.rm = T)))
threshold_snvs$GBH_ref <- with(threshold_snvs, ifelse(control == 0.5, ifelse(GBH_mean >= 0.5, 
                                                                          apply(threshold_snvs[,c('GBH A at T2', 'GBH B at T2', 'GBH C at T2', 'GBH D at T2')], 1, min, na.rm = T), 
                                                                          apply(threshold_snvs[,c('GBH A at T2', 'GBH B at T2', 'GBH C at T2', 'GBH D at T2')], 1, max, na.rm = T)), GBH_ref))
threshold_snvs$pass <- with(threshold_snvs, ifelse(((control == "high" & GBH_ref < control_ref) | (control == "low" & GBH_ref > control_ref) |
                                                    (control == 0.5 & GBH_mean >= 0.5 & GBH_ref > control_ref) | (control == 0.5 & GBH_mean < 0.5 & GBH_ref < control_ref)), "yes", "no"))

write.csv(threshold_snvs, "threshold_snvs_subsamp.csv", row.names = F)

small_mag_hor_pass <- left_join(small_MAG_SNVs_hor, threshold_snvs[,c('groups', 'pass', 'control')], by = c('groups'))
write.csv(small_mag_hor_pass, "small_mag_hor_pass_subsamp.csv", row.names = F)


all_MAG_SNVs$class <- all_MAG_SNVs$class%>%str_sub(-3,-1)

all_finished_SNVs <- data.frame()
for(MAG in mag_list){
  MAG_SNVs <- subset(small_mag_hor_pass, mag == MAG)
  MAG_SNVs <- MAG_SNVs[, c(1:6, which((names(MAG_SNVs) %in% ponds_list[[MAG]]) == T), 21:26)] %>%
  pivot_longer(cols = contains("at"),  names_to = "new_name", values_to = "final_ref_freq", values_drop_na = F) %>%
  left_join(all_MAG_SNVs)
  all_finished_SNVs <- rbind(all_finished_SNVs, MAG_SNVs)
}

write.csv(all_finished_SNVs, "all_MAG_SNVs_med_Apr9.csv", row.names = F)
