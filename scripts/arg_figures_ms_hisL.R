#Script for making the mass-spec and hisL motif figures

#Figures S3, 5, 4, and 7
#Designed to be run sequentially

# Load functions  =================================

#adjust these to proper paths
source('/Users/tessbrewer/Documents/munich/native_cloning/arginine_cloning_june24/arg_project/scripts/arg_functions.R')
setwd('/Users/tessbrewer/Documents/munich/native_cloning/arginine_cloning_june24/arg_project/')

# Massspec, non-rhamnosylated EF-P (Figure S3) =================================
ms_data <- read.delim('data/massspec_data.txt')

#masses for different EF-P
rhamnose = 146.058
methionine = 131.0405
masses <- data.frame(efp = c('heau', 'mepr', 'nico', 'dera', 'gefeR', 'deac', 'eco', 'ppu'), 
                     species = c('Herpetosiphon aurantiacus', 'Mesotoga prima', 'Nitrosomonas communis',
                                 'Deinococcus radiodurans', 'Geoalkalibacter ferrihydriticus', 'Denitrovibrio acetiphilus',
                                 'Escherichia coli', 'Psuedomonas putida'),
                     monoisotopic_mass = c(21756.036, 21935.070, 22356.464 - methionine, 21428.606, 21556.073, 
                                           22119.936, 21544.819 - methionine, 22256.116))

masses$rhamnosylated <- masses$monoisotopic_mass + rhamnose
rm(rhamnose, methionine)


ms_data <- subset(ms_data, abs(Monoisotopic_Mass - target_mass) <= 250)
#subset for relevant peaks -- within 250 Da of each EF-P
ms_data$species <- masses$species[match(ms_data$efp, masses$efp)]
ms_data$Rel_neg <- ifelse(ms_data$earp == 'plus', ms_data$Relative_Abundance, -ms_data$Relative_Abundance)
#this puts the relative abundance upside down for the control plasmid co-expression in figure 5

ms_data$rhamnosylated_mass <- masses$rhamnosylated[match(ms_data$efp, masses$efp)]
ms_data$label <- ifelse(ms_data$earp == 'minus', 'Empty control plasmid', 'EarP')
ms_data$species_ital <- paste('*', ms_data$species, sep = '')
ms_data$species_ital <- paste(ms_data$species_ital, '*', sep = '')
ms_data$peak_status <- ifelse(abs(ms_data$Monoisotopic_Mass - ms_data$target_mass) < 1, 'naked_peak',
                            ifelse(abs(ms_data$Monoisotopic_Mass - ms_data$rhamnosylated_mass) < 1 & ms_data$earp == 'plus', 'rhamnosylated_peak', 'other'))
masses$species_ital <- ms_data$species_ital[match(masses$efp, ms_data$efp)]

#EF-P with no rhamnosylation
no_rhamn <- subset(ms_data, ! efp %in% c('mepr', 'ppu', 'shon') & Relative_Abundance > 10)
no_rhamn$found_label <- paste("MW['Found'] =", round(no_rhamn$Monoisotopic_Mass, 2))
nr_masses <- subset(masses, efp %in% no_rhamn$efp)

figure_s3 <- ggplot(no_rhamn) +
  facet_wrap(~species_ital, ncol = 2, scales = 'free_x', strip.position = 'right') +
  geom_linerange(data = no_rhamn, mapping = aes(x = Monoisotopic_Mass, ymax = Rel_neg, ymin=0, color = label), size = 1.1) + 
  geom_point(data = nr_masses, mapping = aes(x = monoisotopic_mass), y = 0, pch = 19, color = '#666666', size = 2) +
  geom_segment(data = nr_masses, mapping = aes(x = monoisotopic_mass, xend = rhamnosylated), y = 0, color = '#666666', 
               size = 1) +
  geom_point(data = nr_masses, mapping = aes(x = rhamnosylated), y = 0, pch = 4, stroke = 2, color = '#666666', size = 1) +
  geom_richtext(data = subset(no_rhamn, peak_status == 'naked_peak'), mapping = aes(label = paste('MW<sub>Found</sub> =', round(Monoisotopic_Mass, 2)),
                x = Monoisotopic_Mass + 10), y = 90, hjust = 0, label.color = NA, family="Helvetica", size = 5) + #mass closest to calculated mass
  geom_richtext(data = subset(no_rhamn, peak_status == 'naked_peak'), mapping = aes(label = paste('MW<sub>EF-P</sub> =', round(target_mass, 2)), 
                x = Monoisotopic_Mass + 10), y = 75, hjust = 0, label.color = NA, family="Helvetica", size = 5) +
  ylab('Relative Abundance') +
  xlab('Monoisotopic Mass (Da)') +
  scale_color_manual(values = c("#BBECCD"), name = 'Co-expressed with:') +
  scale_y_continuous(breaks = c(-100, -50, 0, 50, 100), labels = c(100, 50, 0, 50, 100)) +
  theme_classic() +
  theme(text = element_text(family="Helvetica", size=16), axis.text = element_text(family="Helvetica", size=16),
        strip.text = element_text(family="Helvetica", size=16)) + 
  theme(strip.text = ggtext::element_markdown(), panel.border = element_rect(colour = "black", fill=NA),
        axis.ticks.x = element_blank(), axis.text.x = element_blank())

# ggsave('../figures/figures3.svg', figure_s3, width = 9, height = 9)
rm(no_rhamn, nr_masses)

# Massspec, rhamnosylated EF-P (Figure 5) =================================

rhamn <- subset(ms_data, efp %in% c('mepr', 'ppu') & Relative_Abundance > 15)
rhamn_mass <- subset(masses, efp %in% rhamn$efp)

figure_5 <- ggplot(rhamn) +
  facet_wrap(~species_ital, ncol = 1, scales = 'free_x', strip.position = 'right') +
  geom_linerange(data = rhamn, mapping = aes(x = Monoisotopic_Mass, ymax = Rel_neg, ymin=0, color = label), size = 1.1) + 
  geom_point(data = rhamn_mass, mapping = aes(x = monoisotopic_mass), y = 0, pch = 19, color = '#666666', size = 2) +
  geom_segment(data = rhamn_mass, mapping = aes(x = monoisotopic_mass, xend = rhamnosylated), y = 0, color = '#666666', 
               arrow = arrow(length=unit(0.20,"cm"), ends="last", type = "closed"), size = 1) +
  geom_segment(data = rhamn_mass, mapping = aes(x = monoisotopic_mass, xend = rhamnosylated), y = -0.5, color = '#666666') +
  #empty plasmid
  geom_richtext(data = subset(rhamn, peak_status == 'naked_peak'), mapping = aes(label = paste('MW<sub>Found</sub> =', round(Monoisotopic_Mass, 2)),
                   x = Monoisotopic_Mass + 10), y = -90, hjust = 0, label.color = NA, family="Helvetica", size = 5) + #mass closest to calculated mass
  geom_richtext(data = subset(rhamn, peak_status == 'naked_peak'), mapping = aes(label = paste('MW<sub>EF-P</sub> =', round(target_mass, 2)),
                   x = Monoisotopic_Mass + 10), y = -65, hjust = 0, label.color = NA, family="Helvetica", size = 5) + 
  #earP
  geom_richtext(data = subset(rhamn, peak_status == 'rhamnosylated_peak'), mapping = aes(label = paste('MW<sub>Found</sub> =', round(Monoisotopic_Mass, 2)),
                   x = Monoisotopic_Mass - 10), y = 90, hjust = 1, label.color = NA, family="Helvetica", size = 5) + #mass closest to rhamnosylated mass
  geom_richtext(data = subset(rhamn, peak_status == 'rhamnosylated_peak'), mapping = aes(label = paste('\\+', round(Monoisotopic_Mass - target_mass, 2)),
                   x = Monoisotopic_Mass - 10), y = 65, hjust = 1, label.color = NA, family="Helvetica", size = 5) + #difference between target and rhamnosylated mass
  ylab('Relative Abundance') +
  xlab('Monoisotopic Mass (Da)') +
  scale_color_manual(values = c("#BBECCD", '#ad85cf'), name = 'Co-expressed with:') +
  scale_y_continuous(breaks = c(-100, -50, 0, 50, 100), labels = c(100, 50, 0, 50, 100)) +
  theme_classic() +
  theme(text = element_text(family="Helvetica", size=16), axis.text = element_text(family="Helvetica", size=16),
        strip.text = element_text(family="Helvetica", size=16)) + 
  theme(strip.text = ggtext::element_markdown(), panel.border = element_rect(colour = "black", fill=NA),
        axis.ticks.x = element_blank(), axis.text.x = element_blank())

# ggsave('../figures/figure5.svg', figure_5, width = 6, height = 6)
rm(masses, ms_data, rhamn, rhamn_mass)


# HisL figure w/out co-expression (Figure 4) =================================
hisL_combined <- read.delim('data/hisL_combined_empty_pbad33.txt')
avgz <- aggregate(hisL_combined, fold_rescue ~ label, median)
avgz <- avgz[order(avgz$fold_rescue, decreasing = T),]
hisL_combined$label <- factor(hisL_combined$label, levels = avgz$label)
#order in terms of average rescue efficiency

hisL_combined$motif_label <- ifelse(grepl('RPDG', hisL_combined$motif), 'Non-stalling motif (RPDG)',
                                    ifelse(grepl('LPPP', hisL_combined$motif), 'Stalling motif (LPPP)', 'Stalling motif (TPPH)'))
hisL_combined$motif_label <- factor(hisL_combined$motif_label, levels = c('Stalling motif (LPPP)', 'Stalling motif (TPPH)', 'Non-stalling motif (RPDG)'))
#order in terms of stalling strength

figure_4 <- ggplot(hisL_combined, aes(x = label, y = fold_rescue)) +
  facet_wrap(~motif_label, ncol = 1) +
  geom_hline(yintercept = 1, color = 'red', linetype = 'dashed') +
  geom_jitter(aes(fill = pbad33), pch = 21, size = 3, position = position_jitterdodge()) +
  geom_boxplot(aes(fill = pbad33), outlier.shape = NA, alpha = 0.6, width = 0.5) +
  ylab("fold change \u0394*efp* / \u0394*efp* + EF-P") +
  xlab("") +
  ggtitle('') +
  theme_bw() +
  theme(legend.position = 'NONE') +
  scale_y_log10() +
  scale_fill_manual(values = c("#884EBA"), name = 'Co-expressed with:') +
  theme(text = element_text(family="Helvetica", size=16, color = "black")) +
  theme(axis.text = element_text(family="Helvetica", size=16, color = "black")) + 
  theme(strip.text = element_text(family="Helvetica", size=16, color = "black"), 
        panel.border = element_rect(colour = "black", fill=NA)) + 
  theme(axis.text.x = ggtext::element_markdown(angle = 45, hjust = 1)) +
  theme(title  = ggtext::element_markdown(family="Helvetica", size=16, color = "black"))

# ggsave('../figures/figure4.svg', figure_4, width = 5.5, height = 10)

# HisL figure with co-expression (Figure 7) =================================
hisL_earp <- read.delim('data/hisL_combined_earp_pbad33_updated.txt')

#label order i want
hisL_earp$label <- factor(hisL_earp$label, 
                          levels = c('Δefp +<br>*E. coli* EF-P', 
                                     'Δefp +<br>*P. putida* EF-P',
                                     'Δefp +<br>*M. prima* EF-P'))

#motif order i want
hisL_earp$motif_label <- factor(hisL_earp$motif_label, levels = c('Stalling motif (LPPP)', 'Stalling motif (TPPH)', 'Non-stalling motif (RPDG)'))

figure_7 <- ggplot(hisL_earp, aes(x = label, y = fold_rescue, fill = pbad33_label)) +
  facet_wrap(~motif_label, ncol = 1) +
  geom_hline(yintercept = 1, color = 'red', linetype = 'dashed') +
  geom_jitter(aes(fill = pbad33_label), pch = 21, size = 3, position = position_jitterdodge()) +
  geom_boxplot(aes(fill = pbad33_label), outlier.shape = NA, alpha = 0.6, width = 0.5) +
  ylab("fold change \u0394*efp* / \u0394*efp* + EF-P") +
  xlab("") +
  ggtitle('') +
  theme_bw() +
  scale_y_log10(limits = c(0.25, 25)) +
  scale_fill_manual(values = c("#BBECCD", '#884EBA'), name = 'Co-expressed with:') +
  theme(text = element_text(family="Helvetica", size=16, color = "black")) +
  theme(axis.text = element_text(family="Helvetica", size=16, color = "black")) + 
  theme(strip.text = element_text(family="Helvetica", size=16, color = "black"), 
        panel.border = element_rect(colour = "black", fill=NA)) + 
  theme(axis.text.x = ggtext::element_markdown(angle = 45, hjust = 1)) +
  theme(title  = ggtext::element_markdown(family="Helvetica", size=16, color = "black"))

#check for significant differences between each sample earP and control
#correct for multiple comparisions
for (i in names(table(hisL_earp$pbad24))){
  print(i)
  p_correct = length(names(table(hisL_earp$pbad24)))
  temp = subset(hisL_earp, pbad24 == i)
  for (motify in names(table(hisL_earp$motif))){
    print(motify)
    temp2 = subset(temp, motif == motify)
    tt = t.test(fold_rescue ~ pbad33, data = temp2)
    print(p.adjust(tt$p.value, method = 'bonferroni', n = 3))
  }
}

#p-value labels set with the above data in inkscape because ggsignif isn't so good with facetting AND dodging 
# ggsave('../figures/figure7.svg', figure_7, width = 7, height = 10)
rm(avgz, hisL_combined, hisL_earp, temp, temp2, tt, i, motify, p_correct)



