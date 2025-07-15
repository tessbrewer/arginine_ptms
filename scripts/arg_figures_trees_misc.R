#Script for making the trees, colony size complementation,
#phylogenetic relatedness, and alignment figures

#Figures 2a, 3, 6, S1, S2, and S5
#Designed to be run sequentially

# Load functions  =================================

#adjust these to proper paths
source('/Users/tessbrewer/Documents/munich/native_cloning/arginine_cloning_june24/arg_project/scripts/arg_functions.R')
setwd('/Users/tessbrewer/Documents/munich/native_cloning/arginine_cloning_june24/arg_project/')

# EF-P tree (Figure 2) =================================
root = '2595564854_2593339269'
efp_tree <- read.tree("data/trimmed_aligned_total_efp.fa.treefile")
efp_tree <- root(efp_tree, root)

plotted_efp <- read.delim('data/mapping_efp.txt')

#same things in both
plotted_efp <- subset(plotted_efp, tree_gene %in% efp_tree$tip.label & conserved_residue != '-')
efp_tree <- drop.tip(efp_tree, efp_tree$tip.label[!efp_tree$tip.label %in% plotted_efp$tree_gene])

#get color of phyla / color palette
plotted_efp$color <- assign_color(plotted_efp, 'phylum', 11)

palettez <- brewer.pal(n = (length(table(plotted_efp$color)) - 1), name = 'Spectral')
palettez <- append(palettez, 'lightgrey', after = (grep('Other', names(table(plotted_efp$color))) - 1)) #i want "Other" to be grey

#target EF-P
targets <- c('Mesotoga.prima.2510065086', 'Deinococcus.radiodurans.2558118174', 'Akkermansia.muciniphila.2854655614',
             'Herpetosiphon.aurantiacus.2508501111', 'Denitrovibrio.acetiphilus.646564527', 'Clostridium.innocuum.2623620482', 
             'Geoalkalibacter.ferrihydriticus.2599185148', 'Nitrosomonas.communis.2675903044') 

target_list <- subset(plotted_efp, species_code %in% targets & 
                      conserved_residue == 'R' & 
                      ! modification_tree %in% c('arginine-YeiP', 'arginine-EarP'))
target_list$full_id <- paste(target_list$gene_oid, target_list$taxon_oid, sep = "_")


#base tree
efp_plot <- ggtree(efp_tree, open.angle=5, color = "#666666", layout = "fan") %<+% plotted_efp +
  geom_tiplab2(mapping=aes(subset = label %in% target_list$full_id,
                           label = species_label), align = T, 
               offset = 0.2, fontface = 3, size = 3) +
  geom_tippoint(aes(fill = color), pch = 21, size = 3, color = "#666666") +
  scale_fill_manual(values = palettez, name = 'Phylum (GTDB)') +
  theme(text = element_text(family="Avenir", size=15, color = "black")) +
  guides(color = guide_legend(override.aes = list(size=5)))


#add modification types heatmap
heatmap_data <- as.data.frame(plotted_efp[, c("modification_tree")])
row.names(heatmap_data) <- plotted_efp$tree_gene

final_plot <- efp_plot + new_scale_fill()
figure_final <- gheatmap(final_plot, heatmap_data, width=0.1, offset = 0, color=NA, colnames = FALSE) +
  scale_fill_manual(values = c('#77529e', '#8dde9f', '#BEBADA', '#2e88a0', 'pink', '#FFFFB3'), name = 'EF-P type') +
  theme(text = element_text(family="Avenir", size=15, color = "black")) + 
  theme(plot.margin=unit(c(0.05,0,0.1,0), "null")) +
  theme(legend.position = 'left')

#add R types heatmap
plotted_efp$conserved_amino <- ifelse(plotted_efp$conserved_residue == 'R', 'Arginine', 'Other')
heatmap_data2 <- as.data.frame(plotted_efp[, c("conserved_amino")])
row.names(heatmap_data2) <- plotted_efp$tree_gene

final_plot2 <- figure_final + new_scale_fill()
figure_2a <- gheatmap(final_plot2, heatmap_data2, width=0.1, offset = 0.52, color=NA, colnames = FALSE) +
  scale_fill_manual(values = c('#f8766D', 'white'), name = 'Conserved residue') +
  theme(text = element_text(family="Avenir", size=15, color = "black")) + 
  theme(legend.position = 'left')

#name labels adjusted in image editor
# ggsave('../figures/efp_tree_new.svg', figure_2a, width = 10, height = 8)
rm(heatmap_data, heatmap_data2, figure_final, final_plot, final_plot2, root, palettez, 
   targets, efp_plot, efp_tree, plotted_efp)

# Complementation assay (Figure 3) =================================
colony_size <- read.delim('data/colony_size.txt') #combined output of imagej with some ids added
colony_size$biggest <- ifelse(colony_size$Width > colony_size$Height, colony_size$Width, colony_size$Height)

#get 95th percentile, 5th percentile
count <- aggregate(biggest ~ sample, colony_size, length)
quant_upper <- aggregate(biggest ~ sample, colony_size, function(x){quantile(x, 0.95)})
quant_lower <- aggregate(biggest ~ sample, colony_size, function(x){quantile(x, 0.05)})
colony_size$count <- count$biggest[match(colony_size$sample, count$sample)]
colony_size$upperq <- quant_upper$biggest[match(colony_size$sample, quant_upper$sample)]
colony_size$lowerq <- quant_lower$biggest[match(colony_size$sample, quant_lower$sample)]
rm(count, quant_upper, quant_lower)

#remove upper and lower quantiles
colony_size <- subset(colony_size, !(biggest < lowerq) & !(biggest > upperq))

#take top 300 largest
colony_size <- colony_size[order(colony_size$biggest, decreasing = TRUE), ]
colony_size <- subsample(colony_size, 'sample', 300)
table(colony_size$sample) #checks out

#get the average for deltaefp with control plasmid & deltaefp with ecoli EF-P complement
deltaefp_size = mean(subset(colony_size, pbad24 == 'empty24' & pbad33 == 'empty33')$biggest)
deltaefp_ecoli_size = mean(subset(colony_size, pbad24 == 'eco24' & pbad33 == 'empty33')$biggest)

#order labels by biggest
avgie <- aggregate(biggest ~ label, colony_size, mean)
avgie <- avgie[order(avgie$biggest),]
colony_size$label <- factor(colony_size$label, levels = avgie$label)

#check for significant differences between each sample and the Eco EF-P complementation
#correct for multiple comparisions
for (i in names(table(subset(colony_size, pbad33 == 'empty33')$sample))){
  print(i)
  p_correct = length(names(table(subset(colony_size, pbad33 == 'empty33' & sample != 'eco24empty33')$sample)))
  if (i != 'eco24empty33'){
    temp = subset(colony_size, sample %in% c('eco24empty33', i))
    tt = wilcox.test(biggest ~ sample, data = temp)
    print(p.adjust(tt$p.value, method = 'bonferroni', n = p_correct))
  }
}

#check for significant differences between each sample and the delta EF-P
#correct for multiple comparisions
for (i in names(table(subset(colony_size, pbad33 == 'empty33')$sample))){
  print(i)
  p_correct = length(names(table(subset(colony_size, pbad33 == 'empty33' & sample != 'empty24empty33')$sample)))
  if (i != 'empty24empty33'){
    temp = subset(colony_size, sample %in% c('empty24empty33', i))
    tt = wilcox.test(biggest ~ sample, data = temp)
    print(p.adjust(tt$p.value, method = 'bonferroni', n = p_correct))
  }
}

colony_size$rescue_fold <- colony_size$biggest / deltaefp_size

#only the colony size data with empty control plasmid, and none of the controls
figure3_data <- subset(colony_size, ! pbad24 %in% c('empty24', 'eco24', 'ppu24', 'shon24') & pbad33 == 'empty33')

figure_3 <- ggplot(figure3_data, aes(y = label, x = biggest)) +
  geom_vline(xintercept = deltaefp_ecoli_size, linetype = 'dashed', color = '#E95F56', size = 1.1) +
  geom_vline(xintercept = deltaefp_size, linetype = 'dashed', color = '#E95F56', size = 1.1) +
  geom_jitter(aes(fill = pbad33), pch = 21, size = 2, position = position_jitterdodge()) +
  geom_boxplot(aes(fill = pbad33), outlier.shape = NA, alpha = 0.6, width = 0.5) +
  xlab("Colony diameter (mm)") +
  ylab("") +
  ggtitle('') +
  theme_bw() +
  theme(legend.position = 'NONE') +
  scale_fill_manual(values = c("#884EBA")) +
  theme(text = element_text(family="Helvetica", size=16, color = "black")) +
  theme(axis.text = element_text(family="Helvetica", size=16, color = "black")) + 
  theme(axis.text.y = ggtext::element_markdown()) +
  theme(plot.margin=unit(c(0.1,0.05,0,0), "null")) +
  theme(title  = ggtext::element_markdown(family="Helvetica", size=16, color = "black"))

# ggsave('../figures/figure3.svg', figure_3, width = 6, height = 4)
rm(figure3_data, avgie, tt, temp, i, p_correct, deltaefp_ecoli_size)


# Complementation assay with EarP co-expression (Figure 6) =================================

#first get y-position for the p-value placement
maxes <- aggregate(subset(colony_size, pbad24 %in% c('eco24', 'ppu24', 'shon24', 'mepr24')), biggest ~ pbad24, max)

#want them in this order
maxes$pbad24 <- factor(maxes$pbad24, levels = c('eco24', 'shon24', 'ppu24', 'mepr24'))
maxes <- maxes[c(order(maxes$pbad24)),]

#run tests and adjust for multiple comparisions
pvaluez <- c(wilcox.test(biggest ~ sample, subset(colony_size, pbad24 == 'eco24'))$p.value,
             wilcox.test(biggest ~ sample, subset(colony_size, pbad24 == 'shon24'))$p.value,
             wilcox.test(biggest ~ sample, subset(colony_size, pbad24 == 'ppu24'))$p.value,
             wilcox.test(biggest ~ sample, subset(colony_size, pbad24 == 'mepr24'))$p.value)
pvaluez <- p.adjust(pvaluez, method = 'bonferroni', n = length(pvaluez))
#*** for all except eco ns

#earp co-expression plot
figure6_data <- subset(colony_size, pbad24 %in% c('eco24', 'ppu24', 'shon24', 'mepr24'))

figure_6 <- ggplot(figure6_data, aes(y = label, x = biggest)) +
  geom_vline(xintercept = deltaefp_size, linetype = 'dashed', color = '#E95F56', size = 1.1) +
  geom_jitter(aes(fill = pbad33), pch = 21, size = 2, position = position_jitterdodge()) +
  geom_boxplot(aes(fill = pbad33), outlier.shape = NA, alpha = 0.6, width = 0.75) +
  geom_signif(annotation = c('***', '***', '***', 'ns'),
              map_signif_level = TRUE,
              y_position = rev(maxes$biggest) + 0.1, 
              xmin = seq(0.8, 0.8*5, by = 1), 
              xmax = seq(1.2, 1.2*4, by = 1),
              tip_length = c(0.01, 0.01)) +
  xlab("Colony diameter (mm)") +
  ylab("") +
  ggtitle('') +
  theme_bw() +
  theme(legend.justification = c(0.975, 0.025), legend.position = c('inside'),
        legend.position.inside = c(0.975, 0.025)) +
  scale_fill_manual(values = c("#BBECCD", "#884EBA"), name = 'Co-expressed with:',
                    labels = c('EarP', 'Control plasmid')) +
  theme(text = element_text(family="Helvetica", size=14, color = "black")) +
  theme(axis.text = element_text(family="Helvetica", size=14, color = "black")) + 
  theme(axis.text.y = ggtext::element_markdown()) +
  theme(plot.margin=unit(c(0.1,0.05,0,0), "null")) +
  theme(title  = ggtext::element_markdown(family="Helvetica", size=14, color = "black"))

# ggsave('../figures/figure6.svg', figure_6, width = 6, height = 4)
rm(figure6_data, pvaluez, maxes)


# Phylogenetic relatedness vs rescue efficiency in E. coli (Supplemental figure 1) =================================

#load tree and mapping file
related <- read.tree('data/trimmed_ecoli_distance.concat.aln.treefile')
related$tip.label <- gsub('.genes', '', related$tip.label)
related <- root(related, '8002550533') #root with firmicute but doesn't really matter for distance between species
tree_map <- read.delim('data/distances_mapping.txt')

#look at the tree for ~fun~
ggtree(related, open.angle=5, color = "#666666") %<+% tree_map +
  geom_tiplab(aes(label = species), fontface = 3, size = 3, offset = 0.03) +
  xlim_tree(2) #mess with this to reduce size of tree and keep labels from being cutoff

#get distances between species on tree
distances <- melt(as.matrix(cophenetic(related)))
distances <- subset(distances, Var1 == 'BW25113') #subset just for e. coli
colnames(distances)[3] <- 'phylogenetic_distance'

#average rescue efficiency for colony size
avg_rescue <- rbind(aggregate(subset(colony_size, pbad33 == 'empty33' & ! pbad24 %in% c('shon24', 'ppu24')), rescue_fold ~ pbad24, mean),
                    aggregate(subset(colony_size, paste(pbad24, pbad33) %in% c('shon24 earp33', 'ppu24 earp33')), rescue_fold ~ pbad24, mean)) 
#use rhamnosylated versions of ppu and shon, they are natively rhamnosylated

#match pbad24 code for each species
distances$pbad24 <- tree_map$pbad24[match(distances$Var2, tree_map$taxon_oid)]
distances$rescue_colony <- avg_rescue$rescue_fold[match(distances$pbad24, avg_rescue$pbad24)]

#tpph and lppp
hisL_rescue <- rbind(read.delim('data/hisL_combined_empty_pbad33.txt'),
                       subset(read.delim('data/hisL_combined_earp_pbad33_updated.txt')[, 1:12], names == 'ppu24_earp33'))
#again, take the rhamnosylated form for p.putida
hisL_avg <- aggregate(hisL_rescue, fold_rescue ~ pbad24 * motif, mean)
distances$rescue_tpph <- subset(hisL_avg, motif == 'TPPH')$fold_rescue[match(distances$pbad24, subset(hisL_avg, motif == 'TPPH')$pbad24)]
distances$rescue_lppp <- subset(hisL_avg, motif == 'LPPP')$fold_rescue[match(distances$pbad24, subset(hisL_avg, motif == 'LPPP')$pbad24)]

distances$species <- hisL_rescue$species[match(distances$pbad24, hisL_rescue$pbad24)]
distances$species[distances$pbad24 == 'shon24'] <- 'S. oneidensis' #wasn't included in hisL so need to manually assign
distances$species <- gsub('\\*', '', distances$species)

#correlations
colony <- cor.test(distances$phylogenetic_distance, distances$rescue_colony)
tpph <- cor.test(distances$phylogenetic_distance, distances$rescue_tpph)
lppp <- cor.test(distances$phylogenetic_distance, distances$rescue_lppp)

correlation_annotation <- data.frame(facet_names = c('Colony size', 'Motif rescue (TPPH)', 'Motif rescue (LPPP)'), 
                                     pvalues = c(colony$p.value, tpph$p.value, lppp$p.value),
                                     corr = c(colony$estimate, tpph$estimate, lppp$estimate))
correlation_annotation$label <- paste(paste('P =', round(correlation_annotation$pvalues, 3)), paste(' corr =', round(correlation_annotation$corr, 3)), sep = ',')


figures1_data <- melt(distances, id.vars = c('phylogenetic_distance', 'species', 'pbad24'), measure.vars = c('rescue_colony', 'rescue_tpph', 'rescue_lppp'))
figures1_data$facet_names <- ifelse(figures1_data$variable == 'rescue_colony', 'Colony size',
                             ifelse(figures1_data$variable == 'rescue_tpph', 'Motif rescue (TPPH)', 'Motif rescue (LPPP)'))
figures1_data <- na.omit(figures1_data)

figure_s1 <- ggplot(figures1_data, aes(y = phylogenetic_distance, x = value)) +
  facet_wrap(~facet_names, ncol = 1, scales = 'free_x') +
  geom_smooth(method = 'lm', fill = 'lightgrey', color = 'purple') +
  geom_point(aes(fill = species), pch = 21, size = 4) +
  geom_text(data = correlation_annotation, aes(label = label), y = 0, x = 0.5, color = 'purple', fontface = 'bold', family= 'Helvetica') +
  geom_text_repel(aes(label = species), family="Helvetica", size = 4, color = "black", point.padding = 2, box.padding = 0.5) +
  xlab("Rescue efficiency<br>(fold change \u0394*efp* / \u0394*efp* + EF-P)") +
  ylab("Evolutionary distance from *E. coli*") +
  scale_x_log10() +
  scale_fill_brewer(palette = 'Set3') +
  ggtitle('') +
  coord_cartesian(clip = "off") +
  theme_bw() +
  theme(legend.position = 'NONE') +
  theme(text = element_text(family="Helvetica", size=18, color = "black")) +
  theme(axis.text = element_text(family="Helvetica", size=18, color = "black")) + 
  theme(axis.text = ggtext::element_markdown()) +
  theme(title  = ggtext::element_markdown(family="Helvetica", size=18, color = "black"))

# ggsave('../figures/figure_s1.svg', figure_s1, width = 4, height = 10)
#adjusted the positions of the correlation data text to be at less insane places in editor
rm(figures1_data, tree_map, tpph, lppp, correlation_annotation, distances, hisL_avg, 
   hisL_rescue, related, colony, colony_size, deltaefp_size, avg_rescue)

# Alignment of EF-P sequences (Supplemental figure 2) =================================
protein_sequences <- readAAStringSet('data/aligned-all_efp.fa')
blastp_efp <- read.delim('data/blastp_efp.txt')
blastp_efp <- blastp_efp[order(blastp_efp$pident, decreasing = TRUE),]

#order them by percent identity except put efpL next EF-P
right_order <- blastp_efp$sseqid[-grep('Escherichia_coli.2619170694', blastp_efp$sseqid)]
right_order <- append(right_order, 'Escherichia_coli.2619170694', after = (grep('Escherichia_coli.2619172497', right_order)))
protein_sequences <- protein_sequences[right_order]

#name EF-P vs EfpL for E. coli
names(protein_sequences) <- gsub('Escherichia_coli.2619172497', 'Escherichia_coli (EF-P)', names(protein_sequences))
names(protein_sequences) <- gsub('Escherichia_coli.2619170694', 'Escherichia_coli (EfpL)', names(protein_sequences))

#clean up names, IMG id not so useful for figure
names(protein_sequences) <- sapply(strsplit(as.character(names(protein_sequences)), '\\.'), '[', 1)
names(protein_sequences) <- gsub('_', ' ', names(protein_sequences))

figure_s2 <- ggmsa(protein_sequences, 30, 71, seq_name = TRUE, char_width = 0.5, border = "white",
                 font = "helvetical", position_highlight = c(38, 54, 67:69), 
                 color = 'Chemistry_AA', show.legend = FALSE) +
          geom_rect(aes(xmin = 34.5, xmax = 35.5, ymin = -Inf, ymax = Inf), color = 'purple', fill = NA, size = .75) +
          geom_rect(aes(ymin = 9.5, ymax = 11.5, xmin = -Inf, xmax = Inf), alpha = 0.2, fill = 'grey', color = NA, size = .75)

#couple annotations added in editor 
# ggsave('../figure_s2.svg', figure_s2, width = 10, height = 4)
rm(right_order, protein_sequences, blastp_efp)

# Gammaproteobacteria tree with EF-P type (Supplemental figure 5) =================================

phylogeny_tree <- read.tree("data/full_species.tre")
phylogeny_tree <- root(phylogeny_tree, "GCF_000009185.1_ASM918v1_genomic") # Haloquadratum walsbyi (archaea)

# unifying species in mapping file and tree
mapping <- read.delim('data/full_species_mapping.txt')
plotted_map <- subset(mapping, phylum_gtdbtk == 'Gammaproteobacteria')
plotted_tree <- ape::drop.tip(phy = phylogeny_tree, 
                              tip = phylogeny_tree$tip.label[ ! phylogeny_tree$tip.label %in% plotted_map$taxon_oid])

# picking what should be indicated on legend
plotted_map$color <- ifelse(plotted_map$order_gtdbtk == 'Betaproteobacteriales', 'Betaproteobacteriales', 'Gammaproteobacteria')

#in this case the distinction between genomes with epmC or not is not improtant
plotted_map$primary_efp <- gsub('lysine-EpmAB$', 'lysine-EpmAB/C', plotted_map$primary_efp)
plotted_map$primary_efp <- gsub('lysine-EpmABC', 'lysine-EpmAB/C', plotted_map$primary_efp)

#find nodes to highlight, non Betaproteobacteriales with earP
look <- subset(plotted_map, grepl('EarP', efp_modification) & order_gtdbtk != 'Betaproteobacteriales')
sort(table(look$family_gtdbtk))

# getting nodes for the many many MANY EarP in Gammas
a_event <- getMRCA(plotted_tree, tip = as.character(subset(plotted_map, family_gtdbtk == 'Aeromonadaceae' & efp_modification == 'EarP')$taxon_oid))
b_event <- getMRCA(plotted_tree, tip = as.character(subset(plotted_map, family_gtdbtk == 'Alteromonadaceae' & efp_modification == 'EarP')$taxon_oid))
c_event <- getMRCA(plotted_tree, tip = as.character(subset(plotted_map, family_gtdbtk == 'Pseudomonadaceae' & efp_modification == 'EarP')$taxon_oid))
d_event <- getMRCA(plotted_tree, tip = as.character(subset(plotted_map, family_gtdbtk == 'Thiomicrospiraceae' & efp_modification == 'EarP')$taxon_oid))
e_event <- getMRCA(plotted_tree, tip = as.character(subset(plotted_map, family_gtdbtk == 'Succinivibrionaceae' & efp_modification == 'EarP')$taxon_oid))
f_event <- getMRCA(plotted_tree, tip = as.character(subset(plotted_map, family_gtdbtk == 'Moraxellaceae' & efp_modification == 'EarP')$taxon_oid))
g_event <- getMRCA(plotted_tree, tip = as.character(subset(plotted_map, family_gtdbtk == 'Shewanellaceae' & efp_modification == 'EarP')$taxon_oid))
h_event <- getMRCA(plotted_tree, tip = as.character(subset(plotted_map, family_gtdbtk == 'Cellvibrionaceae' & efp_modification == 'EarP')$taxon_oid))
i_event <- getMRCA(plotted_tree, tip = as.character(subset(plotted_map, family_gtdbtk == 'Spongiibacteraceae' & grepl('EarP', efp_modification))$taxon_oid))
j_event <- getMRCA(plotted_tree, tip = as.character(subset(plotted_map, grepl('Plasticicumulans', species) & grepl('EarP', efp_modification))$taxon_oid))
k_event <- getParent(plotted_tree, tips = as.character(subset(plotted_map, species == 'Thiofilum flexile' & grepl('EarP', efp_modification))$taxon_oid))
l_event <- getParent(plotted_tree, tips = as.character(subset(plotted_map, species == 'Oceanospirillum multiglobuliferum' & grepl('EarP', efp_modification))$taxon_oid))


highlight_data <- data.frame(node = c(a_event, b_event, c_event, d_event, e_event, f_event, g_event, h_event, i_event, j_event, k_event, l_event), 
                             type = c("Aeromonadaceae", "Alteromonadaceae", 'Pseudomonadaceae', 'Thiomicrospiraceae',
                                      'Succinivibrionaceae', 'Moraxellaceae', 'Shewanellaceae', 'Cellvibrionaceae',
                                      'Spongiibacteraceae', 'Plasticicumulans', 'Thiofilum flexile', 'Oceanospirillum multiglobuliferum'))

tree_plot <- ggtree(plotted_tree, open.angle=20, color = "#666666", layout = "fan") %<+% plotted_map +
  theme(text = element_text(family="Avenir", size=16, color = "black")) +
  geom_highlight(data=highlight_data, mapping=aes(node=node), fill = '#bebada', extendto = 1.1) +
  geom_tippoint(aes(fill = color), pch = 21, size =3, color = "#666666") +
  scale_fill_manual(values = c('#EE5840', '#A6D8DF'), name = 'Group') +
  scale_color_manual(values = brewer.pal(n = 10, name = 'Spectral'), name = 'Phylum (GTDB)') +
  guides(color = guide_legend(override.aes = list(size=5)))

#add heatmap with simplified modifications (epmAB / epmABC / earP / YmfI / unmodified loop / unknown)
heatmap_data3 <- plotted_map[, c("primary_efp", "secondary_efp", "tertiary_efp")]
row.names(heatmap_data3) <- plotted_map$taxon_oid

new_plot <- tree_plot + new_scale_fill()
figure_s5 <- gheatmap(new_plot, heatmap_data3, width=0.21, color=NA, colnames = FALSE) +
  scale_fill_manual(values = c('#77529e', '#8dde9f', '#BEBADA', '#2e88a0'),
                    labels = c('arginine-EarP', 'arginine-EfpL', 'unknown', 'lysine-EpmAB/C', ''), 
                    name = 'EF-P type', na.value = 'white') +
  theme(text = element_text(family="Avenir", size=18, color = "black"))


# ggsave("../figures/figure_s5.svg", plot=figure_s5, width = 14, height = 12)
# labels added in inkscape and some highlights adjusted
rm(list = ls()[grepl('event', ls())], heatmap_data3, highlight_data, look, new_plot, tree_plot, 
   plotted_map, mapping, phylogeny_tree, plotted_tree)



