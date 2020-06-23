library(MASS)
library(nichenetr)
library(Seurat)
library(tidyverse)
library(SingleCellExperiment)
library(ggplot2)
library(cowplot)
library(ggpubr)

ligand_target_matrix = readRDS("/work/shah/ceglian/interactions/ligand_target_matrix.rds")
lr_network = readRDS("/work/shah/ceglian/interactions/lr_network.rds")
weighted_networks = readRDS("/work/shah/ceglian/interactions/weighted_networks.rds")
ligand_tf_matrix = readRDS("/work/shah/ceglian/interactions/ligand_tf_matrix.rds")
sig_network = readRDS("/work/shah/ceglian/interactions/signaling_network.rds")
gr_network = readRDS("/work/shah/ceglian/interactions/gr_network.rds")

weighted_networks_lr = weighted_networks$lr_sig %>% inner_join(lr_network %>% distinct(from,to), by = c("from","to"))

args = commandArgs(trailingOnly=TRUE)
sender <- args[1]
receiver <- args[2]
sender_obj <- readRDS(args[3])
receiver_obj <- readRDS(args[4])
read_genes <- read.csv(args[5])
geneset_oi <- as.character(read_genes$genes)
pathway <- args[6]
png <- args[7]
print(geneset_oi)
combined <- merge(sender_obj, y = receiver_obj, project = "RNASCP")
sce <- as.SingleCellExperiment(combined)
expression <- logcounts(sce)
expression <- t(expression)
geneset_oi <- intersect(geneset_oi, colnames(expression))
#ids sets
sender_ids <- colnames(combined[,combined$cell_type==sender])
receiver_ids <- colnames(combined[,combined$cell_type==receiver])
sample_cell_df <- data.frame(sample=combined$sample, cell=colnames(combined))

#setup
expressed_genes_receiver = get_expressed_genes(receiver, combined, pct = 0.10)
background_expressed_genes = expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]
expressed_genes_sender = get_expressed_genes(sender, combined, pct = 0.05)
ligands = lr_network %>% pull(from) %>% unique()
receptors = lr_network %>% pull(to) %>% unique()
expressed_ligands = intersect(ligands,expressed_genes_sender)
expressed_receptors = intersect(receptors,expressed_genes_receiver)
lr_network_expressed = lr_network %>% filter(from %in% expressed_ligands & to %in% expressed_receptors)
potential_ligands = lr_network_expressed %>% pull(from) %>% unique()

#predict
ligand_activities = predict_ligand_activities(geneset = geneset_oi, background_expressed_genes = background_expressed_genes, ligand_target_matrix = ligand_target_matrix, potential_ligands = potential_ligands)


best_upstream_ligands = ligand_activities %>% top_n(20, pearson) %>% arrange(-pearson) %>% pull(test_ligand)
active_ligand_target_links_df = best_upstream_ligands %>% lapply(get_weighted_ligand_target_links,geneset = geneset_oi, ligand_target_matrix = ligand_target_matrix, n = 250) %>% bind_rows()
active_ligand_target_links = prepare_ligand_target_visualization(ligand_target_df = active_ligand_target_links_df, ligand_target_matrix = ligand_target_matrix, cutoff = 0.25)
order_ligands = intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev()
order_ligands <- intersect(rownames(sce),order_ligands)
order_targets = active_ligand_target_links_df$target %>% unique()
order_targets <- intersect(rownames(sce),order_targets)
vis_ligand_target = active_ligand_target_links[order_targets,order_ligands] %>% t()
p_ligand_target_network = vis_ligand_target %>% make_heatmap_ggplot(paste0("Prioritized ",sender," ligands"),paste0(pathway," genes in ",receiver), color = "purple",legend_position = "top", x_axis_position = "top",legend_title = "Regulatory potential") + scale_fill_gradient2(low = "whitesmoke",  high = "purple", breaks = c(0,0.005,0.01)) + theme(axis.text.x = element_text(face = "italic"))

#top ranked
lr_network_top = lr_network %>% filter(from %in% best_upstream_ligands & to %in% expressed_receptors) %>% distinct(from,to)
best_upstream_receptors = lr_network_top %>% pull(to) %>% unique()
lr_network_top_df = weighted_networks$lr_sig %>% filter(from %in% best_upstream_ligands & to %in% best_upstream_receptors)
lr_network_top_df = lr_network_top_df %>% spread("from","weight",fill = 0)
lr_network_top_matrix = lr_network_top_df %>% select(-to) %>% as.matrix() %>% magrittr::set_rownames(lr_network_top_df$to)
dist_receptors = dist(lr_network_top_matrix, method = "binary")
hclust_receptors = hclust(dist_receptors, method = "ward.D2")
order_receptors = hclust_receptors$labels[hclust_receptors$order]
dist_ligands = dist(lr_network_top_matrix %>% t(), method = "binary")
hclust_ligands = hclust(dist_ligands, method = "ward.D2")
order_ligands_receptor = hclust_ligands$labels[hclust_ligands$order]
vis_ligand_receptor_network = lr_network_top_matrix[order_receptors, order_ligands_receptor]
p_ligand_receptor_network = vis_ligand_receptor_network %>% t() %>% make_heatmap_ggplot(paste0("Prioritized ",sender," ligands"),paste0("Receptors expressed by ", receiver), color = "mediumvioletred", x_axis_position = "top",legend_title = "Prior interaction potential")

ligand_pearson_matrix = ligand_activities %>% select(pearson) %>% as.matrix() %>% magrittr::set_rownames(ligand_activities$test_ligand)
vis_ligand_pearson = ligand_pearson_matrix[order_ligands, ] %>% as.matrix(ncol = 1) %>% magrittr::set_colnames("Pearson")
p_ligand_pearson = vis_ligand_pearson %>% make_heatmap_ggplot(paste0("Prioritized ",sender," ligands"),"Ligand activity", color = "darkorange",legend_position = "top", x_axis_position = "top", legend_title = "Pearson correlation coefficient\n(target gene prediction ability)")

order_samples <- unique(sce[,sender_ids]$sample)

expression_df_sender = expression[sender_ids,order_ligands] %>% data.frame() %>% rownames_to_column("cell") %>% tbl_df() %>% inner_join(sample_cell_df, by =  "cell")
aggregated_expression_sender = expression_df_sender %>% group_by(sample) %>% select(-cell) %>% summarise_all(mean)
aggregated_expression_df_sender = aggregated_expression_sender %>% select(-sample) %>% t() %>% magrittr::set_colnames(aggregated_expression_sender$sample) %>% data.frame() %>% rownames_to_column("ligand") %>% tbl_df()
aggregated_expression_matrix_sender = aggregated_expression_df_sender %>% select(-ligand) %>% as.matrix() %>% magrittr::set_rownames(aggregated_expression_df_sender$ligand)
order_samples <- intersect(colnames(aggregated_expression_matrix_sender),order_samples)
order_ligands <- intersect(rownames(aggregated_expression_matrix_sender),order_ligands)
vis_ligand_tumor_expression = aggregated_expression_matrix_sender[order_ligands,order_samples]
library(RColorBrewer)
color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100)
p_ligand_tumor_expression = vis_ligand_tumor_expression %>% make_heatmap_ggplot(paste0("Prioritized ",sender," ligands"),"Sample", color = color[100],legend_position = "top", x_axis_position = "top", legend_title = "Expression\n(averaged over\nsingle cells)") + theme(axis.text.y = element_text(face = "italic"))


expression_df_target = expression[receiver_ids,geneset_oi] %>% data.frame() %>% rownames_to_column("cell") %>% tbl_df() %>% inner_join(sample_cell_df, by =  "cell")
aggregated_expression_target = expression_df_target %>% group_by(sample) %>% select(-cell) %>% summarise_all(mean)
aggregated_expression_df_target = aggregated_expression_target %>% select(-sample) %>% t() %>% magrittr::set_colnames(aggregated_expression_target$sample) %>% data.frame() %>% rownames_to_column("target") %>% tbl_df()
aggregated_expression_matrix_target = aggregated_expression_df_target %>% select(-target) %>% as.matrix() %>% magrittr::set_rownames(aggregated_expression_df_target$target)
order_samples <- intersect(colnames(aggregated_expression_matrix_target),order_samples)
order_targets <- intersect(rownames(aggregated_expression_matrix_target),order_targets)
vis_target_tumor_expression_scaled = aggregated_expression_matrix_target %>% t() %>% scale_quantile() %>% .[order_samples,order_targets]
p_target_tumor_scaled_expression = vis_target_tumor_expression_scaled  %>% make_threecolor_heatmap_ggplot("Sample","Target", low_color = color[1],mid_color = color[50], mid = 0.5, high_color = color[100], legend_position = "top", x_axis_position = "top" , legend_title = "Scaled expression\n(averaged over\nsingle cells)") + theme(axis.text.x = element_text(face = "italic"))


figures_without_legend = plot_grid(
  p_ligand_pearson + theme(legend.position = "none", axis.ticks = element_blank(),axis.text=element_text(size=6)) + theme(axis.title.x = element_text()),
  p_ligand_tumor_expression + theme(legend.position = "none", axis.ticks = element_blank(),axis.text=element_text(size=6)) + theme(axis.title.x = element_text()) + ylab(""),
  p_ligand_target_network + theme(legend.position = "none", axis.ticks = element_blank(),axis.text=element_text(size=6)) + ylab(""),
  NULL,
  NULL,
  p_target_tumor_scaled_expression + theme(legend.position = "none", axis.ticks = element_blank(),axis.text=element_text(size=6)) + xlab(""),
  align = "hv",
  nrow = 2,
  rel_widths = c(ncol(vis_ligand_pearson)+ 10.5, ncol(vis_ligand_tumor_expression)+10.5, ncol(vis_ligand_target)-12),
  rel_heights = c(1.5, 1))

  legends = plot_grid(
    as_ggplot(get_legend(p_ligand_pearson)),
    as_ggplot(get_legend(p_ligand_tumor_expression)),
    as_ggplot(get_legend(p_ligand_target_network)),
    as_ggplot(get_legend(p_target_tumor_scaled_expression)),
    nrow = 2,
    align = "h")

figure <- plot_grid(figures_without_legend,
          legends,
          rel_heights = c(10,3), nrow = 2, align = "hv")


ligands_all = order_ligands
targets_all = order_targets
active_signaling_network = get_ligand_signaling_path(ligand_tf_matrix = ligand_tf_matrix, ligands_all = ligands_all, targets_all = targets_all, weighted_networks = weighted_networks)
active_signaling_network_min_max = active_signaling_network
active_signaling_network_min_max$sig = active_signaling_network_min_max$sig %>% mutate(weight = ((weight-min(weight))/(max(weight)-min(weight))) + 0.75)
active_signaling_network_min_max$gr = active_signaling_network_min_max$gr %>% mutate(weight = ((weight-min(weight))/(max(weight)-min(weight))) + 0.75)

graph_min_max = diagrammer_format_signaling_graph(signaling_graph_list = active_signaling_network_min_max, ligands_all = ligands_all, targets_all = targets_all, sig_color = "indianred", gr_color = "steelblue")

graph <- DiagrammeR::render_graph(graph_min_max, layout = "tree", as_svg=TRUE)

seurat <- AddModuleScore(object = sender_obj, features = geneset_oi, ctrl = 5, name = "subtype_score")
subfig <- FeaturePlot(seurat,reduction="umapcelltype",features="subtype_score1") + ggtitle(paste0(pathway, " Score")) + theme(plot.title=element_text(face="bold")) + xlab("UMAP-1") + ylab("UMAP-2")

top <- plot_grid(graph,subfig,nrow=1,align="h")

allfigure <- plot_grid(top,figure,nrow=2,ncol=1)

ggsave(png, allfigure, width=14, height=12)
