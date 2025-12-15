library(phangorn)
library(Biostrings)
library(seqinr)
library(msa)
library(ggplot2)
library(ggtree)
library(ggtreeExtra)
library(ggnewscale)
library(paletteer)

raw_file   <- "amblyomma_a3_evasins.fasta"
# alig <- msa(raw_file, method="ClustalOmega")  # Alignment contains substantial N-terminal gapping - trim required
align_file <- "A3eva_trim_align_clustalO_new.fasta"
alig.trim  <- readAAMultipleAlignment(align_file)
aln_matrix <- as.matrix(alig.trim)
aln_seqinr <- msaConvert(alig.trim, type = "seqinr::alignment")

# Distance: D = sqrt(1 - identity)
d.trim  <- dist.alignment(aln_seqinr, "identity")
nj_tree <- nj(d.trim)
nj_tree <- midpoint(nj_tree)

ident_mat  <- 1 - (as.matrix(d.trim)^2)  # convert distance back to identity (0-1 scale) for heatmap
df_heatmap <- as.data.frame(ident_mat)

# D. Extract Metadata (Species)
raw_headers <- names(readAAStringSet(raw_file))
meta_df <- data.frame(
  label = gsub(" \\(.*\\)", "", raw_headers),
  Species = gsub(".*\\((.*)\\).*", "\\1", raw_headers),
  stringsAsFactors = FALSE)

p <- ggtree(nj_tree, layout = "rectangular", ladderize = TRUE, 
            size = 0.2, color = "black", linewidth = 1) +
  geom_tiplab(aes(label = ""), align = TRUE, linetype = "dotted", 
              linesize = 0.5, offset = 0.05, colour = "grey40") +
  geom_fruit(data = meta_df, geom = geom_tile, 
             mapping = aes(y = label, fill = Species),
             width = 0.04, offset = 0.175, pwidth = 0.05, 
             colour = "black", linewidth = 1) + 
  scale_fill_paletteer_d("PNWColors::Sailboat", name = "Species") +
  ggnewscale::new_scale_fill() + 
  geom_fruit(data = meta_df, geom = geom_text, 
             mapping = aes(y = label, label = label),
             hjust = 0, size = 2.5, offset = 0.06, pwidth = 0.15) +
  geom_treescale(x = 0, y = 0, fontsize = 2.5, linesize = 0.5, width = 0.1, 
                 offset = 0.5, color = "grey40")


p_main <- gheatmap(p, df_heatmap, offset = 0.4, width = 4, 
                   colnames = TRUE, color = "grey80",
                   colnames_position = "bottom", colnames_angle = 45,
                   hjust = 1) +
  scale_fill_distiller(palette = "YlOrRd", direction = 1, name = "% Identity",
                       limits = c(0.25, 1), breaks = seq(0.25, 1, 0.25),
                       labels = scales::percent) +
  labs(title = "Main Text: Phylogenetic Divergence & Identity") +
  theme(
    legend.position = "right", legend.box = "vertical",
    legend.spacing.y = unit(0, "pt"),
    legend.title = element_text(size = 9, face = "bold"),
    legend.text = element_text(size = 8)
  )

p_supp <- ggtree(nj_tree, layout = "fan", open.angle = 15, 
                 size = 0.3, color = "grey30", linewidth = 1) +
  geom_fruit(data = meta_df, geom = geom_tile, 
             mapping = aes(y = label, fill = Species),
             width = 0.05, offset = 0.04, colour = "black", linewidth = 1) +
  scale_fill_paletteer_d("PNWColors::Sailboat", name = "Species") +
  geom_tiplab(aes(angle = angle), color = "black", size = 4, 
              offset = 0.08, fontface = "bold") +
  geom_highlight(node = which(nj_tree$tip.label == "EVA-ATL1001"), 
                 fill = "#C464FF", alpha = 0.2, extendto = 0.75) +
  geom_treescale(x = 0, y = 0, width = 0.1, fontsize = 3, 
                 linesize = 0.5, offset = 1) +
  labs(title = "Figure S1: Circular Phylogeny") +
  theme(
    legend.position = "right", legend.box = "vertical",
    legend.spacing.y = unit(0, "pt"),
    legend.title = element_text(size = 9, face = "bold"),
    legend.text = element_text(size = 8)
  )

ggsave("Figure_Main.pdf", p_main, width = 11, height = 6)
ggsave("Figure_Main.png", p_main, width = 11, height = 6, dpi = 300)
ggsave("Figure_Supp.pdf", p_supp, width = 10, height = 10)
ggsave("Figure_Supp.png", p_supp, width = 10, height = 10, dpi = 300)
