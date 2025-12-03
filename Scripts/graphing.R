library(tidyverse)
library(RColorBrewer)
library(ggrepel)
library(biomaRt)

args <- commandArgs(trailingOnly = TRUE)
df <- read.csv(args[1], row.names = 1)
df$geneid <- row.names(df)


ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
ensembl_ids <- df$geneid
cleaned_ids <- sub("\\..*", "", ensembl_ids)
attributes_list <- listAttributes(mart = ensembl)
gene_map <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
                  filters = "ensembl_gene_id",
                  values = cleaned_ids,
                  mart = ensembl)

df$genename <- gene_map$hgnc_symbol[match(cleaned_ids, gene_map$ensembl_gene_id)]
write.csv(df, "Outputs/counts/deseq2results_genenames.csv", row.names = TRUE)

theme_set(theme_classic(base_size = 20) +
            theme(
              axis.title.y = element_text(face = "bold", margin = margin(0,20,0,0), size = rel(1.1), color = 'black'),
              axis.title.x = element_text(hjust = 0.5, face = "bold", margin = margin(20,0,0,0), size = rel(1.1), color = 'black'),
              plot.title = element_text(hjust = 0.5)
            ))

volcano_plot <- function() { 

    df$expression <- "NO"
    df$expression[df$log2FoldChange > 0.6 & df$pvalue < 0.05] <- "UP"
    df$expression[df$log2FoldChange < -0.6 & df$pvalue < 0.05] <- "DOWN"

    Top_Hits = head(arrange(df,pvalue),10) 
    # Add column label, containing the gene name for the top hits or nothing for all others
    df$label = if_else(df$genename %in% Top_Hits$genename,  
                    df$genename, NA)


    output_plot <- ggplot(data = df, aes(x = log2FoldChange, y = -log10(pvalue), col = expression, label = label)) +
    geom_vline(xintercept = c(-0.6, 0.6), col = "gray", linetype = 'dashed') +
    geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') + 
    geom_point(size=1) +

    scale_color_manual(values = c("#EF4538", "grey", "#048ACC"), # to set the colours of our variable<br /><br /><br />
                     labels = c("Downregulated", "Not significant", "Upregulated")) +
    
    coord_cartesian(ylim = c(0, 60), xlim = c(-10, 10)) +
    geom_text_repel(max.overlaps = Inf) +
    ggtitle("Deseq2 Volcano Plot")

        # Display the plot (important for ggsave to capture it by default)
    #print(output_plot) 

    # Save the plot as a PNG file in the current working directory
    ggsave("Outputs/counts/volcano_plot.png", width = 12, height = 10) 
}

volcano_plot()