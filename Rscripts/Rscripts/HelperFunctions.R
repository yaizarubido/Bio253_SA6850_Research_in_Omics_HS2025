# 2025-10-22
# some useful functions for data analysis in R of MtX data



# Function to check if genes overlap with hypo-methylated regions
# this function does not care about the direction or coding strand
find_overlapping_genes <- function(genes_df, hypometh_df) {
    # Create ranges for genes
    gene_ranges <- IRanges(start = genes_df$start,
                           end = genes_df$end,
                           names = genes_df$Acc2)

    # Create ranges for hypo-methylated regions
    # Assuming start column is the beginning of the region
    hypometh_ranges <- IRanges(start = hypometh_df$start,
                               end = hypometh_df$hypoMethRegion_end)

    # Find overlaps
    overlaps <- findOverlaps(gene_ranges, hypometh_ranges)

    # Get the names of overlapping genes
    overlapping_genes <- names(gene_ranges)[queryHits(overlaps)]

    # Return unique gene names
    return(unique(overlapping_genes))
}


# Function to check if a genomic feature overlaps with a bin
# A feature overlaps if its start OR end position falls within the bin range
find_overlapping_features <- function(methylation_bins, gff_data) {

    overlapping_features <- data.frame()

    for (i in 1:nrow(methylation_bins)) {
        bin_start <- methylation_bins$bin[i]
        bin_end <- methylation_bins$endBin[i]
        bin_strand <- methylation_bins$strand[i]

        # Convert strand notation: "-" becomes "-", anything else becomes "+"
        gff_strand_match <- if(bin_strand == "-") "-" else "+"

        # Find features where:
        # 1. The strand matches
        # 2. The feature's start OR end position falls within the bin range
        matches <- gff_data %>%
            filter(
                strand == gff_strand_match,
                (start >= bin_start & start <= bin_end) |
                    (end >= bin_start & end <= bin_end)
            ) %>%
            mutate(
                methylation_bin_start = bin_start,
                methylation_bin_end = bin_end,
                methylation_strand = bin_strand,
                PASN_count = methylation_bins$PASN[i],
                TSB_count = methylation_bins$TSB[i],
            )

        overlapping_features <- rbind(overlapping_features, matches)
    }

    return(overlapping_features)
}


