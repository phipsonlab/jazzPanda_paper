library(matrixStats)
library(reshape2)
library(jazzPanda)
library(tidyr)
library(dplyr)
library(peakRAM)

args <- commandArgs(trailingOnly = TRUE)
task_id <- as.integer(args[1])
repeat_times <- as.integer(args[2])
tile_lens <- as.numeric(strsplit(args[3], ",")[[1]])

# total number of tiles
num_tiles <- length(tile_lens)

# Calculate the tile_len index for the current task
tile_len_index <- ((task_id - 1) %% num_tiles) + 1

# Calculate the repetition number for the current task
repetition_number <- floor((task_id - 1) / num_tiles) + 1

# Extract the current tile length
curr_size <- tile_lens[tile_len_index]
cat("Repeat =", repetition_number, ", curr_size =", curr_size, "\n")

simulate_transcripts <- function(n_transcripts,
                                 n_genes = 500, 
                                 x_range = c(0, 10000), 
                                 y_range = c(0, 10000)) {
    data.frame(
        x = runif(n_transcripts, min = x_range[1], max = x_range[2]),
        y = runif(n_transcripts, min = y_range[1], max = y_range[2]),
        feature_name = sample(paste0("gene", seq_len(n_genes)), n_transcripts, 
                              replace = TRUE)
    )
}


simulate_cluster <- function(n_cells,
                             n_cluster,
                             x_range = c(0, 10000), 
                             y_range = c(0, 10000)) {
    data.frame(
        x = runif(n_cells, min = x_range[1], max = x_range[2]),
        y = runif(n_cells, min = y_range[1], max = y_range[2]),
        cluster = sample(paste0("cluster", seq_len(n_cluster)), n_cells, 
                              replace = TRUE),
        sample = "sample1"
    )
}

seed_number = sample(1:878, size=1)
set.seed(seed_number)
trans_df <- simulate_transcripts(n_transcripts = 1e8, n_genes = 500)
cluster_df <- simulate_cluster(n_cells = 1e5, n_cluster = 10)

transcript_df_size = as.numeric(object.size(trans_df)) / (1024^2 )

test_genes <- unique(trans_df$feature_name)

grid_length = curr_size

# create sv based on n_genes only
usage_sv_df = peakRAM({res = get_vectors(x = list(sample1 = trans_df),
                     cluster_info = cluster_df,
                     sample_names = "sample1",
                     bin_type = "hexagon",
                     bin_param=c(grid_length),
                     test_genes = test_genes,
                     n_cores = 1)})

usage_glm_df = peakRAM({jazzPanda_res_lst = lasso_markers(gene_mt=res$gene_mt,
                      cluster_mt = res$cluster_mt,
                      sample_names=c("sample1"),
                      keep_positive=TRUE,
                      n_fold = 10)})

sv_res_size= as.numeric(object.size(res)) / (1024^2 )

results_df <- data.frame(
    task_id = task_id,
    rd = repetition_number,
    grid_length = grid_length,
    sv_Elapsed_Time_sec = usage_sv_df$Elapsed_Time_sec,
    sv_Peak_RAM_Used_MiB =usage_sv_df$Peak_RAM_Used_MiB,
    glm_Elapsed_Time_sec = usage_glm_df$Elapsed_Time_sec,
    glm_Peak_RAM_Used_MiB =usage_glm_df$Peak_RAM_Used_MiB,
    sv_lst_size=sv_res_size,
    transcript_size=transcript_df_size,
    transcript_n = nrow(trans_df),
    seed_number=seed_number
)


args_all   <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", args_all, value = TRUE)
script_dir <- if (length(script_arg)) {
    dirname(normalizePath(sub("^--file=", "", script_arg[1])))
} else {
    normalizePath(getwd())  # interactive fallback
}



## Output folder name
output_dir_nm <- "ntiles_hexbin_result"

## Path to /scripts/main/cosmx_hlc_simulation_result
out_dir <- file.path(script_dir, output_dir_nm)

setwd(out_dir)

# Output file name uniquely identified by SLURM_ARRAY_TASK_ID
output_file_name <- sprintf("complexity_hexbins_ntiles%d_id%d.csv",grid_length, task_id)
# Write the results to a CSV file
write.csv(results_df, output_file_name, row.names = FALSE)
