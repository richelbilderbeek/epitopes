library(dplyr)
library(epitopes)

# Load epitopes (change folder if needed)
epitopes <- readRDS("./results/epitopes.rds")

# sort(table(epitopes$sourceOrg_id), decreasing = TRUE)[1:5]
# "353153" - T. cruzi (Chagas' disease) ------ 173,689 records
# "9609" --- H. sapiens ---------------------- 52,085 records (self-epitopes?)
# "6282" --- O. volvulus (river blindness) --- 7642 records


# 1) epitope ID 6282: O. volvulus (river blindness)

Org_id <- "6282"
eps_Ov <- epitopes %>%
  filter(sourceOrg_id == Org_id)

uids     <- unique(eps_Ov$molecule_id)
uids     <- uids[!is.na(uids)]
proteins <- get_proteins(uids, save_folder = "./results/Ovolvulus")
proteins <- as_tibble(proteins$proteins)
#proteins <- readRDS("./results/Ovolvulus/proteins.rds")


df_Ov  <- assemble_windowed_dataframe(eps_Ov, proteins,
                                      save_folder  = "./results/Ovolvulus",
                                      min_epit     = 8,
                                      max_epit     = 20,
                                      only_exact   = FALSE,
                                      window_exp   = 0,
                                      step_size    = 1,
                                      ncpus        = parallel::detectCores() - 1)

df_Ov <- df_Ov$windows_df

# Find inconsistent entries (same window_seq with different Class):
to_rm <- df_Ov %>%
  group_by(window_seq) %>%
  summarise(NClasses = length(unique(Class)),
            .groups = "drop") %>%
  filter(NClasses == 2) %>%
  select(window_seq)


# Remove inconsistent entries as well as duplicated windows
df_Ov <- df_Ov %>%
  filter(!window_seq %in% to_rm$window_seq) %>%
  select(-Type) %>%
  group_by(window_seq) %>%
  summarise_all(first)

saveRDS(df_Ov, "./results/Ovolvulus/df_windowed_unique.rds")
