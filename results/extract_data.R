library(dplyr)
library(epitopes)

# Load epitopes (change folder if needed)
epitopes <- readRDS("./results/epitopes.rds")

# sort(table(epitopes$sourceOrg_id), decreasing = TRUE)[1:5]
# "353153" - T. cruzi (Chagas' disease) ------ 173,689 records
# "9609" --- H. sapiens ---------------------- 52,085 records (self-epitopes?)
# "6282" --- O. volvulus (river blindness) --- 7642 records


# 1) epitope ID 353153: Trypanosoma cruzi (Chagas' disease)
# 173k+ records

Org_id <- names(sort(table(epitopes$sourceOrg_id), decreasing = TRUE))[1]
eps_Tc <- epitopes %>%
  filter(sourceOrg_id == Org_id)

# uids     <- unique(eps_Tc$molecule_id)
# uids     <- uids[!is.na(uids)]
# proteins <- get_proteins(uids, save_folder = "./results/Tcruzi")
# proteins <- as_tibble(proteins$proteins)
proteins <- readRDS("./results/Tcruzi/proteins.rds")


df_Tc  <- assemble_windowed_dataframe(eps_Tc, proteins,
                                      save_folder  = "./results/Tcruzi",
                                      min_epit     = 8,
                                      max_epit     = 20,
                                      only_exact   = FALSE,
                                      window_exp   = 0,
                                      step_size    = 1,
                                      ncpus        = parallel::detectCores() - 2)

df_Tc <- df_Tc$windows_df

# Find inconsistent entries (same window_seq with different Class):
to_rm <- df_Tc %>%
  group_by(window_seq) %>%
  summarise(NClasses = length(unique(Class)),
            .groups = "drop") %>%
  ungroup() %>%
  filter(NClasses == 2) %>%
  select(window_seq)


# Remove inconsistent entries as well as duplicated windows
df_Tc2 <- df_Tc %>%
  filter(!window_seq %in% to_rm$window_seq) %>%
  select(-Type) %>%
  group_by(window_seq) %>%
  summarise_all(first)

# saveRDS(df_Tc2, "./results/Tcruzi/df_windowed_unique.rds")


# ==== O. volvulus
# The data below was obtained using the same process used above
df_Ov  <- readRDS("./results/Ovolvulus/df_windowed.rds")

# Find inconsistent entries (same window_seq with different Class):
to_rm <- df_Ov %>%
  group_by(window_seq) %>%
  summarise(NClasses = length(unique(Class)),
            .groups = "drop") %>%
  ungroup() %>%
  filter(NClasses == 2) %>%
  select(window_seq)


# Remove inconsistent entries as well as duplicated windows
df_Ov2 <- df_Ov %>%
  filter(!window_seq %in% to_rm$window_seq) %>%
  select(-Type) %>%
  group_by(window_seq) %>%
  summarise_all(first)

# saveRDS(df_Ov2, "./results/Ovolvulus/df_windowed_unique.rds")




