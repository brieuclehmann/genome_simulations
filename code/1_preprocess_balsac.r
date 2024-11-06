library(readr)
source("misc/pedigree_tools.R")
set.seed(42)

pedigree_df <- read_delim("data/ascendance.csv", delim = ';')
colnames(pedigree_df) <- c('proband_region', 'proband', 'ind', 'sex', 
                           'mother', 'father', 'decade', 'ParentsMarriageLocation', 
                           'ParentsMarriageLocationID', 'Origin', 'OriginID', 'Sosa', 'Generation')

pedigree_df$mother <- as.integer(pedigree_df$mother)
pedigree_df$father <- as.integer(pedigree_df$father)

proband_id <- unique(pedigree_df$proband)
max_generation <- maximum_genealogical_depth(pedigree_df, proband_id)

full_df <- pedigree_df %>%
    left_join(max_generation, by = "ind")
write_csv(full_df, "data/balsac_pedigree.csv")

######################
### Find relatives ###
######################
list_of_probands <- full_df %>% 
  filter(!ind %in% c(full_df$mother, full_df$father)) %>% 
  pull(ind) 
relatives <- identify_related_probands(full_df, list_of_probands, k=3) 
write_csv(relatives, "data/balsac_relatives.csv")

proband_df <- pedigree_df %>%
    distinct(proband_region, proband)

relative_df <- relatives %>%
    left_join(proband_df, by = c("proband1" = "proband")) %>%
    rename(proband_region1 = proband_region) %>%
    left_join(proband_df, by = c("proband2" = "proband")) %>%
    rename(proband_region2 = proband_region) %>%
    select(-ancestor) %>%
    distinct()


##########################
### Select individuals ###
##########################

# Sample probands so that there are at least 5 {siblings, first cousins, second cousins} from each region
relative_sample <- relative_df %>%
    filter(proband_region1 == proband_region2) %>%
    filter(generation <= 3) %>%
    group_by(proband_region1, generation) %>%
    slice_sample(n = 5) %>%
    ungroup() %>%
    select(proband1, proband2) %>%
    unlist() %>%
    unique()

selected_relatives <- full_df %>%
    filter(ind %in% relative_sample) %>%
    left_join(proband_df, by = c("ind" = "proband"))

num_remain <- selected_relatives %>%
    group_by(proband_region) %>%
    summarise(n_sample = 50 - n())

# Continue subsampling to fill 100 individuals from each region
selected_df <- full_df %>%
    left_join(proband_df, by = c("ind" = "proband")) %>%
    filter(generation == 0 & !is.na(proband_region)) %>%
    filter(! ind %in% relative_sample) %>%
    nest(data = -proband_region) %>%
    left_join(num_remain, by = c("proband_region")) %>%
    mutate(Sample = map2(data, n_sample, sample_n)) %>%
    unnest(Sample) %>%
    select(-c(data, n_sample)) %>%
    bind_rows(selected_relatives) %>%
    arrange(proband_region)

# Check (should be 50 distinct probands from each region)
selected_df %>%
    distinct(proband_region, ind) %>%
    group_by(proband_region) %>%
    summarise(n = n())

write_csv(selected_df, "data/balsac_subsample.csv")