library(readr)
source("misc/pedigree_tools.R")
set.seed(42)

################################
## Prepare pedigree dataframe ##
################################
ped_file <- "data/balsac_pedigree.csv"
if (file.exists(ped_file)) {
    full_df <- read_csv(ped_file)
} else {
    pedigree_df <- read_delim("data/ascendance.csv", delim = ';')
    colnames(pedigree_df) <- c('proband_region', 'proband', 'ind', 'sex', 
                            'mother', 'father', 'decade', 'ParentsMarriageLocation', 
                            'ParentsMarriageLocationID', 'Origin', 'OriginID', 'Sosa', 'Generation')

    pedigree_df$mother <- as.integer(pedigree_df$mother)
    pedigree_df$father <- as.integer(pedigree_df$father)

    pedigree_df <- distinct(pedigree_df) # Remove small number of duplicate rows
    proband_id <- unique(pedigree_df$proband)
    max_generation <- maximum_genealogical_depth(pedigree_df, proband_id)

    full_df <- pedigree_df %>%
        left_join(max_generation, by = "ind")
    write_csv(full_df, "data/balsac_pedigree.csv")
}

region_df <- full_df %>%
    distinct(proband_region, proband)

######################
### Find relatives ###
######################

relatives_file <- "data/balsac_relatives.csv"
if (file.exists(relatives_file)) {
    relative_df <- read_csv(relatives_file)
} else {
    list_of_probands <- full_df %>% 
        filter(!ind %in% c(full_df$mother, full_df$father)) %>% 
        pull(ind) 

    relatives <- identify_related_probands(full_df, list_of_probands, k=5) 

    relatives <- relatives %>%
        select(-ancestor) %>%
        group_by(proband1, proband2) %>%
        filter(generation == min(generation)) %>%
        distinct()

    relative_df <- relatives %>%
        left_join(proband_df, by = c("proband1" = "proband")) %>%
        rename(proband_region1 = proband_region) %>%
        left_join(proband_df, by = c("proband2" = "proband")) %>%
        rename(proband_region2 = proband_region) %>%
        distinct()
    
    write_csv(relative_df, relatives_file)
}

##########################
### Select individuals ###
##########################

founder_df <- read_csv("data/balsac_founder_depths.csv")

# Sample probands so that there are at least 5 {siblings, first cousins, second cousins} from each region
relative_sample <- relative_df %>%
    left_join(founder_df, by = c("proband1" = "proband")) %>%
    rename(average_depth1 = average_depth) %>%
    left_join(founder_df, by = c("proband2" = "proband")) %>%
    rename(average_depth2 = average_depth) %>%
    filter(average_depth1 >= 3 & average_depth2 >= 3) %>%
    filter(proband_region1 == proband_region2) %>%
    filter(generation <= 4) %>%
    group_by(proband_region1, generation) %>%
    slice_sample(n = 5) %>%
    ungroup() %>%
    select(proband1, proband2) %>%
    unlist() %>%
    unique()

selected_relatives <- full_df %>%
    filter(ind %in% relative_sample) %>%
    left_join(founder_df, by = join_by(proband)) 

num_remain <- selected_relatives %>%
    group_by(proband_region) %>%
    summarise(n_sample = 50 - n())

# Continue subsampling to fill 100 individuals from each region
selected_df <- full_df %>%
    left_join(founder_df, by = join_by(proband)) %>%
    filter(average_depth >= 3) %>%
#    left_join(proband_df, by = c("ind" = "proband")) %>%
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

selected_df$noid <- 1:nrow(selected_df) # For sharing
write_csv(selected_df, "data/balsac_subsample.csv")

# Save metadata
meta_df <- selected_df %>%
    select(proband, noid, proband_region, average_depth)
write_csv(meta_df, "data/balsac_subsample_meta.csv")
write_csv(subset(meta_df, select = -proband), "data/balsac_subsample_meta_noid.csv")

# Remove ID information
id_map <- meta_df %>%
    select(proband, noid)

relatives_selected_df <- relative_df %>%
    filter((proband1 %in% selected_df$ind) & (proband2 %in% selected_df$ind)) %>%
    left_join(id_map, by = c("proband1" = "proband")) %>%
    select(-proband1) %>%
    rename(proband1 = noid) %>%
    left_join(id_map, by = c("proband2" = "proband")) %>%
    select(-proband2) %>%
    rename(proband2 = noid)
write_csv(relatives_selected_df, "data/balsac_subsample_relatives_noid.csv")
