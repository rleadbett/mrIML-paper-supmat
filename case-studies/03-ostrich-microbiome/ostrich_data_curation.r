library(dplyr)
library(tibble)
library(readr)
library(stringr)
library(here)
library(mrIML)

#----------------------------------------------------
# Prepare ASV table
#----------------------------------------------------
# Read in the otu table
otu_table_ostrich <- read_tsv(
  here(
    "case-studies",
    "03-ostrich-microbiome",
    "data",
    "raw",
    "gut_otu_table.tsv"
  )
)
otu_table_ostrich <- otu_table_ostrich #[1:10]
saveRDS(
  otu_table_ostrich,
  file = here(
    "case-studies",
    "03-ostrich-microbiome",
    "data",
    "raw",
    "asv_data.rds"
  )
)
taxa_table_ostrich <- read_tsv(
  here(
    "case-studies",
    "03-ostrich-microbiome",
    "data",
    "raw",
    "gut_tax_table.tsv"
  )
)
saveRDS(
  taxa_table_ostrich,
  file = here(
    "case-studies",
    "03-ostrich-microbiome",
    "data",
    "raw",
    "asv_taxa_table.rds"
  )
)
# Check data frames match
identical_column <- identical(taxa_table_ostrich$ASV, otu_table_ostrich$ASV)
if (identical_column) {
  print("The 'ColumnName' matches identically.")
} else {
  print("The 'ColumnName' does not match identically.")
}

#ASV column not needed
taxa_table_ostrich$ASV <- NULL
taxa_table_ostrich$Kingdom <- NULL
taxa_table_ostrich$Class <- NULL
taxa_table_ostrich$Phylum <- NULL

# Function to create meaningful ASV names from taxonomic information
create_name <- function(row) {
  non_na_values <- na.omit(row)
  unique_values <- unique(non_na_values)

  if (length(unique_values) > 3) {
    unique_values <- unique_values[1:2]
  }

  # Extract taxonomic levels
  family <- ifelse("Family" %in% names(row), row["Family"], NA)
  genus <- ifelse("Genus" %in% names(row), row["Genus"], NA)
  species <- ifelse("Species" %in% names(row), row["Species"], NA)

  # Combine into meaningful name
  name_parts <- c(family, genus, species, unique_values)
  name_parts <- name_parts[!is.na(name_parts)]
  return(paste(name_parts, collapse = "_"))
}

# Function to handle empty names
replace_empty_names <- function(names_vector) {
  # Find indices where names are empty
  empty_indices <- names_vector == ""

  # Count occurrences of empty names
  empty_counts <- cumsum(empty_indices)

  # Replace empty names with "Taxon" followed by count
  names_vector[empty_indices] <- paste0("Taxon ", empty_counts[empty_indices])

  return(names_vector)
}

# Function to make names unique
make_names_unique <- function(names_vector) {
  # Initialize a counter to keep track of occurrences of each name
  name_counts <- list()

  # Initialize a result vector to store the modified names
  result_names <- character(length(names_vector))

  for (i in seq_along(names_vector)) {
    name <- names_vector[i]

    # Check if the name is empty or has occurred before
    if (name == "" || is.null(name_counts[[name]]) || name_counts[[name]] > 0) {
      # Increment the count for this name
      count <- ifelse(is.null(name_counts[[name]]), 1, name_counts[[name]] + 1)
      name_counts[[name]] <- count

      # Append a number to make the name unique
      result_names[i] <- paste0(name, ".", count)
    } else {
      result_names[i] <- name
      # First occurrence of this name, initialize count to 1
      name_counts[[name]] <- 1
    }
  }

  return(result_names)
}

# Create the 'Name' column
taxa_table_ostrich <- taxa_table_ostrich %>%
  mutate(
    ASV = apply(taxa_table_ostrich, 1, create_name) |>
      replace_empty_names() |>
      make_names_unique()
  )
otu_table_ostrich$ASV <- taxa_table_ostrich$ASV

# Pivote
otu_table_ostrich <- otu_table_ostrich %>%
  column_to_rownames(var = 'ASV')
otu_table_ostrich <- as.data.frame(t(otu_table_ostrich))

# Convert to presence/absence
pa_ASV_table <- otu_table_ostrich %>%
  mutate_all(~ ifelse(. > 0, 1, .))

# Remove rare and common ASVs
Y <- filterRareCommon(pa_ASV_table, lower = 0.2, higher = 0.8) %>%
  dplyr::select(sort(names(.))) %>%
  as.data.frame()
rownames(Y) <- rownames(pa_ASV_table) # Preserve rownames

## Find and fix duplicated column names (should be redundant)
duplicated_cols <- duplicated(colnames(Y))
colnames(Y) <- make.names(colnames(Y), unique = TRUE)
Y <- rename_all(Y, ~ make.names(str_remove_all(., "`")))

# Shorten the part of name to the left of the underscore to 4 characters
abbreviate_names <- function(names_vector) {
  abbreviated_names <- sapply(strsplit(names_vector, "_"), function(parts) {
    first_part <- substr(parts[1], 1, 4)
    second_part <- parts[2]
    return(paste(first_part, second_part, sep = "_"))
  })
  return(abbreviated_names)
}
colnames(Y) <- abbreviate_names(colnames(Y))

# Create co-occurance df for mrIML
X1 <- Y

#----------------------------------------------------
# Prepare metadata
#----------------------------------------------------
# Read in metadata
metadata_ostrich <- read_tsv(
  here(
    "case-studies",
    "03-ostrich-microbiome",
    "data",
    "raw",
    "gut_metadata.tsv"
  )
)
# Remove redundant columns
metadata_ostrich$AgeDied <- NULL # same as age
metadata_ostrich$Method <- NULL # all the same
metadata_ostrich$AgeWeek <- NULL # not needed?
metadata_ostrich$WeightPM <- NULL # highly correlated with weight
metadata_ostrich$TypeShort <- NULL # highly correlated with weigh
metadata_ostrich$TimePM <- NULL # time since post mortem - not interesting
metadata_ostrich$WeightHatch <- NULL # not of interest
# Match with ASVs
identical_column <- identical(row.names(Y), metadata_ostrich$X.SampleID)
if (identical_column) {
  print("The 'ColumnName' matches identically.")
} else {
  print("The 'ColumnName' does not match identically.")
}

# Permute missing values
Xmissing <- metadata_ostrich %>%
  column_to_rownames(var = "X.SampleID")
X_fact <- Xmissing %>%
  select(where(is.character)) %>%
  mutate(across(where(is.character), factor))
Xmissing_num <- Xmissing %>%
  select(where(is.numeric))
X_ <- missForest::missForest(Xmissing_num)
X <- cbind(X_$ximp, X_fact)

# Remove ID
X$ID <- NULL

# Separate by sample type
# Colon, caecum, etc.
data_combined <- cbind(Y, X)
names(data_combined) <- make.names(names(data_combined), unique = TRUE)

#colon
data_combined_colon <- data_combined %>%
  filter(Type == 'Colon')
Y_colon <- data_combined_colon %>%
  select(contains("_"))
Y_colonF <- filterRareCommon(Y_colon, lower = 0.4, higher = 0.7) %>%
  dplyr::select(sort(names(.)))
X_colon <- data_combined_colon %>%
  select(-contains("_"))
X_colon$Type <- NULL

#caecum
data_combined_caecum <- data_combined %>%
  filter(Type == 'Caecum')
Y_caecum <- data_combined_caecum %>%
  select(contains("_"))
Y_caecumF <- filterRareCommon(Y_caecum, lower = 0.3, higher = 0.7) %>%
  dplyr::select(sort(names(.)))
X_caecum <- data_combined_caecum %>%
  select(-contains("_"))
X_caecum$Type <- NULL

#Ileum
data_combined_ileum <- data_combined %>%
  filter(Type == 'Ileum')
Y_ileum <- data_combined_ileum %>%
  select(contains("_"))
Y_ileumF <- filterRareCommon(Y_ileum, lower = 0.2, higher = 0.8) %>%
  dplyr::select(sort(names(.)))
X_ileum <- data_combined_ileum %>%
  select(-contains("_"))
X_ileum$Type <- NULL
