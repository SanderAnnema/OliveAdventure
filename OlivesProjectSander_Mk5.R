#### Introduction: ####
# Master project by Sander Annema. 06/04/2023
# In this project, proteomics data from probiotic bacteria isolated from olives is processed. 
# These bacteria are either resistant, intermediate, or sensitive to bile.
# The expression data of both proteins and peptides was collected.
# The goal of this project is to determine the effects of protein expression on the exposure of the bacterial strains by bile.
# This will be achieved by performing a 2-way ANOVA, then generating a heat-map of the proteins.

#### Preparation ####
# To begin, the session will be cleaned by removing all variables.
rm(list=ls())

## Loading the required libraries
library(ggplot2)
library(reshape2)
library(tidyr)
library(tidyverse)
library(RColorBrewer)
library(broom)
library(multcompView)
library(multcomp)
library(readr)
library(dplyr)
library(ComplexHeatmap)
library(pheatmap)

# Setting the working directory
setwd("C:/Users/Skyar/OneDrive/Documenten/school/Master_BiMoS/2nd_research_project/Project_Olives_New")

# Terms used for streamlining the code
DataFolder = "C:/Users/Skyar/OneDrive/Documenten/school/Master_BiMoS/2nd_research_project/Project_Olives_New/data/"
ProteinFile = "proteins.csv"
PeptideFile = "peptide.csv"
IdxIntensCol = 9:44

#### Loading and formatting data ####
ProtTab_Full = read.csv(paste(DataFolder,ProteinFile, sep = ""), header =  TRUE)
PepTab_Full = read.csv(paste(DataFolder,PeptideFile, sep = ""), header =  TRUE)

# The header names are too long, so 'intensity' will be removed from them for brevity
columnNames = gsub(".Intensity","", colnames(ProtTab_Full) [IdxIntensCol])
colnames(ProtTab_Full) [IdxIntensCol] = columnNames

# The data contains triplicate samples measured in duplicate of 2 factors (control and treated) and 3 groups (resistent, intermediair, sensitive).
# To make handling specific data easier, for example when calculating averages, well-defined label objects are created.
labels = strsplit(columnNames, "_") # The column names are separated based on "_", so a table is formed with the following structure of columns: factor, group, replicate, measurement duplicate.
labelTreatment = unlist(lapply(labels, function(x) x[1])) # The first column, treatment, is isolated
labelStrain = unlist(lapply(labels, function(x) x[2])) # The second column, strain, is isolated
labelReplicate = unlist(lapply(labels, function(x) x[3])) # The third column, replicate, is isolated
labelTreatmentStrainRep = paste(labelTreatment, labelStrain, labelReplicate, sep = "_") # The three columns are reassembled into new labels containing only vital information
labelTreatmentStrainUnique = unique(labelTreatmentStrainRep) # The duplicates are removed

#### Median normalization ####
ProtTab_Norm = ProtTab_Full # The new table for the normalized data
overallMedian = median(as.vector(as.matrix(ProtTab_Full[,IdxIntensCol]))) # The overall median is calculated and used as a point for normalization

ProtTab_Norm[,IdxIntensCol] = apply(ProtTab_Full[,IdxIntensCol], 2, function(x){ # On the intensity columns, a function is applied per column
  x*overallMedian/median(x[which(x>0)])}) # This involves calculating the median for all values above 0, then calculating the ratio between the overall median and this column median, and normalizing each value in the column with it.

#### Visualization of the data ####
# Boxplots will be drawn before- and after normalization, to identify interactions, outliers, and the spread of the data. As well as the assumptions needed for the 2-way ANOVA.

#### Boxplot before normalization ####
# Set up the data, and transform however needed
Data_BNorm = melt(ProtTab_Full[,IdxIntensCol]) # All the data is melted from a many-column table, to a 2 column one
Data_BNormNoZero = subset(Data_BNorm, value !=0) # Remove all values that are 0
Data_BNormNoZero$value = log2(Data_BNormNoZero$value) # Make a log2 of all data

# Add a column containing information based on which the point color is defined
Data_BNormNoZero = cbind(Data_BNormNoZero, apply(as.data.frame(Data_BNormNoZero$variable), 1, function(x){
  unlist(strsplit(as.character(x), "_"))[2] # From the 'variable' column (which contains sample names) the 2nd word of the variable is taken and put in the new column
}))
colnames(Data_BNormNoZero)[3] = "color" # Rename the 3rd column to "color"

# Make a vector of the colors which the datapoints in the plots should have
labelStraincolor = gsub ("Resistent", "green", labelStrain)
labelStraincolor = gsub ("Intermediair", "red", labelStraincolor)
labelStraincolor = gsub ("Sensitive", "blue", labelStraincolor)

# Make the boxplot of the data before normalization
BoxplotFormat1 = theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1), plot.title = element_text(family = "Helvetica", face = "bold", size = (15), hjust = 0.5)) # A format to turn the x-axis 90 degrees and change the title format.
boxPlotBNorm = ggplot(Data_BNormNoZero, aes(x = variable, y = value)) + labs(title = "log2 intensity distribution before normalization", y = "log2(intensity)", x = "samples") +
  geom_violin(aes(col = color)) + geom_boxplot(outlier.color = "black", col = labelStraincolor, width=0.21) + BoxplotFormat1

## Print the plot and save it as a PNG
print(boxPlotBNorm)
ggsave("non_normalized_data.png", plot = boxPlotBNorm,
       scale = 1, width = 25, height = 20, units = "cm", dpi = 600)

#### Boxplot after normalization ####
# Set up the data, and transform however needed (Same as before, just changed the data origin)
Data_ANorm = melt(ProtTab_Norm[,IdxIntensCol]) # All the data is melted from a many-column table, to a 2 column one. Unlike before, the normalization table is used
Data_ANormNoZero = subset(Data_ANorm, value !=0) # Remove all values that are 0
Data_ANormNoZero$value = log2(Data_ANormNoZero$value) # Make a log2 of all data

# Add a column containing information based on which the point color is defined (just changed the data origin)
Data_ANormNoZero = cbind(Data_ANormNoZero, apply(as.data.frame(Data_ANormNoZero$variable), 1, function(x){
  unlist(strsplit(as.character(x), "_"))[2] # From the 'variable' column (which contains sample names) the 2nd word of the variable is taken and put in the new column
}))
colnames(Data_ANormNoZero)[3] = "color" # Rename the 3rd column to "color"

# Make a vector of the colors which the datapoints in the plots should have (same as before)
labelStraincolor = gsub ("Resistent", "green", labelStrain)
labelStraincolor = gsub ("Intermediair", "red", labelStraincolor)
labelStraincolor = gsub ("Sensitive", "blue", labelStraincolor)

# Make the boxplot of the data before normalization
BoxplotFormat1 = theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1), plot.title = element_text(family = "Helvetica", face = "bold", size = (15), hjust = 0.5)) # A format to turn the x-axis 90 degrees and change the title format.
boxPlotANorm = ggplot(Data_ANormNoZero, aes(x = variable, y = value)) + labs(title = "log2 intensity distribution after normalization", y = "log2(intensity)", x = "samples") +
  geom_violin(aes(col = color)) + geom_boxplot(outlier.color = "black", col = labelStraincolor, width=0.21) + BoxplotFormat1

## Print the plot and save it as a PNG
print(boxPlotANorm)
ggsave("normalized_data.png", plot = boxPlotANorm,
       scale = 1, width = 25, height = 20, units = "cm", dpi = 600)

#### Averaging duplicates ####
### The previous 2 plots contain duplicates of the same samples, so here they will be averaged

## Create a dummy frame for the calculations
DummyFrame = Data_ANorm # Create a dummy dataframe to contain the list of normalized intensity per sample.
DummyFrame$variable = gsub("_(1|2)$", "", DummyFrame$variable) # Remove the final '_1' or '_2' so that in the next step the duplicates can be grouped and averaged.
Accession = ProtTab_Norm$Accession # Generate an object containing the Accession numbers.
DummyFrame$accession = rep(Accession, length.out = nrow(Data_ANorm)) # Add a new column to the dummy dataframe containing the accession numbers, repeating the list for every sample.
DummyFrame = DummyFrame[, c("accession", "variable", "value")] # Reorganize the dataframe

# Now each value that has the same Accession number and the same sample duplicate can be averaged.
df_mean = aggregate(value ~ accession + variable, data = DummyFrame, function(x) {
  mean(x[x != 0])
})

# Cast the modified melted dataframe back into a wide format
ProtTab_Nonred = dcast(df_mean, accession ~ variable, value.var = "value")

#### Plot the averaged data in a boxplot ####
# Set up the data, and transform however needed. Different here, though the principle is the same
Data_AMean = df_mean[c("variable", "value")]
Data_AMeanNoZero = subset(Data_AMean, value !=0) # Remove all values that are 0
Data_AMeanNoZero$value = log2(Data_AMeanNoZero$value) # Make a log2 of all data

# Add a column containing information based on which the point color is defined (just changed the data origin)
Data_AMeanNoZero = cbind(Data_AMeanNoZero, apply(as.data.frame(Data_AMeanNoZero$variable), 1, function(x){
  unlist(strsplit(as.character(x), "_"))[2] # From the 'variable' column (which contains sample names) the 2nd word of the variable is taken and put in the new column
}))
colnames(Data_AMeanNoZero)[3] = "color" # Rename the 3rd column to "color"

# Make a vector of the colors which the datapoints in the plots should have (same as before)
labelStraincolor = gsub ("Resistent", "green", labelStrain)
labelStraincolor = gsub ("Intermediair", "red", labelStraincolor)
labelStraincolor = gsub ("Sensitive", "blue", labelStraincolor)

# Make the boxplot of the data before normalization
BoxplotFormat1 = theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1), plot.title = element_text(family = "Helvetica", face = "bold", size = (15), hjust = 0.5)) # A format to turn the x-axis 90 degrees and change the title format.
boxPlotAMean = ggplot(Data_AMeanNoZero, aes(x = variable, y = value)) + labs(title = "log2 intensity distribution after normalization and duplicate averageing", y = "log2(intensity)", x = "samples") +
  geom_violin(aes(col = color)) + geom_boxplot(outlier.color = "black", aes(color = color), width=0.21) + BoxplotFormat1

## Print the plot and save it as a PNG
print(boxPlotAMean)
ggsave("averaged_data.png", plot = boxPlotAMean,
       scale = 1, width = 25, height = 20, units = "cm", dpi = 600)

#### Differential expression analysis with a 2-way ANOVA between controls and treated samples ####
## Filter the data
# Initially, we will remove all proteins where even one value is 0. This is to make sure the ANOVA works. After this, this section will be edited again.
ProtTab_Nonred_complete = ProtTab_Nonred[complete.cases(ProtTab_Nonred), ] # With the 'complete.cases()' function, only those rows where there is a complete dataset are selected.

# Convert the data frame to a long one
Data_long = gather(ProtTab_Nonred_complete, key="variable", value="intensity", -accession) # Uses the 'gather()' function to make turn the data frame into a long format.
Data_long = separate(Data_long, variable, c("variable_type", "sensitivity_type", "replicate")) # Separate the column names from each other.
Data_long = Data_long[, -which(names(Data_long) == "replicate")] # Remove the 'replicate' column.

# Sets the data types
Data_long$variable_type = as.factor(Data_long$variable_type)
Data_long$sensitivity_type = as.factor(Data_long$sensitivity_type)
Data_long$intensity = as.numeric(Data_long$intensity)

## Perform ANOVA for each protein and create the list of tables
list_Tab_ANOVA = Data_long %>%
  group_by(accession) %>%
  summarize(
    mean = mean(intensity),
    sd = sd(intensity),
    anova_result = list(aov(intensity ~ variable_type * sensitivity_type))
  ) %>%
  ungroup()

# Extract information and create tables for each protein
list_of_tables = lapply(list_Tab_ANOVA$anova_result, function(anova) {
  anova_summary = summary(anova)
  p_values = anova_summary[[1]]$`Pr(>F)`
  df = data.frame(
    Variables = row.names(anova_summary[[1]]),
    p_value = p_values,
    stringsAsFactors = FALSE
  )
  return(df)
})

# Assign the table names as the accession numbers
names(list_of_tables) = list_Tab_ANOVA$accession

## Reformat the data in preparation for the heatmap
# Extract the table names into a vector
vector_tableNames = names(list_of_tables)

# Extract the p-values for variable_type
vector_variableP = c()

for (i in 1:length(list_of_tables)) {
  table = list_of_tables[[i]]
  value = table[1, 2]
  vector_variableP = c(vector_variableP, value)
}

# Extract the p-values for sensitivity_type
vector_sensitivityP = c()

for (i in 1:length(list_of_tables)) {
  table = list_of_tables[[i]]
  value = table[2, 2]
  vector_sensitivityP = c(vector_sensitivityP, value)
}

# Extract the p-values for the interaction (variable_type:sensitivity_type)
vector_interactionP = c()

for (i in 1:length(list_of_tables)) {
  table = list_of_tables[[i]]
  value = table[3, 2]
  vector_interactionP = c(vector_interactionP, value)
}

# Reassemble the 4 columns into one table
Data_ANOVA = data.frame(Accession = vector_tableNames, Pvalue_variableType = vector_variableP, Pvalue_sensitivityType = vector_sensitivityP, Pvalue_interaction = vector_interactionP)

## P-value adjustment with the Benjamini & Hochberg method
# Since vectors containing the p-values are already available, the method can be immediately executed
Adj_Pvalue_variableType = p.adjust(vector_variableP, method = "BH")
Adj_Pvalue_sensitivityType = p.adjust(vector_sensitivityP, method = "BH")
Adj_Pvalue_interaction = p.adjust(vector_interactionP, method = "BH")

# Add the adjusted p-values to the heatmap data table
Data_ANOVA$Adj_Pvalue_variableType = Adj_Pvalue_variableType
Data_ANOVA$Adj_Pvalue_sensitivityType = Adj_Pvalue_sensitivityType
Data_ANOVA$Adj_Pvalue_interaction = Adj_Pvalue_interaction

## Data filtering
# Select only those proteins where the p-values of both the adjusted variable_type and sensitivity_type are below 0.01
Selection_rows = Data_ANOVA$Adj_Pvalue_variableType < 0.01 & Data_ANOVA$Adj_Pvalue_sensitivityType < 0.01

# Create the filtered data frame, which will contain the significant proteins
Data_ANOVA_signProteins = Data_ANOVA[Selection_rows, ]

#### Heatmap generation ####
## Data preparation
# Create a vector containing the accession numbers of the significant proteins
Vector_selected_proteins = Data_ANOVA_signProteins$Accession

# Subset 'Data_long to retain only significant proteins
Data_long_sign = subset(Data_long, accession %in% Vector_selected_proteins)
  
# Calculate the overall mean of all the significant protein intensities individually
Data_meanSD = aggregate(intensity ~ accession, data = Data_long_sign, FUN = mean)
colnames(Data_meanSD)[2] = "overall_mean"

# Calculate the standard deviation for each mean
DummyFrame = aggregate(intensity ~ accession, data = Data_long_sign, FUN = sd)
colnames(DummyFrame)[2] = "overall_sd"

# Add the column to the heatmap data frame
Data_meanSD = merge(Data_meanSD, DummyFrame, by = "accession")
rm(DummyFrame)

# Calculate the z-values
Data_heatmap = data.frame(Data_long_sign, z_value = (Data_long_sign$intensity - Data_meanSD$overall_mean[match(Data_long_sign$accession, Data_meanSD$accession)]) / Data_meanSD$overall_sd[match(Data_long_sign$accession, Data_meanSD$accession)])

## Data clustering
# Subset the relevant data
DummyFrame = Data_heatmap[, c("variable_type", "sensitivity_type", "z_value")]

# Convert the data to a matrix
Matrix_heatmap = as.matrix(DummyFrame)

# Perform clustering on rows, here by protein
Data_heatmap_clustering = hclust(dist(Matrix_heatmap))

## Heatmap creation
Heatmap(Matrix_heatmap,
        col = colorRamp2(c(-3, 0, 3), c("blue", "white", "red")),
        cluster_rows = Data_heatmap_clustering,
        show_row_names = FALSE,
        show_column_names = FALSE
)