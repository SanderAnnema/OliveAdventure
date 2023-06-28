#### Introduction: ####
# Master project by Sander Annema. Started at 31/05/2023
# In this project, proteomics data from probiotic bacteria isolated from olives is processed. 
# These bacteria are either resistant, intermediate, or sensitive to human stomach bile.
# The expression data of proteins by these 3 strains was collected with- and without exposure to bile.
# The goal of this project is to determine the effects of protein expression on the exposure of the bacterial strains by bile.
# This will be achieved by performing a 2-way ANOVA, then generating heat-maps of the proteins.




#### Preparation ####
# To begin, the session will be cleaned by removing all variables.
rm(list=ls())

## Loading the required libraries. This part will need to be updated later on, as I still have a bunch of packages loaded that are only used in previous versions of the code.
library(reshape2)
library(ggplot2)
library(dplyr)
library(tidyr)
library(pheatmap)
library(data.table)
library(scales)
library(tidyverse)
library(UniProt.ws)
library(MASS)
library(gProfileR)
library(gprofiler2)

# Setting the working directory
setwd("C:/Users/Skyar/OneDrive/Documenten/school/Master_BiMoS/2nd_research_project/Project_Olives_New")

# Terms used for streamlining the code
DataFolder = "C:/Users/Skyar/OneDrive/Documenten/school/Master_BiMoS/2nd_research_project/Project_Olives_New/data/"
ProteinFile = "proteins.csv"
PeptideFile = "peptide.csv"
IdxIntensCol = 9:44




#### Loading and formatting data ####
ProtTab_Full = read.csv(paste(DataFolder,ProteinFile, sep = ""), header =  TRUE)

# The header names are too long, so 'intensity' will be removed from them for brevity
columnNames = gsub(".Intensity","", colnames(ProtTab_Full) [IdxIntensCol])
colnames(ProtTab_Full) [IdxIntensCol] = columnNames

# Create object defining various labels, to be used in future code.
labels = strsplit(columnNames, "_") # The column names are separated based on "_", so a table is formed with the following structure of columns: factor, group, replicate, measurement duplicate.
labelTreatment = unlist(lapply(labels, function(x) x[1])) # The first column, treatment, is isolated
labelStrain = unlist(lapply(labels, function(x) x[2])) # The second column, strain, is isolated
labelReplicate = unlist(lapply(labels, function(x) x[3])) # The third column, replicate, is isolated
labelTreatmentStrainRep = paste(labelTreatment, labelStrain, labelReplicate, sep = "_") # The three columns are reassembled into new labels containing only vital information
labelTreatmentStrainUnique = unique(labelTreatmentStrainRep) # The duplicates are removed

# Define plot formats
BoxplotFormat1 = theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1), plot.title = element_text(family = "Helvetica", face = "bold", size = (15), hjust = 0.5)) # A format to turn the x-axis 90 degrees and change the title format.




#### Define the functions ####

### Subset data
Function_subset = function(data, index) {
 data_sub = data[, index]
 return(data_sub)
}

### Remove proteins with even one zero from a wide format data frame
Function_NoZero_wide = function(data_wide) {
  row_indices = rowSums(data_wide ==0) == 0
  data_NoZero = data_wide[row_indices, ]
  return(data_NoZero)
}

### Take log2 of values
Function_takeLog2 = function(data) {
  data_log2 = log2(data)
  return(data_log2)
}

### Convert data to long format, with the option of excluding columns.
Function_makeLong = function(data, exclude_columns = NULL) {
  if (is.null(exclude_columns)) {
    data_long = pivot_longer(data, cols = everything(), names_to = "variable", values_to = "value")
  } else {
    data_long = pivot_longer(data, cols = -all_of(exclude_columns), names_to = "variable", values_to = "value")
  }
  return(data_long)
}

### Set the rownames as the accession number of any source file
Function_setAccession = function(data, accession_origin) {
  rownames(data) = accession_origin$Accession
  return(data)
}

### Drawing a histogram, involving only the plotting of the graph.
Function_drawHistogram = function(data, title) {
  # Create a data frame with a single column for the data
  df = data.frame(value = data)
  
  # Create the histogram plot using ggplot2
  ggplot(df, aes(x = value)) +
    geom_histogram(binwidth = 1, fill = "skyblue", color = "black") +
    labs(x = "log2(Intensity)", y = "Frequency", title = title)
}

### Adding a column called 'color' to a dataframe and filling it with the strain types
Function_add_colorColumn = function(data) {
  data$color = sapply(strsplit(as.character(data$variable), "_"), function(x) x[2])
  return(data)
}

### Drawing a box plot, involving the adding of a column defining the color based on the variable, defining the format of the box plot, and the drawing of the plot itself.
Function_drawBoxplot = function(data, title) {
  plot = ggplot(data, aes(x = variable, y = value)) +
    labs(title = title, y = "log2(intensity)", x = "samples") +
    geom_violin(aes(col = color)) +
    geom_boxplot(outlier.color = "black", width = 0.21) +
    BoxplotFormat1
  
  return(plot)
}

### Printing and saving the histogram or boxplot
Function_savePlot = function(plot, filename, plotType) {
  if (plotType == "boxplot") {
    # Print the box plot
    print(plot)
    
    # Save the box plot as a PNG
    ggsave(filename, plot = plot, scale = 1, width = 25, height = 20, units = "cm", dpi = 600)
    
  } else if (plotType == "histogram") {
    # Print the box plot
    print(plot)
    
    # Save the histogram
    ggsave(filename, plot = plot, width = 8, height = 6, dpi = 300)
  } else {
    stop("Invalid plot type. Please specify 'boxplot' or 'histogram'.")
  }
}

### Median normalization
Function_MedianNorm = function(data, index) {
  # Generate a new protein table which will be the output
  ProtTab_ANorm = data
  
  # Calculate the overall median of the intensities within the protein table
  Value_overallMedian = median(as.vector(as.matrix(data[,index])))
  
  # On the intensity columns of the protein table, normalize
  ProtTab_ANorm[, index] = apply(data[, index], 2, function(x) {
    x * Value_overallMedian / median(x[which(x > 0)])
  })
  
  return(ProtTab_ANorm)
}

### Identification of triplicates with too little data
Function_IdentBarrenProt = function(data) {
  # Define a sub-function to count non-zero values in a row
  Function_CountNonZero = function(row) {
    sum(row != 0)
  }
  
  # Apply the sub-function on each row to count non-zero values
  Count_non_zero = apply(data, 1, Function_CountNonZero)
  
  # Find the row names where the non-zero count is zero
  Barren_Proteins = rownames(data)[Count_non_zero == 0]
  
  return(Barren_Proteins)
}

### Mean imputation of a dataset, where each imputation is performed on one triplicate measurement. Those proteins where there's only one value will be removed.
## Note: This function only acts upon a subset, so you need to select the subsets yourself.
Function_TripMeanImput = function(data, min_non_zero = 1) {
  # Define a sub-function for mean imputation over rows
  Function_meanImput = function(row) {
    row_mean = mean(row)
    row[row == 0] = row_mean
    return(row)
  }
  
  # Apply the sub-function on the subset to impute the missing values
  num_non_zero = rowSums(data != 0)
  data[num_non_zero >= min_non_zero, ] = t(apply(data[num_non_zero >= min_non_zero, ], 1, Function_meanImput))

    return(data)

}

### Duplicate averageing
Function_duplicateAverageing = function(data, accession_origin_column) {
  # Convert the data to a long format
  data_long = Function_makeLong(data)
  
  # Remove the duplicate markings within the variable column
  data_long$variable = gsub("_(1|2)$", "", data_long$variable)
  
  # Add another column to the table containing the accession numbers
  data_long$Accession = rep(accession_origin_column, each = 36)
  
  # Reorganize the data frame
  data_long = data_long[, c("Accession", "variable", "value")]
  
  # Average duplicates
  data_mean = aggregate(value ~ Accession + variable, data = data_long, function(x) {
    mean(x[x != 0])
  })
  
  # Cast the modified melted dataframe back into a wide format
  data_mean = data_mean %>%
    pivot_wider(names_from = variable, values_from = value)
  
  # Set the wide data frame as a data frame
  data_mean = as.data.frame(data_mean)
  
  # Make the accession numbers the row names and remove the Accession column
  data_mean = Function_setAccession(data_mean, data_mean)
  data_mean = data_mean[, !(names(data_mean) == "Accession")]
  
  return(data_mean)
}

#### Differential expression analysis with a 2-way ANOVA between controls and treated samples
## Note: Keep in mind to make sure that the format of the 'data' data frame is correct. It should have the accession numbers as the row names, and only the raw intensities in the data frame itself.
Function_performANOVA = function(data, p_value) {  
  #   Create a new column for the accession numbers needed for the analysis
  data$Accession = rownames(data)
  
  # Move the column to the beginning, removing the copy and the row names.
  rownames(data) = NULL
  data = data[, c("Accession", names(data)[-ncol(data)])]
  
  # Convert the data frame to a long format
    Data_long = Function_makeLong(data, exclude_columns = "Accession")
    
    # Separate the variable column into 3
    Data_long = separate(Data_long, variable, c("variable_type", "sensitivity_type", "replicate"))
    
    # Set the data types
    Data_long$variable_type = as.factor(Data_long$variable_type)
    Data_long$sensitivity_type = as.factor(Data_long$sensitivity_type)
    Data_long$intensity = as.numeric(Data_long$value)
    Data_long$replicate = as.factor(Data_long$replicate)
    
    ## Perform ANOVA for each protein and create the list of tables
    list_Tab_ANOVA = Data_long %>%
      group_by(Accession) %>%
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
    names(list_of_tables) = list_Tab_ANOVA$Accession
    
    ## Reformat the data
    # Extract the table names into a vector
    vector_tableNames = names(list_of_tables)
    
    # Extract the p-values for variable_type
    vector_variableP = sapply(list_of_tables, function(table) table[1, 2])
    
    # Extract the p-values for sensitivity_type
    vector_sensitivityP = sapply(list_of_tables, function(table) table[2, 2])
    
    # Extract the p-values for the interaction (variable_type:sensitivity_type)
    vector_interactionP = sapply(list_of_tables, function(table) table[3, 2])
    
    # Reassemble the 4 columns into one table
    Data_ANOVA = data.frame(
      Accession = vector_tableNames,
      Pvalue_variableType = vector_variableP,
      Pvalue_sensitivityType = vector_sensitivityP,
      Pvalue_interaction = vector_interactionP
    )
    
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
    # Select only those proteins where the adjusted p-values of both variable_type and sensitivity_type are below 'p_value'
    Selection_rows = Data_ANOVA$Adj_Pvalue_variableType < p_value & Data_ANOVA$Adj_Pvalue_sensitivityType < p_value & Data_ANOVA$Adj_Pvalue_interaction < p_value
    
    # Create the filtered data frame, which will contain the significant proteins
    Data_ANOVA_signProteins = Data_ANOVA[Selection_rows, ]
    
    output = list(
      list_Tab_ANOVA = list_Tab_ANOVA,
      list_of_tables = list_of_tables,
      Data_ANOVA = Data_ANOVA,
      Data_ANOVA_signProteins = Data_ANOVA_signProteins
    )
    
    return(output)
}
  
### Heatmap generation
Function_drawHeatmap = function(data, ANOVA_output, title = title) {
  ## Data preparation
  # Create a new column for the accession numbers needed for the analysis
  data$Accession = rownames(data)
  
  # Move the column to the beginning, removing the copy and the row names.
  rownames(data) = NULL
  data = data[, c("Accession", names(data)[-ncol(data)])]
  
  # Convert the data frame to a long format
  Data_long = Function_makeLong(data, exclude_columns = "Accession")
  
  # Separate the variable column into 3
  Data_long = separate(Data_long, variable, c("variable_type", "sensitivity_type", "replicate"))
  
  # Set the data types
  Data_long$variable_type = as.factor(Data_long$variable_type)
  Data_long$sensitivity_type = as.factor(Data_long$sensitivity_type)
  Data_long$intensity = as.numeric(Data_long$value)
  Data_long$replicate = as.factor(Data_long$replicate)
  
  ## Heatmap data prep
  # Create a vector containing the accession numbers of the significant proteins
  Sign_Prot = ANOVA_output$Data_ANOVA_signProteins$Accession
  
  # Subset the intensity data frame to retain only significant proteins
  sub_data = subset(Data_long, Accession %in% Sign_Prot)
  
  # Convert the intensities to log2 to decrease variance between samples.
  sub_data = subset(Data_long, Accession %in% Sign_Prot)
  names(sub_data)[5] = "log2_intensity"
  
  # Calculate the overall mean of all the significant protein intensities individually
  Data_meanSD = aggregate(log2_intensity ~ Accession, data = sub_data, FUN = mean)
  colnames(Data_meanSD)[2] = "overall_mean"
  
  # Calculate the standard deviation for each mean and add it to the data frame
  data_sd = aggregate(log2_intensity ~ Accession, data = sub_data, FUN = sd)
  colnames(data_sd)[2] = "overall_sd"
  Data_meanSD = merge(Data_meanSD, data_sd, by = "Accession")
  
  # Calculate the z-values
  Data_heatmap = data.frame(sub_data, z_value = (sub_data$log2_intensity - Data_meanSD$overall_mean[match(sub_data$Accession, Data_meanSD$Accession)]) / Data_meanSD$overall_sd[match(sub_data$Accession, Data_meanSD$Accession)])
  
  # Ensure the z_value column is numeric
  Data_heatmap$z_value = as.numeric(as.character(Data_heatmap$z_value))
  
  # Select the relevant columns from Data_heatmap
  rel_col = Data_heatmap[, c("Accession", "variable_type", "sensitivity_type", "replicate", "z_value")]
  rel_col$z_value = as.numeric(as.character(rel_col$z_value))
  
  # Pivot the data to create a matrix with proteins as rows and combinations as columns
  Matrix_heatmap = reshape2::dcast(rel_col, Accession ~ variable_type + sensitivity_type + replicate, 
                                   value.var = "z_value")
  
  # Convert the 'accession' column into row names and delete the column
  rownames(Matrix_heatmap) = Matrix_heatmap$Accession
  Matrix_heatmap$Accession = NULL
  
  ## Draw the heatmap
  # Specify a randomized number to ensure a reproducible heatmap
  set.seed(12345)
  
  # Create the heatmap based on the matrix
  Plot_heatmap = pheatmap(Matrix_heatmap, 
                          clustering_distance_rows = "euclidean", 
                          clustering_distance_cols = "euclidean",
                          main = title)
  
  return(Plot_heatmap)
}

### Q-Q plot calculation for intensities
# Requires as input a long format and log2(intensity)
Function_calcQQ_int = function(data, title) {
  # Subset the relevant data from the data frame
  sub_data = data[, c("variable", "value")]
  
  # Calculate the Q-Q plot data
  output = qqplot(sub_data$value, ppoints(nrow(sub_data)), main = title)
  
  return(output)
}

### Q-Q plot calculation for p-values
Function_calcQQ_P = function(data) {
  # Subset the data
  sub_data = data[, c("Accession", "Adj_Pvalue_variableType", "Adj_Pvalue_sensitivityType", "Adj_Pvalue_interaction")]
  
  ## Generate the sub-plots
  # Variable type
  QQ_variable = qqplot(sub_data$Adj_Pvalue_variableType, ppoints(nrow(sub_data)), main = "Q-Q Plot: Variable Type")
  
  # Sensitivity type
  QQ_sensitivity = qqplot(sub_data$Adj_Pvalue_sensitivityType, ppoints(nrow(sub_data)), main = "Q-Q Plot: Sensitivity Type")
  
  # Interaction
  QQ_interaction = qqplot(sub_data$Adj_Pvalue_interaction, ppoints(nrow(sub_data)), main = "Q-Q Plot: Interaction")
  
  # Combined factors
  sub_data_full = c(sub_data$Adj_Pvalue_variableType, sub_data$Adj_Pvalue_sensitivityType, sub_data$Adj_Pvalue_interaction)
  QQ_combined = qqplot(sub_data_full, ppoints(length(sub_data_full)), main = "Q-Q Plot: All factors")
  
  # Return the 4 plots as a list
  list(QQ_variable = QQ_variable, QQ_sensitivity = QQ_sensitivity, QQ_interaction = QQ_interaction, QQ_combined = QQ_combined)
}

### Q-Q plot drawing
Function_drawQQ = function (QQPlot_data, plot_type, title) {
  if (plot_type == "intensity") {
    # Convert the plot to a data frame
    QQ_data = data.frame(x = QQPlot_data$x, y = QQPlot_data$y)
    
    # Plot this data frame with ggplot2
    QQ_plot = ggplot(QQ_data, aes(x, y)) +
      geom_point() +
      xlab("Observed Quantiles") +
      ylab("Theoretical Quantiles") +
      ggtitle(title)
    
    return(QQ_plot)
    
  }
  else if (plot_type == "p_values") {
    # Convert the plot to a data frame
    QQ_data = data.frame(x = QQPlot_data$x, y = QQPlot_data$y)
    
    # Plot this data frame with ggplot2
    QQ_plot = ggplot(QQ_data, aes(x, y)) +
      geom_point() +
      xlab("Observed Quantiles") +
      ylab("Theoretical Quantiles") +
      xlim(0, 1) +
      ylim(0, 1) +
      ggtitle(title)
    
    return(QQ_plot)
  }
}

### Calculations for the visualization of the distribution of missing values
Function_calcPercentDistr = function(data, index, title) {
  # Create a subset of the dataset containing only the protein intensities after normalization with a protein per row
  DummyFrame = Function_subset(data, index)
  row.names(DummyFrame) = data$Accession
  
  # Calculate the percentage of missing values per protein
  Data_PercentageZero_ANorm = data.frame(Percentage_Zero = rowMeans(DummyFrame == 0) * 100)
  
  # Create a new variable for the bins (0-10%, 10-20%, ..., 90-100%)
  Data_PercentageZero_ANorm$Bin = cut(Data_PercentageZero_ANorm$Percentage_Zero, breaks = c(seq(0, 100, by = 10), Inf), right = FALSE, labels = FALSE)
  
  # Count the number of proteins in each bin
  Data_PercentageZeroCounts = table(Data_PercentageZero_ANorm$Bin)
  
  # Convert the table to a data frame
  Data_PercentageZeroCounts = as.data.frame(Data_PercentageZeroCounts)
  names(Data_PercentageZeroCounts) = c("Bin", "Count")
  Data_PercentageZeroCounts$Bin = as.numeric(Data_PercentageZeroCounts$Bin) -1
  
  # Generate the labels for the x-axis
  bin_labels = paste0(Data_PercentageZeroCounts$Bin * 10, "%-", Data_PercentageZeroCounts$Bin * 10 + 10, "%")
  
  # Plot the histogram
  Plot_Hist_PercentageZero_ANorm = ggplot(Data_PercentageZeroCounts, aes(x = reorder(bin_labels, -Count), y = Count)) +
    geom_bar(stat = "identity", fill = "skyblue") +
    xlab("Percentage Range") +
    ylab("Number of Proteins") +
    ggtitle(title) +
    scale_x_discrete(labels = bin_labels) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  
  return(Plot_Hist_PercentageZero_ANorm)
}




#### !!! BEGINNING OF OVERALL DATA PROCESSING: NO DIFFERENCES BETWEEN THE 4 METHODS !!! ####
## The data needs to be prepared first, which includes visualization of the data

#### Histogram before imputation and before normalization ####
## Data preparation
Data_BNorm = Function_subset(ProtTab_Full, IdxIntensCol)
Data_BNorm = Function_takeLog2(Data_BNorm)
Data_BNorm = Function_makeLong(Data_BNorm)

## Draw the histogram
Plot_Hist_BNorm = Function_drawHistogram(Data_BNorm$value, "log2(Intensity) distribution before normalization and imputation")

## Save the histogram
Function_savePlot(Plot_Hist_BNorm, filename = "Histogram_Intensity_Distribution_BNorm_Full.png", plotType = "histogram")




#### Boxplot before imputation and before normalization ####
# Add a column to define violin plot color
Data_BNorm = Function_add_colorColumn(Data_BNorm)

# Make the boxplot of the data before normalization
Plot_Box_BNorm = Function_drawBoxplot(Data_BNorm, title = "log2 intensity distribution before normalization")

## Print the plot and save it as a PNG
Function_savePlot(Plot_Box_BNorm, filename = "Boxplot_data_Full.png", plotType = "boxplot")




#### Median normalization of all data ####
# Execute the function, generating the full protein table but normalized
ProtTab_ANorm = Function_MedianNorm(ProtTab_Full, IdxIntensCol)

# Subset relevant data from the full normalized protein table
Data_ANorm = Function_subset(ProtTab_ANorm, IdxIntensCol)

# Add the accession numbers as the row names
Data_ANorm = Function_setAccession(Data_ANorm, ProtTab_Full)




#### Boxplot after normalization ####
## Data preparation
# Take the log2 of the intensities
Data_ANorm_Box = Function_takeLog2(Data_ANorm)

# Replace -Inf values with 0, which were generated by taking the log2 of 0 or negative post-normalization values. Then convert all zeros to NA.
Data_ANorm_Box[Data_ANorm_Box == -Inf] = 0
Data_ANorm_Box[Data_ANorm_Box == 0] = NA

# Make the table long
Data_ANorm_Box = Function_makeLong(Data_ANorm_Box)

# Add a column containing information based on which the point color is defined
Data_ANorm_Box = Function_add_colorColumn(Data_ANorm_Box)

# Make the boxplot of the data after normalization
Plot_Box_ANorm = Function_drawBoxplot(Data_ANorm_Box, title = "log2 intensity distribution after normalization")

# Print and save the boxplot
Function_savePlot(Plot_Box_ANorm, filename = "Boxplot_normalized_data_Full.png", plotType = "boxplot")




#### Q-Q plots of normalized analysis to determine data normality ####
# Calculations
QQData_intensity = Function_calcQQ_int(Data_ANorm_Box, "Q-Q Plot: Normalized Intensities")

# Drawing the plot
Plot_QQ_intensity = Function_drawQQ(QQData_intensity, "intensity", "Q-Q Plot: Normalized Intensities")

## Printing and saving the plot
print(Plot_QQ_intensity)
ggsave("QQPlot_intensity_Full.png", plot = Plot_QQ_intensity, scale = 1, width = 8, height = 6, units = "in", dpi = 300)




#### Histogram showing the degree of missing values ####
# Perform the calculations needed to draw the histogram
Plot_Hist_PercentageZero_ANorm = Function_calcPercentDistr(ProtTab_ANorm, IdxIntensCol, "Distribution of Missing Values per Protein")

# Print the plot
print(Plot_Hist_PercentageZero_ANorm)

# Save the histogram
ggsave("Distribution_MissingValues_PerProtein_ANorm.png", plot = Plot_Hist_PercentageZero_ANorm, width = 10, height = 6, dpi = 300)




#### !!! BEGINNING OF DATA PROCESSING METHOD 1: ONLY PROTEINS WITH FULL DATA !!! ####
## In this method, no imputation is performed. Any proteins that don't have full data will be removed.

#### Remove proteins with missing data ####
# Remove proteins with missing data
ProtTab_NoZero = ProtTab_ANorm[rowSums(ProtTab_ANorm == 0) == 0, ]

# Calculate the percentage of proteins removed
Value_RetainedProt_NoZero = nrow(ProtTab_NoZero) / nrow(ProtTab_ANorm) * 100

# Generate a data file for further processing
Data_NoZero = Function_subset(ProtTab_NoZero, IdxIntensCol)
Data_NoZero = Function_setAccession(Data_NoZero, ProtTab_NoZero)




#### Averaging duplicates ####
### Calculation of averages
Data_NoZero_Mean = Function_duplicateAverageing(Data_NoZero, ProtTab_NoZero$Accession)




#### Visualisation of averaged data ####
### Boxplot of averaged data
# Calculations
Data_NoZero_Mean_Box = Function_takeLog2(Data_NoZero_Mean)
Data_NoZero_Mean_Box = Function_makeLong(Data_NoZero_Mean_Box)
Data_NoZero_Mean_Box = Function_add_colorColumn(Data_NoZero_Mean_Box)

# Drawing and saving the boxplot
Plot_Box_NoZero_Mean = Function_drawBoxplot(Data_NoZero_Mean_Box, title = "log2 intensity distribution after normalization and duplicate averageing of proteins with complete data")
Function_savePlot(Plot_Box_NoZero_Mean, filename = "MeanIntens_Complete.png", plotType = "boxplot")




#### 2-way ANOVA of the NoZero dataset ####
### Differential expression analysis with a 2-way ANOVA between controls and treated samples
## Perform the 2-way ANOVA with p-value adjustment
Output_ANOVA_NoZero = Function_performANOVA(Data_NoZero_Mean, p_value = 0.01)




#### Heatmap generation ####
## Perform the heatmap generation
Plot_Heat_NoZero = Function_drawHeatmap(Data_NoZero_Mean, Output_ANOVA_NoZero, title = "Heatmap depicting z-values of proteins with complete data")

## Save the heatmap as a PNG
ggsave("Heatmap_NoZero.png", plot = Plot_Heat_NoZero,
       scale = 1, width = 25, height = 20, units = "cm", dpi = 600)




#### !!! BEGINNING OF DATA PROCESSING METHOD 2: Imputation method 1, where imputation is performed when there's at least one value per triplicate !!! ####
## In this method, if there's even only one value per triplicate, mean imputation will be performed on the missing values. 

#### Remove proteins lacking too much data (with less than one value per triplicate) ####
# Make 6 subsets for individual triplicates
subset_control_sensitive = subset(Data_ANorm, select = grepl("^Control_Sensitive", colnames(Data_ANorm)))
subset_treated_sensitive = subset(Data_ANorm, select = grepl("^Treated_Sensitive", colnames(Data_ANorm)))
subset_control_intermediair = subset(Data_ANorm, select = grepl("^Control_Intermediair", colnames(Data_ANorm)))
subset_treated_intermediair = subset(Data_ANorm, select = grepl("^Treated_Intermediair", colnames(Data_ANorm)))
subset_control_resistent = subset(Data_ANorm, select = grepl("^Control_Resistent", colnames(Data_ANorm)))
subset_treated_resistent = subset(Data_ANorm, select = grepl("^Treated_Resistent", colnames(Data_ANorm)))

# Determine which proteins have 1 or less values within a row. This allows the future elimination of proteins with too little data.
List_ProtBarren_CSens = Function_IdentBarrenProt(subset_control_sensitive)
List_ProtBarren_TSens = Function_IdentBarrenProt(subset_treated_sensitive)
List_ProtBarren_CInt = Function_IdentBarrenProt(subset_control_intermediair)
List_ProtBarren_TInt = Function_IdentBarrenProt(subset_treated_intermediair)
List_ProtBarren_CRes = Function_IdentBarrenProt(subset_control_resistent)
List_ProtBarren_TRes = Function_IdentBarrenProt(subset_treated_resistent)

# Reassemble the lists into one list of all proteins where even one factor has too much missing data, then remove the subset lists
List_ProtBarren = unique(c(List_ProtBarren_CSens, List_ProtBarren_TSens, List_ProtBarren_CInt, List_ProtBarren_TInt, List_ProtBarren_CRes, List_ProtBarren_TRes))
rm(List_ProtBarren_CSens)
rm(List_ProtBarren_TSens)
rm(List_ProtBarren_CInt)
rm(List_ProtBarren_TInt)
rm(List_ProtBarren_CRes)
rm(List_ProtBarren_TRes)

# Check the removed proteins
Test_CheckProtBarren_Imp1 = Data_ANorm[rownames(Data_ANorm) %in% List_ProtBarren, ]

# Make a subset of Data_ANorm containing only usable proteins, based on the list of barren proteins made above
Data_Imp1 = Data_ANorm[!(rownames(Data_ANorm) %in% List_ProtBarren), ]

# Calculate the percentage of proteins removed
Value_RetainedProt_Imp1 = nrow(Data_Imp1) / nrow(ProtTab_ANorm) * 100




#### Mean imputation upon the remaining proteins ####
# Separate the data file into the 6 factors, and put them in a list
List_subsets = list(
  subset_control_sensitive = subset(Data_Imp1, select = grepl("^Control_Sensitive", colnames(Data_Imp1))),
  subset_treated_sensitive = subset(Data_Imp1, select = grepl("^Treated_Sensitive", colnames(Data_Imp1))),
  subset_control_intermediair = subset(Data_Imp1, select = grepl("^Control_Intermediair", colnames(Data_Imp1))),
  subset_treated_intermediair = subset(Data_Imp1, select = grepl("^Treated_Intermediair", colnames(Data_Imp1))),
  subset_control_resistent = subset(Data_Imp1, select = grepl("^Control_Resistent", colnames(Data_Imp1))),
  subset_treated_resistent = subset(Data_Imp1, select = grepl("^Treated_Resistent", colnames(Data_Imp1)))
)

# On each of the subsets, perform mean imputation within factors with at least two non-zero values
List_subsets = lapply(List_subsets, function(subset) {
  Function_TripMeanImput(subset, min_non_zero = 1)
})

# Reassemble the 6 separate data frames into one full data set. This required some fiddling with the column names
Data_Imp1 = do.call(cbind, lapply(names(List_subsets), function(name) {
  colnames(List_subsets[[name]]) = sub(paste0("^", name, "."), "", colnames(List_subsets[[name]]))
  List_subsets[[name]]
}))

# Remove the subset files
rm(List_subsets)




#### Visualisation of mean imputed data with method 1 ####
### Histogram to show log2(intensity) distribution after normalization and mean imputation, while removing those proteins that had too little data
# Data preparation
Data_Imp1_Hist = Function_takeLog2(Data_Imp1)
Data_Imp1_Hist = Function_makeLong(Data_Imp1_Hist)

# Draw the histogram
Plot_Hist_Imp1 = Function_drawHistogram(Data_Imp1_Hist$value, title = "log2(Intensity) distribution after mean imputation with method 1")

# Save the histogram
Function_savePlot(Plot_Hist_Imp1, filename = "IntensDistrib_Hist_Imp1.png", plotType = "histogram")

### Boxplot of the post imp1 data
# Data preparation
Data_Imp1_Box = Function_takeLog2(Data_Imp1)
Data_Imp1_Box = Function_makeLong(Data_Imp1_Box)
Data_Imp1_Box = Function_add_colorColumn(Data_Imp1_Box)

# Draw the boxplot
Plot_Box_Imp1 = Function_drawBoxplot(Data_Imp1_Box, title = "log2(intensity) distribution after mean imputation with method 1")

# Save the boxplot
Function_savePlot(Plot_Box_Imp1, "IntensDistrib_Box_Imp1.png", plotType = "boxplot")




#### Duplicate averaging after imputation ####
# Perform the calculation
Data_Imp1_Mean = Function_duplicateAverageing(Data_Imp1, rownames(Data_Imp1))

#### Visualization of averaged data ####
### Boxplot visualization
# Data preparation 
Data_Imp1_Mean_Box = Function_takeLog2(Data_Imp1_Mean)
Data_Imp1_Mean_Box = Function_makeLong(Data_Imp1_Mean_Box)
Data_Imp1_Mean_Box = Function_add_colorColumn(Data_Imp1_Mean_Box)

# Draw the boxplot
Plot_Box_Imp1_Mean = Function_drawBoxplot(Data_Imp1_Mean_Box, title = "log2(intensity) distribution of proteins after mean imputation (method 1) and duplicate averaging")

## Print the plot and save it as a PNG
Function_savePlot(Plot_Box_Imp1_Mean, "IntensDistrib_Box_Imp1_Mean.png", plotType = "boxplot")



#### 2-way ANOVA of the Imp1 dataset ####
### Differential expression analysis with a 2-way ANOVA between controls and treated samples
## Perform the 2-way ANOVA with p-value adjustment
Output_ANOVA_Imp1 = Function_performANOVA(Data_Imp1_Mean, p_value = 0.01)




#### Heatmap generation ####
## Perform the heatmap generation
Plot_Heat_Imp1 = Function_drawHeatmap(Data_Imp1_Mean, Output_ANOVA_Imp1, title = "Heatmap depicting z-values of proteins with imputed data (method 1)")

## Save the heatmap as a PNG
ggsave("Heatmap_Imp1.png", plot = Plot_Heat_Imp1,
       scale = 1, width = 25, height = 20, units = "cm", dpi = 600)




#### !!! BEGINNING OF DATA PROCESSING METHOD 3: Imputation method 2, where even if some conditions are missing, smooth histogram imputation is performed !!! ####
## In this method, we seek to include those datasets that show zero expression of a protein under one condition, but then express under another.
## This will require the retention of those conditions which have only zero values, but where other conditions of the same protein have at least 2 >0 values per triplicate.

#### Imputation method 2: Where newly appearing proteins are not discarded
# Separate the data file into the 6 factors, and put them in a list
List_subsets = list(
  subset_control_sensitive = subset(Data_ANorm, select = grepl("^Control_Sensitive", colnames(Data_ANorm))),
  subset_treated_sensitive = subset(Data_ANorm, select = grepl("^Treated_Sensitive", colnames(Data_ANorm))),
  subset_control_intermediair = subset(Data_ANorm, select = grepl("^Control_Intermediair", colnames(Data_ANorm))),
  subset_treated_intermediair = subset(Data_ANorm, select = grepl("^Treated_Intermediair", colnames(Data_ANorm))),
  subset_control_resistent = subset(Data_ANorm, select = grepl("^Control_Resistent", colnames(Data_ANorm))),
  subset_treated_resistent = subset(Data_ANorm, select = grepl("^Treated_Resistent", colnames(Data_ANorm)))
)

# On each of the subsets, perform mean imputation within factors with at least two non-zero values
List_subsets = lapply(List_subsets, function(subset) {
  Function_TripMeanImput(subset, min_non_zero = 2)
  })




#### !!! BEGINNING OF DATA PROCESSING METHOD 4: Imputation method 3, where all proteins are kept, and smooth histogram replacement imputation is performed to fill in missing data !!! ####
## In this method, we seek to keep all the proteins, but drawing a smooth histogram and take the -3SD value to replace all zeros.

#### Replacement imputation using a smooth histogram and fitted Gaussian distributiom ####
## Generate the smooth histogram
# Prepare the data
DummyFrame = Function_makeLong(Data_ANorm)
DummyFrame$value = log10(DummyFrame$value)
DummyFrame = as.data.frame(DummyFrame)

# Since there are many missing values in the dataset, we need to remove these first.
DummyFrame$value[is.infinite(DummyFrame$value)] = 0
DummyFrame = DummyFrame[DummyFrame$value != 0, ]

# Generate a smooth histogram of the data
Plot_SmoothHist_Imp3 = ggplot(DummyFrame, aes(x = value)) +
  geom_density(fill = "skyblue", color = "black", alpha = 0.5, bw = 0.15) +
  labs(x = "log10(intensity)", y = "Density", title = "Smooth Histogram of log10(intensity) for imputation with method 3")

# Print and save the histogram
print(Plot_SmoothHist_Imp3)
ggsave("SmoothHist_Imp3.png", plot = Plot_SmoothHist_Imp3,
       scale = 1, width = 25, height = 20, units = "cm", dpi = 300)

## Fit a Gaussian distribution to the data
# Generate the fit, mean, and standard deviation of the Gaussian distribution
Data_Gaussian_Fit_Imp3 = fitdistr(DummyFrame$value, "normal")
Value_Gaussian_Mu_Imp3 = Data_Gaussian_Fit_Imp3$estimate[1]
Value_Gaussian_Sigma_Imp3 = Data_Gaussian_Fit_Imp3$estimate[2]

# Calculate minus 3 sigma value, which is the log10 intensity point that will be used to replace missing values.
Value_Gaussian_min3Sig_Imp3 = Value_Gaussian_Mu_Imp3 - 3 * Value_Gaussian_Sigma_Imp3

# Convert the log10 intensity back to regular intensity
Value_Rep_min_Imp3 = 10^(Value_Gaussian_min3Sig_Imp3)

## Replace all zeros in the dataset with the found lowest value
Data_Imp3 = Data_ANorm
Data_Imp3[Data_Imp3 == 0] = Value_Rep_min_Imp3

# Remove excessive files
rm(DummyFrame)




#### Visualisation of mean imputed data with method 3 ####
### Histogram to show log2(intensity) distribution after normalization and replacement imputation
# Data preparation
Data_Imp3_Hist = Function_takeLog2(Data_Imp3)
Data_Imp3_Hist = Function_makeLong(Data_Imp3_Hist)

# Draw the histogram
Plot_Hist_Imp3 = Function_drawHistogram(Data_Imp3_Hist$value, title = "log2(Intensity) distribution after smooth histogram imputation with method 3")

# Save the histogram
Function_savePlot(Plot_Hist_Imp3, filename = "IntensDistrib_Hist_Imp3.png", plotType = "histogram")

### Boxplot of the post Imp3 data
# Data preparation
Data_Imp3_Box = Function_takeLog2(Data_Imp3)
Data_Imp3_Box = Function_makeLong(Data_Imp3_Box)
Data_Imp3_Box = Function_add_colorColumn(Data_Imp3_Box)

# Draw the boxplot
Plot_Box_Imp3 = Function_drawBoxplot(Data_Imp3_Box, title = "log2(intensity) distribution after smooth histogram imputation with method 3")

# Save the boxplot
Function_savePlot(Plot_Box_Imp3, "IntensDistrib_Box_Imp3.png", plotType = "boxplot")




#### Duplicate averaging after imputation ####
# Perform the calculation
Data_Imp3_Mean = Function_duplicateAverageing(Data_Imp3, rownames(Data_Imp3))

#### Visualization of averaged data ####
### Boxplot visualization
# Data preparation 
Data_Imp3_Mean_Box = Function_takeLog2(Data_Imp3_Mean)
Data_Imp3_Mean_Box = Function_makeLong(Data_Imp3_Mean_Box)
Data_Imp3_Mean_Box = Function_add_colorColumn(Data_Imp3_Mean_Box)

# Draw the boxplot
Plot_Box_Imp3_Mean = Function_drawBoxplot(Data_Imp3_Mean_Box, title = "log2(intensity) distribution of proteins after smooth histogram imputation (method 3) and duplicate averaging")

## Print the plot and save it as a PNG
Function_savePlot(Plot_Box_Imp3_Mean, "IntensDistrib_Box_Imp3_Mean.png", plotType = "boxplot")



#### 2-way ANOVA of the Imp3 dataset ####
### Differential expression analysis with a 2-way ANOVA between controls and treated samples
## Perform the 2-way ANOVA with p-value adjustment
Output_ANOVA_Imp3 = Function_performANOVA(Data_Imp3_Mean, p_value = 0.01)




#### Heatmap generation ####
## Perform the heatmap generation
Plot_Heat_Imp3 = Function_drawHeatmap(Data_Imp3_Mean, Output_ANOVA_Imp3, title = "Heatmap depicting z-values of proteins with imputed data (method 3)")

## Save the heatmap as a PNG
ggsave("Heatmap_Imp3.png", plot = Plot_Heat_Imp3,
       scale = 1, width = 25, height = 20, units = "cm", dpi = 600)




#### !!!!! The stuff below here is from the previous version, but will be rolled into the code above. !!!!! ####

#### p-value normality analysis with Q-Q plots for the NoZero dataset where all proteins with even one remaining zero after imputation are removed ####
## To determine how the p-values from the entire dataset are distributed, 4 Q-Q plots are made 
## One for each factor, one for the interaction, and one for all the p-values.
# Calculate the basic Q-Q plots
Output_QQ_NoZero = Function_calcQQ_P(Output_ANOVA_NoZero$Data_ANOVA)

## Draw the plots using ggplot2, then save them
# Variable type
Plot_QQ_P_variable_NoZero = Function_drawQQ(Output_QQ_NoZero$QQ_variable, plot_type = "p_values")
print(Plot_QQ_P_variable_NoZero)
ggsave("QQPlot_variable_NoZero.png", plot = Plot_QQ_P_variable_NoZero, scale = 1, width = 8, height = 6, units = "in", dpi = 300)

# Sensitivity type
Plot_QQ_P_sensitivity_NoZero = Function_drawQQ(Output_QQ_NoZero$QQ_sensitivity, plot_type = "p_values")
print(Plot_QQ_P_sensitivity_NoZero)
ggsave("QQPlot_sensitivity_NoZero.png", plot = Plot_QQ_P_sensitivity_NoZero, scale = 1, width = 8, height = 6, units = "in", dpi = 300)

# Interaction
Plot_QQ_P_interaction_NoZero = Function_drawQQ(Output_QQ_NoZero$QQ_interaction, plot_type = "p_values")
print(Plot_QQ_P_interaction_NoZero)
ggsave("QQPlot_interaction_NoZero.png", plot = Plot_QQ_P_interaction_NoZero, scale = 1, width = 8, height = 6, units = "in", dpi = 300)

# Combined factors
Plot_QQ_P_combined_NoZero = Function_drawQQ(Output_QQ_NoZero$QQ_combined, plot_type = "p_values")
print(Plot_QQ_P_combined_NoZero)
ggsave("QQPlot_combined_NoZero.png", plot = Plot_QQ_P_combined_NoZero, scale = 1, width = 8, height = 6, units = "in", dpi = 300)




#### Histogram of the p-value distribution of all proteins of the NoZero dataset ####
## Data preparation
# Subset the relevant data
DummyFrame1 = Output_ANOVA_NoZero$Data_ANOVA[, c("Adj_Pvalue_variableType", "Adj_Pvalue_sensitivityType", "Adj_Pvalue_interaction")]
DummyFrame2 = setNames(data.frame(DummyFrame1$Adj_Pvalue_variableType), "P_values")
DummyFrame3 = setNames(data.frame(DummyFrame1$Adj_Pvalue_sensitivityType), "P_values")
DummyFrame4 = setNames(data.frame(DummyFrame1$Adj_Pvalue_interaction), "P_values")
DummyFrame5 = data.frame(P_values = (unlist(DummyFrame1)))

# Ensure all columns are numeric
DummyFrame1$Adj_Pvalue_variableType = as.numeric(DummyFrame1$Adj_Pvalue_variableType)
DummyFrame1$Adj_Pvalue_sensitivityType = as.numeric(DummyFrame1$Adj_Pvalue_sensitivityType)
DummyFrame1$Adj_Pvalue_interaction = as.numeric(DummyFrame1$Adj_Pvalue_interaction)
DummyFrame2$P_values = as.numeric(DummyFrame2$P_values)
DummyFrame3$P_values = as.numeric(DummyFrame3$P_values)
DummyFrame4$P_values = as.numeric(DummyFrame4$P_values)
DummyFrame5$P_values = as.numeric(DummyFrame5$P_values)

## Create the histograms
# Variable type
PlotHist_variable_NoZero = hist(DummyFrame2$P_values, breaks = 10, col = "skyblue",
                                xlab = "P-values", ylab = "Frequency",
                                main = "P-value distribution of the variable factor of all proteins",
                                xlim = c(0, 1))

# Sensitivity type
PlotHist_sensitivity_NoZero = hist(DummyFrame3$P_values, breaks = 10, col = "skyblue",
                                   xlab = "P-values", ylab = "Frequency",
                                   main = "P-value distribution of the sensitivity factor of all proteins",
                                   xlim = c(0, 1))

# Interaction
PlotHist_interaction_NoZero = hist(DummyFrame4$P_values, breaks = 10, col = "skyblue",
                                   xlab = "P-values", ylab = "Frequency",
                                   main = "P-value distribution of the interaction between factor of all proteins",
                                   xlim = c(0, 1))

# Combined p-values
PlotHist_combined_NoZero = hist(DummyFrame5$P_values, breaks = 10, col = "skyblue",
                                xlab = "P-values", ylab = "Frequency",
                                main = "Total p-value distribution of the all proteins",
                                xlim = c(0, 1))

## Save the plots
dev.copy(png, "PValue_Distribution_variableType_NoZero.png", width = 8, height = 6, units = "in", res = 300)
dev.off()

dev.copy(png, "PValue_Distribution_sensitivityType_NoZero.png", width = 8, height = 6, units = "in", res = 300)
dev.off()

dev.copy(png, "PValue_Distribution_interaction_NoZero.png", width = 8, height = 6, units = "in", res = 300)
dev.off()

dev.copy(png, "PValue_Distribution_combined_NoZero.png", width = 8, height = 6, units = "in", res = 300)
dev.off()

# Remove the dummy frames
rm(DummyFrame1)
rm(DummyFrame2)
rm(DummyFrame3)
rm(DummyFrame4)
rm(DummyFrame5)




#### Protein ID conversion using UniProt ####
## Prepare data
# Make a vector of protein IDs to be converted
Vector_ProteinIDs_Acc_NoZero = rownames(Data_AMean_ImpNoZero)

# Extract only the relevant section of the accession number. So only the part after the second '|' symbol.
Vector_ProteinIDs_Acc_NoZero = sub("^[^|]*\\|[^|]*\\|(.*)", "\\1", Vector_ProteinIDs_Acc_NoZero)

## Convert the accession numbers to ensembl gene IDs
Data_AccMapping = mapUniProt(from = "UniProtKB_AC-ID", to = 'GeneID', query = Vector_ProteinIDs_Acc_NoZero)
Vector_ProteinIDs_NoZero = Data_AccMapping$To




#### !!!! Using gProfiler for gene enrichment ####




#### p-value normality analysis with Q-Q plots for the replacement dataset where all proteins with even one remaining zero after imputation are removed ####
## To determine how the p-values from the entire dataset are distributed, 4 Q-Q plots are made 
## One for each factor, one for the interaction, and one for all the p-values.
# Calculate the basic Q-Q plots
Output_QQ_Rep = Function_calcQQ_P(Output_ANOVA_Rep$Data_ANOVA)

## Draw the plots using ggplot2, then save them
# Variable type
Plot_QQ_P_variable_Rep = Function_drawQQ(Output_QQ_Rep$QQ_variable, plot_type = "p_values")
print(Plot_QQ_P_variable_Rep)
ggsave("QQPlot_variable_Rep.png", plot = Plot_QQ_P_variable_Rep, scale = 1, width = 8, height = 6, units = "in", dpi = 300)

# Sensitivity type
Plot_QQ_P_sensitivity_Rep = Function_drawQQ(Output_QQ_Rep$QQ_sensitivity, plot_type = "p_values")
print(Plot_QQ_P_sensitivity_Rep)
ggsave("QQPlot_sensitivity_Rep.png", plot = Plot_QQ_P_sensitivity_Rep, scale = 1, width = 8, height = 6, units = "in", dpi = 300)

# Interaction
Plot_QQ_P_interaction_Rep = Function_drawQQ(Output_QQ_Rep$QQ_interaction, plot_type = "p_values")
print(Plot_QQ_P_interaction_Rep)
ggsave("QQPlot_interaction_Rep.png", plot = Plot_QQ_P_interaction_Rep, scale = 1, width = 8, height = 6, units = "in", dpi = 300)

# Combined factors
Plot_QQ_P_combined_Rep = Function_drawQQ(Output_QQ_Rep$QQ_combined, plot_type = "p_values")
print(Plot_QQ_P_combined_Rep)
ggsave("QQPlot_combined_Rep.png", plot = Plot_QQ_P_combined_Rep, scale = 1, width = 8, height = 6, units = "in", dpi = 300)




#### Histogram of the p-value distribution of all proteins of the replacement dataset ####
## Data preparation
# Subset the relevant data
DummyFrame1 = Output_ANOVA_Rep$Data_ANOVA[, c("Adj_Pvalue_variableType", "Adj_Pvalue_sensitivityType", "Adj_Pvalue_interaction")]
DummyFrame2 = setNames(data.frame(DummyFrame1$Adj_Pvalue_variableType), "P_values")
DummyFrame3 = setNames(data.frame(DummyFrame1$Adj_Pvalue_sensitivityType), "P_values")
DummyFrame4 = setNames(data.frame(DummyFrame1$Adj_Pvalue_interaction), "P_values")
DummyFrame5 = data.frame(P_values = (unlist(DummyFrame1)))

# Ensure all columns are numeric
DummyFrame1$Adj_Pvalue_variableType = as.numeric(DummyFrame1$Adj_Pvalue_variableType)
DummyFrame1$Adj_Pvalue_sensitivityType = as.numeric(DummyFrame1$Adj_Pvalue_sensitivityType)
DummyFrame1$Adj_Pvalue_interaction = as.numeric(DummyFrame1$Adj_Pvalue_interaction)
DummyFrame2$P_values = as.numeric(DummyFrame2$P_values)
DummyFrame3$P_values = as.numeric(DummyFrame3$P_values)
DummyFrame4$P_values = as.numeric(DummyFrame4$P_values)
DummyFrame5$P_values = as.numeric(DummyFrame5$P_values)

## Create the histograms
# Variable type
PlotHist_variable_Rep = hist(DummyFrame2$P_values, breaks = 10, col = "skyblue",
                         xlab = "P-values", ylab = "Frequency",
                         main = "P-value distribution of the variable factor of all proteins",
                         xlim = c(0, 1))

# Sensitivity type
PlotHist_sensitivity_Rep = hist(DummyFrame3$P_values, breaks = 10, col = "skyblue",
                            xlab = "P-values", ylab = "Frequency",
                            main = "P-value distribution of the sensitivity factor of all proteins",
                            xlim = c(0, 1))

# Interaction
PlotHist_interaction_Rep = hist(DummyFrame4$P_values, breaks = 10, col = "skyblue",
                            xlab = "P-values", ylab = "Frequency",
                            main = "P-value distribution of the interaction between factor of all proteins",
                            xlim = c(0, 1))

# Combined p-values
PlotHist_combined_Rep = hist(DummyFrame5$P_values, breaks = 10, col = "skyblue",
                         xlab = "P-values", ylab = "Frequency",
                         main = "Total p-value distribution of the all proteins",
                         xlim = c(0, 1))

## Save the plots
dev.copy(png, "PValue_Distribution_variableType_Rep.png", width = 8, height = 6, units = "in", res = 300)
dev.off()

dev.copy(png, "PValue_Distribution_sensitivityType_Rep.png", width = 8, height = 6, units = "in", res = 300)
dev.off()

dev.copy(png, "PValue_Distribution_interaction_Rep.png", width = 8, height = 6, units = "in", res = 300)
dev.off()

dev.copy(png, "PValue_Distribution_combined_Rep.png", width = 8, height = 6, units = "in", res = 300)
dev.off()

# Remove the dummy frames
rm(DummyFrame1)
rm(DummyFrame2)
rm(DummyFrame3)
rm(DummyFrame4)
rm(DummyFrame5)