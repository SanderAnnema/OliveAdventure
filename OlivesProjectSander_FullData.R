#### Introduction: ####
# Master project by Sander Annema. 31/05/2023
# In this project, proteomics data from probiotic bacteria isolated from olives is processed. 
# These bacteria are either resistant, intermediate, or sensitive to bile.
# The expression data of both proteins and peptides was collected.
# The goal of this project is to determine the effects of protein expression on the exposure of the bacterial strains by bile.
# This will be achieved by performing a 2-way ANOVA, then generating a heat-map of the proteins.
# Additionally, in this variant of the file, two imputation methods will be used to fill in missing data from the normalized dataset before the log2 is taken. This complete dataset is then processed in the same way as the 'Basic' variant.

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

### Mean imputation
Function_meanImput = function(row) {
  row_mean = mean(row)
  row[row == 0] = row_mean
  return(row)
}

### The data processing of the output of the mean imputation, involving transposing and changing of the table type
Function_ImputProcessing = function(data) {
  # Transpose row- and column names
  data = t(data)
  
  # Revert to data frames
  data = data.frame(data)
  
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
Function_drawHeatmap = function(data, ANOVA_output) {
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
  Plot_heatmap = pheatmap(Matrix_heatmap, clustering_distance_rows = "euclidean", clustering_distance_cols = "euclidean")
  
  return(Plot_heatmap)
}

### Q-Q plot calculation for intensities
Function_calcQQ_int = function(data, title) {
  # Calculate the log2 of the intensities
  data = Function_takeLog2(data)
  
  # Lengthen the data frame
  data_long = Function_makeLong(data)
  
  # Subset the relevant data from the data frame
  sub_data = data_long[, c("variable", "value")]
  
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
Function_drawQQ = function (QQPlot_data, plot_type) {
  if (plot_type == "intensity") {
    # Convert the plot to a data frame
    QQ_data = data.frame(x = QQPlot_data$x, y = QQPlot_data$y)
    
    # Plot this data frame with ggplot2
    QQ_plot = ggplot(QQ_data, aes(x, y)) +
      geom_point() +
      xlab("Observed Quantiles") +
      ylab("Theoretical Quantiles")
    
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
      ylim(0, 1)
    
    return(QQ_plot)
  }
}




#### Histogram before mean imputation and normalization ####
## Data preparation
Data_Hist_BNorm = Function_subset(ProtTab_Full, IdxIntensCol)
Data_Hist_BNorm = Function_takeLog2(Data_Hist_BNorm)
Data_Hist_BNorm = Function_makeLong(Data_Hist_BNorm)

## Draw the histogram
Plot_Hist_intensityBNorm = Function_drawHistogram(Data_Hist_BNorm$value, "log2(Intensity) distribution before normalization and imputation")

## Save the histogram
Function_savePlot(Plot_Hist_intensityBNorm, filename = "Intensity_Distribution_BeforeNormImput_Full.png", plotType = "histogram")




#### Boxplot before normalization, no imputation, remove those proteins where there's even one value that is 0 ####
## Data preparation, subset data, remove those proteins where values are missing, take log2, make long format
Data_Box_BNorm = Function_subset(ProtTab_Full, IdxIntensCol)
Data_Box_BNorm = Function_NoZero_wide(Data_Box_BNorm)
Data_Box_BNorm = Function_takeLog2(Data_Box_BNorm)
Data_Box_BNorm = Function_makeLong(Data_Box_BNorm)

# Add a column containing information based on which the point color is defined
Data_Box_BNorm = Function_add_colorColumn(Data_Box_BNorm)

# Make the boxplot of the data before normalization
Plot_Box_BNorm = Function_drawBoxplot(Data_Box_BNorm, title = "log2 intensity distribution before normalization")

## Print the plot and save it as a PNG
Function_savePlot(Plot_Box_BNorm, filename = "non_normalized_data_Full.png", plotType = "boxplot")




#### Median normalization ####
ProtTab_ANorm = ProtTab_Full # The new table for the normalized data
Value_overallMedian = median(as.vector(as.matrix(ProtTab_Full[,IdxIntensCol]))) # The overall median is calculated and used as a point for normalization

ProtTab_ANorm[,IdxIntensCol] = apply(ProtTab_Full[,IdxIntensCol], 2, function(x){ # On the intensity columns, a function is applied per column
  x*Value_overallMedian/median(x[which(x>0)])}) # This involves calculating the median for all values above 0, then calculating the ratio between the overall median and this column median, and normalizing each value in the column with it.


# Subset relevant data from the full normalized protein table
Data_ANorm = Function_subset(ProtTab_ANorm, IdxIntensCol)

# Add the accession numbers as the row names
Data_ANorm = Function_setAccession(Data_ANorm, ProtTab_Full)

# Remove excess files
rm(Value_overallMedian)




#### Boxplot after normalization, but before imputation ####
## Data preparation
Data_Box_ANorm = Function_takeLog2(Data_ANorm)
Data_Box_ANorm[Data_Box_ANorm == -Inf] = 0 # Replace -Inf values (which were generated by taking the log2 of 0 or negative post-normalization values) with 0
Data_Box_ANorm = Function_makeLong(Data_Box_ANorm)
Data_Box_ANorm = subset(Data_Box_ANorm, value !=0) # Remove all values that are 0

# Add a column containing information based on which the point color is defined
Data_Box_ANorm = Function_add_colorColumn(Data_Box_ANorm)

# Make the boxplot of the data after normalization
Plot_Box_ANorm = Function_drawBoxplot(Data_Box_ANorm, title = "log2 intensity distribution after normalization")

# Print and save the boxplot
Function_savePlot(Plot_Box_ANorm, filename = "normalized_data_Full.png", plotType = "boxplot")




#### Normalized intensity normality analysis with Q-Q plots ####
## Calculations
QQData_intensity_preImp = Function_calcQQ_int(Data_ANorm, "Q-Q Plot: Normalized Intensities")

## Drawing the plot
Plot_QQ_intensity_preImp = Function_drawQQ(QQData_intensity_preImp, "intensity")

## Printing and saving the plot
print(Plot_QQ_intensity_preImp)
ggsave("QQPlot_intensity_Full.png", plot = Plot_QQ_intensity_preImp, scale = 1, width = 8, height = 6, units = "in", dpi = 300)




#### Histogram showing the percentages of zeros per protein before averageing and imputation
# Generate a data frame contain
Data_PercentageZero_ANorm = data.frame(Percentage_Zero = rowMeans(Data_ANorm == 0) * 100)

# Draw the histogram
Plot_Hist_PercentageZero_ANorm = ggplot(Data_PercentageZero_ANorm, aes(x = rownames(Data_PercentageZero_ANorm), y = Percentage_Zero)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  xlab("Protein") +
  ylab("Percentage of Zeros") +
  ggtitle("Percentage of Zeros per protein") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 2))

# Print the plot
print(Plot_Hist_PercentageZero_ANorm)

# Save the histogram
ggsave("Percentage_ZeroPerProtein_ANorm.png", plot = Plot_Hist_PercentageZero_ANorm, width = 20, height = 6, dpi = 600)




#### Averaging duplicates after normalization, but before imputation ####
### The previous plots contain duplicates of the same samples, so here they will be averaged
Data_AMean = Function_duplicateAverageing(Data_ANorm, ProtTab_ANorm$Accession)




#### Plot the averaged pre-imputation data in a boxplot ####
## Set up the data, and transform however needed.
Data_Box_AMean = Function_takeLog2(Data_AMean)
Data_Box_AMean = Function_makeLong(Data_Box_AMean)
Data_Box_AMean = subset(Data_Box_AMean, value !=0) 

# Add a column containing information based on which the point color is defined
Data_Box_AMean = Function_add_colorColumn(Data_Box_AMean)

# Make the boxplot of the data after normalization and averageing
Plot_Box_AMean = Function_drawBoxplot(Data_Box_AMean, title = "log2 intensity distribution after normalization and duplicate averageing")

# Print and save the boxplot
Function_savePlot(Plot_Box_AMean, filename = "averaged_data_Full.png", plotType = "boxplot")




#### Imputation of missing values ####
### Mean imputation of those samples where at least 1 of the triplicates of each sample are >0. This will give a dataset of 6 values from which a mean can be taken, limiting the effects on the overall data.

## Prepare the data
Data_Imput = Data_ANorm

# Make 6 subsets for individual triplicates, so they can be imputed later with greater ease
subset_control_sensitive = subset(Data_Imput, select = grepl("^Control_Sensitive", colnames(Data_Imput)))
subset_treated_sensitive = subset(Data_Imput, select = grepl("^Treated_Sensitive", colnames(Data_Imput)))
subset_control_intermediair = subset(Data_Imput, select = grepl("^Control_Intermediair", colnames(Data_Imput)))
subset_treated_intermediair = subset(Data_Imput, select = grepl("^Treated_Intermediair", colnames(Data_Imput)))
subset_control_resistent = subset(Data_Imput, select = grepl("^Control_Resistent", colnames(Data_Imput)))
subset_treated_resistent = subset(Data_Imput, select = grepl("^Treated_Resistent", colnames(Data_Imput)))

# Apply the function to each subset data frame
subset_control_sensitive = apply(subset_control_sensitive, 1, Function_meanImput)
subset_treated_sensitive = apply(subset_treated_sensitive, 1, Function_meanImput)
subset_control_intermediair = apply(subset_control_intermediair, 1, Function_meanImput)
subset_treated_intermediair = apply(subset_treated_intermediair, 1, Function_meanImput)
subset_control_resistent = apply(subset_control_resistent, 1, Function_meanImput)
subset_treated_resistent = apply(subset_treated_resistent, 1, Function_meanImput)

## Output data management
subset_control_sensitive = Function_ImputProcessing(subset_control_sensitive)
subset_treated_sensitive = Function_ImputProcessing(subset_treated_sensitive)
subset_control_intermediair = Function_ImputProcessing(subset_control_intermediair)
subset_treated_intermediair = Function_ImputProcessing(subset_treated_intermediair)
subset_control_resistent = Function_ImputProcessing(subset_control_resistent)
subset_treated_resistent = Function_ImputProcessing(subset_treated_resistent)

# Reassemble the 6 separate data frames into one full data set
Data_Imput = as.data.frame(cbind(subset_control_sensitive, subset_treated_sensitive, subset_control_intermediair, subset_treated_intermediair, subset_control_resistent, subset_treated_resistent))

# Remove the subset files
rm(subset_control_sensitive)
rm(subset_treated_sensitive)
rm(subset_control_intermediair)
rm(subset_treated_intermediair)
rm(subset_control_resistent)
rm(subset_treated_resistent)




#### Imputation method 1: Remove those proteins which still lack data ####
## All the visualization and extra code will be written immediately below it. Then the replacement imputation method will be below that.
Data_ImpNoZero = Data_Imput[rowSums(Data_Imput == 0) == 0, ]




#### Histogram to show log2(intensity) distribution after normalization and mean imputation, while removing those proteins that had too little data ####
# Prepare the data
Data_Hist_ImpNoZero = Function_takeLog2(Data_ImpNoZero)
Data_Hist_ImpNoZero = Function_makeLong(Data_Hist_ImpNoZero)

## Draw the histogram
Plot_Hist_ImpNoZero = Function_drawHistogram(Data_Hist_ImpNoZero$value, title = "log2(Intensity) distribution after normalization, mean imputation and zero removal")

## Save the histogram
Function_savePlot(Plot_Hist_ImpNoZero, "Intensity_Distribution_AfterNormImput_NoZero_Full.png", plotType = "histogram")




#### Boxplot after normalization, mean imputation and removing proteins lacking data ####
# Data preparation
Data_Box_ImpNoZero = Function_takeLog2(Data_ImpNoZero)
Data_Box_ImpNoZero = Function_makeLong(Data_Box_ImpNoZero)
Data_Box_ImpNoZero = Function_add_colorColumn(Data_Box_ImpNoZero)

# Draw the boxplot
Plot_Box_ImpNoZero = Function_drawBoxplot(Data_Box_ImpNoZero, title = "log2 intensity distribution after normalization and NoZero imputation")

## Print the plot and save it as a PNG
Function_savePlot(Plot_Box_ImpNoZero, "normalized_imput_NoZero_Full.png", plotType = "boxplot")




#### Duplicate averageing after normalization, mean imputation, and removing proteins that still have missing values ####
Accession_ImpNoZero = rownames(Data_ImpNoZero)
Data_AMean_ImpNoZero = Function_duplicateAverageing(Data_ImpNoZero, Accession_ImpNoZero)
rm(Accession_ImpNoZero)



#### Plot the averaged post-imputation (NoZero) data in a boxplot ####
## Data preparation
Data_Box_AMean_ImpNoZero = Function_takeLog2(Data_AMean_ImpNoZero)
Data_Box_AMean_ImpNoZero = Function_makeLong(Data_Box_AMean_ImpNoZero)
Data_Box_AMean_ImpNoZero = Function_add_colorColumn(Data_Box_AMean_ImpNoZero)

# Draw the boxplot
Plot_Box_AMean_ImpNoZero = Function_drawBoxplot(Data_Box_AMean_ImpNoZero, title = "log2 intensity distribution after normalization, mean imputation, zero removal, and duplicate averageing")

## Print the plot and save it as a PNG
Function_savePlot(Plot_Box_AMean_ImpNoZero, "Averaged_ImpNoZero_Full.png", plotType = "boxplot")





#### NoZero: Differential expression analysis with a 2-way ANOVA between controls and treated samples ####
## Perform the 2-way ANOVA with p-value adjustment
Output_ANOVA_NoZero = Function_performANOVA(Data_AMean_ImpNoZero, p_value = 0.01)




#### Heatmap generation ####
## Perform the heatmap generation
Plot_Heat_NoZero = Function_drawHeatmap(Data_AMean_ImpNoZero, Output_ANOVA_NoZero)

## Save the heatmap as a PNG
ggsave("Heatmap_Full_NoZero.png", plot = Plot_Heat_NoZero,
       scale = 1, width = 25, height = 20, units = "cm", dpi = 600)




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




# Instead of removing those proteins with even one value that is 0, here the 0 values are replaced with the lowest non-0 value after normalization

#### Imputation method 2: Replace 0 values with lowest non-0 value ####
## Generate a smooth histogram
# Generate a long variant of the median imputed data
DummyFrame = Function_makeLong(Data_Imput)
DummyFrame$value = log10(DummyFrame$value)
DummyFrame = as.data.frame(DummyFrame)

# Generate and print a smooth histogram of the data
Plot_SmoothHist_ImpRep = ggplot(DummyFrame, aes(x = value)) +
  geom_density(fill = "skyblue", color = "black", alpha = 0.5, bw = 0.15) +
  labs(x = "log10(intensity)", y = "Density", title = "Smooth Histogram of log10 protein intensity")
print(Plot_SmoothHist_ImpRep)

# Generate a data file of the smooth plot to extract the y-values
Data_SmoothHistPlot = ggplot_build(Plot_SmoothHist_ImpRep)$data[[1]]

# Determine the minimum density point  of the smooth plot
Value_Smooth_minDens = which.min(Data_SmoothHistPlot$y)

# Extract the corresponding x-value (log10(intensity)) and convert it back to intensity
Value_Smooth_min = Data_SmoothHistPlot$x[Value_Smooth_minDens]
Value_Smooth_min = 10^(Value_Smooth_min)

## Replace all zeros in the mean imputed dataset with the found lowest value
Data_ImpRep = Data_Imput
Data_ImpRep[Data_ImpRep == 0] = Value_Smooth_min

# Remove excessive files
rm(DummyFrame)
rm(Value_Smooth_min)
rm(Value_Smooth_minDens)



#### Histogram to show log2(intensity) distribution after normalization and mean imputation, and replacing 0-values with the lowest non-0 value ####
# Prepare the data
Data_Hist_ImpRep = Function_takeLog2(Data_ImpRep)
Data_Hist_ImpRep = Function_makeLong(Data_Hist_ImpRep)

## Draw the histogram
Plot_Hist_ImpRep = Function_drawHistogram(Data_Hist_ImpRep$value, title = "log2(Intensity) distribution after normalization, mean imputation and zero replacement")

## Print and save the histogram
print(Plot_Hist_ImpRep)
Function_savePlot(Plot_Hist_ImpRep, "Intensity_Distribution_ANorm_ImpRep_Full.png", plotType = "histogram")




#### Boxplot after normalization, mean imputation and replacing missing values ####
# Data preparation
Data_Box_ImpRep = Function_takeLog2(Data_ImpRep)
Data_Box_ImpRep = Function_makeLong(Data_Box_ImpRep)
Data_Box_ImpRep = Function_add_colorColumn(Data_Box_ImpRep)

# Draw the box plot
Plot_Box_ImpRep = Function_drawBoxplot(Data_Box_ImpRep, title = "log2 intensity distribution after normalization, mean imputation, zero replacement, and duplicate averageing")

## Print the plot and save it as a PNG
Function_savePlot(Plot_Box_ImpRep, "normalized_imput_Replace_Full.png", plotType = "boxplot")




#### Duplicate averageing after normalization, mean imputation, and replacing missing values ####
Accession_ImpRep = rownames(Data_ImpRep)
Data_AMean_ImpRep = Function_duplicateAverageing(Data_ImpRep, Accession_ImpRep)
rm(Accession_ImpRep)




#### Plot the post-imputation, zero replaced, averaged data in a boxplot ####
# Data preparation
Data_Box_AMean_ImpRep = Function_takeLog2(Data_ImpRep)
Data_Box_AMean_ImpRep = Function_makeLong(Data_Box_AMean_ImpRep)
Data_Box_AMean_ImpRep = Function_add_colorColumn(Data_Box_AMean_ImpRep)

# Draw the box plot
Plot_Box_AMean_ImpRep = Function_drawBoxplot(Data_Box_AMean_ImpRep, title = "log2 intensity distribution after normalization, duplicate averageing, mean imputation, and zero replacement")

## Print the plot and save it as a PNG
Function_savePlot(Plot_Box_AMean_ImpRep, "Averaged_ImpRep_Full.png", plotType = "boxplot")




#### Differential expression analysis with a 2-way ANOVA between controls and treated samples ####
# Note: This section tends to give an error if the entire code is launched in one go. Just execute it again and it should work.
Output_ANOVA_Rep = Function_performANOVA(Data_AMean_ImpRep, p_value = 0.01)




#### Heatmap generation ####
## Perform the heatmap generation
Plot_Heat_Rep = Function_drawHeatmap(Data_AMean_ImpRep, Output_ANOVA_Rep)

## Save the heatmap as a PNG
ggsave("Heatmap_Full_Rep.png", plot = Plot_Heat_Rep,
       scale = 1, width = 25, height = 20, units = "cm", dpi = 600)




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