%% FMR1 PAPER CODE SUMMARY 

%% NOTES ON SUPPORTING FUNCTIONS

% In order to produce the figures from the struct, there are very few functions that aren't built 
% into MATLAB that you will need. The list of additional functions and the functions themselves can 
% be found on in a subfolder.

%% GENERATING THE DATA AND REMAPPING STRUCTS

% Rats were added iteratively to the Data and Remapping structs throughout the experiment. The
% following functions have been modified for use if you need to create the data struct completely 
% from scratch.

% dataDir = directory where raw data is kept

Data = a_Put_all_data_in_struct(dataDir); %adds rats, days, unit names/spike times to Data struct and filters LFPs
Data = b_Loop_ratemap_stacks(Data, dataDir); %adds ratemaps to Data struct, makes figures of ratemaps saved in raw data folders
c_Check_arena_coverage(Data); %this function will go through each day and note if the arena coverage threshold is not met
Remapping = d_Remapping_from_struct_byRat(Data, dataDir); %makes a struct with remapping data (rate overlap, spatial correlation, PV corr), separated by rat
% function e makes an altermate Remapping struct that combines rats together based on genotype, not used
Remapping = f_CombineCA2_CA2_3remappingStruct(Remapping); %previous work includes CA2/3 border cells in CA2 experiments
Data = g_add_exploration_time(Data); %this function adds time per bin to the struct. Also makes the heat maps used in Figure 5A-B
[~, olfTest] = h_Put_data_in_olfactory_struct; %makes the olfactory data struct. Data is hard coded.

%% FIGURE 2: EXAMPLE RATEMAPS

 fmr1CA2_example_ratemaps(Data);

%  This function will make ratemaps for all days and all cells.

%% FIGURE 3: SPATIAL CORRELATION COEFFICIENTS ACROSS GENOTYPE AND CONDITION

% Parts B-F of this figure must be made first. This function uses BOTH the Data and Remapping
% structs as inputs.

fmr1CA2_global_remapping(Remapping, Data);

% This function will produce modified dotplots where the values for each rat have different markers.
% Each session combination will be a separate figure. Each figure will have a legend to note which
% markers represent data points from which rats.

% Importantly, this function will keyboard pause before creating the figures, which will allow you
% to create the SPSS Data file. This can be done by copying the pasting the "statAll" variable
% into SPSS. The SPSS Data file has been created with the columns correctly labeled by variable
% name and type and put on the server. 

% The estimated mean values for each genotype and condition have already been hard coded into a
% function to produce Figure 4A. This function does not require any inputs.

fmr1CA2_estMeans_remapping;

%% FIGURE S2: RATE OVERLAP VALUES BY GENOTYPE AND CONDITION

% An important note for the rate overlap analysis: when analyzing data for the CA2 WT paper, we
% decided NOT to use the rate overlap calculation from Colgin et al. 2011, which used the mean 
% firing rate of place cells across the entire ratemap. Instead, we used the average in-field firing
% rate to calculate rate overlap. We also discussed using the peak firing rate of the cell.
% Ultimately, the results did not differ based on which method we used, and we reported data using
% the average in-field firing rate in the paper. The code that follows will create figures for both
% methods; however, only results using the average in-field firing rate method were used.
% 
% The process for making this figure is almost identical to the previous. This function only uses
% the Data struct as an imput. 

fmr1CA2_rate_remapping(Data);

% This will create parts B-F of this figure.
% This function will also create the "statAll"  variable, which can be copied and pasted into SPSS. 
% The estimated means are already hard coded into the function to produce part A of this figure.

fmr1CA2_estMeans_remapping;

%% FIGURE 4: PLACE CELL PROPERTIES BY GENOTYPE

fmr1CA2_placeCellProperties(Data);

% This function makes the figures for place field properties by genotype, condition, and session.
% This function also creates the variables for statistical analysis in SPSS. Note: some measures
% apply to every cell, regardless of if there was a place field identified. For example, spatial
% information was calculated for every cell. However, metrics like in-field firing rate could only
% be calculated wehna place field was identified. That is why there are two separate stats files,
% statCell and statField. As of 4/15/2025, out-of-field firing rate and in-field/out-of-field firing
% rate ratio have been added as figures and to the stats files. 

%% FIGURE 5: EXPLORATION TIME

% The heat maps for parts A-B of this figure are created by function 'g_get_exploration_time' (see
% GENERATING THE DATA AND REMAPPING STRUCTS section above). For Meg's dissertation, the WT heatmaps
% were pulled directly from the CA2 WT paper and the FXS heat maps were created.

fmr1CA2_quanitfyTPB(Data);

% This function will create part C of this figure. This function also creates statRat, which can
% be plugged into SPSS to generate stats for this figure. 


%% FIGURE 6: OLFACTORY TESTING

fmr1CA2_olfactoryTests(OlfTest);

% This function will create plots for both the neutral odor and social odor olfactory tests. This 
% function will make separate figures for each odor and will produce bar plots with two alternate 
% ways of displaying variability. The one used in Meg's dissertation was the social odor figure with
% error bars displaying the 95% confidence intervals. There are also plots with dots representing 
% the values from individual rats. 
% 
% This function also creates the variable "statSocial", which can be inputted into SPSS to run 
% statistics for the social odor data.

%% FIGURE 7: OXTR-2 EXPRESSION

fmr1CA2_estMeans_OXTR2;

% This function uses hard coded values from the stats that Laura ran on the OXTR-2 expression data 
% to make Figure 7C.

%% STATISTICAL ANALYSIS TO ACCOUNT FOR LOWER CELL NUMBERS IN FXS GROUP

%% METHOD 1

% In my dissertation, I described two ways I attempted to account for the lower cell yields in the
% FXS group. I've organized them into separate codes here that hopefully are clear with what I did
% and why. To use these codes, you will need the "statAll" output from the fmr1CA2_global_remapping
% function. Methods 1 requires going back and forth between MATLAB AND SPSS. 

% Method 1 randomly down-samples the data from WT rats to the same cell numbers per condition of FXS
% rats. The goal of this analysis was to determine if the size of the FXS dataset would be
% sufficient to detect differences in spatial correlation coefficient values between conditions.

load('spatialCorrelation_MATLABdata.mat')
statDs = fmr1CA2_downsample_WTdata(statAll); %this function will create 1000 randomly down-sampled datasets

% I then used SPSS to run a generalized linear mixed model analysis on the first 100 datasets,
% exactly as we had done in the original CA2 WT paper (Robson et al. 2025). I saved the p-values
% from the main effect on condition and the post-hoc Bonferroni-corrected comparisons between empty
% and odor and empty and visual + odor in an Excel file. 

tmpP = readmatrix('spatialCorrelation_condnWT_downsample_results.xlsx'); %load the p-values into MATLAB
tmpPow = length(find(tmpP(:,2)<0.05 & tmpP(:,3)<0.05)) ./ length(tmpP);
fprintf('Odor condition power = %s\n', num2str(tmpPow))

tmpPow = length(find(tmpP(:,2)<0.05 & tmpP(:,4)<0.05)) ./ length(tmpP);
fprintf('Visual + Odor condition power = %s\n', num2str(tmpPow))

%% METHOD 2

% Method 2 uses a permutation test. I did this separately for each genotype. This test considers the 
% observed difference of the mean of two samples (for example, the difference in the mean saptial
% correlation coefficient values between the Empty condition and the Odor condition). It then
% creates a number of shuffled distributions where condition is randomly shuffled across the spatial
% correlation coefficient population. For each shuffled dataset, a new mean difference between
% shuffled-Empty and shuffled-Odor is calculated. To determine significance, the p-value is
% calcualted using the formula: (Nsubset + 1) / Nshuffle + 1), where Nsubset is the number of
% shuffles where the mean difference was greater than or eaual to the observed valye, since we are
% using a two-tailed test. Nshuffle is the number of shuffles we are using.

 fmr1CA2_remapping_permutation(statAll);



%% Created 4/15/2025 by MMD
