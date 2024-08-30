
%%========================================================================
 % Implementation of coiTAD algorithm
 %  Language: MATLAB
	
 %	@author: Oluwatosin Oluwadare/Drew Houchens  
 %	UCCS Bioinformatics Lab 
 %	University of Colorado Colorado Springs 
 %	United States 
 %	USA	
 %	email: ooluwada@uccs.edu/dhouchen@uccs.edu
	
 %	Last Update:8/23/2024
	
 %	Instructions: https://github.com/BDM-Lab/coiTAD


  % Set your working directory 
 cd('C:\GitHub Repos\coiTad');
 

clear;
warning('off','all');
%  =============================================================================================
%   coiTAD  Variables
%  =============================================================================================

%	change variable "filepath" to the Input data file path	
%	change variable "name" to the Input data name
%	change variable "chromo" to the chromosome data name
%	change variable "Reso" to the Input data Resolution		
%	specify the maximum TAD size : change variable "Max_TADsize" (optional)

%  ============================================================================================
MaxQuality = 0;			   % Get the Maximum TAD Quality

filepath='C:\GitHub Repos\coiTad\';          % filepath

featurefilepath = 'C:\GitHub Repos\coiTad\featuresGenerated\';  %feature file path 

name= '4noise.hic';          % filename

chromo = 'chr';             % chromosome name

Res = 40000;                  % 40000 = 40KB , 100000 = 100KB

KB = 1000;                    % KB = 1000 :: Constant

Option = 0;                % Option: if option == 0, use rough estimate and if option == 1, use elbow method

Max_TADsize = 800000;          % Maximum TAD size: 800KB

outputfolder_name = 'data_Results';

algorithm = 'HDBSCAN';          % algorithm used code: KM == Kmeans Algorithm
%  ============================================================================================
Resolution=[num2str(Res/KB),'kb'];
filename = strcat(filepath,name);
Chr_Data = dlmread(filename); 

fprintf('Data set loaded.\n');

[path,name,ext] = fileparts(name);
%==========================================================================
% Make directory if it doesn't exist
%--------------------------------------------------------------------------
foldname = [outputfolder_name];
if ~exist(foldname, 'dir')
    % Folder does not exist so create it.
    mkdir(foldname);
end

Resultpath = [foldname,'/'];
N = length(Chr_Data);

% Rough estimate of clusters
M = floor(sqrt(N/2));

Limit = length(Chr_Data);
k_opt =  M;
Start = k_opt- 5;
End = k_opt + 5;
    

row = size(Chr_Data,1); 

leng = End - Start + 1;

clust = zeros(row,leng);


Limit = length(Chr_Data);


%Radius calculations 
radiusOptimalClusters = []; 

min_radius = 2; 
max_radius = (Max_TADsize/Res) + 10;

Feature_Generation; 

cd (featurefilepath); 

  
% Initialize an empty matrix for storing cluster labels
radiusOptimalClusters = [];

% Determine the maximum length of the cluster arrays
maxLength = 0;

% First pass to determine the maximum length
for radius = min_radius:max_radius 

    % Construct the file name based on the current radius value
    file_name = sprintf('feature_radius_%d.txt', radius);
    radius_data = load(file_name);

    % Run HDBSCAN clustering
    currentRadiusObject = HDBSCAN(radius_data); 
    currentRadiusObject.run_hdbscan();
    clusters = currentRadiusObject.labels;

    % Ensure clusters is a column vector
    clusters = clusters(:);

    % Update the maximum length
    if length(clusters) > maxLength
        maxLength = length(clusters);
    end
end

% Initialize the matrix with zeros based on the maximum length
radiusOptimalClusters = zeros(maxLength, max_radius - min_radius + 1);

% Second pass to store the cluster labels in the matrix
for radius = min_radius:max_radius 
    
    % Construct the file name based on the current radius value
    file_name = sprintf('feature_radius_%d.txt', radius);
    radius_data = load(file_name);

    % Run HDBSCAN clustering
    currentRadiusObject = HDBSCAN(radius_data); 
    currentRadiusObject.run_hdbscan();
    clusters = currentRadiusObject.labels;

    % Ensure clusters is a column vector
    clusters = clusters(:);

    % Store clusters in the matrix, padding with zeros if necessary
    radiusOptimalClusters(1:length(clusters), radius - min_radius + 1) = clusters;
end

% Continue with further processing...




cd(filepath);





% gscatter(Chr_Data(:,1),Chr_Data(:,2),graphedClusters,[],"xod")




disp('======================== Quality Assessment==========================');

%% Perform Quality Assessment
Process_Cluster;
Quality_Check;

disp('===================== Quality Assessment Completed==========================');


disp('find the results of the Quality Assessment in the Quality/ directory');
fprintf("Recommended radius = %d\n", bestRadius); 
disp('Find the TADs identified in the TAD/ directory');
disp('===================== coiTAD Completed =============================');



%Best_Radius_Visualization; 
 
