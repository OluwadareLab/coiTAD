%% Find TAD from the result

radius = min_radius;
Len = End - Start + 1;
%Assign_Cluster = zeros(N,Len);
TAD_Q=[]; % Number of TD, Average size of TD

for radiusIndex = 1:size(radiusOptimalClusters, 2)

    fprintf ('Result when radius = %d\n', radius);
    Agn = Order_TADNum(radiusOptimalClusters(:,radiusIndex));
    Assign_Cluster =  Agn';
    tad_assign=strcat(Resultpath,algorithm,'_',name,'_',num2str(radius),'.mat');

    %==========================================================================
    % Quality Assessment
    %--------------------------------------------------------------------------    
    %Run the TAD Assessment      
    Extract_TAD; % Extract TAD from clusters identified   
    TAD_Q = [TAD_Q; Q];    
    radius = radius + 1;
     

end
%--------------------------------------------------------------------------
% Write results to file
%--------------------------------------------------------------------------
foldname = [Resultpath,'/Quality'];
if ~exist(foldname, 'dir')
    %Folder does not exist so create it.
    mkdir(foldname);
end
out_path_1 = [foldname,'/'];


TAD_name = strcat(out_path_1,'TADQ_',name,'.csv');




function [ Assign]  = Order_TADNum(Found_TAD)
%% This Function number the TAD clusters found.
%--------------------------------------------------------------------------
len = length(Found_TAD); % Length of the current TAD
count = 1;
Assign(1) = count; 
for i = 2:len
    if(Found_TAD(i-1) == Found_TAD(i))
        count = count + 0;
    else
        count = count + 1;
    end
        Assign(i) = count;
end

end