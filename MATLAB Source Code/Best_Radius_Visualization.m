%Script to visualize the best radius generated from coiTAD 

radius = min_radius; 

for radiusIndex = 1:size(radiusOptimalClusters, 2) 
        
    fprintf("Current radius is %d\n", radius)
    
    if radius == bestRadius 

        fprintf ('Result when radius = %d\n', radius);
        Agn = Order_TADNum(radiusOptimalClusters(:,radiusIndex));
        Assign_Cluster =  Agn';
        tad_assign=strcat(Resultpath,algorithm,'_',name,'_',num2str(radius),'.mat');
    
        %==========================================================================
        % Quality Assessment
        %--------------------------------------------------------------------------    
        %Run the TAD Assessment      
        nn =  ['_',num2str(radius)];
        foldname_1 = [Resultpath,'TADs'];
        if ~exist(foldname_1, 'dir')
            % Folder does not exist so create it.
            mkdir(foldname_1);
        end

        out_path = [foldname_1,'/'];
        filename = [out_path,algorithm,nn,'TADReport_.txt'];
        fprintf('TAD for %s CLUSTERING\n',algorithm);
        n = size(Chr_Data,1); % specify the size of matrix to display :::default
        % n = 100;   
        Visualize(log(Chr_Data(1:n,1:n)));
        drawnow limitrate       
        newB = F_TAD(Chr_Data(1:n,1:n), Assign_Cluster, algorithm, nn,out_path,Res );
        if ~isempty(newB)
            avg_size = (sum(newB(:,2))- sum(newB(:,1)))/size(newB,1);
            fprintf('Average size = %f\n',avg_size);         
            Q = [length(newB(:,1)) avg_size]; % Number of TD, Average size of TD
            title_text = sprintf('TD for %s Implementation',algorithm);
            title(title_text)
       
        end
                    
    break
    end

    radius = radius + 1; 
    

end



function [ output_args ] = Visualize( Data )
% This function allows for the visualization of data in Heat Map
colormap('hot');
imagesc(Data);
colorbar;
xlabel('Genomic bin (resolution: 40kb)');
ylabel('Genomic bin (resolution: 40kb)');
end


function [newB] = F_TAD(Chr_Data, Assign_Cluster, algorithm,nn,Dir,Res )
count = 1;
Border = [];
start = 1;
Limit = length(Assign_Cluster);
for i = 2:Limit
    if(Assign_Cluster(i) ~= Assign_Cluster(start))
        Border = [Border;[start,i-1]];
        start = i;
        count = count + 1;
    end
end
Border = [Border;[start,Limit]];

%when the number of bins×resolution>100kb.  
%Research has shown that the length of TAD is between 100 kb and 5 Mb 
%citation. Rocha P.P., Raviram R., Bonneau R., et al. Breaking TADs: insights into hierarchical genome organization[J] Epigenomics. 2015;7(4):523–526.
min_num_tad = round(100000/Res);
 newB = [];
    for j = 1:length(Border(:,1))
        if ((Border(j,2)- Border(j,1) + 1) > min_num_tad)
          newB = [newB ; Border(j,:) ] ;    
        end      
    end

disp('List the border');
%===========================================================================
% Find the gaps classified as domain
% load('Chr_Data.mat')
zeroRows = [];
for i = 1:size(Chr_Data,1)
    if isequal(Chr_Data(i,:),zeros(1,size(Chr_Data,2)))
        zeroRows = [zeroRows; i];    
    end
end
%===========================================================================

if (isempty(newB) == 0)
     
    Redefine = [];
    newB_domain = [];
    B_end = newB(:,2);
    for i = 1:length(B_end)
     if(sum(zeroRows==B_end(i))==0)
         X = [newB(i,:)  1];   % Recognized as domain = 1
         newB_domain = [newB_domain ; newB(i,:)];
     else
         X = [newB(i,:)  0];   %  Recognized as void  = 0
     end
         Redefine = [Redefine; X];
    end
    
    newB = newB_domain;
    disp(Redefine);
    name = strcat(Dir,algorithm,nn,'_TAD_BinID.txt');
    

    dlmwrite(name,newB); 
    
    fname = strcat(Dir,algorithm,nn,'_domain.bed');
    fileID = fopen(fname,'w');
    
    List = [];
    for i = 1:size(newB, 1)
        %List =[List , newB(i,1)*Res, newB(i,2)*Res ];
        List = [List, newB(i,1), newB(i,2)]; 

    end
        
    
    fprintf(fileID,'%12d %12d\n',  List);
    fclose(fileID);
    
    for i = 1:length(newB(:,1))
         hold on;
         Start = newB(i,1);
         Last = newB(i,2);
         for j = Start:Last
                 plot(Start,j,'b.','MarkerSize',10);
                 plot(j,Start,'b.','MarkerSize',10);
                 plot(j,Last,'b.','MarkerSize',10);
                 plot(Last,j,'b.','MarkerSize',10);
                 drawnow
         end        
    end
    fprintf('The Number of TAD = %d\n', length(newB(:,1)));
 end
end





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
