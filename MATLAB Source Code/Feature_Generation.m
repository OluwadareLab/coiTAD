filename = 'C:\GitHub Repos\coiTAD';
cd(filename);

contactMatrix = dlmread("4noise.hic");
contactMatrix; 


higherRadius = (Max_TADsize/Res) + 10; 
lowerRadius = 2; 

folder = 'featuresGenerated'; 
if ~exist(folder, 'dir')
    % Folder does not exist so create it.
    mkdir(folder);
end


for radius = lowerRadius:higherRadius 

    clusteringInput = createEntireFeature(contactMatrix, radius);  

    %for pca feature analysis
    %[coeff, score, ~, ~, explained] = pca(clusteringInput);
    %num_components_to_keep = find(cumsum(explained) >= 85, 1, 'first');
    %reduced_features = score(:, 1:num_components_to_keep);
    
    
    %for raw feature analysis 
    reduced_features = clusteringInput; 
   

    % Create a filename for the current matrix based on the radius
    filename = sprintf('feature_radius_%d.txt', radius);

    filePath = fullfile(folder, filename);
    
    % Open the file for writing
    fileID = fopen(filePath, 'w');
    
    % Write the matrix data to the file
    for row = 1:size(reduced_features, 1)
        fprintf(fileID, '%f ', reduced_features(row, :));
        fprintf(fileID, '\n');
    end
    
    % Close the file
    fclose(fileID);

end 

function appendedMatrix = fillFinalMatrix(currentFinalMatrix, contactMatrix, startingPoint, radius)

    matrixSize = size(contactMatrix);
    numRows = matrixSize(1);
    numCols = matrixSize(2);

    operator = 2;

    c = contactMatrix(startingPoint(1), startingPoint(2));

    if (startingPoint(1) - 1) <= 0

        tC = [0];

    else 

        tC = [contactMatrix(startingPoint(1) - 1, startingPoint(2))];

    end 


    if (startingPoint(1) - 1) <= 0 || (startingPoint(2) + 1) > numCols

        tR = [0];

    else

        tR = [contactMatrix(startingPoint(1) - 1, startingPoint(2) + 1)];

    end


    if (startingPoint(2) + 1) > numCols

        r = [0];

    else

    r = [contactMatrix(startingPoint(1), startingPoint(2) + 1)];

    end

    if (startingPoint(1) + 1) > numRows || (startingPoint(2) + 1) > numCols

        bR = [0]; 

    else

        bR = [contactMatrix(startingPoint(1) + 1, startingPoint(2) + 1)];

    end

    if (startingPoint(1) + 1) > numRows

        b = [0];
   
    else

        b = [contactMatrix(startingPoint(1) + 1, startingPoint(2))];

    end

    if (startingPoint(1) + 1) > numRows || (startingPoint(2) - 1) <= 0

        bL = [0];

    else

        bL = [contactMatrix(startingPoint(1) + 1, startingPoint(2) - 1)];

    end

    if (startingPoint(2) - 1) <= 0

        l = [0];
    else

        l = [contactMatrix(startingPoint(1), startingPoint(2) - 1)];

    end

    if (startingPoint(1) - 1) <= 0 || (startingPoint(2) - 1) <= 0

        tL = [0];

    else

        tL = [contactMatrix(startingPoint(1) - 1, startingPoint(2) - 1)];

    end


    for index = 1:(radius-2)
        
        if (startingPoint(1) - operator) <= 0

            tC = [tC 0;];

        else 

            tC = [tC contactMatrix(startingPoint(1) - operator, startingPoint(2));];

        end 
    
        if (startingPoint(1) - operator) <= 0 || (startingPoint(2) + operator) > numCols

            tR = [tR 0;];

        else

            tR = [tR contactMatrix(startingPoint(1) - operator, startingPoint(2) + operator);];

        end
    

        if (startingPoint(2) + operator) > numCols

            r = [r 0];

        else

            r = [r contactMatrix(startingPoint(1), startingPoint(2) + operator)];

        end

        if (startingPoint(1) + operator) > numRows || (startingPoint(2) + operator) > numCols

            bR = [bR 0]; 

        else

            bR = [bR contactMatrix(startingPoint(1) + operator, startingPoint(2) + operator)];

        end

        if (startingPoint(1) + operator) > numRows

            b = [b 0];
   
        else

            b = [b contactMatrix(startingPoint(1) + operator, startingPoint(2))];

        end

        if (startingPoint(1) + operator) > numRows || (startingPoint(2) - operator) <= 0

            bL = [bL 0];

        else

            bL = [bL contactMatrix(startingPoint(1) + operator, startingPoint(2) - operator)];

        end
    
        if (startingPoint(2) - operator) <= 0

            l = [l 0];

        else

            l = [l contactMatrix(startingPoint(1), startingPoint(2) - operator)];

        end

        if (startingPoint(1) - operator) <= 0 || (startingPoint(2) - operator) <= 0

            tL = [tL 0];

        else

            tL = [tL contactMatrix(startingPoint(1) - operator, startingPoint(2) - operator)];

        end
    
        operator = operator + 1;

    end

    %for full circle feature
    appendedMatrix = [tL tC tR r c l bL b bR];

    %for semi circle feature
    %appendedMatrix = [tL tC tR r c bR];

    
    appendedMatrix = [currentFinalMatrix; appendedMatrix];

end

%for seperate feature order
% function appendedMatrix = fillFinalMatrix(currentFinalMatrix, contactMatrix, startingPoint, radius)
% 
%     matrixSize = size(contactMatrix);
%     numRows = matrixSize(1);
%     numCols = matrixSize(2);
% 
%     operator = 2;
% 
%     % Initialize vectors
%     center = contactMatrix(startingPoint(1), startingPoint(2));
%     topLeft = [];
%     topCenter = [];
%     topRight = [];
%     right = [];
%     left = [];
%     bottomLeft = [];
%     bottomCenter = [];
%     bottomRight = [];
% 
%     % Extract values for each direction based on the radius
%     for r = 1:radius
%         % Top left
%         if (startingPoint(1) - r) > 0 && (startingPoint(2) - r) > 0
%             topLeft = [topLeft, contactMatrix(startingPoint(1) - r, startingPoint(2) - r)];
%         else
%             topLeft = [topLeft, 0];
%         end
% 
%         % Top center
%         if (startingPoint(1) - r) > 0
%             topCenter = [topCenter, contactMatrix(startingPoint(1) - r, startingPoint(2))];
%         else
%             topCenter = [topCenter, 0];
%         end
% 
%         % Top right
%         if (startingPoint(1) - r) > 0 && (startingPoint(2) + r) <= numCols
%             topRight = [topRight, contactMatrix(startingPoint(1) - r, startingPoint(2) + r)];
%         else
%             topRight = [topRight, 0];
%         end
% 
%         % Right
%         if (startingPoint(2) + r) <= numCols
%             right = [right, contactMatrix(startingPoint(1), startingPoint(2) + r)];
%         else
%             right = [right, 0];
%         end
% 
%         % Left
%         if (startingPoint(2) - r) > 0
%             left = [left, contactMatrix(startingPoint(1), startingPoint(2) - r)];
%         else
%             left = [left, 0];
%         end
% 
%         % Bottom left
%         if (startingPoint(1) + r) <= numRows && (startingPoint(2) - r) > 0
%             bottomLeft = [bottomLeft, contactMatrix(startingPoint(1) + r, startingPoint(2) - r)];
%         else
%             bottomLeft = [bottomLeft, 0];
%         end
% 
%         % Bottom center
%         if (startingPoint(1) + r) <= numRows
%             bottomCenter = [bottomCenter, contactMatrix(startingPoint(1) + r, startingPoint(2))];
%         else
%             bottomCenter = [bottomCenter, 0];
%         end
% 
%         % Bottom right
%         if (startingPoint(1) + r) <= numRows && (startingPoint(2) + r) <= numCols
%             bottomRight = [bottomRight, contactMatrix(startingPoint(1) + r, startingPoint(2) + r)];
%         else
%             bottomRight = [bottomRight, 0];
%         end
%     end
% 
%     % full circle feature
%     %appendedMatrix = [center, topLeft, topCenter, topRight, right, bottomRight, bottomCenter, bottomLeft, left];
% 
%     %semi-circle option
%     appendedMatrix = [center, topLeft, topCenter, topRight, right, bottomRight]
% 
%     % Append to the current matrix of features
%     appendedMatrix = [currentFinalMatrix; appendedMatrix];
% 
% end



function clusteringInput = createEntireFeature(contactMatrix, radius)

    matrixSize  = size(contactMatrix);
    numRows = matrixSize(1);
    numCols = matrixSize(2); 

    startingPoint = [1 1];

    clusteringInput = [];
 
    if numCols > numRows 
        for index = 1: numRows 

            if startingPoint(1) <= numRows && startingPoint(2) <= numRows

                clusteringInput = fillFinalMatrix(clusteringInput, contactMatrix, startingPoint, radius);
                startingPoint(1) = startingPoint(1) + 1; 
                startingPoint(2) = startingPoint(2) + 1;
                

            end

        end 

    elseif  numRows > numCols 

        for index = 1: numCols 
            
            if startingPoint(1) <= numCols && startingPoint(2) <= numCols

            clusteringInput = fillFinalMatrix(clusteringInput, contactMatrix, startingPoint, radius);
            startingPoint(1) = startingPoint(1) + 1; 
            startingPoint(2) = startingPoint(2) + 1;
            

            end

        end

    elseif numCols == numRows 
        for index = 1: numRows
            
            if startingPoint(1) <= numRows && startingPoint(2) <= numRows

                clusteringInput = fillFinalMatrix(clusteringInput, contactMatrix, startingPoint, radius);
                startingPoint(1) = startingPoint(1) + 1; 
                startingPoint(2) = startingPoint(2) + 1;
                 

            end

       end 
    end

end
