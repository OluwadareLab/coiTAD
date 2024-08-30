    % Citation: Rocha P.P., Raviram R., Bonneau R., et al. Breaking TADs: insights into hierarchical genome organization [J] Epigenomics. 2015; 7(4): 523–526.
diary off;
nn = ['_', num2str(radius)];
foldname_1 = [Resultpath, 'TADs'];
if ~exist(foldname_1, 'dir')
    mkdir(foldname_1);
end
out_path = [foldname_1, '/'];
filename = [out_path, algorithm, nn, 'TADReport_.txt'];
diary(filename)
diary on;
fprintf('TAD for %s CLUSTERING\n', algorithm);

n = size(Chr_Data, 1); % specify the size of the matrix to display
max_tad_size = 2000000; % Set the maximum TAD size

try
    [newB, tad_sizes] = F_TAD(Chr_Data(1:n, 1:n), Assign_Cluster, algorithm, nn, out_path, Res, max_tad_size);
catch ME
    fprintf('Error occurred: %s\n', ME.message);
    diary off;
    return;
end

if ~isempty(newB)
    avg_size = (sum(newB(:, 2)) - sum(newB(:, 1))) / size(newB, 1);
    fprintf('Average size = %f\n', avg_size);
    Q = [length(newB(:, 1)), avg_size]; % Number of TD, Average size of TD
else
    Q = [0 0]; % Number of TD, Average size of TD
end
diary off;
close;

function [newB, tad_sizes] = F_TAD(Chr_Data, Assign_Cluster, algorithm, nn, Dir, Res, max_tad_size)
    count = 1;
    Border = [];
    start = 1;
    Limit = length(Assign_Cluster);
    
    for i = 2:Limit
        if (Assign_Cluster(i) ~= Assign_Cluster(start))
            Border = [Border; [start, i - 1]];
            start = i;
            count = count + 1;
        end
    end
    Border = [Border; [start, Limit]];
    
    min_num_tad = round(100000 / Res);
    newB = [];
    tad_sizes = []; % Array to store TAD sizes
    
    for j = 1:size(Border, 1)
        % Ensure indices are within bounds
        if Border(j, 2) <= size(Chr_Data, 1) && Border(j, 1) <= size(Chr_Data, 1)
            if ((Border(j, 2) - Border(j, 1) + 1) >= min_num_tad)
                tad_size = (Border(j, 2) - Border(j, 1) + 1) * Res; % Calculate TAD size
                if tad_size > max_tad_size
                    sub_tads = break_down_TAD(Chr_Data(Border(j, 1):Border(j, 2), Border(j, 1):Border(j, 2)), Res, max_tad_size);
                    for k = 1:size(sub_tads, 1)
                        new_tads = sub_tads(k, :) + Border(j, 1) - 1;
                        if new_tads(1) ~= new_tads(2)
                            newB = [newB; new_tads];
                            tad_sizes = [tad_sizes; (sub_tads(k, 2) - sub_tads(k, 1) + 1) * Res];
                        end
                    end
                else
                    if Border(j, 1) ~= Border(j, 2)
                        newB = [newB; Border(j, :)];
                        tad_sizes = [tad_sizes; tad_size]; % Store TAD size
                    end
                end
            end
        else
            fprintf('Skipping Border index out of bounds: [%d, %d]\n', Border(j, 1), Border(j, 2));
        end
    end
    
    disp('List the border');
    
    zeroRows = [];
    for i = 1:size(Chr_Data, 1)
        if isequal(Chr_Data(i, :), zeros(1, size(Chr_Data, 2)))
            zeroRows = [zeroRows; i];
        end
    end
    
    if ~isempty(newB)
        Redefine = [];
        newB_domain = [];
        B_end = newB(:, 2);
        
        for i = 1:length(B_end)
            if (sum(zeroRows == B_end(i)) == 0)
                X = [newB(i, :) 1]; % Recognized as domain = 1
                newB_domain = [newB_domain; newB(i, :)];
            else
                X = [newB(i, :) 0]; % Recognized as void = 0
            end
            Redefine = [Redefine; X];
        end
        
        newB = newB_domain;
        disp(Redefine);
        name = strcat(Dir, algorithm, nn, '_TAD_BinID.txt');
        dlmwrite(name, newB);

        OutputTAD; % Write domain to file

        fprintf('The Number of TAD = %d\n', length(newB(:, 1)));
    end
end

function new_tads = break_down_TAD(TAD_data, Res, max_tad_size)
    MIN_TAD_SIZE = 100000;
    MAX_TAD_SIZE = max_tad_size;
    MIN_BINS = ceil(MIN_TAD_SIZE / Res);
    MAX_BINS = floor(MAX_TAD_SIZE / Res);
    
    new_tads = [];
    tad_size = size(TAD_data, 1) * Res;
    
    if tad_size <= max_tad_size
        new_tads = [1, size(TAD_data, 1)];
        return;
    end
    
    num_bins = size(TAD_data, 1);
    max_bins = floor(max_tad_size / Res);
    
    start_bin = 1;
    while start_bin <= num_bins
        end_bin = min(start_bin + max_bins - 1, num_bins);
        if start_bin ~= end_bin
            new_tads = [new_tads; start_bin, end_bin];
        end
        start_bin = end_bin + 1;
    end
end
