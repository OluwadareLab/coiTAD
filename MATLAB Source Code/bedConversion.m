% Define input and output file paths
inputFilePath = 'C:\GitHub Repos\coiTAD\data_Results\TADs\HDBSCAN_9_domain.txt'; % Change this to your input file path
outputFilePath = 'C:\GitHub Repos\coiTad\e.bed'; % Change this to your desired output file path

firstNumber = 3; 
secondNumber = 5; 

% Open input file for reading
fid = fopen(inputFilePath, 'r');
if fid == -1
    error('Unable to open input file');
end

% Open output file for writing
fout = fopen(outputFilePath, 'w');
if fout == -1
    error('Unable to open output file');
end

% Read lines from input file
line = fgetl(fid); % Read the first line
if ischar(line) && any(isletter(line)) % Check if the line contains any letters
    % Skip the line
else
    % Split line by spaces
    values = str2double(strsplit(line));
    
    % Check if line is valid (contains at least 4 numbers)
    if numel(values) >= 4
        % Extract 2nd and 4th numbers
        start = values(firstNumber);
        stop = values(secondNumber);
        
        % Write to output file in BED format (0-based, half-open)
        fprintf(fout, '%d\t%d\n', start, stop);
    end
end

% Read remaining lines from input file
while ~feof(fid)
    line = fgetl(fid);
    
    % Split line by spaces
    values = str2double(strsplit(line));
    
    % Check if line is valid (contains at least 4 numbers)
    if numel(values) >= 4
        % Extract 2nd and 4th numbers
        start = values(firstNumber);
        stop = values(secondNumber);
        
        % Write to output file in BED format (0-based, half-open)
        fprintf(fout, '%d\t%d\n', start, stop);
    end
end

% Close files
fclose(fid);
fclose(fout);

disp('Conversion completed successfully.');
