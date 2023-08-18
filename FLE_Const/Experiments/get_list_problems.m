function [problems, sizes] = get_list_problems(files)
% Retrieves a list of problems and sizes from input files.

% Inputs:
%   - files: A cell array of file names.

% Outputs:
%   - problems: A column vector containing the list of problems.
%   - sizes: A column vector containing the corresponding sizes.

n = length(files);
problems = [];
sizes = [];

for i = 1:n
    % Read the table from the file
    table = readtable(files(i));
    
    % Extract the problem names and sizes from the table
    list_problems = table2array(table(:, 1));
    problem_sizes = table2array(table(:, 2));
    
    % Append the problem names and sizes to the respective lists
    problems = [problems; list_problems];
    sizes = [sizes; problem_sizes];
end

