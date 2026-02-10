%% FMM test data 


clc
clear

% Define the sample sizes
sample_sizes = [1e4, 2e4, 3e4, 4e4, 5e4];

% Preallocate a cell array to hold the test data
bb_data_1 = cell(length(sample_sizes), 1);
source_data_1 = cell(length(sample_sizes), 1);

% Generate random numbers in (0.5, 1) for each sample size
for k = 1:length(sample_sizes)
    n = sample_sizes(k);
    bb_data_1{k} = 0.5 + 0.5 * rand(n, 1);  % n×1 vector of random numbers in (0.5, 1)
    source_data_1{k} = sort(rand(n,1));
end
















