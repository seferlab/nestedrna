clear all;
close all;

% Add helper functions
addpath('HelperFunctions');
addpath('gL_partitions');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CHOOSE CHROMOSOME DATA
chr1 = 10; %choose the first chromosome
chr2 = 10; %choose the second chromosome (for interchromsomal Hi-C) or repeat the first one (for intrachromosomal Hi-C)

chromosomes = chr1:-1:chr2;
chr_list = [int2str(chromosomes(1))];

for chr = chromosomes(2:end)
  chr_list = [chr_list '-' int2str(chr)];
end


% Load chrosomoes' Hi-C matrix KR normalized
data_path = ['/MATLAB Drive/mapping2023bernenko/HiC_data/chr' num2str(chr1) '_dna_00per.mat']
info_path = ['/MATLAB Drive/mapping2023bernenko/HiC_data/chr_' num2str(chr1) '.data.info']
save_coms_path = ['/MATLAB Drive/mapping2023bernenko/gL_partitions/']
data = load(data_path);
data = data.data;
%data = data.hicmap;
%data = double(data);
% Load chromosomes' length
info = load(info_path);

fprintf('Data has been loaded for chr%s\n', num2str(chr1))

% PLOT the HI-C map 
FIG = 1;
figure(FIG);
imshow(log(data), [min(min(log(data))) max(max(log(data)))]);
colormap(hot);
hold on;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GenLouvain community partition

% we set three parameters for the null-model:
% ---> alpha - the exponent of the decay of contact frequencies between node i and node j
% ---> gamma - the resolution paramater
% ---> distance matrix D


alpha = 1; % FIX one parameter

for gamma = 0.6% [0.4, 0.5, 0.6, 0.7, 0.8, 0.81, 0.82, 0.83, 0.84, 0.85, 0.86, 0.87, 0.88, 0.89, 0.9]        
    S = []; % to save community assignemnts for nodes

    for k = 1:size(info, 1)          

      chromosome = info(k, 3);
      
      if (k == 1) % if gL will be done on a single chromosome Hi-C data (the CASE)
        data_part = data(1:info(k, 1), 1:info(k, 2)); % gets a slice of the Hi-C matrix accordingly to the INFO datafile
        
        % calculate the distance matrix between node i and node j;
        % Dij = |i-j| -- this is how D matrix is calculated
        D = toeplitz(1-1:info(k,1)-1); % It makes the D matrix to have zeros on the main diagonal, all following diagonals are +1
      else
        data_part = data((info(k-1, 1) + 1):info(k, 1), (info(k-1, 2) + 1):info(k, 2));
        D = toeplitz(1:(info(k,1) - info(k-1, 1)));
      end

      % Calculate the modularity matrix
      B = fractalglobulepoly_f(data_part, D, gamma, alpha);
      
      % Get communities with genlouvain
      [S_p, Q] = iterated_genlouvain(B);


      % Save results into S matrix
      S = [S; S_p ones(length(S_p), 1) * chromosome ones(length(S_p), 1) * Q];          

      gamma_n = num2str(gamma *100);
      alpha_n = num2str(alpha *100);
      save([save_coms_path 'chr' num2str(chr1) '_gamma' gamma_n '_alpha' alpha_n '.mat'], 'S')
    end
end
