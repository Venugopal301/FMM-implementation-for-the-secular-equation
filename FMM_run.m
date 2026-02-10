%% FMM implementation for secular equation


addpath(genpath('./'))


% load   FMM_test_data.mat file     from test_data file here
load("D:\Venugopal_files\FMM_upload\test_data\FMM_test_data.mat")

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Initialization

s = 6.311887342690112;
iter = 10;  % newton iterations
elapsed_time = zeros(iter,2,5);
new_sum_fmm = cell(5,1);
new_sum_dir = cell(5,1);

targets_fmm_conv =  cell(5,1);
targets_dir_conv = cell(5,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

for jj = 1:5
jj

bb = bb_data_1{jj};
sources = source_data_1{jj};   


% predefine no of levels and number of expansion coefficients
J = 4;  % predefine no of levels
pp = 30; % no of expansion coefficients  % Please vary number of expansion coefficients and check the results

%%%%%%%%%%%%%%%%% build tree structure with prescribed levels
tree = struct('level', {}, 'interval', {}, 'center', {}, 'parent', {}, 'child', {}, 'neighbor', {},...
    'far_ilist', {}, 'local_ilist', {}, 'sources', {}, 'targets',{}, 'x_coeff', {},'x2_coeff', {},...
    'y_coeff', {}, 'blk_2',{}, 'blk', {}, 'der_blk',{}, 'der_blk_2', {},  'near_sum', {}, 'der_near_sum', {});

%%%%%%%%%%%%%%%%%  updating root node 1 at level 0
a = floor(sources(1));  b = ceil(sources(end));
tree(1).interval = [a,b];
tree(1).center = mean([a,b]);
tree(1).parent = [];
tree(1).child = [2,3];  % child nodes for the root node
tree(1).sources = sources;  %
% tree(1).targets = targets; %
n = length(sources);

% recursively build the basic tree structure with update of levels, interval, center, sources, parent, child
for j = 1:J
    for i = 1:2^j
        k = 2^j+(i-1);
        tree(k) = tree_update_1(tree(floor(0.5*k)),k,pp);
    end
end

for j = 2:J
    for k = 2^j:2^(j+1)-1
        tree(k).neighbor = assignneighbor(k);
    end
end
%     assigning interaction list for level 2 nodes
tree(4).far_ilist = [6,7];
tree(5).far_ilist = 7;
tree(6).far_ilist = 4;
tree(7).far_ilist = [4,5];
%

for i = J:-1:3  % assigning interaction list for nodes from bottom up
    tree = update_ilist(tree,i);
end
%
for j = J       % updating multipole expansion coefficients for bottom level
    for k = 2^j:2^(j+1)-1
        tree(k).x_coeff = update_x_coeff(tree,k,pp,sources,bb);
        tree(k).x2_coeff = update_x2_coeff(tree,k,pp,sources,bb);
    end
end

for j = J-1:-1:2   % M2M  Multipole-to-Multipole translation bottom to top
    for k = 2^j:2^(j+1)-1
        tree(k).x_coeff = M2M(tree,k,pp);

    end
end

for j = J-1:-1:2   % M2M2  Multipole-to-Multipole translation bottom to top for the derivative
    for k = 2^j:2^(j+1)-1
        tree(k).x2_coeff = M2M2(tree,k,pp);
    end
end


for k =  2^J:2^(J+1)-1
     tree(k).blk = M2L_L2L(tree,k,pp);  % % Multipole-to-Local and local to local transfer of far field sources
     tree(k).blk_2 = M2L_2(tree,k,pp);  % multipole to local transfer for near field sources at same level
end

  sources_ = [sources(2:end); sources(end)+10];
  targets =  sources + 0.5*(sources_ - sources);  % selecting the mid point of each interval for testing function evaluation
  targets_ = targets;
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

 for it = 1:iter
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  FMM sum
  fmm_sum =  zeros(n,1); 
  der_fmm_sum = zeros(n,1);

  tree =  update_target(tree,targets,pp,J); % update target expansions for bottom level
tic
    for k = 2^J:2^(J+1)-1
        node_target = tree(k).targets;  % node targets
        [~, idx] = ismember(node_target,targets); % index of node targets
        tree(k).near_sum = near_sum(tree,k,sources,bb);  % direct sum for neighbor sources
        tree(k).der_near_sum = der_near_sum(tree,k,sources,bb);  % direct sum for neighbor sources
        y = tree(k).y_coeff; 
        fmm_sum(idx) = tree(k).near_sum + ((tree(k).blk)'*y)' + ((tree(k).blk_2)'*y)';
        der_fmm_sum(idx) = tree(k).der_near_sum + der_fmm_sum_1(tree,k,pp);
    end

    fmm_sum = 1 + s*fmm_sum;
    der_fmm_sum = s*der_fmm_sum;

elapsed_time(it,1,jj) = toc;

%%%%%%%%%%%%%%%%%%%% direct sum

dir_sum = zeros(n,1);
der_dir_sum = zeros(n,1);
nr = bb.^2;

tic
for i = 1:length(targets_)
        dn = sources-targets_(i);  % (x-y)
        der_dn = (sources-targets_(i)).^2;  % (x-y)^2
        dir_sum(i) = sum(nr./dn); % function estimation
        der_dir_sum(i) = sum(nr./der_dn); 
end
dir_sum = 1 + s*dir_sum;
der_dir_sum = s*der_dir_sum;

elapsed_time(it,2,jj) = toc;


%%%%%%%%%%%%%%%%%%% 

targets_new = targets - fmm_sum./der_fmm_sum;   % Newton approximation
targets_dir = targets_ - dir_sum./der_dir_sum;   % Newton approximation

targets = targets_new;
targets_ = targets_dir;

end  % iter loop

%%%%%%%%%%%%%%%%%%%% function evaluation at the convergence

targets_fmm_conv{jj} = targets;
targets_dir_conv{jj} = targets_; 

new_sum = zeros(n,1);
new_sum_ = zeros(n,1);


 for i = 1:length(targets)
        dn = sources-targets_new(i);  % (x-y)
        dn_ = sources-targets_dir(i);
        new_sum(i) = sum(nr./dn); % function estimation
        new_sum_(i) = sum(nr./dn_);
end
new_sum_fmm{jj} = 1+s*new_sum ;
new_sum_dir{jj} = 1+s*new_sum_ ;

end  % test data end loop

 


% Result analysis
rel_err = zeros(5,1);

for i9 = 1:5
    fmm_t = targets_fmm_conv{i9};
    dir_t = targets_dir_conv{i9};
    zz = abs(fmm_t- dir_t);
    idx = find(~isinf(zz)); length(idx);

    rel_err(i9) = norm(zz(idx),2)/norm(dir_t(idx),2);
end

fprintf('Relative error of FMM computed roots for %d expansion coefficients \n \n', pp);
disp(rel_err)

fprintf('Elapsed time in seconds for FMM vs direct computation for %d expansion coefficients for 10 iterations \n \n',pp)

for i9 = 1:5
fprintf('particle size = %d \n', i9*1e4);
fprintf('%-20s %-20s\n', 'FMM evaluation', 'Direct evaluation');   
disp(elapsed_time(:,:,i9))
end













