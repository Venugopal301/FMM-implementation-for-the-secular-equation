%% build tree from top down, assign sources, interval, center, parent cell id and cell id

function tree1 = tree_update_1(tree,k,pp)  % parent cell, cell index of the current cell

tree1 = struct('level', {}, 'interval', {}, 'center', {}, 'parent', {}, 'child', {}, 'neighbor', {},...
    'far_ilist', {}, 'local_ilist', {}, 'sources', {}, 'targets',{}, 'x_coeff', {}, 'x2_coeff', {}, ...
    'y_coeff', {}, 'blk_2',{}, 'blk', {}, 'der_blk',{}, 'der_blk_2', {}, 'near_sum', {}, 'der_near_sum', {} );

a = tree.interval(1);
b = tree.interval(2);

center = a + 0.5*(b-a);

if mod(k,2) == 0
    tree1(1).sources = tree.sources(tree.sources <= center);
    tree1(1).interval = [a center];
    tree1(1).center = mean(tree1.interval);
    tree1(1).parent = 0.5*k;
    tree1(1).child =  [2*k, 2*k+1]; % update child nodes

else
    tree1(1).sources = tree.sources(tree.sources > center);
    tree1(1).interval = [center b];
    tree1(1).center = mean(tree1.interval);
    tree1(1).parent = floor(0.5*k);
    tree1(1).child =  [2*k, 2*k+1]; % update child nodes
end

tree1(1).blk = zeros(pp+1,1);
tree1(1).der_blk = zeros(pp+1,1);
tree1(1).blk_2 = zeros(pp+1,1);
tree1(1).der_blk_2 = zeros(pp+1,1);

[k1,~] = decomposeNumber(k);  % k = 2^k1 + m, k1 is the largest power of 2 contained in k
tree1.level = k1;

    function [k1, remainder] = decomposeNumber(N)  %
        if N <= 0
            error('Input must be a natural number (positive integer).');
        end
        k1 = floor(log2(N));  % largest power of 2 contained in the number
        powerOf2 = 2^k1; % Largest power of 2 ≤ N
        remainder = N-powerOf2;    % N = powerOf2 + remainder
    end


end % function end


