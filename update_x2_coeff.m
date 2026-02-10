%% function to update multipole expansion coefficients for derivative

function M = update_x2_coeff(tree,k,pp,sources,bb)

c = tree(k).center;
x = tree(k).sources.';
[~, idx] = ismember(x,sources);
bb = bb(idx);
bb = bb.^2; bb = bb';  % 
i = (0:pp)';               % column vector of powers (pp x 1)
M = (x - c) .^ i;  
b1 = ones(pp+1,1)*bb;
M = b1.*M;
end