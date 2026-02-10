%% evaluating the derivative sum by FMM

function sum1 = der_fmm_sum_check7(tree,k,pp)  % input child node

sum1 = 0; tt = tree(k).targets;

l1 = tree(k).far_ilist;

for k1 = 1:numel(l1)
    n1 = size(tree(l1(k1)).x2_coeff);
    blk = ((1:pp+1).')*ones(1,n1(2)).*tree(l1(k1)).x2_coeff;
    blk = sum(blk,2);
    Y = 1./(tt-tree(l1(k1)).center);
    Y = Y.^(2:pp+2);
    sum1 = sum1 + (Y*blk)';
end

l2 = tree(k).local_ilist;
if ~isempty(l2)
    n2 = size(tree(l2).x2_coeff);
    blk_2 = ((1:pp+1).')*ones(1,n2(2)).*tree(l2).x2_coeff; 
    blk_2 = sum(blk_2,2);
    Y = 1./(tt-tree(l2).center);
    Y = Y.^(2:pp+2);
    sum2 = (Y*blk_2)';
  else
    sum2  = zeros(1,size(tt,1));
end

sum1 = sum1' + sum2';
end















