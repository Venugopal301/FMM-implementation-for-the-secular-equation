%% function to multiply each multipole coefficient with b_j

function tree = update_bj(tree,bb)
% Given
A1 = tree(4).multipole_coeff;   
A2 = tree(5).multipole_coeff;
A3 = tree(6).multipole_coeff;
A4 = tree(7).multipole_coeff;
bb  = -bb.^2;     % negative sign due to series formulation 1/(x-y)

cols = [size(A1,2), size(A2,2), size(A3,2), size(A4,2)];

% compute index ranges
idx1 = 1:cols(1);
idx2 = (idx1(end)+1) : (idx1(end)+cols(2));
idx3 = (idx2(end)+1) : (idx2(end)+cols(3));
idx4 = (idx3(end)+1) : (idx3(end)+cols(4));

% scale each matrix: each column j is multiplied by scalar b(j)
A1s = A1 .* (ones(size(A1,1),1) * bb(idx1)');   % 
A2s = A2 .* (ones(size(A2,1),1) * bb(idx2)');   % 
A3s = A3 .* (ones(size(A3,1),1) * bb(idx3)');   % 
A4s = A4 .* (ones(size(A4,1),1) * bb(idx4)');   % 

tree(4).multipole_coeff = A1s;
tree(5).multipole_coeff = A2s;
tree(6).multipole_coeff = A3s;
tree(7).multipole_coeff = A4s;

end  % end of function

















