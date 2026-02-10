%% function to translate multipole-to-multipole coefficients from bottom to top (M2M) for derivatives
function M = M2M3(tree,k,pp)

M = [];
for cc = 1:2  % loop over each child
    k1 = 2*k+cc-1;  % k1 is the child node index
    x = tree(k1).x2_coeff;
    h = tree(k1).center - tree(k).center;    % interval length h
    R = zeros(pp+1);  % lower triangular matrix R

    for i = 0:pp  % row index = power
        for j = 0:i
            coeff = nchoosek(i, j);
            R(i+1,j+1) = (coeff*(h^(i-j)))*(i+1); %  h^{i-j} * binomial coefficient
        end
    end
    M = [M, R*x];
end

end
