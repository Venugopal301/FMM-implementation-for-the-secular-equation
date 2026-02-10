%% Local-to-Local translation from parent to child of the derivative

function tree = L2L_der(tree,k,pp)  % input child node
k1 = floor(0.5*k);  % parent node
h = tree(k).center - tree(k1).center;
R = L2L_matrix(pp+1,h); 
R = R'; % prolongation matrix
der_blk = tree(k1).der_blk;
tree(k).der_blk = tree(k).der_blk + R*der_blk;

    function A = L2L_matrix(n,h)
        A = zeros(n);               % preallocate
        for i = 1:n
            for j = 1:i
                A(i,j) = nchoosek(i-1, j-1) * h^(i-j);
            end
        end
    end

 
end  % function end







