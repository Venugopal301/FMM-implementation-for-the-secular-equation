%% Local-to-Local translation from parent to child

function tree = L2L(tree,k,pp)  % input child node
k1 = floor(0.5*k);  % parent node
h = tree(k).center - tree(k1).center;
R = L2L_matrix(pp+1,h); 
R = R'; % prolongation matrix
blk = tree(k1).blk;
tree(k).blk = tree(k).blk +  R*blk;

    function A = L2L_matrix(n,h)
        A = zeros(n);               % preallocate
        for i = 1:n
            for j = 1:i
                A(i,j) = nchoosek(i-1, j-1) * h^(i-j);
            end
        end
    end

 
end  % function end







