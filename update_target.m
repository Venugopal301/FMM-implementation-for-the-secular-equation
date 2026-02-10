%% update the targets for the top level

function tree = update_target(tree,targets,pp,J)  % parent cell, cell index of the current cell

for k = 2^J:2^(J+1)-1 
    a = tree(k).interval(1);
    b = tree(k).interval(2);
    y = targets(targets > a & targets <= b);
    if ~isempty(y)
    tree(k).targets = y; y = y.';
    c = tree(k).center;
    i = (0:pp)';               % column vector of powers (ppx1)
    tree(k).y_coeff = (y - c).^i;  
    else
        break
    end
end

end