%% near field summation from neighbors and cell sources

function sum2 = near_sum(tree,k,sources,bb)

nb = tree(k).neighbor;
ss = [];

for i = 1:numel(nb)
    ss = [ss;tree(nb(i)).sources];
end
ss = [ss; tree(k).sources];


[~, idx] = ismember(ss,sources);
nr = bb(idx).^2;
tt = tree(k).targets;
sum2 = zeros(length(tt),1);

for i = 1:length(tt)
    dn = ss - tt(i);
    sum2(i) = sum(nr./dn);
end

end  % function end