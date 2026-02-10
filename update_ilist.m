%% function to update the interaction list for each node from level 4
function tree = update_ilist(tree,J)

X = 2^J: 2^(J+1)-1;

c = 2^J/4;  % no of columns
r = 2^J/c;  % no of rows

A = reshape(X,c,r); A = A.';
B = A(:,1)/2^(J-2);

%%%%%%%%%%%%%%%%%%%%%%%%%%% updating extreme leaf nodes
x = A(1,:);
x(1:2) = 0;  % left leaf
tree(A(1,1)).far_ilist = [B(2:end).', int_list(x)];
x = A(end,:);
x(end-1:end) = 0;  % right leaf
tree(A(end,end)).far_ilist = [B(1:end-1).', int_list(x)];
%%%%%%%%%%%%%%%%%%%%%%%% updating interior border nodes
for t = 1:3
    e = c;
    x = [A(t,:),A(t+1,:)];  % consider first column
    x(e-1:e+1)=0;   %
    [list1, list2] = int_list(x);
    idx = setdiff(1:4, [t, t+1]);  % Indices other than t and t+1
    picked = B(idx);               % Pick elements other than B(t), B(t+1)
    tree(A(t,e)).far_ilist = [picked.', list1];
    tree(A(t,e)).local_ilist = list2;

    e = e+1;
    x = [A(t,:),A(t+1,:)];  % consider first column
    x(e-1:e+1)=0;    %
    [list1, list2] = int_list(x);
    tree(A(t+1,1)).far_ilist = [picked.', list1];
    tree(A(t+1,1)).local_ilist = list2;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% updating interior nodes

for t = 1:r
    for e = 2:c-1
        x = A(t,:);  % consider first column
        x(e-1:e+1)=0;   %
        [list1, list2] = int_list(x);
        tree(A(t,e)).far_ilist = [B([1:t-1, t+1:end]).', list1];
        tree(A(t,e)).local_ilist = list2;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function [far_list,local_list] = int_list(x)

        x_grouped = reshape(x, 2, []); x_grouped = x_grouped.'; %

         far_list = [];    % For storing half of first element when both are nonzero

        for i1 = 1:size(x_grouped, 1)
            row = x_grouped(i1, :);
            nz = row(row ~= 0);  % Nonzero elements in the row

            if isscalar(nz)
                local_list = nz;  % Append single nonzero element
            elseif length(nz) == 2
                far_list(end+1) = row(1)/2;  % Both nonzero: take half of first
            end
            % If both zero: do nothing
        end
        far_list =    reduce_consecutive_even_pairs(far_list);
    end


    function Y = reduce_consecutive_even_pairs(X)
        while true
            i = 1;
            changed = false;
            Y = [];

            while i <= length(X)
                if i < length(X) && mod(X(i), 2) == 0 && X(i+1) == X(i) + 1
                    Y(end+1) = X(i) / 2;
                    i = i + 2;      % skip next (merged)
                    changed = true; % mark that a change occurred
                else
                    Y(end+1) = X(i);
                    i = i + 1;
                end
            end

            if ~changed
                break; % stop if no replacements were made in this pass
            end

            X = Y; % update input for next pass
        end
    end
end  % function end






