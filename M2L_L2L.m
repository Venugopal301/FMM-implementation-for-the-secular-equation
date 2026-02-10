%% function to compute the Multipole_to_Local coefficient expansion

function blk = M2L_L2L_check5(tree,k,pp)

list1 = tree(k).far_ilist;
list2 = groupInteractionByLevel(list1,k);
list1 = list2(1,:);
list2 = list2(2,:);

bin_coeff = zeros(pp+1,pp+1); % matrix to hold binomial expansion coefficients
  for b = 1:pp+1
      bin_coeff(:,b) = ((-1)^(b))*negbin_coeffs(b,pp+1); % (x-y)-n coefficients
  end

for i1 =   1:numel(list1)
 k2 = list2(i1);    
 kk = list1(i1);    
    d = 1/(tree(kk).center - tree(k2).center);  % 1/(x* - y*)
    d_matrix = zeros(pp+1,pp+1);   % matrix to hold d = 1/(x* - y*)^b

    for b = 1:pp+1
        d_matrix(:,b) = d.^(b:b+pp);
    end
    
    M2L_matrix = bin_coeff.*d_matrix;  % multipole to local translation matrix
    
    farfield_coeff = tree(kk).x_coeff;
    farfield_sum   = sum(farfield_coeff,2);

    tree(k2).blk = tree(k2).blk + M2L_matrix*farfield_sum;
end

path_all = findPathsBottomUp(list2,k);  % local to local translations from top to bottom
path_1 = path_all{1};
    for  i2 =  1:length(path_1)-1
        k4 = i2+1;
        tree = L2L(tree,path_1(k4),pp); % L2L translation till target node
    end
blk = tree(k).blk;
    function c = negbin_coeffs(n,K)
        %NEGBIN_COEFFS Computes coefficients for (1 - x)^(-n)
        % n  - positive real (can be non-integer too)
        % K  - number of terms required (starting from k=0)
        %
        % Returns c(1) = coefficient of x^0, c(2) = coefficient of x^1, etc.
        c = zeros(1, K);
        c(1) = 1; % k = 0 term

        for k1 = 1:K-1
            c(k1+1) = c(k1) * (n + k1 - 1) / k1;
        end
    end

function M = groupInteractionByLevel(v, targetNode)
%GROUPINTERACTIONBYLEVEL Return interaction list with matching ancestor per level
%   v           - interaction list (vector of node numbers)
%   targetNode  - node whose ancestors will be attached
%
%   Output: 2 x N matrix:
%           Row 1 = interaction list sorted by level
%           Row 2 = ancestor for each entry in Row 1

    % Sort interaction list
    v = sort(v(:))';

    % Step 1: Compute exponents (levels) for interaction list
    exponentV = floor(log2(v));

    % Step 2: Find all ancestors of targetNode
    ancestors = [];
    cur = targetNode;
    while cur > 1
        cur = floor(cur / 2);
        ancestors(end+1) = cur; %
    end

    % Step 3: Exponents for ancestors
    ancestorExp = floor(log2(ancestors));

    % Step 4: Build output
    M = zeros(2, numel(v));
    M(1, :) = v;

    for i = 1:numel(v)
        level = exponentV(i);
        idx = find(ancestorExp == level, 1);
        if ~isempty(idx)
            M(2, i) = ancestors(idx);
        else
            M(2, i) = NaN; % if no matching ancestor exists
        end
    end
end

function paths = findPathsBottomUp(startNodes, targetNode)
    % startNodes: vector of start nodes
    % targetNode: single target node
    startNodes = unique(startNodes, 'stable'); % preserve order
    paths = cell(1, numel(startNodes));
    
    for i = 1:numel(startNodes)
        startNode = startNodes(i);
        path = targetNode;
        currentNode = targetNode;
        
        while true
            % Determine level of current node
            level = floor(log2(currentNode));
            currLevelStart = 2^level;
            indexInLevel = currentNode - currLevelStart + 1;
            
            % Previous level
            prevLevel = level - 1;
            if prevLevel < 0
                break; % reached top
            end
            prevLevelStart = 2^prevLevel;
            prevLevelEnd   = 2^(level) - 1;
            prevNodes = prevLevelStart:prevLevelEnd;
            
            % Find parent using logic
            parentIndex = ceil(indexInLevel/2);
            parentNode = prevNodes(parentIndex);
            
            % Prepend parent to path
            path = [parentNode path]; %
            
            % Stop if we reached the start node
            if parentNode == startNode
                break;
            end
            
            % Move up one level
            currentNode = parentNode;
        end
        
        paths{i} = path;
    end
end


end % end of function