%% assign neighbors

function out = assignneighbor(k)

[power, ~] = decomposeNumber(k);

if k == 2^power
    out = k+1;
elseif k == 2^(power+1) -1
      out = k-1;
else
    out = [k-1 k+1];
end



   function [kk, remainder] = decomposeNumber(N)
        if N <= 0
            error('Input must be a natural number (positive integer).');
        end
        kk = floor(log2(N));   
        powerOf2 = 2^kk; % Largest power of 2 ≤ N
        
        remainder = N - powerOf2;    % Compute remainder
    end

end


