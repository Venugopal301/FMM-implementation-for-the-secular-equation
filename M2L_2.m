%% multipole to local summation at same level

function blk_2 = M2L_2(tree,k,pp)  % 

k2 = tree(k).local_ilist;
if isempty(k2)
     blk_2 = zeros(pp+1,1);
        return
end

d = 1/(tree(k2).center - tree(k).center);     % 1/(x*-y*)
    bin1 = zeros(pp+1,pp+1); % matrix to hold binomial expansion coefficients
    bin2 = zeros(pp+1,pp+1);   % matrix to hold d = 1/(x* - y*)^(-b)
    bin2(1,:) = d.^(1:pp+1);  %    
    for b = 1:pp+1
        c = negbin_coeffs(b,pp+1);
        bin1(b,:) = ((-1)^b)*c;
        bin2(b,:) = bin2(1,:)*d^(b-1);
    end
    bin2 = bin2.*(bin1.');  % multipole to local translation matrix
    farfield_coeff = tree(k2).x_coeff;
    farfield_sum   = sum(farfield_coeff,2);

    blk_2 = bin2*farfield_sum;

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

end















