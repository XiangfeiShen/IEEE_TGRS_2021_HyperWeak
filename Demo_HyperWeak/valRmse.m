function [error_min,r] = valRmse(H, Hest)
[P, K] = size(H);
error_all = zeros(P);
for i = 1:P
    for j = 1:P
        error_all(i, j) = norm(H(i, :) - Hest(j, :))/sqrt(K);
    end
end
error_min = min(error_all, [], 2);
r = mean(error_min);
fprintf('final_RMSE estimated:%.4f \n',r);