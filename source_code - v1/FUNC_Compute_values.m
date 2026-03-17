% NumUE_in_SAT: 1xM
% IdClusUE: 1xU (the id cluster of each UE)
% Rho: MxU (numSAT X numUE)
% Pow: MxU (numSAT X numUE)
function [sumrate, sumpow, Rate] = FUNC_Compute_values(NumUE_in_SAT, IdClusUE, Rho, sigma_u, Pow)
numSAT = size(NumUE_in_SAT,2); numUE = sum(NumUE_in_SAT);
Rate = zeros(1, numUE);
sumpow = sum(sum(Pow));
for m = 1:numSAT
    % users in SAT m
    Um = find(IdClusUE==m);
    for iu = 1:length(Um)
        u = Um(iu);
        Rate(u) = log2(1 + Pow(m,u)*Rho(m,u)/sigma_u);
    end
end
sumrate = sum(Rate);
end