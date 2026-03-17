% H = SAT x UE x N
% N is the number of antennas per SAT
function H = FUNC_PATH_SAT(MODEL, N, Bw, f)
numSAT = size(MODEL.SAT,2); numUE = size(MODEL.UE,2);
% Distance from SATs to users, row-SAT, column-UE
Distance = sqrt( (MODEL.SAT(1,:)'-MODEL.UE(1,:)).^2 + ...
                 (MODEL.SAT(2,:)'-MODEL.UE(2,:)).^2 + ...
                 (MODEL.SAT(3,:)'-MODEL.UE(3,:)).^2 );
H = zeros(numSAT, numUE, N);
for m = 1:numSAT
    for u = 1:numUE
        H(m,u,:) = func_path_1sat_1ss(N, Distance(m,u), f, Bw);
    end
end
end


% N: the number of radiation elements
% d_mu: the distance from SAT to sensor u
function H_mu = func_path_1sat_1ss(N, d_mu, f, Bw)
% lamda: wavelength
lambda = 3e8 / f;
% kB: the Boltzmann constant
kB = 1.38*1e-23;
% TR: the receiver noise temperature
TR = 235.3; % Kelvin xxxx
% Bw: the carrier bandwidth
%Bw = 500*1e6; % 500MHz xxxx
% GR: Receibed antennas power gain at the sensor/ ground terminal
GR = 10^(40.7/10); % 40.7 dBi
% % Channel 1: PHI*C
% Phi = 2*pi*rand(N,1);
% PHI = exp(1i*Phi);
% G_mu = ones(N,1); % normalize
% C_mu = sqrt(GR*G_mu) / (4*pi*d_mu/lambda*sqrt(kB*TR*Bw));
% H_mu = PHI.*C_mu;

% Channel 2: A_mu*G_mu
% Asmopheric fading
A_mu = 1; % xxxx
% Fading
Phi = 2*pi*rand(N,1);
G_mu = exp(1i*Phi);
H_mu = 1/sqrt(A_mu) .* (sqrt(GR)*G_mu) / (4*pi*d_mu/lambda*sqrt(kB*TR*Bw));

end
