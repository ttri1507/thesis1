clc
clear
close all

% Scenario
N = 100; % The number of radiation feeds/antennas
B = 25; % The maximum beams per SAT
numSAT = 4; % The number of SATs M = 4
numUE = 80;
Pmax = 10^(20/10); % The maximum power of each SAT 20 dBW
P0 = 10^(10/10); % Pcircuit = 10 dBW
Bw = 500*1e6; % 500MHz
% noiseVariancedBm = -174 + 10*log10(Bw) + 9;
% Noise_var=db2pow(noiseVariancedBm-30);
Noise_var = 1;
f = 20e9; % 20GHz
Area = [2000e3 2000e3]; % 2000kmx2000km
zm = 1000e3; % The altitude of SATs
t_thres = 2000e3/3e8; % The distance must be less than 1546km
Rmin = 0.2; % xxxx the minimum data rate each UE


%% Create model
MODEL = struct();
[SAT, UE] = FUNC_Create_SAT_UE(Area, numSAT, numUE, zm);
MODEL.SAT = SAT; MODEL.UE = UE;
%FUNC_Plot_model_2D(MODEL, Area, "normal");
MODEL.Channel = FUNC_PATH_SAT(MODEL, N, Bw, f);

%% Clustering based Coalition game
theta = 0.5; % Trade-off between latency and channel gain
% start_RC = tic;
[MODEL_RC] = FUNC_randomClustering(MODEL,B,t_thres,theta,Area);
% time_RC = toc(start_RC);
% figure(1),FUNC_Plot_model_2D(MODEL_RC, Area, "cluster")

IdClusUE_RC = MODEL_RC.IdClusUE; NumUE_in_SAT_RC = MODEL_RC.NumUE_in_SAT;
start_GTC = tic;
[MODEL_GTC] = FUNC_CoalitionGame(MODEL,B,t_thres,theta,Area,IdClusUE_RC,NumUE_in_SAT_RC);
time_GTC = toc(start_GTC);
% figure(2),FUNC_Plot_model_2D(MODEL_GTC, Area, "cluster")
% save_to_pdf(2)

%% Precoding
H_norm2 = sqrt(abs(sum( conj(MODEL.Channel).*MODEL.Channel, 3)));
% precoding for all UEs
Precoding = conj(MODEL.Channel)./ repmat(H_norm2,1,1,N);

%% Power allocation
% Rho square
Rho = abs(sum(MODEL.Channel.*Precoding,3)).^2;
sigma_u = Noise_var; % sigma square
Pow_overline = sigma_u*(2^Rmin - 1)./Rho;
Beta = 1 + Pow_overline.*Rho/sigma_u;

% %%%%%%%%%%%% GTC-PA
NumUE_in_SAT = MODEL_GTC.NumUE_in_SAT;
IdClusUE = MODEL_GTC.IdClusUE;
lambda_ref = 1; tau_ref2 = 1;
tau_underline2 = 1e-10; tau_overline2 = 1000;
SumP_overline = zeros(1, numSAT);
for m = 1:numSAT
    Um = find(IdClusUE==m);
    SumP_overline(m) = sum(Pow_overline(m,Um));
end
Pmax_overline = Pmax - SumP_overline;
Pow_cir_overline = P0 + SumP_overline;
[sumrate_GTC, sumpow_GTC, obj_GTC, Pow_GTC, Rate_GTC, LoopEE_GTC] = FUNC_multiSAT_OPT(NumUE_in_SAT, IdClusUE, lambda_ref, Beta, ...
            Rho, sigma_u, Pmax_overline, Pow_cir_overline, tau_ref2, tau_underline2, tau_overline2);
EE_GTC = sumrate_GTC/(sum(Pow_cir_overline)+sumpow_GTC);

% %%%%%%%%%%%% RC-PA
NumUE_in_SAT = MODEL_RC.NumUE_in_SAT;
IdClusUE = MODEL_RC.IdClusUE;
lambda_ref = 1; tau_ref2 = 1;
tau_underline2 = 1e-10; tau_overline2 = 1000;
SumP_overline = zeros(1, numSAT);
for m = 1:numSAT
    Um = find(IdClusUE==m);
    SumP_overline(m) = sum(Pow_overline(m,Um));
end
Pmax_overline = Pmax - SumP_overline;
Pow_cir_overline = P0 + SumP_overline;
[sumrate_RC, sumpow_RC, obj_RC, Pow_RC, Rate_RC, LoopEE_RC] = FUNC_multiSAT_OPT(NumUE_in_SAT, IdClusUE, lambda_ref, Beta, ...
            Rho, sigma_u, Pmax_overline, Pow_cir_overline, tau_ref2, tau_underline2, tau_overline2);
EE_RC = sumrate_RC/(sum(Pow_cir_overline)+sumpow_RC);

% Show the result
X = 1:20; n = size(X,2);
semilogy(X,LoopEE_GTC(1:n),'r-*','markersize',4,'Linewidth',1);
hold on; grid on;
semilogy(X,LoopEE_RC(1:n),'b--^','markersize',4,'Linewidth',1);

legend('EE-GTC','EE-RC')
xlabel('Number of iterations')
ylabel('EE (bits/Joule/Hz)')

% %% Save EE to txt file
% Data = [X', LoopEE_GTC(1:n)', LoopEE_RC(1:n)'];
% fileID = fopen('EE-convergence.txt','w');
% fprintf(fileID,'%20s %20s %20s\n','X','EE_GTC','EE_RC');
% for i = 1:n
%     fprintf(fileID,'%20.5f %20.5f %20.5f \n',Data(i,:));
% end
% fclose(fileID);

%% Save to pdf
function save_to_pdf(num)
%Position plot at left hand corner with width 5 and height 5.
set(figure(num), 'PaperPosition', [0 0 5 5]);
%Set the paper to have width 5 and height 5.
set(figure(num), 'PaperSize', [5 5]); 
saveas(figure(num), 'GT_Clustering', 'pdf')
end
