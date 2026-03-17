function [SAT, UE] = FUNC_Create_SAT_UE(Area, numSAT, numUE, zm)
x_area = Area(1); y_area = Area(2);

%% Create SATs
% Firstly, assume numSAT = 4
SAT = zeros(3, 4);
SAT(1,:) = [x_area/4 x_area*3/4 x_area*3/4 x_area/4];
SAT(2,:) = [y_area/4 y_area/4 y_area*3/4 y_area*3/4];
SAT(3,:) = [zm zm zm zm];

%% Create UEs
UE = zeros(3, numUE);
UE(1,:) = x_area*rand(1,numUE);
UE(2,:) = y_area*rand(1,numUE);

end