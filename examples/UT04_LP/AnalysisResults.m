% UT04 with LP with python library OSPA
% Cleaning
clear all
close all
clc

% Definition of the results file to load
filename = ["LatDispl.txt";...
        "ColSpringB.txt";...
        "ColSpringT.txt";...
        "PZDef.txt";...
        "BeamSpringE.txt";...
        "BeamSpringW.txt";...
        "TipColumnForce.txt"];
    
newPath_and_filename_outside = ["/UT04_Histories/", "Test_F.txt";
    "/UT04_Histories/", "Test_u.txt";
    "/UT04_Histories/", "Test_PZ_V_gamma.txt";
    "/UT04_Histories/", "Test_theta_pz.txt";
    "/UT04_Histories/", "Test_theta_b.txt";
    "/UT04_Histories/", "Test_theta_c.txt"];

% Initialize class for loading data
Load = LoadData;

% Loading the results
data = Load.Results(filename);

% Selecting the data to be analysed
delta = Load.ChooseResult(data, 1, 2);
theta_c_T = Load.ChooseResult(data, 2);
theta_c_B = Load.ChooseResult(data, 3);
gamma = Load.ChooseResult(data, 4);
theta_b_W = Load.ChooseResult(data, 5);
theta_b_E = Load.ChooseResult(data, 6);
LoadTip = -Load.ChooseResult(data, 7, 2);

% Load data from sections and material models
tlines = Load.InitilizeLoading("SavedInfos.txt");
[ColB, next_line] = Load.LoadNextInfo(tlines, 1);
[ColT, next_line] = Load.LoadNextInfo(tlines, next_line);
[BeamE, next_line] = Load.LoadNextInfo(tlines, next_line);
[BeamW, next_line] = Load.LoadNextInfo(tlines, next_line);
[IMKColB, next_line] = Load.LoadNextInfo(tlines, next_line);
[IMKColT, next_line] = Load.LoadNextInfo(tlines, next_line);
[IMKBeamE, next_line] = Load.LoadNextInfo(tlines, next_line);
[IMKBeamW, next_line] = Load.LoadNextInfo(tlines, next_line);
[PZMat, ~] = Load.LoadNextInfo(tlines, next_line);

L = BeamE.L*2 + ColB.d; % m
H = ColT.L + ColB.L + BeamE.d; % m
n = 10;
E = ColB.E; % N/m²
tp = ColB.tw; % m
kips_to_kN = 4.4482216;

% Import data from ouside
Test_F = Load.ImportDataFromOutside(newPath_and_filename_outside(1, :));
Test_u = Load.ImportDataFromOutside(newPath_and_filename_outside(2, :));
Test_PZ = Load.ImportDataFromOutside(newPath_and_filename_outside(3, :), [1, 2], '\t');
Test_gamma = Test_PZ(:, 1);
% Test_Vpz = Test_PZ(:, 2);
Test_theta_pz = Load.ImportDataFromOutside(newPath_and_filename_outside(4, :));
Test_theta_b = Load.ImportDataFromOutside(newPath_and_filename_outside(5, :));
Test_theta_c = Load.ImportDataFromOutside(newPath_and_filename_outside(6, :));


% Compute checking parameters

% K
K_c = (ColB.L^3/(3*E*ColB.Iy) + ColT.L^3/(3*E*ColB.Iy))^(-1)*H;

K_b = (3*E*BeamE.Iy/BeamE.L + 3*E*BeamW.Iy/BeamE.L)/H * L/(2*BeamE.L);

Ke_pz = E/2.6*0.95*tp*ColB.d;
K_pz = Ke_pz /(1-ColB.d/L-BeamE.d/H) /(H*2*BeamE.L/(L*BeamE.d)-1);

K = (1/K_b + 1/K_c + 1/K_pz)^(-1);

% Theta
SDR = delta/H;
theta_c = (-theta_c_T*ColT.L+theta_c_B*ColB.L)/H;
theta_pz = gamma*(1 - ColB.d/L - BeamE.d/H);
theta_b = -theta_b_W*(L - ColB.d)/L;


% Plots
x_lim = 7; % [% rad]

% Col
PlotSetup(1, 'Rotation [% rad]','Tip Load [kN]', x_lim);
plot(theta_c*100, LoadTip, '-k', 'DisplayName', '\theta_c');
plot(LoadTip/K_c*100, LoadTip, '--r', 'DisplayName', 'Check');
title('SDR Col');

% SDR PZ
PlotSetup(2, 'Rotation [% rad]','Tip Load [kN]', x_lim);
plot(theta_pz*100, LoadTip, '-k', 'DisplayName', '\theta_{pz}');
plot(LoadTip/K_pz*100, LoadTip, '--r', 'DisplayName', 'Check');
ylabel('Tip Load [kN]');
title('SDR PZ');

% SDR beam
PlotSetup(3, 'Rotation [% rad]','Tip Load [kN]', x_lim);
plot(theta_b*100, LoadTip, '-k', 'DisplayName', '\theta_b');
plot(LoadTip/K_b*100, LoadTip, '--r', 'DisplayName', 'Check');
title('SDR Beam');

% SDR and SUM SDR
PlotSetup(4, 'Rotation [% rad]','Tip Load [kN]', x_lim);
plot(SDR*100, LoadTip, '-k', 'DisplayName', 'SDR', 'LineWidth', 1);
plot(LoadTip/K*100, LoadTip, '--r', 'DisplayName', 'Check');
title('SDR');

%% Compare results with Shin

% Compare F - delta
PlotSetup(1, 'Story drift ratio \theta [% rad]', 'Applied Force \itF\rm [kN]', x_lim);
plot(-Test_u*100, -Test_F*kips_to_kN, 'LineStyle', '-', 'color', 'black', 'LineWidth', 1.0, 'DisplayName', 'Test data');
plot(SDR*100, LoadTip, '-b', 'DisplayName', 'Model', 'LineWidth', 1);
plot(LoadTip/K*100, LoadTip, '--r', 'DisplayName', 'Hand calculation', 'LineWidth', 1);

% Panel Zone Contribution
PlotSetup(2, 'Story drift angle due to panel zone \theta_{pz} [% rad]', 'Applied Force \itF\rm [kN]', x_lim);
plot(Test_theta_pz*100,Test_F*kips_to_kN,'LineStyle','-','color','black','LineWidth',1.0);
plot(theta_pz*100, LoadTip, '-r', 'DisplayName', '\theta_{pz}', 'LineWidth', 1);
plot(LoadTip/K_pz*100, LoadTip, '--r', 'DisplayName', 'Check', 'LineWidth', 1);

% Beam Contribution
PlotSetup(3, 'Story drift angle due to beam \theta_{b} [% rad]', 'Applied Force \itF\rm [kN]', x_lim);
plot(Test_theta_b*100,Test_F*kips_to_kN,'LineStyle','-','color','black','LineWidth',1.0);
plot(theta_b*100, LoadTip, '-r', 'DisplayName', '\theta_b', 'LineWidth', 1);
plot(LoadTip/K_b*100, LoadTip, '--r', 'DisplayName', 'Check', 'LineWidth', 1);

% Column Contribution
PlotSetup(4, 'Story drift angle due to column \theta_{c} [% rad]', 'Applied Force \itF\rm [kN]', x_lim);
plot(Test_theta_c*100,Test_F*kips_to_kN,'LineStyle','-','color','black','LineWidth',1.0);
plot(theta_c*100, LoadTip, '-r', 'DisplayName', '\theta_c', 'LineWidth', 1);
plot(LoadTip/K_c*100, LoadTip, '--r', 'DisplayName', 'Check', 'LineWidth', 1);

%% Save data for sharing
save("BeamW.mat", "BeamW");
save("BeamE.mat", "BeamE");
save("ColB.mat", "ColB");
save("ColT.mat", "ColT");
save("IMKBeamE.mat", "IMKBeamE");
save("IMKBeamW.mat", "IMKBeamW");
save("IMKColT.mat", "IMKColT");
save("IMKColB.mat", "IMKColB");

















