%% Altman Z-Score Calculation
% Specify the path to the .mat file for the scoring models dataset
filename_Z_O = 'C:\\Users\\Ildar Shakirzianov\\Documents\\ID_11355000_CRM_Data_1.mat';

% Load the .mat file
loadedData = load(filename_Z_O);
dataTable_Z_and_O_Scores = array2table(loadedData.data);  % Replace 'data' with your actual matrix variable name

% Assign variable names to the columns in the table
dataTable_Z_and_O_Scores.Properties.VariableNames = {'TA', 'CL', 'CA', 'WC', 'RE','EBIT', 'TL', 'MVE', 'REV', 'GNP', 'OPCF', 'NI','INTWO', 'OENEG'};

% Extract the first and second rows of data
dataRows = dataTable_Z_and_O_Scores(1:2, :);

% Create an array to store the Z-scores
Z_scores = zeros(1, height(dataRows));

% Calculate the Z-score for each row
for i = 1:height(dataRows)
    % Calculate the ratios
    WCTA = (dataRows.WC(i) / dataRows.TA(i));
    RETA = (dataRows.RE(i) / dataRows.TA(i));
    EBITTA = (dataRows.EBIT(i) / dataRows.TA(i));
    MVETL = (dataRows.MVE(i) / dataRows.TL(i));
    REVTA = (dataRows.REV(i) / dataRows.TA(i));
    
    % Calculate the Altman Z-score
    Z_scores(i) = 1.2 * WCTA + 1.4 * RETA + 3.3 * EBITTA + 0.6 * MVETL + 0.999 * REVTA;
end

% Display the Z-scores
fprintf('2020 Z-Score = %.4f\n', Z_scores(1));
fprintf('2022 Z-Score = %.4f\n', Z_scores(2));

%% O-Score Calculation & O-Score Probability of Default

% Create an array to store the O-scores
O_scores = zeros(1, height(dataRows));

% Define the values for NI for 2019 and 2021
NI_past_values = [-137.22, -150.77]; % You can add more values if needed

% Calculate the O-score for each row
for i = 1:height(dataRows)
    % Calculate the ratios
    SIZE = dataRows.TA(i) / dataRows.GNP(i);
    TLTA = dataRows.TL(i) / dataRows.TA(i);
    WCTA = dataRows.WC(i) / dataRows.TA(i);
    CLCA = dataRows.CL(i) / dataRows.CA(i);
    OENEG = dataRows.OENEG(i);
    NITA = dataRows.NI(i) / dataRows.TA(i);
    FUTL = dataRows.OPCF(i) / dataRows.TL(i);
    INTWO = dataRows.INTWO(i);
    
    % Use different NI_past values for 2020 and 2022
    if i == 1
        CHIN = (dataRows.NI(i) - NI_past_values(1)) / (abs(dataRows.NI(i)) + abs(NI_past_values(1)));
    elseif i == 2
        CHIN = (dataRows.NI(i) - NI_past_values(2)) / (abs(dataRows.NI(i)) + abs(NI_past_values(2)));
    else
        % Handle other iterations as needed
        % ...
    end
    
    % Calculate the O-scores for 2020 and 2022
    O_scores(i) = -1.32 - 0.407*log(SIZE) + 6.03*TLTA - 1.43*WCTA + 0.0757*CLCA - 1.72*OENEG - 2.37*NITA - 1.83*FUTL - 0.285*INTWO - 0.521*CHIN;
    PD(i) = exp(O_scores(i))/(1+exp(O_scores(i)));
end

% Display the O-scores and probabilities of default
fprintf('2020 O-Score = %.4f\n', O_scores(1));
fprintf('2022 O-Score = %.4f\n', O_scores(2));
fprintf('2020 O-Score Probability of default = %.4f\n', PD(1));
fprintf('2022 O-Score Probability of default = %.4f\n', PD(2));

%% Merton Model
clearvars
% Set the paths for 2020 and 2022 .mat data files


Data_2020_Merton = 'C:\\Users\\Ildar Shakirzianov\\Documents\\ID_11355000_CRM_Data_2.mat';
Data_2022_Merton = 'C:\\Users\\Ildar Shakirzianov\\Documents\\ID_11355000_CRM_Data_3.mat';

% Load the 2020 .mat file
loadedData = load(Data_2020_Merton);

% Convert the data into table
dataTable_Z_and_O_Scores = array2table(loadedData.data);  

% Assign variable names to the columns in the table
dataTable_Z_and_O_Scores.Properties.VariableNames = {'P_2020', 'N_2020', 'S_2020', 'DTB3_2020', 'R_2020'};

% Calculating the debt value
D_2020 = (3559.05/2 + 2466.53) * 1000000;

T = 1;

% Introducing variables
dataTable_Z_and_O_Scores.V_2020 = dataTable_Z_and_O_Scores.S_2020 + D_2020;
V_0_2020 = dataTable_Z_and_O_Scores.V_2020(1);
r_2020 = dataTable_Z_and_O_Scores.DTB3_2020/100;

% Calculating initial values for sigma_v
sigma_e_2020 = std(dataTable_Z_and_O_Scores.R_2020)*sqrt(252);
sigma_v_0_2020 = sigma_e_2020 * dataTable_Z_and_O_Scores.S_2020(1) / (dataTable_Z_and_O_Scores.V_2020(1));

% Initialize sigma_old and sigma_new
sigma_old_2020 = sigma_v_0_2020;
sigma_new_2020 = 0;

tolerance = 0.0001;

% Iterate until the change in sigma is less than the tolerance
while abs(sigma_new_2020 - sigma_old_2020) > tolerance
    % Store the old value of sigma
    sigma_old_2020 = sigma_new_2020;

    % Calculate the Black-Scholes call option price for each day
    Vt_estimation_2020 = zeros(length(dataTable_Z_and_O_Scores.S_2020), 1);
    C_2020 = zeros(length(dataTable_Z_and_O_Scores.S_2020), 1);
    for i = 1:length(dataTable_Z_and_O_Scores.S_2020)
        d1_2020 = @(Vt) (log(Vt/D_2020)+r_2020(i).*T+0.5*sigma_old_2020^2.*T)./(sigma_old_2020.*sqrt(T));
        d2_2020 = @(Vt) d1_2020(Vt)-sigma_old_2020*sqrt(T);
        BS_2020 = @(Vt) Vt.*normcdf(d1_2020(Vt))-D_2020*exp(-r_2020(i).*T).*normcdf(d2_2020(Vt));
        f_2020 = @(Vt) BS_2020(Vt) - dataTable_Z_and_O_Scores.S_2020(i);
        Vt_estimation_2020(i) = fzero(f_2020, dataTable_Z_and_O_Scores.V_2020(i));
        C_2020(i) = BS_2020(Vt_estimation_2020(i));
    end

    % Store the solutions in the Vt column
    dataTable_Z_and_O_Scores.Vt_2020 = Vt_estimation_2020;

    % Calculate mu_daily
    mu_daily_2020 = (log(dataTable_Z_and_O_Scores.("Vt_2020")(end)) - log(dataTable_Z_and_O_Scores.("Vt_2020")(1))) / 252;
    delta_t = 1/252;

    N = 252;

    % Calculating the new sigma
    sigma_new_2020 = sqrt((1/(N*delta_t)) * sum((log(dataTable_Z_and_O_Scores.Vt_2020(2:end)./dataTable_Z_and_O_Scores.Vt_2020(1:end-1))-mu_daily_2020).^2));
end

%% PD Estimation 2020
% Introducing a brownian motion 
dW_2020 = normrnd(0,sqrt(delta_t),[1000,T*252]); 


V_2020 = NaN(1000,1+T*252);
V_2020(:,1) = V_0_2020;

for i = 2:1+T*252
    V_2020(:,i) = V_2020(:,i-1)+V_2020(:,i-1)*mu_daily_2020*delta_t+V_2020(:,i-1).*sigma_new_2020.*dW_2020(:,i-1);
end

PD_2020 = sum(V_2020(:,end)<D_2020)/1000;

% Calculating the PD for 2020 and storing it
PD_Merton_2020 = normcdf(-(log(V_0_2020/D_2020)+(mu_daily_2020-0.5*sigma_new_2020^2)*T)/(sigma_new_2020*sqrt(T)));


%% 2022 Estimation
% Replication the solution above for year 2022

% Load the 2022 .mat file
loadedData = load(Data_2022_Merton);

% Convert the data into table
dataTable_2022 = array2table(loadedData.data);  % Replace 'data' with your actual matrix variable name

% Assign variable names to the columns in the table
dataTable_2022.Properties.VariableNames = {'P_2022', 'N_2022', 'S_2022', 'DTB3_2022', 'R_2022'};

% Calculating the debt value
D_2022 = (3559.05/2 + 2466.53) * 1000000;

T = 1;

% Introducing variables
dataTable_2022.V_2022 = dataTable_2022.S_2022 + D_2022;
V_0_2022 = dataTable_2022.V_2022(1);
r_2022 = dataTable_2022.DTB3_2022/100;

% Calculating initial values for sigma_v
sigma_e_2022 = std(dataTable_2022.R_2022)*sqrt(252);
sigma_v_0_2022 = sigma_e_2022 * dataTable_2022.S_2022(1) / (dataTable_2022.V_2022(1));

% Initialize sigma_old and sigma_new
sigma_old_2022 = sigma_v_0_2022;
sigma_new_2022 = 0;

tolerance = 0.0001;

% Iterate until the change in sigma is less than the tolerance
while abs(sigma_new_2022 - sigma_old_2022) > tolerance
    % Store the old value of sigma
    sigma_old_2022 = sigma_new_2022;

    % Calculate the Black-Scholes call option price for each day
    Vt_estimation_2022 = zeros(length(dataTable_2022.S_2022), 1);
    C_2022 = zeros(length(dataTable_2022.S_2022), 1);
    for i = 1:length(dataTable_2022.S_2022)
        d1_2022 = @(Vt) (log(Vt/D_2022)+r_2022(i).*T+0.5*sigma_old_2022^2.*T)./(sigma_old_2022.*sqrt(T));
        d2_2022 = @(Vt) d1_2022(Vt)-sigma_old_2022*sqrt(T);
        BS_2022 = @(Vt) Vt.*normcdf(d1_2022(Vt))-D_2022*exp(-r_2022(i).*T).*normcdf(d2_2022(Vt));
        f_2022 = @(Vt) BS_2022(Vt) - dataTable_2022.S_2022(i);
        Vt_estimation_2022(i) = fzero(f_2022, dataTable_2022.V_2022(i));
        C_2022(i) = BS_2022(Vt_estimation_2022(i));
    end

    % Store the solutions in the Vt column
    dataTable_2022.Vt_2022 = Vt_estimation_2022;

    % Calculate mu_daily
    mu_daily_2022 = (log(dataTable_2022.("Vt_2022")(end)) - log(dataTable_2022.("Vt_2022")(1))) / 252;
    delta_t = 1/252;

    N = 252;

    % Calculating the new sigma
    sigma_new_2022 = sqrt((1/(N*delta_t)) * sum((log(dataTable_2022.Vt_2022(2:end)./dataTable_2022.Vt_2022(1:end-1))-mu_daily_2022).^2));
end


%% PD Estimation 2022

% Introducing a brownian motion 
dW_2022 = normrnd(0,sqrt(delta_t),[1000,T*252]); 


V_2022 = NaN(1000,1+T*252);
V_2022(:,1) = V_0_2022;

for i = 2:1+T*252
    V_2022(:,i) = V_2022(:,i-1)+V_2022(:,i-1)*mu_daily_2022*delta_t+V_2022(:,i-1).*sigma_new_2022.*dW_2022(:,i-1);
end



PD_2022 = sum(V_2022(:,end)<D_2022)/1000;

% Calculating the PD for 2022 and storing it
PD_Merton_2022 = normcdf(-(log(V_0_2022/D_2022)+(mu_daily_2022-0.5*sigma_new_2022^2)*T)/(sigma_new_2022*sqrt(T)));

%% Results
% Printing out the probabilities of default for 2020 and 2022

fprintf('2020 Merton Probability of default = %.4f\n', PD_Merton_2020);
fprintf('2022 Merton Probability of default = %.4f\n', PD_Merton_2022);

% Some additional results
fprintf('2020 sigma after calibrating = %.4f\n', sigma_new_2020);
fprintf('2022 sigma after calibrating = %.4f\n', sigma_new_2022);


% Graph for V_2020 simulations
figure(2)
plot(V_2020')
line([0, T*252], [D_2020, D_2020], 'Color', 'g','LineWidth', 2);
fontsize(gca, 15, "points")
grid on;


% Graph for V_2022 simulations
figure(3)
plot(V_2022')
line([0, T*252], [D_2022, D_2022], 'Color', 'g','LineWidth', 2);
fontsize(gca, 15, "points")
grid on;




