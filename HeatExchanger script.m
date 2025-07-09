% HeatExchanger.m
% Analyzes a counter-flow heat exchanger to calculate heat transfer rate,
% outlet temperatures, effectiveness, fluid flow, and optimize performance for automotive cooling.
% Created by Uday Kaurav, May 2025, for Volvo Eicher interview preparation.

% Clear workspace and command window
clear all;
clc;

% Step 1: Define Parameters
T_h_in = 80;       % Hot fluid inlet temperature (°C)
T_c_in = 20;       % Cold fluid inlet temperature (°C)
m_h = 0.5;         % Hot fluid mass flow rate (kg/s)
m_c = 1.0;         % Cold fluid mass flow rate (kg/s)
c_ph = 4180;       % Hot fluid specific heat (J/kg·K)
c_pc = 1005;       % Cold fluid specific heat (J/kg·K)
U = 500;           % Overall heat transfer coefficient (W/m²·K)
A = 2.0;           % Initial heat transfer area (m²)

% Step 2: Initial Calculations (Baseline)
C_h = m_h * c_ph;  % Heat capacity rate of hot fluid (W/K)
C_c = m_c * c_pc;  % Heat capacity rate of cold fluid (W/K)
C_min = min(C_h, C_c);
C_max = max(C_h, C_c);
C_r = C_min / C_max;  % Capacity ratio

NTU = U * A / C_min;  % Number of Transfer Units
epsilon_guess = (1 - exp(-NTU * (1 - C_r))) / (1 - C_r * exp(-NTU * (1 - C_r)));
Q_max = C_min * (T_h_in - T_c_in);
Q = epsilon_guess * Q_max;

% Step 3: Calculate Outlet Temperatures
T_h_out = T_h_in - Q / C_h;
T_c_out = T_c_in + Q / C_c;

% Step 4: Calculate LMTD
DT1 = T_h_in - T_c_out;
DT2 = T_h_out - T_c_in;
if abs(DT1 - DT2) < 1e-6
    LMTD = DT1;
else
    LMTD = (DT1 - DT2) / log(DT1 / DT2);
end

% Step 5: Recalculate Q using LMTD
Q_LMTD = U * A * LMTD;

% Step 6: Adjust and Iterate
for i = 1:5
    Q = Q_LMTD;
    T_h_out = T_h_in - Q / C_h;
    T_c_out = T_c_in + Q / C_c;
    DT1 = T_h_in - T_c_out;
    DT2 = T_h_out - T_c_in;
    if abs(DT1 - DT2) < 1e-6
        LMTD = DT1;
    else
        LMTD = (DT1 - DT2) / log(DT1 / DT2);
    end
    Q_LMTD = U * A * LMTD;
    if abs(Q - Q_LMTD) < 1e-3
        break;
    end
end

% Step 7: Calculate Effectiveness
epsilon = Q / Q_max;

% Step 8: Display Baseline Results
fprintf('Baseline Heat Exchanger Analysis (A = %.1f m²):\n', A);
fprintf('--------------------------------\n');
fprintf('Heat Transfer Rate (Q): %.2f W\n', Q);
fprintf('Hot Fluid Outlet Temperature: %.2f °C\n', T_h_out);
fprintf('Cold Fluid Outlet Temperature: %.2f °C\n', T_c_out);
fprintf('Log Mean Temperature Difference (LMTD): %.2f °C\n', LMTD);
fprintf('Effectiveness (epsilon): %.3f\n', epsilon);

% Step 9: Optimization - Vary Heat Transfer Area (A) to Maximize Effectiveness
A_values = [2.0, 2.5, 3.0];  % Test different areas (m²)
epsilon_values = zeros(size(A_values));
Q_values = zeros(size(A_values));
T_h_out_values = zeros(size(A_values));
T_c_out_values = zeros(size(A_values));

for i = 1:length(A_values)
    A = A_values(i);
    NTU = U * A / C_min;
    epsilon_guess = (1 - exp(-NTU * (1 - C_r))) / (1 - C_r * exp(-NTU * (1 - C_r)));
    Q = epsilon_guess * Q_max;
    T_h_out = T_h_in - Q / C_h;
    T_c_out = T_c_in + Q / C_c;
    DT1 = T_h_in - T_c_out;
    DT2 = T_h_out - T_c_in;
    if abs(DT1 - DT2) < 1e-6
        LMTD = DT1;
    else
        LMTD = (DT1 - DT2) / log(DT1 / DT2);
    end
    Q_LMTD = U * A * LMTD;
    for j = 1:5
        Q = Q_LMTD;
        T_h_out = T_h_in - Q / C_h;
        T_c_out = T_c_in + Q / C_c;
        DT1 = T_h_in - T_c_out;
        DT2 = T_h_out - T_c_in;
        if abs(DT1 - DT2) < 1e-6
            LMTD = DT1;
        else
            LMTD = (DT1 - DT2) / log(DT1 / DT2);
        end
        Q_LMTD = U * A * LMTD;
        if abs(Q - Q_LMTD) < 1e-3
            break;
        end
    end
    epsilon = Q / Q_max;
    epsilon_values(i) = epsilon;
    Q_values(i) = Q;
    T_h_out_values(i) = T_h_out;
    T_c_out_values(i) = T_c_out;
end

% Select the best A (maximum epsilon)
[best_epsilon, best_idx] = max(epsilon_values);
best_A = A_values(best_idx);
best_Q = Q_values(best_idx);
best_T_h_out = T_h_out_values(best_idx);
best_T_c_out = T_c_out_values(best_idx);

% Recalculate LMTD for best A
A = best_A;
DT1 = T_h_in - best_T_c_out;
DT2 = best_T_h_out - T_c_in;
if abs(DT1 - DT2) < 1e-6
    LMTD = DT1;
else
    LMTD = (DT1 - DT2) / log(DT1 / DT2);
end

% Display Optimization Results
fprintf('\nOptimization Results:\n');
fprintf('--------------------------------\n');
for i = 1:length(A_values)
    fprintf('Area = %.1f m², Q = %.2f W, Effectiveness = %.3f, T_h_out = %.2f °C, T_c_out = %.2f °C\n', ...
        A_values(i), Q_values(i), epsilon_values(i), T_h_out_values(i), T_c_out_values(i));
end
fprintf('\nBest Configuration:\n');
fprintf('Area = %.1f m², Q = %.2f W, Effectiveness = %.3f, T_h_out = %.2f °C, T_c_out = %.2f °C\n', ...
    best_A, best_Q, best_epsilon, best_T_h_out, best_T_c_out);

% Plot Effectiveness vs. Area
figure;
plot(A_values, epsilon_values, 'b-o', 'DisplayName', 'Effectiveness');
xlabel('Heat Transfer Area (m²)');
ylabel('Effectiveness');
title('Effectiveness vs. Heat Transfer Area');
legend('Location', 'best');
grid on;
saveas(gcf, 'optimization_plot.png');

% Step 10: Fluid Flow Analysis - Calculate Pressure Drop for Air
rho_c = 1.2;         % Density of air (kg/m³)
A_cross = 0.02;      % Cross-sectional area of air flow path (m²)
V_c = m_c / (rho_c * A_cross);  % Velocity of air (m/s)
f = 0.02;            % Friction factor (assumed)
L = 1.0;             % Length of exchanger (m)
D = 0.01;            % Hydraulic diameter (m)
delta_P = f * (L / D) * (rho_c * V_c^2 / 2);  % Pressure drop (Pa)

fprintf('\nFluid Flow Analysis:\n');
fprintf('--------------------------------\n');
fprintf('Air Velocity: %.2f m/s\n', V_c);
fprintf('Pressure Drop (Air): %.2f Pa (%.2f kPa)\n', delta_P, delta_P/1000);

% Step 11: Plot Temperature Profiles (Using Optimized Parameters)
x = 0:0.01:1;  % Normalized length (0 to 1)
T_h = T_h_in - (T_h_in - best_T_h_out) * x;  % Hot fluid: T_h_in to best_T_h_out
T_c = T_c_in + (best_T_c_out - T_c_in) * (1 - x);  % Cold fluid: best_T_c_out at x=0 to T_c_in at x=1

figure;
plot(x, T_h, 'r-', 'DisplayName', 'Hot Fluid (Water)'); hold on;
plot(x, T_c, 'b-', 'DisplayName', 'Cold Fluid (Air)');
xlabel('Normalized Length');
ylabel('Temperature (°C)');
title('Temperature Profiles in Counter-Flow Heat Exchanger (Optimized)');
legend('Location', 'best');
grid on;
hold off;

% Save the updated plot
saveas(gcf, 'temperature_profile_optimized.png');