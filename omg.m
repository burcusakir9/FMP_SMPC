%% Otter Square Wave Input Simulation
clear; clc; close all;

% Zaman ayarları
dt = 0.01;
T_total = 20; % Toplam 20 saniye
t_vec = 0:dt:T_total;
steps = length(t_vec);

%% 1. Çalışma Noktası ve Model Kurulumu
% Modeli durur vaziyette (u=0) ve n=0 etrafında lineerleştiriyoruz
x_eq = [0; 0; 0; 0; 0; 0]; 
u_eq = [0; 0];

[A_c, B_c] = linearize_otter3d(x_eq, u_eq);
sys_d = c2d(ss(A_c, B_c, eye(6), 0), dt);
A = sys_d.A;
B = sys_d.B;

%% 2. Kare Dalga Giriş Tanımlama
% İlk 10 sn: n1=100, n2=0
% Son 10 sn: n1=0, n2=100
u_input = zeros(2, steps);
for k = 1:steps
    if t_vec(k) <= 10
        u_input(:, k) = [100; 0];
    else
        u_input(:, k) = [0; 100];
    end
end

%% 3. Simülasyon Döngüsü
x_lin = x_eq;
x_nonlin = x_eq;

history_lin = zeros(6, steps);
history_nonlin = zeros(6, steps);

for k = 1:steps
    history_lin(:, k) = x_lin;
    history_nonlin(:, k) = x_nonlin;
    
    u_k = u_input(:, k);
    
    % Lineer Model
    x_lin = x_eq + A * (x_lin - x_eq) + B * (u_k - u_eq);
    
    % Nonlinear Model (Gerçek)
    [xdot, ~] = otter3d(x_nonlin, u_k);
    x_nonlin = x_nonlin + dt * xdot;
    x_nonlin(3) = atan2(sin(x_nonlin(3)), cos(x_nonlin(3))); 
end

%% 4. Grafikleme
figure('Color', 'w', 'Position', [100, 100, 1200, 800]);

% --- Kartezyen Pozisyon (X-Y) ---
subplot(3, 2, 1);
plot(history_nonlin(1,:), history_nonlin(2,:), 'b', 'LineWidth', 2); hold on;
plot(history_lin(1,:), history_lin(2,:), 'r--', 'LineWidth', 1.5);
xlabel('X [m]'); ylabel('Y [m]'); title('Cartesian Position (X-Y)');
legend('Nonlinear', 'Linear'); grid on; axis equal;

% --- Girişler (n1, n2) ---
subplot(3, 2, 2);
plot(t_vec, u_input(1,:), 'k', 'LineWidth', 1.5); hold on;
plot(t_vec, u_input(2,:), 'r', 'LineWidth', 1.5);
ylabel('Propeller Speed [rad/s]'); title('Inputs');
legend('n1 (Left)', 'n2 (Right)'); grid on;

% --- Hızlar (u, v) ---
subplot(3, 2, 3);
plot(t_vec, history_nonlin(4,:), 'b'); hold on;
plot(t_vec, history_nonlin(5,:), 'g');
plot(t_vec, history_lin(4,:), 'b--', 'LineWidth', 0.5);
plot(t_vec, history_lin(5,:), 'g--', 'LineWidth', 0.5);
ylabel('Velocity [m/s]'); title('Surge (u) and Sway (v)');
legend('u (Nonlin)', 'v (Nonlin)', 'u (Lin)', 'v (Lin)'); grid on;

% --- Yaw ve Yaw Rate (psi, r) ---
subplot(3, 2, 4);
plot(t_vec, history_nonlin(3,:), 'k'); hold on;
plot(t_vec, history_lin(3,:), 'r--');
ylabel('Yaw [rad]'); title('Heading (psi)'); grid on;

subplot(3, 2, 5);
plot(t_vec, history_nonlin(6,:), 'm'); hold on;
plot(t_vec, history_lin(6,:), 'k--');
ylabel('r [rad/s]'); xlabel('Time [s]'); title('Yaw Rate (r)'); grid on;

% --- Hata Analizi (Distance Error) ---
subplot(3, 2, 6);
err = sqrt((history_nonlin(1,:)-history_lin(1,:)).^2 + (history_nonlin(2,:)-history_lin(2,:)).^2);
plot(t_vec, err, 'r', 'LineWidth', 1.5);
ylabel('Error [m]'); xlabel('Time [s]'); title('X-Y Position Error'); grid on;