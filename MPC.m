%% MPC
%% Parameters
dt = 0.05; 
N  = 30;  
T  = 100; 
nx = 3; % [px, py, theta]
nu = 2; % [v, w]

% Limits
v_max = 10.0; 
w_max = 10.0; 

% Waypoints
% Waypoints extraction from pathIds
num_wp = length(pathIds);
waypoints = zeros(num_wp, 2);

% İlk hedef başlangıç noktanız
waypoints(1, :) = q_start;

for i = 1:(num_wp - 1)
    curr_node = nodes(pathIds(i)).poly;
    next_node = nodes(pathIds(i+1)).poly;
    
    % Hedef = İki düğümün kesişiminin merkezi (Her iki poligonun da içindedir!)
    intersection_poly = intersect(curr_node, next_node);
    [cx, cy] = centroid(intersection_poly);
    waypoints(i+1, :) = [cx, cy];
end

% Son hedef varış noktanız
waypoints(end+1, :) = q_goal; 
% --------------------------------------------------

wp_tol = 0.5; 
current_wp_idx = 2; 
target = waypoints(current_wp_idx, :)'; 
x_true = [waypoints(1,1); waypoints(1,2); 0];

% Storage
history_x = zeros(nx, T);
history_u = zeros(nu, T);

options = optimoptions('quadprog','Display','off');

%% Simulation Loop
for t = 1:T
    % --- 1. Waypoint Logic (UPDATED for Stopping) ---
    dist_to_target = norm(x_true(1:2) - target(1:2));
    
    % Son waypoint'te miyiz kontrolü
    is_last_waypoint = (current_wp_idx == size(waypoints, 1));
    
    if dist_to_target < wp_tol
        if is_last_waypoint
            % Eğer son noktadaysak: BREAK YAPMA. Sadece bekle.
            disp(['Final destination reached. Station Keeping... (Step: ' num2str(t) ')']);
            
            % Daha hassas duruş için toleransı daraltabiliriz (Opsiyonel)
            wp_tol = 0.5; 
            
            % Eğer tamamen durduysa (Hız ~ 0 ve Konum < 0.1) simülasyonu bitir
            if dist_to_target < 0.1 && norm(history_u(:, t-1)) < 0.05
                 disp('Robot has fully stopped. Simulation End.');
                 break;
            end
        else
            % Son nokta değilse, bir sonrakine geç
            disp(['Waypoint ' num2str(current_wp_idx) ' reached!']);
            current_wp_idx = current_wp_idx + 1;
            target = waypoints(current_wp_idx, :)';
        end
    end

    % --- 2. Linearization (LTV) ---
    th = x_true(3);
    v_ref = 1.0; 
    
    A = [1 0 -v_ref*sin(th)*dt;
         0 1  v_ref*cos(th)*dt;
         0 0  1];
    B = [cos(th)*dt  0;
         sin(th)*dt  0;
         0           dt];

   % --- 3. Local Coordinate Shift (Karagoz Method) ---
    xr = target(1);
    yr = target(2);
    
    dx = xr - x_true(1);
    dy = yr - x_true(2);
    raw_target_theta = atan2(dy, dx);
    
    % Yerel duruma (Local State) geçiş: x_bar = x - x_ref
    x_bar = zeros(3, 1);
    x_bar(1) = x_true(1) - xr;
    x_bar(2) = x_true(2) - yr;
    
    if dist_to_target < 0.5
        x_bar(3) = 0; 
    else
        delta_theta = raw_target_theta - x_true(3);
        x_bar(3) = -atan2(sin(delta_theta), cos(delta_theta)); 
    end
    
    % --- 4. MPC Matrices ---
    Q = diag([10, 10, 0.5]); 
    R = diag([0.1, 0.1]);
    
    Phi = zeros(nx*N, nx);
    Gamma = zeros(nx*N, nu*N);
    temp_A = A;
    for i = 1:N
        Phi((i-1)*nx+1:i*nx, :) = temp_A;
        for j = 1:i
             if j == i, mult = B; else, mult = (A^(i-j))*B; end
             Gamma((i-1)*nx+1:i*nx, (j-1)*nu+1:j*nu) = mult;
        end
        temp_A = temp_A * A;
    end
    
    Q_bar = kron(eye(N), Q);
    R_bar = kron(eye(N), R);
    
    % Yerel çerçevede referans noktası orijindir (0,0,0)
    Ref = zeros(nx*N, 1);
    
    H_qp = Gamma' * Q_bar * Gamma + R_bar;
    f_qp = (x_bar' * Phi' * Q_bar * Gamma - Ref' * Q_bar * Gamma)';
    
    % --- 5. Funnel (Polygon) Constraints ---
    % Aracın içinde bulunduğu node'un ID'sini bul
    current_node_id = pathIds(current_wp_idx - 1);
    curr_poly = nodes(current_node_id).poly;
    
    % Poligon sınırlarını [A_poly * x <= b_poly] formatına çevir
    [vx, vy] = boundary(curr_poly);
    if isnan(vx(end)), vx(end)=[]; vy(end)=[]; end
    if vx(1)==vx(end) && vy(1)==vy(end), vx(end)=[]; vy(end)=[]; end
    
    num_edges = length(vx);
    A_poly = zeros(num_edges, 2);
    b_poly = zeros(num_edges, 1);
    node_centroid = nodes(current_node_id).c;
    
    for k = 1:num_edges
        p1 = [vx(k), vy(k)];
        p2 = [vx(mod(k, num_edges)+1), vy(mod(k, num_edges)+1)];
        edge_vec = p2 - p1;
        normal = [edge_vec(2), -edge_vec(1)]; % Normal vektör
        
        % Dışa doğru baktığından emin ol
        if dot(normal, node_centroid - p1) > 0
            normal = -normal; 
        end
        normal = normal / norm(normal);
        
        A_poly(k, :) = normal;
        b_poly(k) = dot(normal, p1);
    end
    
    % Kısıtlamaları N adımlık MPC ufku için genişlet
    % Yerel eksen kayması: A_poly * (x_bar + x_r) <= b_poly => A_poly * x_bar <= b_poly - A_poly * x_r
    A_rep = kron(eye(N), A_poly);
    b_rep = repmat(b_poly - A_poly * [xr; yr], N, 1);
    
    % Phi ve Gamma'dan sadece (x,y) satırlarını çıkaran filtre matrisi
    C_pos = kron(eye(N), [1 0 0; 0 1 0]); 
    
    % Eşitsizlik Kısıtlamaları: A_ineq * U <= b_ineq
    A_ineq = A_rep * C_pos * Gamma;
    b_ineq = b_rep - A_rep * C_pos * Phi * x_bar;
    
    % Fiziksel Hız Kısıtlamaları
    lb = repmat([-v_max; -w_max], N, 1);
    ub = repmat([ v_max;  w_max], N, 1);
    
    % --- 6. Solve ---
    [U_opt, ~, exitflag] = quadprog(H_qp, f_qp, A_ineq, b_ineq, [], [], lb, ub, [], options);
    
    if exitflag ~= 1
        disp(['QP Infeasible at step ' num2str(t) ' - applying brakes.']);
        u_apply = [0;0];
    else
        u_apply = U_opt(1:2);
    end
    
    % !!! BU KISIM EKSİKTİ - GERİ EKLENDİ !!!
    % --- 7. Apply to Real Plant ---
    x_true(1) = x_true(1) + u_apply(1) * cos(x_true(3)) * dt;
    x_true(2) = x_true(2) + u_apply(1) * sin(x_true(3)) * dt;
    x_true(3) = x_true(3) + u_apply(2) * dt;
    x_true(3) = atan2(sin(x_true(3)), cos(x_true(3)));
    
    history_x(:, t) = x_true;
    history_u(:, t) = u_apply;
    % !!! ------------------------------- !!!
    
end % <-- For döngüsü burada bitmeli

% Truncate
if t < T
    final_step = t - 1;
else
    final_step = T;
end

history_x = history_x(:, 1:final_step);
history_u = history_u(:, 1:final_step);

%% Plotting

% --- FİGÜR 2: Dijkstra Çıktısı (Sadece kullanılan düğümler ve Waypoint'ler) ---
figure('WindowState','maximized', 'Color','w'); hold on; axis equal;
xlim([W(1) W(2)]); ylim([W(3) W(4)]);
title('Figür 2: Dijkstra - Seçilen Yol ve Kesişim Waypointleri');

% Engelleri çiz
for i=1:numel(obs)
    plot(obs{i}, 'FaceColor',[0 0 0], 'FaceAlpha',0.6, 'EdgeColor','none');
end

if ~isempty(pathIds)
    % SADECE yoldaki düğümleri belirgin bir renkle (ör: turuncu) çiz
    for k = 1:length(pathIds)
        node_idx = pathIds(k);
        plot(nodes(node_idx).poly, 'FaceColor',[1.0 0.85 0.7], 'FaceAlpha',0.40, 'EdgeColor',[1.0 0.5 0.0], 'LineWidth',1.5);
    end
    
    num_nodes = length(pathIds);
    route_points = zeros(num_nodes + 1, 2); 
    route_points(1, :) = q_start; 
    
    for k = 1:(num_nodes - 1)
        curr_node = nodes(pathIds(k)).poly;
        next_node = nodes(pathIds(k+1)).poly;
        
        intersection_poly = intersect(curr_node, next_node);
        [cx, cy] = centroid(intersection_poly);
        route_points(k+1, :) = [cx, cy];
    end
    
    route_points(end, :) = q_goal; 
    
    % Kesişim noktalarını ve aralarındaki yolu çiz
    plot(route_points(:,1), route_points(:,2), 'bo', 'MarkerSize',6, 'LineWidth',1.5);
end

% Başlangıç ve bitiş
plot(q_start(1), q_start(2), 'go', 'MarkerSize',9, 'LineWidth',2);
plot(q_goal(1),  q_goal(2),  'ro', 'MarkerSize',9, 'LineWidth',2);
grid on;

hold on
% plot(waypoints(:,1), waypoints(:,2), 'k--s', 'LineWidth', 1.5, 'MarkerSize', 8); hold on;
plot(history_x(1,:), history_x(2,:), 'r-', 'LineWidth', 2);
title('MPC Waypoint Tracking (Fixed)');
legend('Waypoints', 'Robot Path');
grid on; axis equal;
xlabel('X [m]'); ylabel('Y [m]');


t = (0:size(history_u,2)-1) * dt;

figure('Color','w');
subplot(2,1,1);
plot(t, history_u(1,:), 'b', 'LineWidth',1.5);
grid on;
ylabel('u');
title('Control inputs');

subplot(2,1,2);
plot(t, history_u(2,:), 'r', 'LineWidth',1.5);
grid on;
xlabel('Time [s]');
ylabel('w');