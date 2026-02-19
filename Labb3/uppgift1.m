clear; clc; close all;

a = 0.5; 
t_end = 100;
% Begynnelsevillkor: q1, q2, p1, p2 
y0 = [1-a; 0; 0; sqrt((1+a)/(1-a))]; 

% Tidssteg och motsvarande antal steg
h_euler = 0.0005; 
N_euler = round(t_end / h_euler);
h_mid   = 0.025;  
N_mid   = round(t_end / h_mid);

y_fe = zeros(4, N_euler+1); y_fe(:,1) = y0; 
y_be = zeros(4, N_euler+1); y_be(:,1) = y0; 
y_se = zeros(4, N_euler+1); y_se(:,1) = y0; 
y_mp = zeros(4, N_mid+1);   y_mp(:,1) = y0; 

% Definition av första ordningens vektorfält f(y)
% y = (q1,q2,p1,p2) = (q1,q2, q1_p, q2_p)
f = @(y) [y(3); 
          y(4); 
          -y(1)/norm(y(1:2))^3; 
          -y(2)/norm(y(1:2))^3]; 

% Definition av Hamiltonfunktionen
energy = @(y) 0.5*(y(3)^2 + y(4)^2) - 1/norm(y(1:2));

% Jacobianen fr fr
Jf = @(y) [0, 0, 1, 0; 
           0, 0, 0, 1; 
           (2*y(1)^2 - y(2)^2)/norm(y(1:2))^5, 3*y(1)*y(2)/norm(y(1:2))^5, 0, 0; 
           3*y(1)*y(2)/norm(y(1:2))^5, (2*y(2)^2 - y(1)^2)/norm(y(1:2))^5, 0, 0];

for n = 1:N_euler
    % 1. Framåt Euler (Explicit)
    y_fe(:, n+1) = y_fe(:, n) + h_euler * f(y_fe(:, n));
    
    % 2. Symplektisk Euler (Semi-implicit)
    % Uppdaterar först rörelsemängden p, och använder sedan detta nya p 
    % för att uppdatera positionen q.
    q_n = y_se(1:2, n); 
    p_n = y_se(3:4, n);
    p_next = p_n - h_euler * (q_n / norm(q_n)^3);
    q_next = q_n + h_euler * p_next;
    y_se(:, n+1) = [q_next; p_next];
    
    % 3. Bakåt Euler (Implicit)
    y_guess = y_be(:, n) + h_euler * f(y_be(:, n)); % Explicit Euler som startgissning
    for iter = 1:10
        F_val = y_guess - y_be(:, n) - h_euler * f(y_guess);
        J_val = eye(4) - h_euler * Jf(y_guess);
        delta = J_val \ F_val;
        y_guess = y_guess - delta;
        if norm(delta) < 1e-8, break; end
    end
    y_be(:, n+1) = y_guess;
end

for n = 1:N_mid
    % 4. Implicit Mittpunktsmetod
    y_guess = y_mp(:, n) + h_mid * f(y_mp(:, n));
    for iter = 1:10
        y_mid = 0.5 * (y_mp(:, n) + y_guess);
        F_val = y_guess - y_mp(:, n) - h_mid * f(y_mid);
        J_val = eye(4) - 0.5 * h_mid * Jf(y_mid);
        delta = J_val \ F_val;
        y_guess = y_guess - delta;
        if norm(delta) < 1e-8, break; end
    end
    y_mp(:, n+1) = y_guess;
end

figure('Name', 'Planetbanor i q1-q2 planet');
subplot(2,2,1); plot(y_fe(1,:), y_fe(2,:)); title('Framåt Euler'); axis equal; grid on;
subplot(2,2,2); plot(y_be(1,:), y_be(2,:)); title('Bakåt Euler'); axis equal; grid on;
subplot(2,2,3); plot(y_se(1,:), y_se(2,:)); title('Symplektisk Euler'); axis equal; grid on;
subplot(2,2,4); plot(y_mp(1,:), y_mp(2,:)); title('Mittpunktsmetoden'); axis equal; grid on;

% Beräkning av energi
E_fe = arrayfun(@(i) energy(y_fe(:,i)), 1:N_euler+1);
E_be = arrayfun(@(i) energy(y_be(:,i)), 1:N_euler+1);
E_se = arrayfun(@(i) energy(y_se(:,i)), 1:N_euler+1);
E_mp = arrayfun(@(i) energy(y_mp(:,i)), 1:N_mid+1);

figure('Name', 'Energins tidsutveckling'); hold on; grid on;
plot(linspace(0, t_end, N_euler+1), E_fe, 'r', 'DisplayName', 'Framåt Euler');
plot(linspace(0, t_end, N_euler+1), E_be, 'b', 'DisplayName', 'Bakåt Euler');
plot(linspace(0, t_end, N_euler+1), E_se, 'g', 'DisplayName', 'Symplektisk Euler');
plot(linspace(0, t_end, N_mid+1),   E_mp, 'k', 'DisplayName', 'Mittpunktsmetoden');
xlabel('Tid t'); ylabel('Energi H(p,q)');
legend('show', 'Location', 'best');
title('Total energi över tid för Keplerproblemet');