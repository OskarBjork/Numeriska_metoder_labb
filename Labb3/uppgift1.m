% Laboration 3 - Keplerproblemet och Numeriska Metoder för Hamiltonska System
% Denna kod implementerar fyra numeriska metoder: Framåt Euler, Bakåt Euler, 
% Symplektisk Euler samt den Implicita Mittpunktsmetoden.
clear; clc; close all;

% Fysikaliska och numeriska parametrar
a = 0.5; % Excentricitet som styr banans form
t_end = 100; % Total simuleringstid
% Begynnelsevillkor: q1, q2, p1, p2 (kroppen befinner sig vid perihelium)
y0 = [1-a; 0; 0; sqrt((1+a)/(1-a))]; 

% Tidssteg och motsvarande antal steg (enligt rekommendation för stabilitet)
h_euler = 0.0005; 
N_euler = round(t_end / h_euler);
h_mid   = 0.025;  
N_mid   = round(t_end / h_mid);

% Allokering av minne för de fyra metoderna för att optimera exekveringstid
y_fe = zeros(4, N_euler+1); y_fe(:,1) = y0; % Framåt Euler
y_be = zeros(4, N_euler+1); y_be(:,1) = y0; % Bakåt Euler
y_se = zeros(4, N_euler+1); y_se(:,1) = y0; % Symplektisk Euler
y_mp = zeros(4, N_mid+1);   y_mp(:,1) = y0; % Mittpunktsmetoden

% Definition av det första ordningens vektorfält f(y)
% norm(y(1:2)) beräknar |q| = sqrt(q1^2 + q2^2)
% y = (q1,q2,p1,p2) = (q1,q2, q1_p, q2_p)
f = @(y) [y(3); % p1 = q1_p
          y(4); % p2 = q2_p
          -y(1)/norm(y(1:2))^3; % p1_p
          -y(2)/norm(y(1:2))^3]; % p2_p

% Definition av Hamiltonfunktionen H(p,q) för beräkning av energi
energy = @(y) 0.5*(y(3)^2 + y(4)^2) - 1/norm(y(1:2));

% Analytisk Jacobian av f(y) för användning i Newton-Raphson.
% Detta är kritiskt för att den implicita lösaren ska konvergera snabbt.
Jf = @(y) [0, 0, 1, 0; 
           0, 0, 0, 1; 
           (2*y(1)^2 - y(2)^2)/norm(y(1:2))^5, 3*y(1)*y(2)/norm(y(1:2))^5, 0, 0; 
           3*y(1)*y(2)/norm(y(1:2))^5, (2*y(2)^2 - y(1)^2)/norm(y(1:2))^5, 0, 0];

% =========================================================================
% Tidsstegning för de första ordningens metoderna (Euler)
% =========================================================================
for n = 1:N_euler
    % 1. Framåt Euler (Explicit)
    % Värdet vid n+1 bestäms enbart av värdet och lutningen vid n.
    y_fe(:, n+1) = y_fe(:, n) + h_euler * f(y_fe(:, n));
    
    % 2. Symplektisk Euler (Semi-implicit)
    % Uppdaterar först rörelsemängden p, och använder sedan detta nya p 
    % för att uppdatera positionen q. Detta bevarar den symplektiska formen.
    q_n = y_se(1:2, n); 
    p_n = y_se(3:4, n);
    p_next = p_n - h_euler * (q_n / norm(q_n)^3);
    q_next = q_n + h_euler * p_next;
    y_se(:, n+1) = [q_next; p_next];
    
    % 3. Bakåt Euler (Implicit)
    % Använder Newtons metod för att hitta roten till F(y^{n+1}) = 0
    y_guess = y_be(:, n) + h_euler * f(y_be(:, n)); % Explicit Euler som startgissning
    for iter = 1:10
        % Beräkna residualen F
        F_val = y_guess - y_be(:, n) - h_euler * f(y_guess);
        % Beräkna Jacobianen J_F
        J_val = eye(4) - h_euler * Jf(y_guess);
        % Lös linjärt system för uppdateringssteget delta
        delta = J_val \ F_val;
        y_guess = y_guess - delta;
        % Avbryt iterationen om toleransen är uppnådd
        if norm(delta) < 1e-8, break; end
    end
    y_be(:, n+1) = y_guess;
end

% =========================================================================
% Tidsstegning för Mittpunktsmetoden (Andra ordningen)
% =========================================================================
for n = 1:N_mid
    % 4. Implicit Mittpunktsmetod
    % Integrerar via vektorfältet utvärderat vid medelvärdet av y^n och y^{n+1}
    y_guess = y_mp(:, n) + h_mid * f(y_mp(:, n)); % Startgissning
    for iter = 1:10
        y_mid = 0.5 * (y_mp(:, n) + y_guess);
        F_val = y_guess - y_mp(:, n) - h_mid * f(y_mid);
        % Notera att h multipliceras med 0.5 på grund av kedjeregeln 
        % från medelvärdesbildningen inuti f(...)
        J_val = eye(4) - 0.5 * h_mid * Jf(y_mid);
        delta = J_val \ F_val;
        y_guess = y_guess - delta;
        if norm(delta) < 1e-8, break; end
    end
    y_mp(:, n+1) = y_guess;
end

% =========================================================================
% Visualisering av fasrumsbanor och energidynamik
% =========================================================================
figure('Name', 'Planetbanor i q1-q2 planet');
subplot(2,2,1); plot(y_fe(1,:), y_fe(2,:)); title('Framåt Euler'); axis equal; grid on;
subplot(2,2,2); plot(y_be(1,:), y_be(2,:)); title('Bakåt Euler'); axis equal; grid on;
subplot(2,2,3); plot(y_se(1,:), y_se(2,:)); title('Symplektisk Euler'); axis equal; grid on;
subplot(2,2,4); plot(y_mp(1,:), y_mp(2,:)); title('Mittpunktsmetoden'); axis equal; grid on;

% Beräkning av energi för alla tidssteg
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