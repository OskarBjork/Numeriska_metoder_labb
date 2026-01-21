% --- Initialisering ---
clc; clear; close all;

% Startvärden (vila) och tidsintervall
v0 = [0; 0; 0; 0]; % [z1; z2; z1_dot; z2_dot] [cite: 73]
T_end = 2.0; % Välj en tid lång nog för att se svängningarna klinga av [cite: 74]

% --- Lösning med ode45 (Uppgift 2a & 2b) ---
% Inställningar: RelTol 1e-6 och Refine 1 för att se faktiska tidssteg [cite: 75, 79]
options = odeset('RelTol', 1e-6, 'Refine', 1); 

[t_ode, y_ode] = ode45(@quartercar, [0 T_end], v0, options);

z1_ode = y_ode(:, 1);
z2_ode = y_ode(:, 2);

% Plotta ode45 resultat (Uppgift 2a)
figure(1);
plot(t_ode, z1_ode, 'b', 'LineWidth', 1.5); hold on;
plot(t_ode, z2_ode, 'r', 'LineWidth', 1.5);
legend('z_1 (Chassi)', 'z_2 (Hjul)');
title('Kvartsmodell simulering (ode45)');
xlabel('Tid [s]'); ylabel('Position [m]');
grid on;

% Visualisera tidssteg (Uppgift 2b)
dt_ode = diff(t_ode); % Skillnad mellan tidpunkter
figure(2);
plot(t_ode(1:end-1), dt_ode, '.-');
title('Tidsstegsstorlek för ode45 över tid');
xlabel('Tid [s]'); ylabel('\Delta t [s]');
grid on;

% --- Lösning med Euler Framåt (Uppgift 2c) ---
dt_vals = [5e-3, 5e-4]; % Tidssteg att testa [cite: 84]
colors = {'g--', 'k:'};

figure(3);
plot(t_ode, z2_ode, 'r', 'LineWidth', 1.5); hold on; % Referens från ode45
legend_str = {'ode45 (z_2)'};

for i = 1:length(dt_vals)
    dt = dt_vals(i);
    t_eu = 0:dt:T_end;
    n_steps = length(t_eu);
    
    % Allokera minne
    v_eu = zeros(4, n_steps);
    v_eu(:, 1) = v0;
    
    % Euler-loopen
    for k = 1:n_steps-1
        dvdt = quartercar(t_eu(k), v_eu(:, k));
        v_eu(:, k+1) = v_eu(:, k) + dt * dvdt;
    end
    
    % Plotta z2 för jämförelse
    plot(t_eu, v_eu(2, :), colors{i}, 'LineWidth', 1.2);
    legend_str{end+1} = sprintf('Euler \\Delta t = %.0e', dt);
end

title('Jämförelse: Euler vs ode45 (Hjulposition z_2)');
xlabel('Tid [s]'); ylabel('Position z_2 [m]');
legend(legend_str);
grid on;
ylim([-0.1 0.4]); % Begränsa y-axel om Euler spårar ur