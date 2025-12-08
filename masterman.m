function dvdt = quartercar(t, v_state)
    % v_state är vektorn [z1; z2; z1_prick; z2_prick]

    % Parametrar från Tabell 1 
    m1 = 475;
    m2 = 53;
    k1 = 5400;
    k2 = 135000;
    c1 = 310;
    c2 = 1200;
    vel = 65/3.6; 
    H = 0.24;
    L = 1;

    z1 = v_state(1);
    z2 = v_state(2);
    z1_dot = v_state(3);
    z2_dot = v_state(4);

    % Beräkna vägprofilen h(t) och h_prick(t) enligt Ekv (2) 
    if t <= L/vel
        % disp("Hej")
        arg = (2*pi*vel*t) / L;
        h = (H/2)*(1 - cos(arg));
        h_dot = (H/2) * (2*pi*vel/L) * sin(arg);
    else
        h = 0;
        h_dot = 0;
    end

    % Rörelseekvationer baserade på Ekv (39) och (40) 
    % Vi löser ut accelerationerna (z1_dubbelprick och z2_dubbelprick)

    size(c1);
    % Ekv för m1: m1z1_dd + c1(z1_d - z2_d) + k1(z1 - z2) = 0
    z1_ddot = (-c1*(z1_dot - z2_dot) - k1*(z1 - z2)) / m1;


    % Ekv för m2: m2z2_dd + c1(z2_d - z1_d) + c2(z2_d - h_d) + k1(z2 - z1) + k2(z2 - h) = 0
    z2_ddot = (-c1*(z2_dot - z1_dot) - c2*(z2_dot - h_dot) - k1*(z2 - z1) - k2*(z2 - h)) / m2;

    % Returnera derivatorna [z1_prick; z2_prick; z1_dubbelprick; z2_dubbelprick]
    dvdt = [z1_dot; z2_dot; z1_ddot; z2_ddot];
end

%dvdt = quartercar(0,[10;20;30;40])

% Startvärden (vila) och tidsintervall
v0 = [0;0;0;0]; % [z1; z2; z1_dot; z2_dot] [cite: 73]
T_end = 10; % Välj en tid lång nog för att se svängningarna klinga av [cite: 74]

dt = 5e-3;
t_eu = 0:dt:T_end;
n_steps = length(t_eu);

% Allokera minne
v_eu = zeros(4, n_steps);
v_eu(:, 1) = v0;

% Euler-loopen
for k = 1:n_steps-1;
    dvdt = quartercar(t_eu(k), v_eu(:, k));
    v_eu(:, k+1) = v_eu(:, k) + dt * dvdt;
end
v_eu
plot(t_eu,v_eu(2,:),'o-')