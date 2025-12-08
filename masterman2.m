b = 0.0107807258388
c = 0.000111810524139

function dvdt = quartercar(t, v_state)
    % v_state är vektorn [z1; z2; z1_prick; z2_prick]
 
    m1 = 475;
    m2 = 53;
    k1 = 5400*0+690;
    k2_ref = 135000;
    k2 = 100*k2_ref*0 + 182120;
    c1 = 310;
    c2 = 1200;
    vel = 65/3.6; 
    H = 0.24;
    L = 1;

    z1 = v_state(1);
    z2 = v_state(2);
    z1_dot = v_state(3);
    z2_dot = v_state(4);
  
    if t <= L/vel
        arg = (2*pi*vel*t) / L;
        h = (H/2)*(1 - cos(arg));
        h_dot = (H/2) * (2*pi*vel/L) * sin(arg);
    else
        h = 0;
        h_dot = 0;
    end


    z1_ddot = (-c1*(z1_dot - z2_dot) - k1*(z1 - z2)) / m1;

    z2_ddot = (-c1*(z2_dot - z1_dot) - c2*(z2_dot - h_dot) - k1*(z2 - z1) - k2*(z2 - h)) / m2;

    dvdt = [z1_dot; z2_dot; z1_ddot; z2_ddot];
end

%tspan = linspace(0,1,100)
tspan = [0 10]
v0 = [0 0 0 0]

options = odeset("RelTol", 1e-6, 'Refine',1, 'InitialStep',0.001)

[t,v] = ode45(@quartercar,tspan,v0, options)
%plot(t,v(:,1),'o-')
%plot(t,v(:,1),'o-',t,v(:,2),'-o')
%plot(t,v(:,1),'o-',t,v(:,2),'-o',t,v(:,3),'o-',t,v(:,4),'-o',LineWidth=.1)

t
v_m = [max(v(:,1)), max(v(:,2))]

v0 = [0;0;0;0]; % [z1; z2; z1_dot; z2_dot] [cite: 73]
T_end = 10; % Välj en tid lång nog för att se svängningarna klinga av [cite: 74]

%dt = 5e-3;
dt = .1*c
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
%plot(t_eu,v_eu(2,:),'o-',t,v(:,2),'o-',)
plot(t_eu,v_eu(2,:),'o-',t_eu,v_eu(1,:),'o-')
legend('z2','z1')

