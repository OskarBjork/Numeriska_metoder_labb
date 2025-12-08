function trapezoid2()
    
    p.m1 = 475; p.m2 = 53; p.k1 = 5400; p.c1 = 310; p.c2 = 1200;
    p.k2 = 100 * 135000; 
    p.v = 65/3.6; p.H = 0.24; p.L = 1;
    
    A = [0, 0, 1, 0;
         0, 0, 0, 1;
         -p.k1/p.m1, p.k1/p.m1, -p.c1/p.m1, p.c1/p.m1;
         p.k1/p.m2, -(p.k1+p.k2)/p.m2, p.c1/p.m2, -(p.c1+p.c2)/p.m2];
    
    v0 = [0;0;0;0];
    T_end = 0.05; 

    opts_ref = odeset('RelTol', 1e-10, 'AbsTol', 1e-11);
    ode_fun = @(t, v) system_dynamics_explicit(t, v, A, p);
    [t_ref, y_ref] = ode45(ode_fun, [0 T_end], v0, opts_ref);
    
    dt0 = 1e-4; 
    alphas = [1, 0.5, 0.25, 0.125];
    errors = zeros(size(alphas));

    [t,v] = trapezoidal_solver(A, p, v0, dt0, T_end);

    plot(t,v(2,:))
    
    % fprintf('\n--- Uppgift 4b: Noggranhet trapetsmetoden ---\n');
    % for i = 1:length(alphas)
    %     dt = dt0 * alphas(i);
    %     [t_trap, y_trap] = trapezoidal_solver(A, p, v0, dt, T_end);
        
    %     % Beräkna fel (jämfört med referens vid sluttidpunkten)
    %     % Interpolera referenslösning till sluttidpunkten exakt
    %     v_ref_end = interp1(t_ref, y_ref, T_end, 'spline');
    %     v_trap_end = y_trap(:,end)';
        
    %     % Felnorm (max absolut fel i z2)
    %     errors(i) = abs(v_trap_end(2) - v_ref_end(2));
        
    %     fprintf('dt = %.2e s, Fel = %.4e\n', dt, errors(i));
    % end
    
    % slopes = diff(log(errors))./ diff(log(dt0 * alphas));
    % fprintf('Estimerad noggranhetssordning: %.2f\n', mean(slopes));
end

function [t, y] = trapezoidal_solver(A, p, v0, dt, T_end)
    t = 0:dt:T_end;
    N = length(t);
    n_states = length(v0);
    y = zeros(n_states, N);
    y(:,1) = v0;
    
    I = eye(n_states);
    M_lhs = I - 0.5 * dt * A;
    M_rhs = I + 0.5 * dt * A;
    
    for n = 1:N-1
        tn = t(n);
        tnp1 = t(n+1);
        
        g_n = evaluate_g(tn, p);
        g_np1 = evaluate_g(tnp1, p);
        
        rhs = M_rhs * y(:,n) + 0.5 * dt * (g_n + g_np1);
        
        y(:,n+1) = M_lhs \ rhs;
    end
    plot(t,y(2,:))
end

function g = evaluate_g(t, p)
    if t <= p.L/p.v
        arg = 2*pi*p.v*t/p.L;
        h = (p.H/2)*(1-cos(arg)); h_dot = (p.H*pi*p.v/p.L)*sin(arg);
    else
        h=0; h_dot=0;
    end
    g = [0;0;0; (p.k2*h + p.c2*h_dot)/p.m2];
end

function dv = system_dynamics_explicit(t, v, A, p)
    g = evaluate_g(t, p);
    dv = A*v + g;
end