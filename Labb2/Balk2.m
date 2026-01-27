function Balk()
    A = 9; alpha = 0.5; omega = 4;
    
    y = @(t) A * exp(-alpha*t) .* cos(omega*t);
    y_prim = @(t) A * (-alpha * exp(-alpha*t) .* cos(omega*t) - omega * exp(-alpha*t) .* sin(omega*t));
    
    tol = 1e-8;
    
    %Uppgift 2a
    H1 = 0.5;
    f1 = @(t) y(t) - H1;
    df1 = y_prim;
    
    fprintf('Uppgift 2a\n');
    t = 1.5; 
    for k = 1:100
        val = f1(t);
        fprintf('Iter %d: t = %.10f, fel = %.2e\n', k, t, abs(val));
        if abs(val) < tol, break; end
        t = t - val / df1(t);
    end
    
    %Uppgift 2b
    fprintf('\nUppgift 2b \n');
    t0 = 1.4; t1 = 1.6; 
    for k = 1:100
        f0_val = f1(t0);
        f1_val = f1(t1);
        if abs(f1_val) < tol, break; end
        
        t_new = t1 - f1_val * (t1 - t0) / (f1_val - f0_val);
        t0 = t1;
        t1 = t_new;
        fprintf('Iter %d: t = %.10f, fel = %.2e\n', k, t1, abs(f1_val));
    end
    
    %Uppgift 2c
    H2 = 4.1355432362;
    f2 = @(t) y(t) - H2;
    fprintf('Uppgift 2c\n');
    t = 0.5; 
    for k = 1:100
        val = f2(t);
        deriv = df1(t); 
        fprintf('Iter %d: t = %.10f, derivata = %.4e\n', k, t, deriv);
        if abs(val) < tol, break; end
        t = t - val / deriv;
    end
end