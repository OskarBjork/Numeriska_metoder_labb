function TrapzMain()
    % Uppgift 4
    
    %Uppgift 4a
    x_vals = [0.11, 0.32, 1.14];
    N_vals = round(logspace(1.3, 2.7, 10));
    
    figure(5); clf;
    colors = {'b', 'r', 'y'};
    
    for i = 1:length(x_vals)
        x = x_vals(i);
        col = colors{i};
        errors = [];
        
        for N = N_vals
            integral_val = my_trapz(x, N);
            exact_val = erf(x); 
            errors = [errors, abs(integral_val - exact_val)];
        end
        
        loglog(N_vals, errors, ['-o' col], 'LineWidth', 1.5, ...
               'DisplayName', sprintf('x=%.2f', x)); hold on;
    end
    title('Fel i trapetsregeln vs N');
    xlabel('N'); ylabel('Absolut fel');
    legend('Location', 'best'); grid on;

    %Uppgift 4b
    N_fixed = [50, 120, 400];
    x_range = linspace(0, 6, 100);
    
    figure(6); clf;

    ylim([1e-16, 5e-5]); 

    for N = N_fixed
        errors_x = [];
        for x = x_range
            val = my_trapz(x, N);
            errors_x = [errors_x, abs(val - erf(x))];
        end
        semilogy(x_range, errors_x, '-', 'DisplayName', sprintf('N=%d', N)); hold on;
    end
    title('Fel i trapetsregeln vs x (Fix N)');
    xlabel('x'); ylabel('Absolut fel');
    legend; grid on;
    
end

function val = my_trapz(x, N)
    a = 0; b = x;
    h = (b-a)/N;
    t = linspace(a, b, N+1);
    
    g = @(v) (2/sqrt(pi)) * exp(-v.^2);
    y = g(t);
    
    T = sum(y) - 0.5*y(1) - 0.5*y(end);
    val = h * T;
end


% Felet minskar igen ty andra derivatan är liten för stora x