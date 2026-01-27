function Road()
    % --- Uppgift 3a: Hitta koordinater ---
    % Data för P1, P2, P3 från tabell
    % Format: [xA, yA, xB, yB, LA, LB]
    data = [
        170, 950, 160, 1008, 60, 45;   % P1
        420, 2400, 370, 2500, 75, 88;  % P2
        670, 1730, 640, 1760, 42, 57   % P3
    ];
    
    P_found = zeros(3, 2);
    
    fprintf('--- Uppgift 3a: Newtons metod för system ---\n');
    for i = 1:3
        xA = data(i,1); yA = data(i,2);
        xB = data(i,3); yB = data(i,4);
        LA = data(i,5); LB = data(i,6);
        
        % Startgissning: Medelvärdet av A och B (eller grafisk gissning)
        x = [(xA+xB)/2 + 50; (yA+yB)/2]; 
        
        for k = 1:10
            % Funktion F
            F = [ (x(1)-xA)^2 + (x(2)-yA)^2 - LA^2;
                  (x(1)-xB)^2 + (x(2)-yB)^2 - LB^2 ];
              
            % Jacobian J
            J = [ 2*(x(1)-xA), 2*(x(2)-yA);
                  2*(x(1)-xB), 2*(x(2)-yB) ];
              
            h = -J \ F;
            x = x + h;
            if norm(h) < 1e-6, break; end
        end
        P_found(i,:) = x';
        fprintf('Punkt P%d: (%.2f, %.2f)\n', i, x(1), x(2));
    end
    
    % --- Uppgift 3b: Polynominterpolation ---
    P0 = [0, 0]; P4 = [1020, 0];
    All_P = [P0; P_found; P4]; % 5 punkter
    
    x_nodes = All_P(:,1);
    y_nodes = All_P(:,2);
    
    % Anpassa polynom grad 4 (5 punkter)
    coeffs = polyfit(x_nodes, y_nodes, 4);
    
    xx = linspace(0, 1020, 200);
    yy = polyval(coeffs, xx);
    
    figure(3);
    plot(x_nodes, y_nodes, 'o', xx, yy, '-');
    title('3b: Polynominterpolation av väg');
    axis equal;
    
    % --- Uppgift 3c & 3d: roadcoord filer ---
    % OBS: Detta är pseudokod då filerna saknas. 
    % Ladda filerna: load roadcoord.mat; load roadcoord2.mat;
    load roadcoord.mat; 
    coeffs_road = polyfit(x, y, length(x)-1);
    xx_road = linspace(min(x), max(x), 300);
    yy_road = polyval(coeffs_road, xx_road);
    figure(4);
    plot(x, y, 'o', xx_road, yy_road, '-');
    %hold on


    title('3c: Polynominterpolation av väg från roadcoord.mat');
    axis equal;

    xx_interp = linspace(min(x), max(x), 1000);
    v = interp1(x, y, xx_interp);
    figure(5);
    plot(x, y, 'o', xx_interp, v, '-');
    title('3c: Linjär interpolation av väg från roadcoord.mat');

    load roadcoord2.mat;
    coeffs_road2 = polyfit(x2, y2, length(x2)-1);
    xx_road2 = linspace(min(x2), max(x2), 300);
    yy_road2 = polyval(coeffs_road2, xx_road2);
    figure(6);
    plot(x2, y2, 'o', xx_road2, yy_road2, '-');
    title('3d: Polynominterpolation av väg från roadcoord2.mat');
    %axis equal;

    % Grafen för roadcoord2.mat ser bätre ut eftersom punkterna ej har samma avstånd till omgivande punkter, 
    % vilket skulle ge upphov till Runges fenomen.


    % t = 1:length(x); % Parameter (index)
    % tt = linspace(1, length(x), 1000);
    % xx_int = spline(t, x, tt);
    % yy_int = spline(t, y, tt);
    % figure(4);
    % plot(x, y, 'o', xx_int, yy_int);
    % title('3c: Spline-interpolation av väg från roadcoord.mat');

    
    % Generell princip för parameter-interpolation (3d):
    % Antag att du har vektorer x2 och y2 från roadcoord2.mat
    % t = 1:length(x2); % Parameter (index)
    % tt = linspace(1, length(x2), 1000);
    % xx_int = spline(t, x2, tt);
    % yy_int = spline(t, y2, tt);
    % figure(4); plot(x2, y2, 'o', xx_int, yy_int);
end
