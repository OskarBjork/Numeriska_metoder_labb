function main()
    k1 = 5400
    k2 = 135000
    ck = 0.71
    cs = 1
    x = [k1 k2]'; % Startgissning
    tol = 1e-6; 
    F = transfer_functions(x(1), x(2), ck, cs);
    iter = 0; maxiter = 10000;
    diffx = norm(F);
    while diffx > tol && iter < maxiter
        
        iter = iter + 1;
        J = Jacobian_transfer_functions(x(1), x(2));
        inv_J = inv(J);
        xnew = x - mtimes(inv_J, F) 
        x = xnew;
        F = transfer_functions(x(1), x(2), ck, cs);
        diffx = norm(F); 
    disp([iter diffx])
    disp("Hej")
    disp(xnew)
    end
    if iter >= maxiter
        disp("Varning Maximalt antal iterationer natt")
    end
end
