function S = switch_func_eps1(lambdaV_mag, lambdam, m, c)

    % Computes switching function for epsilon=1 case (minimum energy, L2
    % norm)
    
    S = -lambdaV_mag*c/m - lambdam;
    
end