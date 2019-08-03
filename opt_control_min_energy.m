function u = opt_control_min_energy(chi,params)

    % chi is combined state and costate vector
    % c is exhaust velocity
    
    c = params.c; % nondim
    lambdaV = chi(11:13); % velocity costate
    lambdaV_mag = norm(lambdaV);
    m = chi(7); % already nondim
    lambdam = chi(14);

    S = -lambdaV_mag*c/m - lambdam + 1; % switching function for min energy case
    %fprintf("S = %d\n",S);

    if S < -1
        u = 1;
%     elseif abs(S) <= 1
%         u = -S/2;
    elseif S >= -1 && S <= 1
        u = (1-S)/2;
    elseif S > 1
        u = 0;
    else
        error("S is not real, S=%d\n",S)
    end

end