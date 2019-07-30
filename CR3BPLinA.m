function A = CR3BPLinA(X,mu)
    % A matrix for linearized CR3BP equations of motion
    % A = [0_3 I_3
    %      U_ii T], T = [0 2 0; -2 0 0; 0 0 0]
    % U_ii is hessian of potential wrt to x,y,z
    
    x = X(1);
    y = X(2);
    z = X(3);
    
    r13 = sqrt((x+mu)^2 + y^2 + z^2);
    r23 = sqrt((x-1+mu)^2 + y^2 + z^2);
    
    A = zeros(6,6);
    A(1:3,4:6) = eye(3);
    
    Uii = zeros(3,3);
    Uii(1,1) = 1 - (1-mu)/r13^3 - mu/r23^3 + 3*(1-mu)*(x+mu)^2/r13^5 + 3*mu*(x-1+mu)^2/r23^5;
    Uii(1,2) = 3*(1-mu)*(x+mu)*y/r13^5 + 3*mu*(x-1+mu)*y/r23^5;
    Uii(1,3) = 3*(1-mu)*(x+mu)*z/r13^5 + 3*mu*(x-1+mu)*z/r23^5;
    Uii(2,1) = Uii(1,2);
    Uii(2,2) = 1 - (1-mu)/r13^3 - mu/r23^3 + 3*(1-mu)*y^2/r13^5 + 3*mu*y^2/r23^5;
    Uii(2,3) = 3*(1-mu)*y*z/r13^5 + 3*mu*y*z/r23^5;
    Uii(3,1) = Uii(1,3);
    Uii(3,2) = Uii(2,3);
    Uii(3,3) = -(1-mu)/r13^3 - mu/r23^3 + 3*(1-mu)*z^2/r23^5 + 3*mu*z^2/r23^5;
    
    A(4:6,1:3) = Uii;
    
    A(4,5) = 2;
    A(5,4) = -2;
end