function g = gfunc(R,mu)

    x = R(1);
    y = R(2);
    z = R(3);
    
    r1 = sqrt((x+mu)^2 + y^2 + z^2);
    r2 = sqrt((x+mu-1)^2 + y^2 + z^2);

    g = [x - (1-mu)*(x+mu)/r1^3 - mu*(x+mu-1)/r2^3;
        y - (1-mu)*y/r1^3 - mu*y/r2^3;
        -(1-mu)*z/r1^3 - mu*z/r2^3];
    
end