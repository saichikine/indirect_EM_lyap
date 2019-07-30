% Sara Powell, ASEN 5050, Project

% This function calculates A using the jacobian. It only needs to be run 
% once to write the A-matrix to memory. 
% initialize symbols used in the program.
syms x y z dotx doty dotz dotdotx dotdoty dotdotz mu r1 r2 

% DistAnce from the third body to the larger and smaller body respectively. Parker, 2014 page 36.
r1 = sqrt((x+mu)^2 + y^2 + z^2); 
r2 = sqrt((x-1+mu)^2 + y^2 + z^2); 
% Accelerations
dotdotx = 2*doty + x -(1-mu)*((x+mu)/(r1^3)) - mu*(x-1+mu)/(r2^3);
dotdoty = -2*dotx + y - (1-mu)*(y/(r1^3)) - mu*(y)/(r2^3);
dotdotz = -(1-mu)*((z)/(r1^3)) - mu*(z)/(r2^3);
x_state = [x,y,z,dotx,doty,dotz];

MatrixAinit = [dotx doty dotz dotdotx dotdoty dotdotz];

% Get the partial derivatives. 
A = jacobian(MatrixAinit, x_state);

fid = fopen('AfromSymbolic.m','w');
fprintf(fid,'function A = AfromSymbolic(x,y,z,dotx,doty,dotz, mu)\n');
for i=1:6
   for j=1:6
      fprintf(fid,'A(%i,%i) = %s;\n',i,j,char(A(i,j)));
   end
   fprintf(fid,'\n\n');
end

fclose(fid);