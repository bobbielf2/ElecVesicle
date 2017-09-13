% Parameter values
mu=1;
shear = 0.0;
T=1000;
dt=1e-2;
kb=0.05;
N = 128;

% Creat initial vesicle positions
t=linspace(0,2*pi,N+1)';

M = 10;
for k=1:M
    x{k}=3.8*cos(2*pi*k/M+pi/2/M)+.5*cos(t);
    y{k}=4.8*sin(2*pi*k/M+pi/2/M)+.7*sin(t);
end

% Electrical parameters
% conduct_rat = 0.8;
beta = 40;
E_infty = 0+sqrt(beta/5)*1i;
tol = 1e-13;
Gm = 5;

% Name of output file
Name ='Test';