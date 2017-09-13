% Parameter values
mu=1;
shear = 0;
T=2000;
dt=1e-2;
a = 5/2/pi;
% chi = 7.9e-1;
% kb=chi*mu*a^2;
kb = 5e-2;
N = 64;

% Creat initial vesicle positions
t=linspace(0,2*pi,N+1)';
% x{1}=pi+2*0.1*cos(t)/2-2;
% y{1}=2*0.1*(sin(t)/2+cos(2*t)/4);
% for k=2:2
%     x{k}=x{k-1}+0.21;
%     y{k}=y{1};
% end
M = 1;
for k=1:M
% 2-1 vesicle, param = 5;
    x{k} = 0.51607852*cos(t)-1.1*(k-1);
    y{k} = 1.03215704*sin(t)+1.3*(k-1);

% Reduce area = 0.99, perimeter = 5
%     x{k}=0.7295*cos(t);
%     y{k}=0.8595*sin(t);

% Delta = .95
% M = 1;
% for k=1:M
%     x{k}=0.6440*cos(t);
%     y{k}=0.9343*sin(t);
% end

% Delta = .9
%     x{k}=0.5775*cos(t); 
%     y{k}=0.9872*sin(t);

% Delta = .8
%     x{k}=0.479*cos(t); 
%     y{k}=1.0576*sin(t);

% Delta = .7
%     x{k}=0.4004*cos(t);
%     y{k}=1.1070*sin(t);
    
% Delta = .6    
%     x{k}=0.332*cos(t);
%     y{k}=1.145*sin(t);

% multiple
%     x{k}=2*cos(2*pi*k/M+pi/2/M)+.5*cos(t);
%     y{k}=2*sin(2*pi*k/M+pi/2/M)+.7*sin(t);
end

% Electrical parameters
% beta = 80; %field strength compared to shear
% E_infty = 0+sqrt(beta/a*mu)*1i;
Mn = 100; %field strength compared to 
E_infty = 0+sqrt(Mn*kb/a^3)*1i;
G = 0;
Gm = G/a;
conduct_rat = .2;

tol = 1e-13;

% Name of output file
Name='Test';