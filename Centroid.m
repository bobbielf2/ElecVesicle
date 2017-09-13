function val = Centroid(X,Y)
% find the centroid of a planar region enclosed by (X(t),Y(t))

if nargin == 0
    testCentroid;
    return
end
n= length(X);
if length(Y)~=n
    warning('X and Y must be of same dimension')
    return;
end
c0 = fftshift(1i*(-n/2:1:n/2-1)');
D1 = @(f) real(ifft(c0.*fft(f)));
Xa = D1(X); Ya = D1(Y);

A = pi/n*sum(X.*Ya - Y.*Xa); %area
Mx = -pi/n*sum(Y.^2.*Xa); %moment about x-axis
My = pi/n*sum(X.^2.*Ya); %moment about y-axis

val = [My, Mx]/A; %centroid

function testCentroid
n = 32;
t = linspace(0,2*pi,n+1).';
t = t(2:end);

% X = 2*cos(t)+4;
% Y = sin(t)+5;

% equilateral triangle (non smooth)

X = linspace(-1,1,n+1).';
X(end) = [];
X = [X;linspace(1,0,n+1).'];
X(end) = [];
X = [X;linspace(0,-1,n+1).'];
X(end) = [];

Y = linspace(0,0,n+1).';
Y(end) = [];
Y = [Y;linspace(0,sqrt(3),n+1).'];
Y(end) = [];
Y = [Y;linspace(sqrt(3),0,n+1).'];
Y(end) = [];

plot(X,Y)
axis equal
val = Centroid(X,Y);

disp('equilateral triangle test, error: '),disp(val(2)*3-sqrt(3))