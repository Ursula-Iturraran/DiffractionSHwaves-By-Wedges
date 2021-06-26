function [ X,Y ] = genMeshWedge(m,n,a,b,c,e,nu)
%genMeshWedge Generates an algebraic mesh ( lineal interpolation)
%   Ejemplo: [X,Y] = genMeshWedge(35,35,0,1,1,0,1/3)
%
%  Details
%               D
%       E               C
%       |               |
%       A---------------B
%
%   A = (a,e)
%   B = (b,e)
%   D = ( (a+b)/2, c)
%
%   ang(EDC) = nu*pi
%
%  m, n are the dimension of the resulting mesh m x n
%

alpha = (pi/2)*(1-nu);
f = (a+b)/2;
dist = abs(a-b)/2;
hipo = dist/cos(alpha);
x = hipo*sin(alpha);

d = c-x;

A = [a,e];
B = [b,e];
C = [b,d];
D = [f,c];
E = [a,d];

X = zeros(m,n);
Y = X;

% Firstly, we fill the edges 
m1 = floor(m/2);
m2 = m - m1;


t = linspace(0,1,m1);
for i=1:m1
    T = E + t(i)*(D-E);
    X(i,1)=T(1);
    Y(i,1)=T(2);
end

t = linspace(0,1,m2+1);
for i=1:m2+1
	T = D + t(i)*(C-D);
    X(i+(m1-1),1)=T(1);
    Y(i+(m1-1),1)=T(2);
end

t1 = linspace(0,0.5,m1);
t2 = linspace(0.5,1,m2+1);
t = [t1, t2(2:m2+1)];
for i=1:m
	T = A + t(i)*(B-A);
    X(i,n)=T(1);
    Y(i,n)=T(2);
end

% Fill the mesh interpolating 
t = linspace(0,1,n);
for i=1:m
	for j=2:n-1
		P1 = [X(i,1),Y(i,1)];
		P2 = [X(i,n),Y(i,n)];
		
		T = P1 + t(j)*(P2-P1);
		X(i,j) = T(1);
		Y(i,j) = T(2);
	end
end

