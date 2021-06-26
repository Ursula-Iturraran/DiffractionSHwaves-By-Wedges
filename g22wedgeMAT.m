function g22w=g22wedgeMAT(nu,phi,phi0,r,r0,k,n)
%with vectorization makes computation faster
rmin=min(r,r0);
rmax=max(r,r0);
ui=complex(0.0,1.0);
j=0:n-1;
en=2*ones(1,n);
en(1)=1;
%normailized
%g22=en.*cos((j/nu)*(phi+nu*(pi/2))).*cos((j/nu)*(phi0+nu*(pi/2))).* besselj(j/nu,k*rmin).*besselh(j/nu,2,k*rmax)/besselh(0,2,k*r0);
% not normailized
g22=en.*cos((j/nu)*(phi+nu*(pi/2))).*cos((j/nu)*(phi0+nu*(pi/2))).* besselj(j/nu,k*rmin).*besselh(j/nu,2,k*rmax);
%nu, phi, pi, phi0, k, rmin, rmax, r0 needs to be a scalar.
% normalized
g22w=sum(g22(1:1:n));
% not normalized
%g22w=(ui/4)*sum(g22(1:1:n));

