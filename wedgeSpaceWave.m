function [ g22fEx ] = wedgeSpaceWave( nb, nu,  nt,  dt,  ts,  tp, bet, xa, za,   xi,   zi, nx,   nz)
%wedgeSpaceWave  Compute the Green's function in a mediumum with a wedge 
%   nb  = Number of terms to include in the Bessel series
%   nu  = Argument that defines the wedge
%   nt  = Number of time steps
%   dt  = Time step increment
%   ts  = Ricker parameter 
%   tp  = Ricker parameter (width)
%   bet = Velocity of S-wave
%
%   (xa,za) = Location of the source 
%   (xi,zi) = Initial x and z for the mesh grid
%   (nx,nz) = Number of grid points for x and z
%
%   (ni,zi)------------------(-ni,zi)
%      |                        |
%      |                        |
%       \                      /
%           \               /
%                \     /
%                 (0,0)
%
%   EXAMPLES
%
%        wedgeSpaceWave(40, 1/2, 256, 0.04, 0.5,  0.4,  1, 0, 1.0,  -2.0,  4, 101, 101)
%        wedgeSpaceWave(40, 1/2, 256, 0.04, 0.5,  0.4,  1, 0, 1.0,  -2.0,  4, 41,  41)
%        wedgeSpaceWave(40, 1/2, 256, 0.04, 0.5,  0.4,  1, 0, 1.0,  -1.5,  4, 7,   7)
%        wedgeSpaceWave(40, 1/2, 256, 0.04, 0.5,  0.4,  1, 0, 3.0,  -4,    4, 41,  41)
%        wedgeSpaceWave(40, 1/2, 256, 0.04, 0.2,  0.4,  1, 0, 3.0,  -4,    4, 17,  17)
%        wedgeSpaceWave(40, 1/2, 256, 0.04, 1.0,  0.4,  1, 0, 3.0,  -4,    4, 17,  17)
%        wedgeSpaceWave(40, 1/2, 256, 0.01, 0.25, 0.02, 2, 0, 3.0,  -4.5,  4, 11,  11)
%        wedgeSpaceWave(40, 1/2, 256, 0.04, 1.0,  0.08, 2, 0, 2.25, -4.5,  4, 101, 101)
%        wedgeSpaceWave(40, 1/2, 256, 0.04, 1.0,  0.08, 2, 0, 2.25, -4.5,  4, 11,  11)
%        wedgeSpaceWave(40, 1/2, 256, 0.02, 0.1,  0.15, 2, 0, 2.25, -4.5,  4, 101, 101)
%        wedgeSpaceWave(40, 1/2, 256, 0.03, 0.5,  0.2,  2, 0, 3,    -4.5,  4, 31,  315)
%        wedgeSpaceWave(40, 1/2, 256, 0.01, 0.2,  0.1,  2, 0, 3,    -4.5,  4, 31,  31)
%        wedgeSpaceWave(30, 1/2, 256, 0.02, 1.0,  0.3,  1, 0, 1,    -2,    4, 21,  21)
%        wedgeSpaceWave(30, 1/2, 256, 0.02, 2.0,  0.1,  2, 0, 0,    -2,    4, 11,  11)
%

%-- Initialization of variables
disp('Initialization of variables');
%-Constants

%complex zero
zer=complex(0,0);

%-Variables

%delta frequency
df=1/(nt*dt);
% shear modulus
rho=1;
mu=bet^2/rho;
ra=sqrt(xa^2+za^2);
phia = pi/2 - atan2(za,xa);
%Vectors

%frequencies
fr=0:df:(nt/2)*df;
%omega
om=2.0*pi*fr;
%complex k
ck=complex(om/bet,0.0);
%time
t=0:dt:(nt-1)*dt;
% Ricker wavelet
cr=complex(ricker1(t,tp,ts),0.0);
% Ricker in the frequency domain 
F_w = fft(cr,nt);

%Work space
g22fEx = zeros(nx,nz,nt/2+1,'double');

%Generating the mesh 
disp('Generating the mesh');
[X,Z] = genMeshWedge(nx,nz,xi,-xi,0,zi,-nu);

%Starting the computation
for j=1:nz
   for i=1:nx
       ri=sqrt(X(i,j)^2+Z(i,j)^2);
       phib = pi/2 - atan2(Z(i,j),X(i,j));             
       for f=1:nt/2+1 %number of frequencies    
          % inisde the wedge
          g22fEx(i,j,f)=(2/nu)*g22wedgeMAT(nu,phib,phia,ri,ra,ck(f),nb);
       end
   end
   disp([ num2str((i*j)/(nx*nz)*100,4) '% Computing']);
end

g22fEx=(1/mu)*g22fEx;
prodRick(:,:,nt/2)=F_w(nt/2)*g22fEx(:,:,nt/2);        

disp('Convolucionamos');
for j=1:nz
    for i=1:nx
        for k=1:nt/2-1
            prodRick(i,j,k)=F_w(k)*g22fEx(i,j,k);      %Convolution with a Ricker's
            prodRick(i,j,nt-k+1)=conj(prodRick(i,j,k+1));      %La Crepe..!
        end
    end
end
prodRick(:,:,1)=zer;
disp('Back to time domain');
for ii=1:nx
    for jj=1:nz
        g22tEx(ii,jj,1:nt)=ifft(prodRick(ii,jj,1:nt),nt);
    end    
end  

disp('Dibujamos');
mi=min(min(min(real(g22tEx))));
ma=max(max(max(real(g22tEx))));

%Animation
nombre = '';
while isempty(nombre)
    nombre = input('Name for the video? (i.e. ''wedge2'') ');
end
writerObj = VideoWriter(nombre);
open(writerObj);

for k=1:1:nt-1
    surf(X,Z,real(g22tEx(:,:,k)));
    axis([-abs(xi) abs(xi) 0 zi mi ma]);
    view(0,90);
    xlabel('x');
    ylabel('z');
    shading interp;
    F = getframe;
    writeVideo(writerObj,F);   
end 
close(writerObj);


end

