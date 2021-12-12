function bpm2d

xmin=-50.0;
xmax=50.0;
Nx=200;

ymin=-50.0;
ymax=50.0;
Ny=100;

Dz=1.0;
amp=1.0;
wd=10;

%filename='mode.dat';
% -----------------------

Dx=(xmax-xmin)/(Nx-1);
x=(xmin:Dx:xmax);
Dy=(ymax-ymin)/(Ny-1);
y=(ymin:Dy:ymax);
Dx2=Dx^2;
Dy2=Dy^2;
alpha=i*Dz/2;
N=Nx*Ny;

[Y,X]=meshgrid(y,x);

% index profile
%n2=n0^2-a^2*(X.^2+Y.^2);
%n2=rect_wg(X,Y,n,ns,w,d);
%surf(x,y,n2); shading interp; pause;

psi=amp*exp((-X.^2-Y.^2)/wd^2);
%surf(x,y,psi); shading interp; pause;

b=1.0+alpha*(-2.0/Dx2+0)*ones(N,1);
c=(alpha/Dx^2)*ones(N,1);
c(Nx:Nx:N,1)=0.0;
Mx=sparse(1:N,1:N,b,N,N)+sparse(2:N,1:N-1,c(1:N-1),N,N)+sparse(1:N-1,2:N,c(1:N-1,1),N,N);

b=1.0-alpha*(-2.0/Dy2+0)*ones(N,1);
c=(-alpha/Dy^2)*ones(N,1);
c(Ny:Ny:N,1)=0.0;
My=sparse(1:N,1:N,b,N,N)+sparse(2:N,1:N-1,c(1:N-1),N,N)+sparse(1:N-1,2:N,c(1:N-1,1),N,N);

ID=sparse(1:N,1:N,1);
cnt=0;


while(1)

psi(:)=Mx*psi(:);
psi=psi';

for k=1:Ny:N
    if abs(psi(k+2))>0
       t1=psi(k+1)/psi(k+2);
    else t1=0; end
    if imag(t1)<0
        t1=1; end
    
    if abs(psi(k+Ny-3))>0
        t2=psi(k+Ny-2)/psi(k+Ny-3);
     else t2=0; end
    if imag(t2)<0
        t2=1; end
    
    My(k,k)=1.0;
    My(k,k+1)=-t1;
    psi(k)=0.0;
    My(k+Ny-1,k+Ny-1)=1.0;
    My(k+Ny-1,k+Ny-2)=-t2;
    psi(k+Ny-1)=0.0;
end
psi(:)=My\psi(:);

Mx=2*ID-Mx;
My=2*ID-My;

psi(:)=My*psi(:);
psi=psi';

for k=1:Nx:N
    if abs(psi(k+2))>0
       t1=psi(k+1)/psi(k+2);
    else t1=0; end
    if imag(t1)<0
        t1=1; end
    
    if abs(psi(k+Nx-3))>0
        t2=psi(k+Nx-2)/psi(k+Nx-3);
     else t2=0; end
    if imag(t2)<0
        t2=1; end
    
    Mx(k,k)=1.0;
    Mx(k,k+1)=-t1;
    psi(k)=0.0;
    Mx(k+Nx-1,k+Nx-1)=1.0;
    Mx(k+Nx-1,k+Nx-2)=-t2;
    psi(k+Nx-1)=0.0;
end

psi(:)=Mx\psi(:);

Mx=2*ID-Mx;
My=2*ID-My;

cnt=cnt+1











surf(y,x,abs(psi)), view(2); shading interp; 
pause;
end


function n2=rect_wg(X,Y,n,ns,w,d)
    n2=ones(size(X))*(ns^2);
    I1=find(abs(X)<w/2 & abs(Y)<d/2);
    n2(I1)=n^2;
  

    
    
        
