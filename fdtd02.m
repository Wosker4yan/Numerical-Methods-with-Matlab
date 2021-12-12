function fdtd01

xmin=-20.0;
xmax=20.0;
Nx=100;

ymin=-20.0;
ymax=20.0;
Ny=100;

nsteps=100;
w=1.0;

% PML
Nc=10;
m=3.5;


%filename='mode.dat';
% -----------------------

Dx=(xmax-xmin)/(Nx-1);
x=(xmin:Dx:xmax);
Dy=(ymax-ymin)/(Ny-1);
y=(ymin:Dy:ymax);

N1x=Nc;
N2x=Nx-Nc+1;
N1y=Nc;
N2y=Ny-Nc+1;

Hx=zeros(Ny,Nx);
Hy=zeros(Ny,Nx);
Ez=zeros(Ny,Nx);

eps=ones(Ny,Nx);

%eps(:,find(abs(x)<2.5))=1.5;

surf(x,y,eps),view(2); shading interp;
%pause

Dt=1.0/sqrt(1.0/Dx^2+1.0/Dy^2);

eps1=eps(Ny/2,N1x+1);
smax=0.8*(m+1)/(Dx*sqrt(eps1));
d=x(N1x)-x(1);
s1m=smax*((d-x(1:N1x)+x(1))/d).^m;
s1ms=s1m/eps1;
a1=1.0./(1+Dt*s1m/2);
b1=1-Dt*s1m/2;
p1=1.0./(1+Dt*s1m/(2*eps1));
q1=1-Dt*s1m/(2*eps1);

N3y=N2y-N1y-1;
E1zx=zeros(N3y,N1x);
E1zy=zeros(N3y,N1x);

%plot(x(1:N1x),[s1m;s1ms]);
%pause;

for n=1:nsteps
    Ez(Ny/2,Nx/2)=sin(w*n*Dt);
    
    Hx(N1y+1:N2y-1,N1x+1:N2x-1)=Hx(N1y+1:N2y-1,N1x+1:N2x-1)-(Dt/Dy)*(Ez(N1y+1:N2y-1,N1x+1:N2x-1)-Ez(N1y:N2y-2,N1x+1:N2x-1));
    Hy(N1y+1:N2y-1,N1x+1:N2x-1)=Hy(N1y+1:N2y-1,N1x+1:N2x-1)+(Dt/Dx)*(Ez(N1y+1:N2y-1,N1x+1:N2x-1)-Ez(N1y+1:N2y-1,N1x:N2x-2)); 
    Ez(N1y+1:N2y-1,N1x+1:N2x-1)=Ez(N1y+1:N2y-1,N1x+1:N2x-1)+ ( (Dt/Dx)*(Hy(N1y+1:N2y-1,N1x+2:N2x)-Hy(N1y+1:N2y-1,N1x+1:N2x-1)) - (Dt/Dy)*(Hx(N1y+2:N2y,N1x+1:N2x-1)-Hx(N1y+1:N2y-1,N1x+1:N2x-1)) )./eps(N1y+1:N2y-1,N1x+1:N2x-1);
    
    % left boundary 
    Hx(N1y+1:N2y-1,1:N1x) = Hx(N1y+1:N2y-1,1:N1x) -(Dt/Dy)*(Ez(N1y+1:N2y-1,1:N1x)-Ez(N1y:N2y-2,1:N1x));
    Hy(N1y+1:N2y-1,2:N1x) = a1(2:N1x).*( b1(2:N1x).*Hy(N1y+1:N2y-1,2:N1x) + (Dt/Dx)*( Ez(N1y+1:N2y-1,2:N1x)- Ez(N1y+1:N2y-1,1:N1x-1))  );
    
    size(p1)
    size(q1)
    size(E1zx(1:N3y,1:N1x))
    
    E1zx(1:N3y,1:N1x)= p1(:)'.*( q1(:)'.*E1zx(1:N3y,1:N1x) + (Dt/(Dx*eps1))*(Hy(N1y+1:N2y-1,2:N1x+1)-Hy(N1y+1:N2y-1,1:N1x) ) );
    E1zy(1:N3y,1:N1x)= E1zy(1:N3y,1:N1x) - (Dt/(Dy*eps1))*(Hx(N1y+2:N2y,1:N1x)-Hx(N1y+1:N2y-1,1:N1x) );
    Ez(N1y+1:N2y-1,1:N1x)=E1zx(1:N3y,1:N1x)+E1zy(1:N3y,1:N1x);
    
    
    surf(x,y,Ez), view(2);  shading interp;
    caxis([-.3,.3]);
    disp(n);
    pause(0.1);
end

