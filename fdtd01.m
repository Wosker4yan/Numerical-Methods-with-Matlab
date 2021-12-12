function fdtd01

xmin=-20.0;
xmax=20.0;
Nx=100;

ymin=-20.0;
ymax=20.0;
Ny=100;

nsteps=100;

w=1.0;

%filename='mode.dat';
% -----------------------

Dx=(xmax-xmin)/(Nx-1);
x=(xmin:Dx:xmax);
Dy=(ymax-ymin)/(Ny-1);
y=(ymin:Dy:ymax);

Hx=zeros(Ny,Nx);
Hy=zeros(Ny,Nx);
Ez=zeros(Ny,Nx);

eps=ones(Ny,Nx);

Dt=1.0/sqrt(1.0/Dx^2+1.0/Dy^2);

for n=1:nsteps
    Ez(Ny/2,Nx/2)=sin(w*n*Dt);
    
    Hx(2:Ny,:)=Hx(2:Ny,:)-(Dt/Dy)*(Ez(2:Ny,:)-Ez(1:Ny-1,:));
    Hy(:,2:Nx)=Hy(:,2:Nx)+(Dt/Dx)*(Ez(:,2:Nx)-Ez(:,1:Nx-1));
    Ez(1:Ny-1,1:Nx-1)=Ez(1:Ny-1,1:Nx-1) + ( (Dt/Dx)*(Hy(1:Ny-1,2:Nx)-Hy(1:Ny-1,1:Nx-1)) - (Dt/Dy)*(Hx(2:Ny,1:Nx-1)-Hx(1:Ny-1,1:Nx-1)) )./eps(1:Ny-1,1:Nx-1);
    
    surf(x,y,Hy), view(2);  shading interp;
    caxis([-.3,.3]);
    disp(n);
    pause(0.1);
end

