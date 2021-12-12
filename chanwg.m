function chanwg

n0=1.5;
ns=1.5;
a=1.0;

n=1.55;
w=30;
d=30;


xmin=-50.0;
xmax=50.0;
Nx=1000;

ymin=-50.0;
ymax=50.0;
Ny=1000;

p=1.55;
err=0.00001;

filename='mode.dat';
% -----------------------

Nef1=0;
Nef=1;

Dx=(xmax-xmin)/(Nx-1);
x=(xmin:Dx:xmax);
Dy=(ymax-ymin)/(Ny-1);
y=(ymin:Dy:ymax);
Dx2=Dx^2;
g=(Dx/Dy)^2;

N=Nx*Ny;

p1=p^2*Dx2;

% index profile
%nx2=step_index(x,n,ns,nc,d);
%nx2=parabolic_index(x,1.5,1.0);
[Y,X]=meshgrid(y,x);

%n2=n0^2-a^2*(X.^2+Y.^2);
n2=rect_wg(X,Y,n,ns,w,d);
%surf(x,y,n2); shading interp; pause;

b=-2.0*(1.0+g)+Dx2*n2-p1;
%c=ones(Nx,Ny);
%d=ones(Nx,Ny)*g;

j=1:N;
M=sparse(j,j,b(:))+sparse(j(2:N),j(1:N-1),1,N,N)+sparse(j(1:N-1),j(2:N),1,N,N)+sparse(j(Nx+1:N),j(1:N-Nx),g,N,N)+sparse(j(1:N-Nx),j(Nx+1:N),g,N,N);
%disp(full(M));

psi=ones(Nx,Ny);
%psi=Y.*exp((-X.^2-Y.^2)/30^2);

cnt=0;
while(abs(Nef-Nef1)>=err)
    phi=M\psi(:);
    [m j]=max(abs(phi));
    r=phi(j);
    psi(:)=phi/r;
    ev=p1+1.0/r;
    
    Nef1=Nef;
    Nef=sqrt(ev/Dx^2);
    %p1=ev;

    disp(Nef);

    surf(x,y,psi); shading interp;
    cnt=cnt+1;
    %pause;
end

surf(x,y,psi); shading interp;

str=sprintf('Nef = %f\nno. iterations = %d\n',Nef,cnt);
disp(str);

fi=fopen(filename,'wt');
fprintf(fi,'%f %f\n',[x; psi';]);
fclose(fi);


function n2=rect_wg(X,Y,n,ns,w,d)
    n2=ones(size(X))*(ns^2);
    I1=find(abs(X)<w/2 & abs(Y)<d/2);
    n2(I1)=n^2;
  
    
    
function nx2=parabolic_index(x,n0,a)
    nx2=n0^2-a*x.^2;
    
    
    
        
