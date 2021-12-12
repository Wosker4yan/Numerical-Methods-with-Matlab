% 1D BMP algorithm based on Crank Nicolson scheme

dz = 1.0;
Nz = 200;

xmin=-50;
xmax=50;
N=1000;

% initial field
A=1.0;
w=10.0;

% waveguide
wg=5;
n=1.55;
ns=1.5;
nc=1.5;

% ---------

D=(xmax-xmin)/(N-1);
D2=D^2;
x=xmin:D:xmax;

z=0:dz:(Nz-1)*dz;
phi=zeros(Nz,N);

phi(1,:)=A*exp(-x.^2/w^2);
u=phi(1,:)'
plot(x,phi(1,:));

% refractive index

dn=ones(1,N)*ns;
I=find(abs(x)<wg/2);
I1=find(x<-w/2);
%dn(I)=n;
%dn(I1)=nc;
plot(x,dn);
%pause;

alpha=-i*dz*0.5;
a=(1.0 - alpha*(-2.0/D2 + 0))*ones(N,1);
b=-alpha*ones(N-1,1)/D2;

M=diag(a,0)+diag(b,1)+diag(b,-1);


for j=2:Nz
     u=M*u;
     M=2*eye(N)-M;
      
     if abs(u(3))>0
       t1=u(2)/u(3);
     else t1=0; end
         
     if abs(u(N-2))>0
        t2=u(N-1)/u(N-2);
     else t2=0; end
      
     if imag(t1)<0
        t1=1; end
     if imag(t2)<0
        t2=1; end      
    M(1,1)=1;
    M(1,2)=-t1;
    u(1)=0;
    M(N,N)=1;
    M(N,N-1)=-t2;
    u(N)=0;
      
    u=M\u;
    M=2*eye(N)-M;
    phi(j,:)=u;

    disp(j);
end



disp(size(x))
disp(size(z))
disp(size(phi))
surf(x,z,abs(phi));
shading interp;


