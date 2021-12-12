function grinwg

n=1.8;
ns=1.5;
nc=1.0;
d=10.0;

xmin=-50.0;
xmax=50.0;
N=1000;

p=2.0;
err=0.00001;
% -----------------------
Nef1=0;
Nef=1;

D=(xmax-xmin)/(N-1);
x=(xmin:D:xmax);

p=p^2*D^2;

nx2=step_index(x,n,ns,nc,d);
psi=ones(N,1);
%plot(x,nx2);

while(abs(Nef-Nef1)>=err)
b=-2.0+D^2*nx2-p;
a=ones(N-1,1);
A=diag(b)+diag(a,-1)+diag(a,1);

phi=A\psi;
[m j]=max(abs(phi));
r=phi(j);
psi=phi/r;
Nef1=Nef;
Nef=sqrt((p+1.0/r)/D^2);
p=p+1.0/r;

disp(Nef);

plot(x,psi);
pause;
end
%disp(A);

function nx2=step_index(x,n,ns,nc,d)
    nx2=x;
    I1=find(x<0);
    I2=find(x>=0 & x<=d);
    I3=find(x>d);
    nx2(I1)=nc^2;
    nx2(I2)=n^2;
    nx2(I3)=ns^2;
