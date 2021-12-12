% Modes of a step-index planar waveguide 

function planarwg

n=1.55;
ns=1.5;
nc=1.0;
d=5.0e-6;
lambda=1.5e-6;
nmodes=5;
nu=0;

xmin=-5e-6;
xmax=15e-6;
N=100;

% -------------------

D=(xmax-xmin)/(N-1);
x=(xmin:D:xmax)';

k0=2*pi/lambda;
V=k0*d*sqrt(n^2-ns^2);
a=(n^2-nc^2)/(n^2-ns^2);

nmodes=ceil(dispersion(0,V,a,0.0)/pi);
display(nmodes);
hold off;

%for nu = 0:nmodes-1
    b=fzero(@(b)dispersion(b,V,a,nu),0.5);
    nef2=ns^2+b*(n^2-ns^2);
    bt2=k0^2*nef2;
    nef=sqrt(nef2);
    bt=sqrt(bt2);

    disp(nef);

    al=k0*sqrt(n^2-nef2);
    gm=k0*sqrt(nef2-ns^2);
    dl=k0*sqrt(nef2-nc^2);

    E = electric_component(x,d,al,gm,dl);
    plot(x,E);
    hold on;
%end
   

function f=dispersion(b,V,a,nu)
   f=V*sqrt(1-b)-atan(sqrt((b+a)/(1-b)))-atan(sqrt(b/(1-b)))-nu*pi;
   
function E = electric_component(x,d,al,gm,dl)
   E=x;
   I1=find(x<0);
   I2=find(x>=0 & x<=d);
   I3=find(x>d);
   E(I1)=exp(dl*x(I1));
   E(I2)=cos(al*x(I2))+dl*sin(al*x(I2))/al;
   E(I3)=(cos(al*d)+dl*sin(al*d)/al)*exp(-gm*(x(I3)-d));
   