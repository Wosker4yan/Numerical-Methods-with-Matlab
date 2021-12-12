function planarwg

n=1.55;
ns=1.5;
nc=1.0;
d=30.0;

xmin=-10.0;
xmax=100.0;
N=1000;

p=1.51;
err=0.00001;

filename='mode.dat';
% -----------------------

Nef1=0;
Nef=1;

D=(xmax-xmin)/(N-1);
x=(xmin:D:xmax);

p1=p^2*D^2;

% index profile
nx2=step_index(x,n,ns,nc,d);
%nx2=parabolic_index(x,1.5,1.0);

%plot(x,nx2); pause;

psi=ones(N,1);

cnt=0;
while(abs(Nef-Nef1)>=err)
    b=-2.0+D^2*nx2-p1;
    a=ones(N-1,1);
    A=diag(b)+diag(a,-1)+diag(a,1);

    phi=A\psi;
    [m j]=max(abs(phi));
    r=phi(j);
    psi=phi/r;
    ev=p1+1.0/r;
    
    Nef1=Nef;
    Nef=sqrt(ev/D^2);
    %p1=ev;

    disp(Nef);

    plot(x,psi);
    cnt=cnt+1;
    %pause;
end

str=sprintf('Nef = %f\nno. iterations = %d\n',Nef,cnt);
disp(str);

fi=fopen(filename,'wt');
fprintf(fi,'%f %f\n',[x; psi';]);
fclose(fi);



function nx2=step_index(x,n,ns,nc,d)
    nx2=x;
    I1=find(x<0);
    I2=find(x>=0 & x<=d);
    I3=find(x>d);
    nx2(I1)=nc^2;
    nx2(I2)=n^2;
    nx2(I3)=ns^2;
    
    
function nx2=parabolic_index(x,n0,a)
    nx2=n0^2-a*x.^2;
    
    
    
        
