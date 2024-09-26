load cvpv_sparse.dat
load cfv_sparse.dat
v=cvpv_sparse;
y=cfv_sparse;
y(1)=0;
fm=zeros(100,100);
dt=0.025;
for k=2:1:100
fm(k,1)=0.5*dt*v(k);
fm(k,k)=0.5*dt*v(1);
for i=2:1:k-1
fm(k,i)=dt*v(k+1-i);
end
end
size(y)
nx=fm'*fm;
nf=fm'*y;
reg=zeros(100,100);
for i=1:1:100
reg(i,i)=0.00001;%0.0001;
end

rr=-pinv(nx+reg)*nf;

