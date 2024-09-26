function value = frictions(a)
load friction.dat
pp=friction;
f = zeros(101,1);
for s=0:1:100
    f(s+1)=-pp(s+1);
    for k=0:1:(100-s)
        f(s+1) = f(s+1) + a(k+1)*a(k+1+s);
    end
end
value = f(3:end)'*f(3:end)+10*f(1:2)'*f(1:2);
