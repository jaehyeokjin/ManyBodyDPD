intvalue=[];
for i=1:500
   int_comput = cumtrapz(realval(1:i));
   intvalue = [intvalue int_comput];
end
