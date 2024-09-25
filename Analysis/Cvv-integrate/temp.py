import sys
import math
from scipy.integrate import simps

fname=sys.argv[1]
time_val = int(sys.argv[2])
outname="./integrate/"+fname
f=open(fname,'r')
out=open(outname,'w')

cvv=[]
t=[]
value=[]
if time_val == 100:
	rescale_factor = 5.0
else:
	rescale_factor = 1.0
print(rescale_factor)
i=0
for line in f:
	line_e=line.split()
	t.append(i*rescale_factor)
	cvv.append(float(line_e[1]))
	area = simps(cvv, dx=rescale_factor)/3.0
	value.append(area)
	i+=1
for i in range(0,len(t)):
	str_out = "%.6f %.6f\n" % (t[i],value[i])
	out.write(str_out)
f.close()
out.close()
