import sys
import math
from scipy.integrate import simps

fname=sys.argv[1]
time_val = int(sys.argv[2])
outname="./integrate/num_"+fname
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

def cvv_value(x):
	fit_val = 0.304*math.exp(-(0.1931*x)**0.8243)
	return (fit_val)
for line in f:
	line_e=line.split()
	t.append(i*rescale_factor)
	cvv.append(cvv_value(i*rescale_factor))
	area = simps(cvv, dx=rescale_factor)/3.0
	value.append(area)
	i+=1
for i in range(0,len(t)):
	str_out = "%.6f %.6f\n" % (t[i],value[i])
	out.write(str_out)
f.close()
out.close()
