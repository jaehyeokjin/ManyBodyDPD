import math
f=open("limit_coarse.out",'w')
#xx=[320+8*i for i in range(0,24)]
xx=[200+10*i for i in range(0,31)]

def y1(x):
	return (0.909206/(1.0+math.exp(-0.00877338*(x+46.6606))))
def y2(x):
	return (0.315139/(1.0+math.exp(-0.217467*(x-5.94482))))


for i in range(0,len(xx)):
	str_out = "%.6f %.6f\n" % (xx[i], y2(xx[i]))
	f.write(str_out)
f.close()

