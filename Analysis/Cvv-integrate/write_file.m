fp=fopen("real_coarse_dpd.out",'w');
for i=1:500
    fprintf(fp,'%.6f %.6f\n',xval(i),yval(i));
end
fclose(fp);