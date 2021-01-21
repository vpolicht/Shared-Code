function [ output ] = read_bin(bin_file)
fid=fopen(bin_file,'r','ieee-le');
dim1=fread(fid,1,'*uint32');
dim2=fread(fid,1,'*uint32');
output=fread(fid,[dim2 dim1],'float32=>float64')';
fclose(fid);



end

