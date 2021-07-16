clear; clc; close all;

fname = 'rptr.dat';
fid = fopen(fname,'r');
A = textscan(fid,'%d','HeaderLines',1);
fclose(fid);
rptr = int32(A{1});

fname = 'cptr.dat';
fid = fopen(fname,'r');
A = textscan(fid,'%d','HeaderLines',1);
fclose(fid);
cptr = int32(A{1});

nNo = size(rptr,1) - 1;
A = int32(zeros(nNo));
for i=1:nNo
    for k=rptr(i):rptr(i+1)-1
        j = cptr(k);
        A(i,j) = 1;
    end
end

spy(A);

return;