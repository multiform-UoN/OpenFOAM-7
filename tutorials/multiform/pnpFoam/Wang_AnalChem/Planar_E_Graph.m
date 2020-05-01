
%Open file to read
fid = fopen('Electric_Field.dat','r');
%Number of lines to read (excluding time and header lines)
numLines = 60;
%Dump first two lines
dump = fgetl(fid);
dump = fgetl(fid);
A = textscan(fid,'%f %f',numLines);
plot(A{2},A{1});
xlabel('x');
ylabel('grad(V)')