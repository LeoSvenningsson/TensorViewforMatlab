function [Tensors] = LoadTensor(fileName)

%Tensors = readmatrix(fileName,'Delimiter', ' ', 'ConsecutiveDelimitersRule', 'join') % old style of reading matrix  
fid = fopen(fileName,'r');

if fid ~= -1
    n = 0;
    tline = fgetl(fid);
    tline = fgetl(fid);
    while ischar(tline)
        tline = fgetl(fid);
        n = n+1;
    end
end
fclose(fid);
Tensors  = NaN(n,16);
 fileid = fopen(fileName,'r');
 if fileid ~= -1
     line = fgetl(fileid);
     line = fgetl(fileid);
         i=1;
         while ischar(line)
                 linelength = length(sscanf(line(1:end),'%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f')');
                 Tensors(i,1:linelength) =  sscanf(line(1:end),'%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f')';
                 i=i+1;
                 line = fgetl(fileid);
         end
end

end