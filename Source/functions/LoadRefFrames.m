function [Refframes] = LoadRefFrames(fileName)
  
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
Refframes  = cell(n,9);
 fileid = fopen(fileName,'r');
 if fileid ~= -1
     line = fgetl(fileid);
     line = fgetl(fileid);
         i=1;
         while ischar(line)
                 linelength = length(sscanf(line(1:end),'%f %f %f %f %f %f %f %f')');
                 pat = ("AZYZ"|"PZYZ"|"AZXZ"|"PZXZ");
                 modeind = strfind(sscanf(line(1:end),'%c'),pat);
                 mode = sscanf(line(modeind:(modeind+3)),'%c');
                 Refframes(i,1:linelength) = num2cell(sscanf(line(1:end),'%f %f %f %f %f %f %f %f')');
                 Refframes(i,9) = {mode};
                 %Refframes{i,1:linelength} =  sscanf(line(1:end),'%f %f %f %f %f %f %f %f %4c')';
                 i=i+1;
                 line = fgetl(fileid);
         end
 end
end