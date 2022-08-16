function SaveTensor(fileName,Tensor)
TensorNames = 'space delimited xx xy xz yx yy yz zx zy zz X(Å) Y(Å) Z(Å) Scale R G B';
fileID = fopen(fileName,'w');
fprintf(fileID,[TensorNames '\n']);%,Tensor
fprintf(fileID,'%.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f %.5f\n',Tensor);
fclose(fileID);
end