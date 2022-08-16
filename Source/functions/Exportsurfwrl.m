function Exportsurfwrl(fig,filename,path)
fullfilename = fullfile(path,filename);

if contains(filename(end-4+1:end) ,'.wrl')
fileID = fopen(fullfilename,'w');
else
   fileID = fopen([fullfilename '.wrl'],'w'); 
end

for j = 1:length(fig.Children.Children)
    if fig.Children.Children(j).Type == "surface"
        text3 = '    Transform {\n      translation 0 0 0\n      rotation 0 0 1 0\n      scale 1 1 1\n      children [\n        Shape {\n          appearance Appearance {\n            material Material {\n              ambientIntensity 0\n              diffuseColor 1 1 1\n              specularColor 0 0 0\n              shininess 0.0078125\n              transparency %1.4f\n              }\n            }\n          geometry IndexedFaceSet {\n            solid FALSE\n            coord DEF VTKcoordinates Coordinate {\n              point [\n';
        
        fprintf(fileID,text3,1-fig.Children.Children(j).AlphaData(1));
        
        
        
        
        
        surf = [fig.Children.Children(j).XData(:) fig.Children.Children(j).YData(:) fig.Children.Children(j).ZData(:)];
        formatSpec = '              %4.7f %4.7f %4.7f,\n';
        fprintf(fileID,formatSpec,surf');
        
        text4 = '              ]\n            }\n            normal DEF VTKnormals Normal {\n              vector [\n';
        fprintf(fileID,text4);
        
        normsX = fig.Children.Children(j).VertexNormals(:,:,1);
        normsX = normsX(:);
        normsY = fig.Children.Children(j).VertexNormals(:,:,2);
        normsY = normsY(:);
        normsZ = fig.Children.Children(j).VertexNormals(:,:,3);
        normsZ = normsZ(:);
        norms = [normsX normsY normsZ];
        formatSpec = '           %1.6f %1.6f %1.6f,\n';
        fprintf(fileID,formatSpec,norms');
        
        text5 = '            ]\n          }\n            color DEF VTKcolors Color {\n              color [\n';
        fprintf(fileID,text5);
        if isequal(fig.Children.Children(j).FaceColor,'flat')
        color = reshape(fig.Children.Children(j).CData,length(normsX),3);
        else
        [color,~]=meshgrid([fig.Children.Children(j).FaceColor],1:length(normsX)); 
        end
        formatSpec = '           %1.6f %1.6f %1.6f,\n';
        fprintf(fileID,formatSpec,color');
        
        text6 = '            ]\n          }\n            coordIndex  [\n';
        fprintf(fileID,text6);
        
        
        l1 = length(fig.Children.Children(j).XData(:,1));
        l2 = length(fig.Children.Children(j).XData(2,:));
        Conmat=zeros((l1),(l2));
        Con1=zeros((l1-1),(l2-1));
        Con2=zeros((l1-1),(l2-1));
        Con3=zeros((l1-1),(l2-1));
        Con4=zeros((l1-1),(l2-1));
        for p = 0:(l1-1)
            for q = 0:(l2-1)
                Conmat(p+1,q+1) = p+q*l1;
            end
        end
        
        for p = 1:(l1-1)
            for q = 1:(l2-1)
                Con1(p,q) = Conmat(p,q);
                Con2(p,q) = Conmat(p+1,q);
                Con3(p,q) = Conmat(p+1,q+1);
                Con4(p,q) = Conmat(p,q+1);
            end
        end
        Conflat1 = Con1(:);
        Conflat2 = Con2(:);
        Conflat3 = Con3(:);
        Conflat4 = Con4(:);
        
        Conflat = [Conflat1 Conflat2 Conflat3 Conflat4];
        
        formatSpec = '              %1d, %1d, %1d, %1d, -1,\n';
        fprintf(fileID,formatSpec,Conflat');
        text7 = '            ]\n          }\n        }\n      ]\n    }\n';
        fprintf(fileID,text7);
    end
end

fclose(fileID);

end