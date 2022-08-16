function Exportsurfglb(fig,filename,path)
fullfilename = fullfile(path,filename);
objcount = 1; % used to avoid non-surface objects
Nobjects = 0;
for Nchild = 1:length(fig.Children.Children)
    if fig.Children.Children(Nchild).Type == "surface"
        Nobjects = Nobjects + 1;
    end
end
collection = cell(1,Nobjects);
for j = 1:length(fig.Children.Children)
    if fig.Children.Children(j).Type == "surface"
        clear Model
        fvc = surf2patch(fig.Children.Children(j),'triangles');
        if isequal(fig.Children.Children(j).FaceColor,'flat')
        Model.COLOR_0 = fvc.facevertexcdata;
        else
        [Model.COLOR_0,~] = meshgrid([fig.Children.Children(j).FaceColor],1:length(fvc.vertices(:,1))); 
        end
        d = Model.COLOR_0; % precise sRGB convert
        s = d <= 0.04045;
        d(s) = d(s)/12.92;   
        d(~s) = ((d(~s) + 0.055)/1.055).^2.4;
        Model.COLOR_0 = d;
            %Model.COLOR_0 = (Model.COLOR_0).^(2.2); % aprox sRGB convert
        Model.POSITION = fvc.vertices;
        Model.indices = fvc.faces;        
        Model.NORMAL = reshape(fig.Children.Children(j).VertexNormals, [], 3);
        Model.prop.material.pbrMetallicRoughness.metallicFactor = 0;
        Model.prop.material.pbrMetallicRoughness.roughnessFactor = 0.5;
        Model.prop.material.alphaMode = 'OPAQUE'; %'BLEND''OPAQUE'
        Model.prop.material.doubleSided = true;
        collection{objcount} = Model; 
        objcount = objcount + 1;
    end
end

write_glb(fullfilename, collection)

end