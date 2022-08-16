function TensorViewApp(Fname,Tensor,x,y,z,Scale,Trans,R,G,B,OvEl,MolTen,Bond,MoreT)
%%
% TensorView for Matlab is a GUI and script based tool to visualize tensors in a molecular context. 
% TensorView for Matlab includes: Reading arbitrary .pdb and .xyz files for molecular visualization. 
% Ovaloid and ellipsoid tensor visualisation. 3D model file conversion to .glb and .wrl. Euler angle and relative angle decoder.    
%
% TensorView for Matlab is licenced by creative commons CC BY. https://creativecommons.org/licenses/
% Free to share and adapt. Give appropriate credits to authors:
% ------ Credits
% github.com/LeoSvenningsson/TensorViewforMatlab
% 
% mathworks.com/matlabcentral/fileexchange/55231-molecoule3d
% 
% onlinelibrary.wiley.com/doi/full/10.1002/mrc.4793
% 
% https://se.mathworks.com/matlabcentral/fileexchange/109264-matlab2glb
% 
%    Version: 1.15
% 
% Authors: 
% Dr. Leo Svenningsson (leo.svenningsson@fkem1.lu.se) 
% Prof. Leonard Mueller (leonard.mueller@ucr.edu)
% Contributors:
% Dr. André Ludwig (aludwig@alumni.ethz.ch)
% Dmitri Sastin
% ------ 
% TensorView for Matlab is a collaboration with works derrived from molecule3D.m (André Ludwig: mathworks.com/matlabcentral/fileexchange/55231-molecule3d),
% matlab2glb ( Dmitri Sastin: https://se.mathworks.com/matlabcentral/fileexchange/109264-matlab2glb),
% 
%%% User input starts here
FileName = Fname; % the molecular file in either .pdb or .xyz format

atomcoord = [x, y, z]; % Atom coordinate for the CSA tensor

OvaloidEllipsoid = OvEl; % "ovaloid" for ovaloid;  "ellipsoid" for elipsoid tensor %

TensorScale = Scale; % Tensor scaling

PlotTensor = MolTen; % 0 Plots the tensor with the molecule; 1 plots only the molecule; 2 only the tensor

styles = 'ballstick'; % 'ballstick','licorice','large','superlarge'

Transparency = Trans; % Sets the transperancy of the tensor

Bondlimit = Bond; % 1.6Å This is only used for the .xyz and connectionless .pdb format

Brightnes = 70; % Increases brightness for negative component
BrR = min([255 (R*255+Brightnes)])/255;
BrG = min([255 (G*255+Brightnes)])/255;
BrB = min([255 (B*255+Brightnes)])/255;

%%% User input ends here
fileid = fopen(FileName,'r');
if fileid ~= -1
    line = fgetl(fileid);
    
    if contains(FileName,".pdb") % reads .pdb files
        i=1;
        j=1;
        Conlist = 0;
        while ischar(line)
            if(strncmp('HETATM',line,6) || strncmp('ATOM',line,4))
                xyz(i,:) =  sscanf(line(31:54),'%f %f %f')';
                labels(i) = {sscanf(line(76:78),'%s')};
                i=i+1;
            end
            if(strncmp('CONECT',line,6))
                while j ~= sscanf(line(9:11),'%f')
                    Conlist(j,1:5) = [j,0,0,0,0];
                    j = j+1;
                end
                Conlist(j,1:length(sscanf(line(9:end),'%f')')) =  sscanf(line(9:end),'%f')';
                j = j+1;
            end
            line = fgetl(fileid);
        end
        fclose(fileid);
        filetype = ".pdb";
    elseif contains(FileName,".xyz") % reads .xyz files
        NrOfAtoms=str2double(line);
        line = fgetl(fileid); % Skip a line
        for i = 1:NrOfAtoms
            line = fgetl(fileid);
            xyz(i,:) =  sscanf(line(10:48),'%f %f %f')';
            labels(i) = {sscanf(line(1:3),'%s')};
        end
        fclose(fileid);
        filetype = ".xyz";
        Conlist = 0;
    else
        error('fileformat .pdb or .xyz was not found')
    end
    
end

[X,Y,Z,AlphaMap]=CreateTensor(Tensor,atomcoord,OvaloidEllipsoid,Transparency,TensorScale);

ColorMap = zeros([size(AlphaMap),3]);
Alphaind = AlphaMap<Transparency; % inderectly finds negative tensor blob
ColorMapR = BrR*Alphaind;
ColorMapR(ColorMapR==0) = R;
ColorMapG = BrG*Alphaind;
ColorMapG(ColorMapG==0) = G;
ColorMapB = BrB*Alphaind;
ColorMapB(ColorMapB==0) = B;
ColorMap(:,:,1) = ColorMapR;
ColorMap(:,:,2) = ColorMapG;
ColorMap(:,:,3) = ColorMapB;

if fileid ~= -1
    %%%% Adjust the camera zoom to the size of the molecule
    Xbox=max(xyz(:,1))-min(xyz(:,1));
    Ybox=max(xyz(:,2))-min(xyz(:,2));
    Zbox=max(xyz(:,3))-min(xyz(:,3));
    zoomFactor = min([Xbox Ybox Zbox])/max([Xbox Ybox Zbox]);
    %%%%
end


if PlotTensor==0 % plots the tensor with the molecule
    fig = figure;
    surface(X,Y,Z,ColorMap,'EdgeColor','none','FaceLighting','gouraud','AlphaData',AlphaMap,'AlphaDataMapping','none','FaceAlpha','flat','AmbientStrength',0.7);
    hold on
    set(gcf,'Color','w')
    molecule3D(xyz,labels,styles,filetype,Conlist,Bondlimit) % plots the molecule
    axis equal
    %cameratoolbar % Use the camera control, not the "zoom options/magnifying glass"
    %camzoom(zoomFactor)
    light('Position',[-2 -2 -4]);
    ax = gca;
    ax.Interactions = [rotateInteraction];
elseif PlotTensor==1
    fig = figure;
    set(gcf,'Color','w')
    molecule3D(xyz,labels,styles,filetype,Conlist,Bondlimit) % plots the molecule
    axis equal
    %cameratoolbar % Use the camera control, not the "zoom options/magnifying glass"
    %camzoom(zoomFactor)
    light('Position',[-2 -2 -4]);
    ax = gca;
    ax.Interactions = [rotateInteraction];
elseif PlotTensor==2
    if MoreT == false
        fig = figure;
    end
    g = groot;
    if isempty(g.Children)
        light('Position',[2 2 4]);
        light('Position',[-2 -2 -4]);
    end
    if MoreT == false
        light('Position',[2 2 4]);
        light('Position',[-2 -2 -4]);
    end
    set(gcf,'Color','w')
    hold on
    surf(X,Y,Z,ColorMap,'EdgeColor','none','FaceLighting','gouraud','AlphaData',AlphaMap,'AlphaDataMapping','none','FaceAlpha','flat','AmbientStrength',0.7); % 'FaceAlpha',Transparency
    axis equal
    set(gca,'visible','off')
    %cameratoolbar
    camtarget('auto')
    %material shiny
    ax = gca;
    ax.Interactions = [rotateInteraction];
    hold on
else
    error("PlotTensor must be either 0 or 1 or 2")
end

% Active rotation matrices for reference
%
%   rx =    [1, 0, 0;
%            0, Cos(rad), -Sin(rad);
%            0, Sin(rad), Cos(rad)];
%
%   ry =    [Cos(rad), 0, Sin(rad);
%             0, 1, 0;
%             -Sin(rad), 0, Cos(rad)];
%
%   rz =   [Cos(rad), -Sin(rad), 0;
%             Sin(rad), Cos(rad), 0;
%             0, 0, 1];
end