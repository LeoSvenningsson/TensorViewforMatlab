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
% and TensorView (Prof. Leonard Mueller: https://onlinelibrary.wiley.com/doi/full/10.1002/mrc.4793 and https://sites.google.com/ucr.edu/Mueller/home/tensorview)
%%
clear all

funfolder = fullfile(pwd,'functions'); % adds the path to the functions folder
addpath(funfolder);


%%% User input starts here
FileName ='ALA-ALA-ALA_NMR.pdb'; % the molecular file in either .pdb or .xyz and connectionless .pdb format

Tensor =  [-5.9766, -60.3020, -10.8928;
    -65.5206, -23.0881, -25.2372;
    -9.5073, -28.2399, 56.2779]; % The molecular frame tensor in the either non-symmetric or symmetric form. The script will transform to the symetrical version regardles of input.

%Tensor = RotateTensor(Alpha,Beta,Gamma,PAS,RotationMode); % RotateTensor(Alpha,Beta,Gamma,PAS,RotationMode) uncomment this line if you have the PAS CSA and the euler angles; 0 = ZYZactive; 1 = ZYZpassive; 2 = ZXZactive; 3 = ZXZpassive; PAS = [Sxx 0 0; 0 Syy 0; 0 0 Szz]

atomcoord = [-4.299848 -0.400088 -0.033384]; % Atom coordinate for the CSA tensor

OvaloidEllipsoid = "ovaloid"; % "ovaloid" for ovaloid;  "ellipsoid" for elipsoid tensor %

TensorScale = 1.4; % Tensor scaling

R = 1;
G = 0.5;
B = 0;
Color = [R,G,B]; % color of the CSA tensor

PlotTensor = 0; % 0 Plots the tensor with the molecule; 1 plots only the molecule; 2 only the tensor

styles = 'ballstick'; % 'ballstick','licorice','large','superlarge'

Transparency = 0.7; % Sets the transperancy of the tensor

Bondlimit = 1.6; % 1.6Å This is only used for the .xyz format and connectionless .pdb format

Brightnes = 70; % Increases brightness fornegative compinent
BrR = min([255 (R*255+Brightnes)])/255;
BrG = min([255 (G*255+Brightnes)])/255;
BrB = min([255 (B*255+Brightnes)])/255;

%%%% Add extra tensors
NR = 0; % Number of extra tensors
% T1 = [ 50.675  61.810 -29.154;
%            53.606 -50.064  30.894;
%           -39.382  60.395  4.832  ]; % Tensor naming convention T1, T2, T...
%
% A1 = [-4.745, 10.679, 6.273];% atomcoord naming convention A1, A2, A...
% C = {T1;A1}; % cell structure {T1,T2,T3;A1,A2,A3}
%%%%

% %%% dipolar tensors
% NRDP = 1;% Number of extra dipolar tensors
% DPcoorda1 = [-7.745, 10.679, 6.273];
% DPcoordb1 = [-6.745, 10.679, 6.273];
% DP1 = CreateDipolar(DPcoorda1,DPcoordb1);
% CDP = {DP1;DPcoorda1}; % cell structure {DT1,DT2,DT3;DPcoorda1,DPcoorda2,DPcoorda3}
% %%%

%%% User input ends here


fileid = fopen(FileName,'r');
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

%%%% Adjust the camera zoom to the size of the molecule
Xbox=max(xyz(:,1))-min(xyz(:,1));
Ybox=max(xyz(:,2))-min(xyz(:,2));
Zbox=max(xyz(:,3))-min(xyz(:,3));
zoomFactor = min([Xbox Ybox Zbox])/max([Xbox Ybox Zbox]);
%%%%



if PlotTensor==0 % plots the tensor with the molecule
    fig = figure
    surface(X,Y,Z,ColorMap,'EdgeColor','none','FaceLighting','gouraud','AlphaData',AlphaMap,'AlphaDataMapping','none','FaceAlpha','flat','AmbientStrength',0.7);
    hold on
    if NR ~= 0
        for i = 1:NR
            ExtraTensor = cell2mat(C(1,i));
            Extraatomcoord = cell2mat(C(2,i));
            [X,Y,Z,AlphaMap]=CreateTensor(ExtraTensor,Extraatomcoord,OvaloidEllipsoid,Transparency,TensorScale);
            surface(X,Y,Z,ColorMap,'EdgeColor','none','FaceLighting','gouraud','AlphaData',AlphaMap,'AlphaDataMapping','none','FaceAlpha','flat','AmbientStrength',0.7);
        end
    end
    
    set(gcf,'Color','w')
    molecule3D(xyz,labels,styles,filetype,Conlist,Bondlimit) % plots the molecule
    axis equal
    cameratoolbar % Use the camera control, not the "zoom options/magnifying glass"
    %camzoom(zoomFactor)
elseif PlotTensor==1
    fig = figure
    set(gcf,'Color','w')
    molecule3D(xyz,labels,styles,filetype,Conlist,Bondlimit) % plots the molecule
    axis equal
    cameratoolbar % Use the camera control, not the "zoom options/magnifying glass"
    %camzoom(zoomFactor)
elseif PlotTensor==2
    fig = figure
    set(gcf,'Color','w')
    surf(X,Y,Z,ColorMap,'EdgeColor','none','FaceLighting','gouraud','AlphaData',AlphaMap,'AlphaDataMapping','none','FaceAlpha','flat','AmbientStrength',0.7); % 'FaceAlpha',Transparency
    axis equal
    set(gca,'visible','off')
    light('Position',[1 1 2]);
    cameratoolbar
    camtarget('auto')
    %material shiny
else
    error("PlotTensor must be either 0 or 1 or 2")
end


% Active rotation matrices for reference
%
%   rx =    [1, 0, 0;
%            0, cos(rad), -sin(rad);
%            0, sin(rad), cos(rad)];
%
%   ry =    [cos(rad), 0, sin(rad);
%             0, 1, 0;
%             -sin(rad), 0, cos(rad)];
%
%   rz =   [cos(rad), -sin(rad), 0;
%             sin(rad), cos(rad), 0;
%             0, 0, 1];
