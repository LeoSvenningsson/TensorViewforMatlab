function TensorViewApp(Fname,s11,s21,s31,s12,s22,s32,s13,s23,s33,x,y,z,CSARef,Scale,Trans,R,G,B,shieldShift,OvEl,MolTen,Bond,MoreT)
%
%   TensorView for Matlab is a tool to visualize chemical shift/shielding tensors in a
%   molecular context. TensorView for Matlab can read arbitrary .pdb and .xyz files
%   for molecular visualization. Any 3D tensor can be used for visualisation though
%   the chemical shift/shielding tensor is used as an input in script.
%
%   TensorView for Matlab is licenced with creative commons CC BY.
%   https://creativecommons.org/licenses/ Free to share and adapt.
%   Give appropriate credits to authors:
%   github.com/LeoSvenningsson/TensorViewforMatlab
%   mathworks.com/matlabcentral/fileexchange/55231-molecule3d
%   onlinelibrary.wiley.com/doi/full/10.1002/mrc.4793
%
%   Version: 1.12
%
%   Authors: Dr. Leo Svenningsson (leo.svenningsson@chalmers.se)
%            Dr. André Ludwig (aludwig@alumni.ethz.ch)
%            Prof. Leonard Mueller (leonard.mueller@ucr.edu)
%
%   TensorView for Matlab is a collaboration with works derrived from
%   molecule3D.m (André Ludwig: mathworks.com/matlabcentral/fileexchange/55231-molecule3d)
%   and TensorView (Prof. Leonard Mueller: https://onlinelibrary.wiley.com/doi/full/10.1002/mrc.4793  and https://sites.google.com/ucr.edu/Mueller/home/tensorview)
%
%%% User input starts here
FileName = Fname; % the molecular file in either .pdb or .xyz format

Tensor = [ s11  s12  s13;
    s21  s22  s23;
    s31  s32  s33  ]; % The molecular frame tensor in the either non-symmetric or symmetric form. The script will transform to the symetrical version regardles of input.

atomcoord = [x, y, z]; % Atom coordinate for the CSA tensor

CSAref = CSARef; % reference shift to go from "chemical shielding" to "chemical shift"


ShieldingShift = shieldShift; % 0 for Shielding;  1 for Shift % The app is always in shielding mode
OvaloidEllipsoid = OvEl; % "ovaloid" for ovaloid;  "ellipsoid" for elipsoid tensor %

TensorScale = Scale; % Tensor scaling

Color = [R,G,B]; % color of the CSA tensor

PlotTensor = MolTen; % 0 Plots the tensor with the molecule; 1 plots only the molecule; 2 only the tensor

styles = 'ballstick'; % 'ballstick','licorice','large','superlarge'

Transparency = Trans; % Sets the transperancy of the tensor

Bondlimit = Bond; % 1.6Å This is only used for the .xyz and connectionless .pdb format


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

[X,Y,Z,AlphaMap]=CreateTensor(Tensor,atomcoord,CSAref,OvaloidEllipsoid,ShieldingShift,Transparency,TensorScale);


if fileid ~= -1
    %%%% Adjust the camera zoom to the size of the molecule
    Xbox=max(xyz(:,1))-min(xyz(:,1));
    Ybox=max(xyz(:,2))-min(xyz(:,2));
    Zbox=max(xyz(:,3))-min(xyz(:,3));
    zoomFactor = min([Xbox Ybox Zbox])/max([Xbox Ybox Zbox]);
    %%%%
end


if PlotTensor==0 % plots the tensor with the molecule
    figure
    surface(X,Y,Z,'FaceColor',Color,'EdgeColor','none','FaceLighting','gouraud','AlphaData',AlphaMap,'AlphaDataMapping','none','FaceAlpha','interp','AmbientStrength',0.7);
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
    figure
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
        figure
    end
    set(gcf,'Color','w')
    surf(X,Y,Z,'FaceColor',Color,'EdgeColor','none','FaceLighting','gouraud','AlphaData',AlphaMap,'AlphaDataMapping','none','FaceAlpha','interp','AmbientStrength',0.7); % 'FaceAlpha',Transparency
    axis equal
    set(gca,'visible','off')
    if MoreT == false
        light('Position',[2 2 4]);
        light('Position',[-2 -2 -4]);
    end
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