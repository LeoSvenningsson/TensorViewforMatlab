function PlotAllRef(Tensor,x,y,z,Order,Scale,Brightness)
SymTensor = (Tensor+transpose(Tensor))/2; % symmetrizes the tensor
[U, D] = eig(SymTensor, 'vector'); % Rotation matrix is gained from the Eigenvectors
Sym1 = length(D)-length(unique(round(D,7)));

if Order=="Ascending"
    [D, ind] = sort(D,'ascend');
    U = U(:, ind);
end
if Order=="Descending"
    [D, ind] = sort(D,'descend');
    U = U(:, ind);
end
if Order=="Absascending"
    [~, ind] = sort(abs(D),'ascend');
    U = U(:, ind);
    D = D(ind);
end

R = U;

Xref = R*[1;0;0]*Scale;
Yref = R*[0;1;0]*Scale;
Zref = R*[0;0;1]*Scale;

BrR1 = min([1 (1+Brightness)]);
BrG1 = min([1 (Brightness)]);
BrB1 = min([1 (Brightness)]);

BrG1 = max([0 (Brightness)]);
BrB1 = max([0 (Brightness)]);

BrR2 = min([1 (Brightness)]);
BrG2 = min([1 (1+Brightness)]);
BrB2 = min([1 (Brightness)]);

BrR2 = max([0 (Brightness)]);
BrB2 = max([0 (Brightness)]);

BrR3 = min([1 (Brightness)]);
BrG3 = min([1 (Brightness)]);
BrB3 = min([1 (1+Brightness)]);

BrR3 = max([0 (Brightness)]);
BrG3 = max([0 (Brightness)]);

hold on
quiver3(x,y,z,Xref(1),Xref(2),Xref(3),'Color',[BrR1,BrG1,BrB1],'LineWidth',3,'AutoScale','off')
axis equal
set(gcf,'Color','w')
set(gca,'visible','off')
quiver3(x,y,z,Yref(1),Yref(2),Yref(3),'Color',[BrR2,BrG2,BrB2],'LineWidth',3,'AutoScale','off')
quiver3(x,y,z,Zref(1),Zref(2),Zref(3),'Color',[BrR3,BrG3,BrB3],'LineWidth',3,'AutoScale','off')


end