function [Y,d] = matrixY(d,Node,LineB,LineI,LineJ,LineY,TransformerI,TransformerJ,TransformerK,TransformerY,GroundI,GroundB)
%% Line parameter
B0 = sparse([LineI LineJ],[LineJ LineI],[LineB LineB],Node,Node); 
Y0 = sparse([LineI LineJ],[LineJ LineI],-[LineY LineY],Node,Node); 
Y0 = Y0 + sparse(1:Node,1:Node,sum(-Y0 + 1i.*B0,1),Node,Node);             
%% Transformer parameter
Y1 = sparse([TransformerI TransformerJ],[TransformerJ TransformerI],-[TransformerY./TransformerK TransformerY./TransformerK],Node,Node); 
Y1 = Y1 + sparse([TransformerI TransformerJ],[TransformerI TransformerJ],[TransformerY./TransformerK.^2 TransformerY],Node,Node);
                                                   
Y2 = sparse(GroundI,GroundI,1i.*GroundB,Node,Node);
Y = Y0 + Y1 + Y2; 
d.Ybus = Y ;
[IndexI,IndexJ,nzeroY] = find(Y);
d.IndexI = IndexI;
d.IndexJ = IndexJ;
d.ImpAngle = angle(nzeroY);
d.Yabs = abs(nzeroY);
end

