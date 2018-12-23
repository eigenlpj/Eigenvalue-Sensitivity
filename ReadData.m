%% ***************************************************************
function [Node,MaxK,IterationPrecision,BalanceNode,PVNode,VoltageAmplitude,P,Q,VoltageAngle,d,LineB,LineI,LineJ,LineY,TransformerI,TransformerJ,TransformerK,TransformerY,GroundI,GroundB,PGIndex] = ReadData(FileName,PathName)
%% 
cd(PathName); 
Data = dlmread(FileName);
Node = Data(1,1);
Row0 = find(Data(:,1)==0);
BalanceNode = Data(3,2);
d.BalanceNode = BalanceNode;
d.Refs = Data(1,3);
PVNode = Data(Row0(5) + 1:Row0(6)-1,1);
MaxK = Data(1,4);
IterationPrecision = Data(2,1);    
%% Network parameter
Line = sparse(Data(Row0(1) + 1:Row0(2)-1,:));                              %line data
Ground = sparse(Data(Row0(2) + 1:Row0(3)-1,:));                            %ground data
Transformer = sparse(Data(Row0(3) + 1:Row0(4)-1,:));                       %transformer data
LineR = Line(:,4);
LineX = Line(:,5);
LineB = Line(:,6);
LineI = Line(:,2);
LineJ = Line(:,3);
LineY = 1./(LineR + 1j*LineX);
TransformerR = Transformer(:,4);
TransformerX = Transformer(:,5);
TransformerI = Transformer(:,2);
TransformerJ = Transformer(:,3);
TransformerK = Transformer(:,6);
TransformerY = 1./(TransformerR + 1j*TransformerX);
GroundI = Ground(:,1);
GroundB = Ground(:,2);
PQ = Data(Row0(4) + 1:Row0(5) - 1,1);                                      
%% Gnerator consumption characteristics
GenMatrix = Data(Row0(6) + 1:Row0(7)-1,:); 
d.gen_idx = find(GenMatrix(:,7) == 1);
PIndex  = GenMatrix(:,1);
PGIndex = PIndex(d.gen_idx); 
m = size(PIndex ,1);
d.C0 = GenMatrix(:,2);
d.C1 = GenMatrix(:,3);
d.C2 = GenMatrix(:,4);
d.Pmin = GenMatrix(:,5)/d.Refs;
d.Pmax = GenMatrix(:,6)/d.Refs; 
d.fr = GenMatrix(:,1);                                  
d.to = GenMatrix(:,2);
d.NodeNum = Node;                                                          
d.PNum = m ;                                                               
d.PGenNumber = size(PGIndex,1);                                            
d.PIndex = GenMatrix(:,1);                                                 
d.PGIndex = PGIndex;                                                       

%% Generator parameter
MachineMatrix = Data(Row0(7) + 1:Row0(8)-1,:); 
d.Xd = MachineMatrix(:,2);
d.Xdp = MachineMatrix(:,3);
d.Xq = MachineMatrix(:,4);
d.Xqp = MachineMatrix(:,5);
d.Rs = MachineMatrix(:,6);
d.D = MachineMatrix(:,7);
d.Tdop = MachineMatrix(:,8);
d.Tqop = MachineMatrix(:,9);
d.H = MachineMatrix(:,10);
%% Exciter parameter
ExciterMatrix = Data(Row0(8) + 1:Row0(9)-1,:);                             
d.KA = ExciterMatrix(:,2);
d.TA = ExciterMatrix(:,3);
d.KE = ExciterMatrix(:,4);
d.TE = ExciterMatrix(:,5);
d.KF = ExciterMatrix(:,6);
d.TF = ExciterMatrix(:,7);
d.Ae = ExciterMatrix(:,8);
d.Be = ExciterMatrix(:,9);
%% Initial value
VoltageAmplitude = ones(Node,1);
VoltageAmplitude(PVNode) = Data(Row0(5) + 1:Row0(6) - 1,2);
P = sparse(PQ,1,Data(Row0(4) + 1:Row0(5) - 1,2) - Data(Row0(4) + 1:Row0(5) - 1,4),Node,1)/d.Refs;  
Q = sparse(PQ,1,Data(Row0(4) + 1:Row0(5) - 1,3) - Data(Row0(4) + 1:Row0(5) - 1,5),Node,1)/d.Refs;  
d.genP = sparse(PQ,1,Data(Row0(4) + 1:Row0(5) - 1,2),Node,1)/d.Refs;   
d.genQ = sparse(PQ,1,Data(Row0(4) + 1:Row0(5) - 1,3),Node,1)/d.Refs;
d.PL = sparse(PQ,1,Data(Row0(4) + 1:Row0(5) - 1,4),Node,1)/d.Refs;
d.QL = sparse(PQ,1,Data(Row0(4) + 1:Row0(5) - 1,5),Node,1)/d.Refs;
VoltageAngle = sparse(Node,1);
d.Qmin = Data(Row0(5) + 1:Row0(6) - 1,3)/d.Refs;
d.Qmax = Data(Row0(5) + 1:Row0(6) - 1,4)/d.Refs;

end



