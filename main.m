function  main
%% Date selection dialog box
clear;
clc;                                                              
format short e                                                 
wd = cd ;                                                                  
cd([wd '\Data']);   
[FileName,PathName]= uigetfile('*.txt','CHOOSE DATA');  
if isequal(FileName,0) || isequal(PathName,0)                              %If"CANCEL"
    cd(wd);
    return;
end
cd(wd);
tic
%% Read data
[Node,MaxK,IterationPrecision,BalanceNode,PVNode,V0,P0,Q0,...  
theta0,d,LineB,LineI,LineJ,LineY,TransformerI,TransformerJ,...
TransformerK,TransformerY,GroundI,GroundB,PGIndex] = ReadData(FileName,PathName);         
cd(wd);
%% Powerflow
[V1,theta1,d,P,Q,U] = pf(Node,MaxK,IterationPrecision,BalanceNode,PVNode,...
V0,P0,Q0,theta0,d,LineB,LineI,LineJ,LineY,TransformerI,TransformerJ,...
TransformerK,TransformerY,GroundI,GroundB);
%% Eigenvalue calculation
[data,var] = sss(d,V1,U,theta1,P,Q);
%% Sensitivity calculation
a = asysderiv_debug(d,data,var);
y = eigderiv(a,var,data,V1,theta1,BalanceNode,PVNode);
ym = eigderiv_m(a,var,data,V1,theta1);
%% Closed-form spectral abscissa sensitivity
sensitivity_p = y.deig.p;
sensitivity_p = real(sensitivity_p);
sensitivity_p = sensitivity_p(d.gen_idx);

sensitivity_q = y.deig.q;
sensitivity_q = real(sensitivity_q);
sensitivity_q = sensitivity_q(d.gen_idx);
%% Mathematical spectral abscissa sensitivity
sensitivity_p1 = ym.deig.p;
sensitivity_p1 = real(sensitivity_p1);
sensitivity_p1 = sensitivity_p1(d.gen_idx);

sensitivity_q1 = ym.deig.q;
sensitivity_q1 = real(sensitivity_q1);
sensitivity_q1 = sensitivity_q1(d.gen_idx);
%% Order   
o_p = [PGIndex sensitivity_p sensitivity_p1];
o_p = sortrows(o_p,1);
o_q = [PGIndex sensitivity_q sensitivity_q1];
o_q = sortrows(o_q,1);
%% Display
disp('Active power:CFSAS vs MSAD');
disp(o_p(:,2:3));
disp('Reactive power:CFSAS vs MSAD');
disp(o_q(:,2:3));
end

