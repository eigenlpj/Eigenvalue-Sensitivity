function a=asysderiv_debug(d,data,var)
%*****************************************************************
%var.eigenvalue derivatives
%*****************************************************************

%%----------------some common value---------------------------------------
D1 = var.eigen.D1;
D3 = var.eigen.D3;
D5 = var.eigen.D5;
D6 = var.eigen.D6;
C1 = var.eigen.C1;
C2 = var.eigen.C2;
B2 = var.eigen.B2; 
B = var.eigen.B;
C = var.eigen.C; 
D = var.eigen.D; 
zero1 = var.eigen.zero1;
theta = var.theta;
delta = var.delta;
Iq = var.Iq;
GenNumber = data.ngen;
GenNode = data.Pgenbuses;
GenNodeOff = data.nongen_ind;
Node = data.nbuses;

M = data.M;
Ae = data.machine.Ae;
Be = data.machine.Be;
Efd = var.Efd;
TE = data.machine.TE;
Xdp = data.machine.Xdp;
Xqp = data.machine.Xqp;

ejTheta = exp(1i*theta);
ejdelta_c = exp(-1i*delta);
Vp = var.V.* ejTheta;
Ig_c = var.Id - 1i*Iq;
psi = var.eigen.psi.' ;

B_D = B/D;
psi_B_D = psi.'*B_D;
psi_B_D = psi_B_D.';

Vdelta = Vp(GenNode).*ejdelta_c;
ejTheta_delta = ejTheta(GenNode).*ejdelta_c;
Vtheta_delta = Vdelta;
Ig_Vdelta = Ig_c.*Vdelta;
Ig_Theta_delta = Ig_c.*ejTheta_delta;

%% ---------------------A---------------------------------------------------
A1Id_row_ind = 2:7:7*GenNumber;
A1Id_col_ind = 3:7:7*GenNumber;
A1Efd_row_ind = 5:7:7*GenNumber;
A1Efd_col_ind = 5:7:7*GenNumber;

A1_Iq_vec = -psi(A1Id_row_ind)./M;
a.dA1.dIq = sparse(1:GenNumber, A1Id_col_ind,A1_Iq_vec,GenNumber,7*GenNumber);

A1_Id_vec = A1_Iq_vec;
a.dA1.dId = sparse(1:GenNumber, A1Id_col_ind + 1,A1_Id_vec,GenNumber,7*GenNumber);

f = Ae.*Be.*exp(Be.*Efd);
f1 = Ae.*Be.^2.*exp(Be.*Efd);
A1_Efd_vec = -psi(A1Efd_row_ind)./TE.*(2*f + Efd.*f1);
a.dA1.dEfd = sparse(1:GenNumber,A1Efd_col_ind,A1_Efd_vec,GenNumber,7*GenNumber);

a.dA.dId = a.dA1.dId;
a.dA.dIq = a.dA1.dIq;
a.dA.dEfd = a.dA1.dEfd;
%% ---------------------B---------------------------------------------------
B1row_ind = 2:7:7*GenNumber;
B1col_ind = 1:2:2*GenNumber;
B1_Id_vec = psi(B1row_ind).*(Xdp - Xqp)./M;
a.dB1.dId =  sparse(1:GenNumber,B1col_ind + 1,B1_Id_vec,GenNumber,2*GenNumber);
B1_Iq_vec = B1_Id_vec;
a.dB1.dIq = sparse(1:GenNumber,B1col_ind,B1_Iq_vec,GenNumber,2*GenNumber);
B1_Edp_vec = -psi(B1row_ind)./M;
a.dB1.dEdp = sparse(1:GenNumber,B1col_ind,B1_Edp_vec,GenNumber,2*GenNumber);
a.dB1.dEqp = sparse(1:GenNumber,B1col_ind + 1,B1_Edp_vec,GenNumber,2*GenNumber);

[rowdB1,~] = size(a.dB1.dId); 
[~,columnB2] = size(B2);
[~,columnzero1] = size(zero1);
zeroB2 = sparse(rowdB1,columnB2);
zero5 = sparse(rowdB1,columnzero1);

a.dB.dId = [a.dB1.dId,zeroB2,zero5];
a.dB.dIq = [a.dB1.dIq,zeroB2,zero5];
a.dB.dEdp = [a.dB1.dEdp,zeroB2,zero5];
a.dB.dEqp = [a.dB1.dEqp,zeroB2,zero5];
%% ----------------------C--------------------------------------------------
[rowC1,~] = size(C1); 
[rowC2,~] = size(C2); 
psi_B_D_C1 = psi_B_D(1:rowC1);
psi_B_D_C2 = psi_B_D(rowC1+1:rowC1+rowC2);
C1row_ind = 1:2:2*GenNumber;
C1col_ind =1:7:7*GenNumber;  
delta_ind =1:7:GenNumber*7;  
real_ind = 1:GenNumber;
imag_ind = real_ind + GenNumber; 
            
C1_V_real_vec =  psi_B_D_C1(C1row_ind).*real(-ejTheta_delta);
C1_V_imag_vec =  psi_B_D_C1(C1row_ind+1).*imag(-ejTheta_delta);             
a.dC1.dV = sparse(GenNode,C1col_ind,C1_V_real_vec + C1_V_imag_vec,Node,7*GenNumber);  
C2_V_vec = psi_B_D_C2(real_ind).*real(Ig_Theta_delta)+ psi_B_D_C2(imag_ind).*imag(Ig_Theta_delta);
a.dC2.dV  =  sparse(GenNode,delta_ind,C2_V_vec,Node,7*GenNumber);  
a.dC.dV = a.dC1.dV + a.dC2.dV;

C1_theta_real_vec =  psi_B_D_C1(C1row_ind).* real(-1i.*Vtheta_delta);
C1_theta_imag_vec =  psi_B_D_C1(C1row_ind+1).*imag(-1i.*Vtheta_delta);
a.dC1.dtheta = sparse(GenNode,C1col_ind,C1_theta_real_vec + C1_theta_imag_vec,Node,7*GenNumber);  
C2_theta_vec = psi_B_D_C2(real_ind).*real(1i.*Ig_Vdelta) + psi_B_D_C2(imag_ind).*imag(1i.*Ig_Vdelta);  
a.dC2.dtheta  = sparse(GenNode,delta_ind,C2_theta_vec,Node,7*GenNumber);  
a.dC.dtheta = a.dC1.dtheta + a.dC2.dtheta;

C1_delta_real_vec = psi_B_D_C1(C1row_ind).*real(1i.*Vtheta_delta);
C1_delta_imag_vec = psi_B_D_C1(C1row_ind+1).*imag(1i.*Vtheta_delta);      
a.dC1.ddelta = sparse(1:GenNumber,C1col_ind, C1_delta_real_vec,GenNumber,7*GenNumber)...
             + sparse(1:GenNumber,C1col_ind, C1_delta_imag_vec,GenNumber,7*GenNumber);
C2_delta_real_vec = psi_B_D_C2(real_ind).*real(-1i.*Ig_Vdelta);
C2_delta_imag_vec = psi_B_D_C2(imag_ind).*imag(-1i.*Ig_Vdelta);       
a.dC2.ddelta =  sparse(1:GenNumber,delta_ind,C2_delta_real_vec,GenNumber,7*GenNumber)...
              + sparse(1:GenNumber,delta_ind,C2_delta_imag_vec,GenNumber,7*GenNumber);  
a.dC.ddelta = a.dC1.ddelta + a.dC2.ddelta;

C2_Id_real_vec = psi_B_D_C2(real_ind).*real(Vtheta_delta);
C2_Id_imag_vec = psi_B_D_C2(imag_ind).*imag(Vtheta_delta);
a.dC2.dId =sparse(1:GenNumber,delta_ind,C2_Id_real_vec,GenNumber,7*GenNumber)...
         + sparse(1:GenNumber,delta_ind,C2_Id_imag_vec,GenNumber,7*GenNumber);
a.dC.dId = a.dC2.dId;
        
C2_Iq_real_vec = psi_B_D_C2(real_ind).*real(-1i.*Vtheta_delta);
C2_Iq_imag_vec = psi_B_D_C2(imag_ind).*imag(-1i.*Vtheta_delta);     
a.dC2.dIq =sparse(1:GenNumber,delta_ind,C2_Iq_real_vec,GenNumber,7*GenNumber)...
         + sparse(1:GenNumber,delta_ind,C2_Iq_imag_vec,GenNumber,7*GenNumber);  
a.dC.dIq = a.dC2.dIq;   
%% ---------------------D---------------------------------------------------
[rowD1,~] = size(D1); 
[rowD3,~] = size(D3); 
[rowD6,~] = size(D6); 
psi_B_D_D1 = psi_B_D(1:rowD1);
psi_B_D_D3 = psi_B_D(rowD1+1:rowD1+rowD3);
psi_B_D_D6 = psi_B_D(rowD1+rowD3+1:rowD1+rowD3+rowD6);
Iq_ind = 2:2:2*GenNumber;
Id_ind = Iq_ind - 1;
%---------------------D3---------------------------------------------------       
D3_V_real_vec = psi_B_D_D3(real_ind).*real(1i.*ejTheta_delta) + psi_B_D_D3(imag_ind).*imag(1i.*ejTheta_delta);
D3_V_imag_vec = psi_B_D_D3(real_ind).*real(ejTheta_delta) + psi_B_D_D3(imag_ind).*imag(ejTheta_delta);
a.dD3.dV = sparse(GenNode,Id_ind,D3_V_real_vec,Node,2*GenNumber)+ sparse(GenNode,Iq_ind,D3_V_imag_vec,Node,2*GenNumber);      
          
D3_theta_real_vec = psi_B_D_D3(real_ind).*real(-Vtheta_delta)+ psi_B_D_D3(imag_ind).*imag(-Vtheta_delta);
D3_theta_imag_vec = psi_B_D_D3(real_ind).*real(1i*Vtheta_delta)+ psi_B_D_D3(imag_ind).*imag(1i*Vtheta_delta);
a.dD3.dtheta = sparse(GenNode,Id_ind,D3_theta_real_vec,Node,2*GenNumber)+ sparse(GenNode,Iq_ind,D3_theta_imag_vec,Node,2*GenNumber);     
             
D3_delta_real_vec = psi_B_D_D3(real_ind).*real(Vtheta_delta) + psi_B_D_D3(imag_ind).*imag(Vtheta_delta);
D3_delta_imag_vec = psi_B_D_D3(real_ind).*real(-1i.*Vtheta_delta) + psi_B_D_D3(imag_ind).*imag(-1i.*Vtheta_delta);             
a.dD3.ddelta = sparse(1:GenNumber,Id_ind,D3_delta_real_vec,GenNumber,2*GenNumber)+ sparse(1:GenNumber,Iq_ind,D3_delta_imag_vec,GenNumber,2*GenNumber); 
% ---------------------D2---------------------------------------------------             
D2_real_ind = 1:2:2*GenNumber;
D2_imag_ind = D2_real_ind +1;

D2_V_vec = psi_B_D_D1(D2_real_ind).*real(ejTheta_delta) + psi_B_D_D1(D2_imag_ind).*imag(ejTheta_delta);         
a.dD2.dV = sparse(GenNode,1 + GenNumber:2*GenNumber,D2_V_vec,Node,2*GenNumber);

D2_dtheta_vec = psi_B_D_D1(D2_real_ind).*real(1i.*Vtheta_delta) + psi_B_D_D1(D2_imag_ind).*imag(1i.*Vtheta_delta);         
a.dD2.dtheta = sparse(GenNode,1 + GenNumber:2*GenNumber,D2_dtheta_vec,Node,2*GenNumber);
D2_dtheta_vec = psi_B_D_D1(D2_real_ind).*real(ejTheta_delta) + psi_B_D_D1(D2_imag_ind).*imag(ejTheta_delta); 
a.dD2.dtheta = a.dD2.dtheta + sparse(GenNode,1:GenNumber,D2_dtheta_vec,Node,2*GenNumber);
           
D2_delta_imag_vec = psi_B_D_D1(D2_real_ind).*real(-ejTheta_delta)+ psi_B_D_D1(D2_imag_ind).*imag(-ejTheta_delta);              
D2_delta_real_vec = psi_B_D_D1(D2_real_ind).*real(-1i*Vtheta_delta)+ psi_B_D_D1(D2_imag_ind).*imag(-1i*Vtheta_delta);                        
a.dD2.ddelta = sparse(1:GenNumber,1:GenNumber,D2_delta_imag_vec,GenNumber,2*GenNumber)...
             + sparse(1:GenNumber,1 + GenNumber:2*GenNumber,D2_delta_real_vec,GenNumber,2*GenNumber);
% ---------------------D4---------------------------------------------------                     
Yp = psi_B_D_D3(1:GenNumber);
Yq = psi_B_D_D3(GenNumber + 1:2*GenNumber);   
D4_V_vec = psi_B_D_D3(real_ind).*real(-Ig_Theta_delta) + psi_B_D_D3(imag_ind).*imag(-Ig_Theta_delta);   
D4_theta_real_vec = psi_B_D_D3(real_ind).*real(-Ig_Theta_delta) + psi_B_D_D3(imag_ind).*imag(-Ig_Theta_delta);
D4_theta_imag_vec = psi_B_D_D3(real_ind).*real(-1i.*Ig_Vdelta) + psi_B_D_D3(imag_ind).*imag(-1i.*Ig_Vdelta);  
[Pvv_gen,Qvv_gen,Pva_gen,Qva_gen,Pav_gen,Qav_gen,Paa_gen,Qaa_gen] = cal_Hessian(var,data,GenNode,Yp,Yq);
a.dD4.dV = -[Pvv_gen(:,GenNode) + Qvv_gen(:,GenNode)  Pav_gen(:,GenNode) + Qav_gen(:,GenNode)];  
a.dD4.dV = a.dD4.dV +  sparse(GenNode ,1 + GenNumber:2*GenNumber,D4_V_vec,Node,2*GenNumber); 
a.dD4.dtheta = -[Pva_gen(:,GenNode) + Qva_gen(:,GenNode)  Paa_gen(:,GenNode) + Qaa_gen(:,GenNode)];            
a.dD4.dtheta = a.dD4.dtheta + sparse([GenNode  GenNode],[1:GenNumber  1+GenNumber:2*GenNumber ],...
            [D4_theta_real_vec D4_theta_imag_vec],Node,2*GenNumber); 
        
D4_Id_real_vec = psi_B_D_D3(real_ind).*real(1i.*ejTheta_delta) + psi_B_D_D3(imag_ind).*imag(1i.*ejTheta_delta);
D4_Id_imag_vec = psi_B_D_D3(real_ind).*real(-Vtheta_delta) + psi_B_D_D3(imag_ind).*imag(-Vtheta_delta);
a.dD4.dId =sparse([1:GenNumber 1:GenNumber],[1:GenNumber  1 + GenNumber:2*GenNumber],...
           [D4_Id_real_vec  D4_Id_imag_vec], GenNumber,2*GenNumber) ;
             
D4_Iq_real_vec = psi_B_D_D3(real_ind).*real(ejTheta_delta) + psi_B_D_D3(imag_ind).*imag(ejTheta_delta);
D4_Iq_imag_vec = psi_B_D_D3(real_ind).*real(1i.*Vtheta_delta) + psi_B_D_D3(imag_ind).*imag(1i.*Vtheta_delta);        
a.dD4.dIq =sparse([1:GenNumber 1:GenNumber],[1:GenNumber  1 + GenNumber:2*GenNumber],...
           [D4_Iq_real_vec  D4_Iq_imag_vec],GenNumber,2*GenNumber) ;        

D4_delta_real_vec = psi_B_D_D3(real_ind).*real(Ig_Theta_delta) + psi_B_D_D3(imag_ind).*imag(Ig_Theta_delta);
D4_delta_imag_vec = psi_B_D_D3(real_ind).*real(1i.*Ig_Vdelta) + psi_B_D_D3(imag_ind).*imag(1i.*Ig_Vdelta);    
a.dD4.ddelta =sparse([1:GenNumber 1:GenNumber],[1:GenNumber  1 + GenNumber:2*GenNumber],...
           [D4_delta_real_vec  D4_delta_imag_vec],GenNumber,2*GenNumber) ;  
%---------------------D5---------------------------------------------------     
a.dD5.dV = -[Pvv_gen(:,GenNodeOff) + Qvv_gen(:,GenNodeOff)   Pav_gen(:,GenNodeOff) + Qav_gen(:,GenNodeOff)];  
a.dD5.dtheta = -[Pva_gen(:,GenNodeOff) + Qva_gen(:,GenNodeOff)   Paa_gen(:,GenNodeOff) + Qaa_gen(:,GenNodeOff)];

%---------------------D6---------------------------------------------------  
Yp = psi_B_D_D6(1:Node-GenNumber);
Yq = psi_B_D_D6(Node-GenNumber+1:2*(Node-GenNumber));   
[Pvv_load,Qvv_load,Pva_load,Qva_load,Pav_load,Qav_load,Paa_load,Qaa_load] = cal_Hessian(var,data,GenNodeOff,Yp,Yq);
a.dD6.dV = -[(Pvv_load(:,GenNode) + Qvv_load(:,GenNode))  (Pav_load(:,GenNode) + Qav_load(:,GenNode))];    
a.dD6.dtheta =  -[(Pva_load(:,GenNode) + Qva_load(:,GenNode))  (Paa_load(:,GenNode) + Qaa_load(:,GenNode))];   

%---------------------D7---------------------------------------------------        
a.dD7.dV = -[(Pvv_load(:,GenNodeOff) + Qvv_load(:,GenNodeOff))  (Pav_load(:,GenNodeOff) + Qav_load(:,GenNodeOff))];             
a.dD7.dtheta = -[(Pva_load(:,GenNodeOff) + Qva_load(:,GenNodeOff))  (Paa_load(:,GenNodeOff) + Qaa_load(:,GenNodeOff))];

[~,columnD1] = size(D1);
[~,columnD5] = size(D5);
zero6 = sparse(GenNumber,columnD1);
zero7 = sparse(GenNumber,columnD5);
%---------------------D--------------------------------------------------- 
a.dD.dId = [zero6,a.dD4.dId,zero7];
a.dD.dIq = [zero6,a.dD4.dIq,zero7];
a.dD.dV = [a.dD3.dV,a.dD2.dV+a.dD4.dV+a.dD6.dV,a.dD5.dV+a.dD7.dV];
a.dD.dtheta = [a.dD3.dtheta,a.dD2.dtheta+a.dD4.dtheta+a.dD6.dtheta,a.dD5.dtheta+a.dD7.dtheta];
a.dD.ddelta = [a.dD3.ddelta,a.dD2.ddelta+a.dD4.ddelta,zero7];
%% Asys derivatives
D_C = D\C;   
a.dAsys.dV =  - a.dC.dV + a.dD.dV*D_C;   
a.dAsys.dtheta = - a.dC.dtheta + a.dD.dtheta*D_C; 
a.dAsys.ddelta = - a.dC.ddelta + a.dD.ddelta*D_C;  
a.dAsys.dId = a.dA.dId - a.dB.dId*D_C - a.dC.dId + a.dD.dId*D_C;        
a.dAsys.dIq = a.dA.dIq - a.dB.dIq*D_C - a.dC.dIq + a.dD.dIq*D_C;  
a.dAsys.dEfd = a.dA.dEfd;
a.dAsys.dEdp = -a.dB.dEdp*D_C;
a.dAsys.dEqp = -a.dB.dEqp*D_C;
end
