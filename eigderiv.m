function y = eigderiv(a,var,data,V,theta,BalanceNode,PVNode) 
%First-order eigenvalue derivatives 
%-----------------------------------------------------------------
%% State matrix sensitivity
phi = var.eigen.phi; 
y.deig.dV     = a.dAsys.dV*phi; 
y.deig.dtheta = a.dAsys.dtheta*phi; 
y.deig.ddelta = a.dAsys.ddelta*phi; 
y.deig.dId    = a.dAsys.dId*phi; 
y.deig.dIq    = a.dAsys.dIq*phi; 
y.deig.dEdp   = a.dAsys.dEdp*phi; 
y.deig.dEqp   = a.dAsys.dEqp*phi; 
y.deig.dEfd   = a.dAsys.dEfd*phi; 
%% System parameter
GenNumber = data.ngen;
GenNode = data.Pgenbuses;
Node = data.nbuses;
%% Active sensitivity calculation
Y = data.Ybus;
U = V.*exp(1i.*theta);                                                     
s = U.*conj(Y*U);                                                          
y1 = real(s);                                                                   
y2 = imag(s);  
dU = sparse(diag(U));                                                      
HIJ = -dU*conj(Y*dU);  
ImagHIJ = imag(HIJ);
RealHIJ = real(HIJ);
y4 = sparse(diag(y1));
y5 = sparse(diag(y2));
Pq0 = ImagHIJ + y5;
Pv0 = (RealHIJ - y4)/diag(V);
Qq0 = -RealHIJ - y4;
Qv0 = (ImagHIJ - y5)/diag(V);

Pq = -Pq0;
Pv = -Pv0 ;
Qq = -Qq0;
Qv = -Qv0 ;
%% PV node
Qv(:,PVNode) = 0;
Qv(PVNode,:) = 0;
Qq(PVNode,:) = 0;
Pv(:,PVNode) = 0;

Pq(:,BalanceNode) = 0;
Pq(BalanceNode,:) = 0;
Pv(:,BalanceNode) = 0;
Pv(BalanceNode,:) = 0;
Qq(:,BalanceNode) = 0;
Qq(BalanceNode,:) = 0;
Qv(:,BalanceNode) = 0;
Qv(BalanceNode,:) = 0;

data1 = find(Pq(:)~=0);
data2 = find(Pv(:)~=0);
data3 = find(Qq(:)~=0);
data4 = find(Qv(:)~=0);

Pq_T = sparse(Node,Node);
Pv_T = sparse(Node,Node);
Qq_T = sparse(Node,Node);
Qv_T = sparse(Node,Node);

Pq_T(data1(:)) = 1./Pq(data1(:));
Pv_T(data2(:)) = 1./Pv(data2(:));
Qq_T(data3(:)) = 1./Qq(data3(:));
Qv_T(data4(:)) = 1./Qv(data4(:));


x0 = Pv_T(GenNode,GenNode);
y0 = Pq_T(GenNode,GenNode);
x1 = Qv_T(GenNode,GenNode);
y1 = Qq_T(GenNode,GenNode);

data.x = diag(x0);
data.y = diag(y0);
data.x1 = diag(x1);
data.y1 = diag(y1);

y.dV_Pg = Pv_T(GenNode,:);
y.dtheta_Pg = Pq_T(GenNode,:);
y.dV_Qg = Qv_T(GenNode,:);
y.dtheta_Qg = Qq_T(GenNode,:);
%% Active sensitivity calculation
%  Solve the equation
b1 = data.b1 - (data.a11*data.x + data.a12*data.y);
b2 = data.b2 - (data.a21*data.x + data.a22*data.y);
b3 = data.b3 - (data.a31*data.x + data.a32*data.y);
b4 = data.b4 - (data.a41*data.x + data.a42*data.y);
A = [data.a13,data.a14,data.a15,sparse(GenNumber,3*GenNumber);...
     data.a23,data.a24,data.a25,sparse(GenNumber,3*GenNumber);...
     data.a33,data.a34,data.a35,data.a36,sparse(GenNumber,2*GenNumber);...
     data.a43,data.a44,data.a45,sparse(GenNumber,GenNumber),data.a47,sparse(GenNumber,GenNumber);...
     sparse(GenNumber,2*GenNumber),data.a55,data.a56,sparse(GenNumber,2*GenNumber);...
     sparse(GenNumber,GenNumber),data.a64,sparse(GenNumber,2*GenNumber),data.a67,data.a68];
b = [b1;b2;b3;b4;data.b5;data.b6];
x1 = A\b;
 
y.ddelta_Pg = x1(1:GenNumber);
y.Id_Pg = x1(GenNumber+1:2*GenNumber);
y.Iq_Pg = x1(2*GenNumber+1:3*GenNumber);
y.Edp_Pg = x1(3*GenNumber+1:4*GenNumber);
y.Eqp_Pg = x1(4*GenNumber+1:5*GenNumber);
y.Efd_Pg = x1(5*GenNumber+1:6*GenNumber);
%% Reactive sensitivity calculation
%  Solve the equation
b01 = data.b01 - (data.a11*data.x1 + data.a12*data.y1);
b02 = data.b02 - (data.a21*data.x1 + data.a22*data.y1);
b03 = data.b03 - (data.a31*data.x1 + data.a32*data.y1);
b04 = data.b04 - (data.a41*data.x1 + data.a42*data.y1);
B = [b01;b02;b03;b04;data.b5;data.b6];
x = A\B;
 
y.ddelta_Qg = x(1:GenNumber);
y.Id_Qg = x(GenNumber+1:2*GenNumber);
y.Iq_Qg = x(2*GenNumber+1:3*GenNumber);
y.Edp_Qg = x(3*GenNumber+1:4*GenNumber);
y.Eqp_Qg = x(4*GenNumber+1:5*GenNumber);
y.Efd_Qg = x(5*GenNumber+1:6*GenNumber);
%% Active power:Closed-form spectral abscissa abscissa derivatives
y.deig.p = y.dV_Pg*y.deig.dV + y.dtheta_Pg*y.deig.dtheta + y.deig.ddelta.*y.ddelta_Pg...
         + y.deig.dId.*y.Id_Pg + y.deig.dIq.*y.Iq_Pg + y.deig.dEdp.*y.Edp_Pg...
         + y.deig.dEqp.*y.Eqp_Pg + y.deig.dEfd.*y.Efd_Pg;

%% Reactive power:Closed-form spectral abscissa abscissa derivatives
y.deig.q = y.dV_Qg*y.deig.dV + y.dtheta_Qg*y.deig.dtheta + y.deig.ddelta.*y.ddelta_Qg...
         + y.deig.dId.*y.Id_Qg + y.deig.dIq.*y.Iq_Qg + y.deig.dEdp.*y.Edp_Qg...
         + y.deig.dEqp.*y.Eqp_Qg + y.deig.dEfd.*y.Efd_Qg;
end

 

