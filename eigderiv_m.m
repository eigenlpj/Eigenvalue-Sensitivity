function ym = eigderiv_m(a,var,data,V,theta) 
%First-order eigenvalue derivatives 
%-----------------------------------------------------------------
%% State matrix sensitivity
phi = var.eigen.phi; 
ym.deig.dV     = a.dAsys.dV*phi; 
ym.deig.dtheta = a.dAsys.dtheta*phi; 
ym.deig.ddelta = a.dAsys.ddelta*phi; 
ym.deig.dId    = a.dAsys.dId*phi; 
ym.deig.dIq    = a.dAsys.dIq*phi; 
ym.deig.dEdp   = a.dAsys.dEdp*phi; 
ym.deig.dEqp   = a.dAsys.dEqp*phi; 
ym.deig.dEfd   = a.dAsys.dEfd*phi; 
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
%%
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

ym.dV_Pg = Pv_T(GenNode,:);
ym.dtheta_Pg = Pq_T(GenNode,:);
ym.dV_Qg = Qv_T(GenNode,:);
ym.dtheta_Qg = Qq_T(GenNode,:);
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
 
ym.ddelta_Pg = x1(1:GenNumber);
ym.Id_Pg = x1(GenNumber+1:2*GenNumber);
ym.Iq_Pg = x1(2*GenNumber+1:3*GenNumber);
ym.Edp_Pg = x1(3*GenNumber+1:4*GenNumber);
ym.Eqp_Pg = x1(4*GenNumber+1:5*GenNumber);
ym.Efd_Pg = x1(5*GenNumber+1:6*GenNumber);
%% Reactive sensitivity calculation
%  Solve the equation
b01 = data.b01 - (data.a11*data.x1 + data.a12*data.y1);
b02 = data.b02 - (data.a21*data.x1 + data.a22*data.y1);
b03 = data.b03 - (data.a31*data.x1 + data.a32*data.y1);
b04 = data.b04 - (data.a41*data.x1 + data.a42*data.y1);
B = [b01;b02;b03;b04;data.b5;data.b6];
x = A\B;
 
ym.ddelta_Qg = x(1:GenNumber);
ym.Id_Qg = x(GenNumber+1:2*GenNumber);
ym.Iq_Qg = x(2*GenNumber+1:3*GenNumber);
ym.Edp_Qg = x(3*GenNumber+1:4*GenNumber);
ym.Eqp_Qg = x(4*GenNumber+1:5*GenNumber);
ym.Efd_Qg = x(5*GenNumber+1:6*GenNumber);
%% Active power:Mathematical spectral abscissa derivatives
ym.deig.p = ym.dV_Pg*ym.deig.dV + ym.dtheta_Pg*ym.deig.dtheta + ym.deig.ddelta.*ym.ddelta_Pg...
         + ym.deig.dId.*ym.Id_Pg + ym.deig.dIq.*ym.Iq_Pg + ym.deig.dEdp.*ym.Edp_Pg...
         + ym.deig.dEqp.*ym.Eqp_Pg + ym.deig.dEfd.*ym.Efd_Pg;

%% Reactive power:Mathematical spectral abscissa derivatives
ym.deig.q = ym.dV_Qg*ym.deig.dV + ym.dtheta_Qg*ym.deig.dtheta + ym.deig.ddelta.*ym.ddelta_Qg...
         + ym.deig.dId.*ym.Id_Qg + ym.deig.dIq.*ym.Iq_Qg + ym.deig.dEdp.*ym.Edp_Qg...
         + ym.deig.dEqp.*ym.Eqp_Qg + ym.deig.dEfd.*ym.Efd_Qg;

end

 

