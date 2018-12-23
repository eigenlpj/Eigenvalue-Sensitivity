function [data,var] =sss(d,V,U,theta,P,Q)  
f = 60;
GenNumber = d.PNum;
Node = d.NodeNum;
GenNode = d.PIndex;
Ybus = d.Ybus;
Pg = d.genP(GenNode);
Qg = Q(GenNode) + d.QL(GenNode);
%-----------------------------------------------------------------------
% Determine inital operating point 
%-----------------------------------------------------------------------
Ix = (Pg-1i*Qg)./conj(U(GenNode));
delta = angle(U(GenNode) + 1i*d.Xq.*Ix);
Im = Ix.*exp(-1i*(delta - 0.5*pi));
Id = real(Im);
Iq = imag(Im);
Etp = U(GenNode); 
Et = V(GenNode); 
nongen_ind = 1:Node;
nongen_ind(GenNode) = [];  
%-----------------------------------------------------------------------
%Calculate Edp,Eqp,Efd,VR,Rf,Vref,TM 
%-----------------------------------------------------------------------
Vq = imag(V(GenNode).*exp(1i*theta(GenNode)).*exp(-1i*(delta-pi/2)));
Ig = zeros(Node,1); 
Ig(GenNode) = Id + 1i*Iq; 
Edp = (d.Xq - d.Xqp).*Iq;
Eqp = Vq + d.Xdp.*Id;
Efd = Eqp +(d.Xd - d.Xdp).*Id;
dS_dEfd = d.Be.*d.Ae.*exp(d.Be.*full(Efd)); 
SE = d.Ae.*exp(d.Be.*Efd); 
fsi = (d.KE + Efd.*dS_dEfd + SE)./d.TE; 
M = 2*d.H/2/pi/f; 
d.D = d.D./2./pi./f;
%% Forming A1
A1 = sparse(1:7:7*GenNumber,2:7:7*GenNumber,ones(GenNumber,1),7*GenNumber,7*GenNumber);
A1 = A1 + sparse(2:7:7*GenNumber,2:7:7*GenNumber, -d.D./M,7*GenNumber,7*GenNumber);
A1 = A1 + sparse(2:7:7*GenNumber,3:7:7*GenNumber,-Iq./M,7*GenNumber,7*GenNumber); 
A1 = A1 + sparse(2:7:7*GenNumber,4:7:7*GenNumber,-Id./M,7*GenNumber,7*GenNumber); 
A1 = A1 + sparse(3:7:7*GenNumber,3:7:7*GenNumber, -1./d.Tdop,7*GenNumber,7*GenNumber);
A1 = A1 + sparse(3:7:7*GenNumber,5:7:7*GenNumber, 1./d.Tdop,7*GenNumber,7*GenNumber);
A1 = A1 + sparse(4:7:7*GenNumber,4:7:7*GenNumber, -1./d.Tqop,7*GenNumber,7*GenNumber);
A1 = A1 + sparse(5:7:7*GenNumber,5:7:7*GenNumber, -fsi,7*GenNumber,7*GenNumber);
A1 = A1 + sparse(5:7:7*GenNumber,6:7:7*GenNumber,1./d.TE,7*GenNumber,7*GenNumber);
A1 = A1 + sparse(6:7:7*GenNumber,5:7:7*GenNumber, -d.KA.*d.KF./d.TA./d.TF,7*GenNumber,7*GenNumber);
A1 = A1 + sparse(6:7:7*GenNumber,6:7:7*GenNumber, -1./d.TA,7*GenNumber,7*GenNumber);
A1 = A1 + sparse(6:7:7*GenNumber,7:7:7*GenNumber, d.KA./d.TA,7*GenNumber,7*GenNumber);
A1 = A1 + sparse(7:7:7*GenNumber,5:7:7*GenNumber,d.KF./d.TF.^2,7*GenNumber,7*GenNumber);
A1 = A1 + sparse(7:7:7*GenNumber,7:7:7*GenNumber,-1./d.TF,7*GenNumber,7*GenNumber);
%% Forming B1、B2
B1 = sparse(2:7:7*GenNumber,1:2:2*GenNumber,(Iq.*(d.Xdp - d.Xqp)-Edp)./M,7*GenNumber,2*GenNumber);
B1 = B1 + sparse(2:7:7*GenNumber,2:2:2*GenNumber,(Id.*(d.Xdp - d.Xqp)-Eqp)./M,7*GenNumber,2*GenNumber);
B1 = B1 + sparse(3:7:7*GenNumber,1:2:2*GenNumber,-(d.Xd - d.Xdp)./d.Tdop,7*GenNumber,2*GenNumber);
B1 = B1 + sparse(4:7:7*GenNumber,2:2:2*GenNumber,(d.Xq - d.Xqp)./d.Tqop,7*GenNumber,2*GenNumber);

B2 = sparse(6:7:7*GenNumber,1:GenNumber,-d.KA./d.TA,7*GenNumber,2*GenNumber);
%% Forming C1
cos_delta_theta = cos(delta-theta(GenNode));
sin_delta_theta = sin(delta-theta(GenNode));
C1 = sparse(1:2:2*GenNumber,1:7:7*GenNumber,-abs(Etp).*cos_delta_theta,2*GenNumber,7*GenNumber);
C1 = C1 + sparse(2:2:2*GenNumber,1:7:7*GenNumber,abs(Etp).*sin_delta_theta,2*GenNumber,7*GenNumber);
C1 = C1 + sparse(1:2:2*GenNumber,4:7:7*GenNumber,ones(GenNumber,1),2*GenNumber,7*GenNumber);
C1 = C1 + sparse(2:2:2*GenNumber,3:7:7*GenNumber,ones(GenNumber,1),2*GenNumber,7*GenNumber);
%% Forming D1
D1 = sparse(1:2:2*GenNumber,1:2:2*GenNumber,- d.Rs,2*GenNumber,2*GenNumber);
D1 = D1 + sparse(1:2:2*GenNumber,2:2:2*GenNumber,d.Xqp,2*GenNumber,2*GenNumber);
D1 = D1 + sparse(2:2:2*GenNumber,1:2:2*GenNumber,- d.Xdp,2*GenNumber,2*GenNumber);
D1 = D1 + sparse(2:2:2*GenNumber,2:2:2*GenNumber,- d.Rs,2*GenNumber,2*GenNumber);
%% Forming D2
D2 = sparse(1:2:2*GenNumber,1:GenNumber,-sin_delta_theta,2*GenNumber,2*GenNumber);
D2 = D2 + sparse(1:2:2*GenNumber,1+GenNumber:GenNumber*2,Et.*cos_delta_theta,2*GenNumber,2*GenNumber);
D2 = D2 + sparse(2:2:2*GenNumber,1:GenNumber,-cos_delta_theta ,2*GenNumber,2*GenNumber);                             
D2 = D2 + sparse(2:2:2*GenNumber,1+GenNumber:GenNumber*2,  -Et.*sin_delta_theta,2*GenNumber,2*GenNumber);

%% Forming C2
delta_ind = 1:7:((GenNumber-1)*7+1); 
real_ind2 = 1:GenNumber; 
imag_ind2 = real_ind2+GenNumber; 
Vcos = V(GenNode).*cos(delta -theta(GenNode));
Vsin = V(GenNode).*sin(delta -theta(GenNode));
C2 = sparse(real_ind2,delta_ind,Id.*Vcos - Iq.*Vsin,2*GenNumber,7*GenNumber); 
C2 = C2 + sparse(imag_ind2,delta_ind,-Id.*Vsin - Iq.*Vcos,2*GenNumber,7*GenNumber);
%% Forming D3
imag_ind = 2:2:2*GenNumber; 
real_ind = 1:2:2*GenNumber; 
D3 = sparse(real_ind2,real_ind,Vsin,2*GenNumber,2*GenNumber);
D3 = D3 + sparse(real_ind2,imag_ind,Vcos,2*GenNumber,2*GenNumber);
D3 = D3 + sparse(imag_ind2,real_ind,Vcos,2*GenNumber,2*GenNumber);
D3 = D3 + sparse(imag_ind2,imag_ind,-Vsin,2*GenNumber,2*GenNumber);
%% Forming D4、D5、D6、D7
ejTheta = exp(1i*theta); 
ejDelta = sparse(Node,1); 
ejDelta(GenNode) = exp(1i*(delta-pi/2)); 
Idelta = conj(Ig.*ejDelta);                                                      
DiagPcos = sparse(diag(P));
DiagPsin = sparse(diag(Q));
DiagVp = sparse(diag(U));                                                    
HIJ = -DiagVp*conj(Ybus*DiagVp);  
ImagHIJ = imag(HIJ);
RealHIJ = real(HIJ);
DiagV = sparse(diag(V));
Pa = ImagHIJ + DiagPsin;                                                   %Jacobian element H
Pv = (RealHIJ - DiagPcos)/DiagV;                                           %Jacobian element N
Qa = -RealHIJ - DiagPcos;                                                  %Jacobian element J
Qv = (ImagHIJ - DiagPsin)/DiagV;                                           %Jacobian element L
f = sparse(1:Node,1:Node,ejTheta.*Idelta,Node,Node);
g = sparse(1:Node,1:Node,1i*U.*Idelta,Node,Node);
D4_gen = [real(f(GenNode,GenNode)),real(g(GenNode,GenNode));imag(f(GenNode,GenNode)),imag(g(GenNode,GenNode))];
D4 = [Pv(GenNode,GenNode),Pa(GenNode,GenNode);Qv(GenNode,GenNode),Qa(GenNode,GenNode)] + D4_gen;
D5 = [Pv(GenNode,nongen_ind),Pa(GenNode,nongen_ind);Qv(GenNode,nongen_ind),Qa(GenNode,nongen_ind)];
D6 = [Pv(nongen_ind,GenNode),Pa(nongen_ind,GenNode);Qv(nongen_ind,GenNode),Qa(nongen_ind,GenNode)];
D7 = [Pv(nongen_ind,nongen_ind),Pa(nongen_ind,nongen_ind);Qv(nongen_ind,nongen_ind),Qa(nongen_ind,nongen_ind)];

%% Forming state matrix
[rowB2,~] = size(B2); 
[rowD2,~] = size(D2); 
[~,columnD5] = size(D5); 
[rowD6,~] = size(D6); 
[~,columnC2] = size(C2); 
[~,columnD3] = size(D3); 

zero1 = sparse(rowB2,columnD5);
zero2 = sparse(rowD2,columnD5); 
zero3 = sparse(rowD6,columnC2); 
zero4 = sparse(rowD6,columnD3); 

A = A1;
B = [B1,B2,zero1];
C = [C1;C2;zero3];
D = [D1,D2,zero2;D3,D4,D5;zero4,D6,D7];

Asys = A - B/D*C;

%% right eigenvector
[phi,Lambda] = eig(full(Asys));
lambda = diag(Lambda); 
lambda2 = real(lambda); 
[~,Ind] = sort(lambda2,'descend'); 
phi = phi(:,Ind);
%% left eigenvector
[psi,Lambda] = eig(full(Asys.'));
lambda = diag(Lambda); 
lambda2 = real(lambda); 
[~,Ind] = sort(lambda2,'descend'); 
lambda = lambda(Ind);
psi = psi(:,Ind);
%% Obtain spectral abscissa and the eigenvector of it 
nonzero_idx = find(d.D ~= 0, 1);
if abs(real(lambda(1)))>= 1e-6
    var.eigen.lambda = lambda(1);
    psi = psi(:,1);
    var.eigen.psi = psi.'; 
    var.eigen.phi = phi(:,1);
else
        if ~isempty(nonzero_idx)
            var.eigen.lambda = lambda(2);
            psi = psi(:,2);
            var.eigen.psi = psi.'; 
            var.eigen.phi = phi(:,2);
        else
            var.eigen.lambda = lambda(3);
            psi = psi(:,3);
            var.eigen.psi = psi.'; 
            var.eigen.phi = phi(:,3);
        end
end
%% Verification of spectral abscissa
RetValue = Asys*var.eigen.phi - var.eigen.lambda*var.eigen.phi;
RetValue1 = var.eigen.psi*Asys - var.eigen.lambda*var.eigen.psi;
if abs(RetValue)>=1.0e-08 
    msgbox('ERROR','IMFORMATION')
elseif abs(RetValue1)>=1.0e-8
    msgbox('ERROR','IMFORMATION')
end
%% Normalization processing of left/right eigenvector
addvar = sqrt(var.eigen.psi*var.eigen.phi);
var.eigen.psi = var.eigen.psi/addvar;
var.eigen.phi = var.eigen.phi/addvar;
%% Return parameter
data.ngen = d.PNum; 
data.nbuses = d.NodeNum;
data.gen_bus = d.PIndex ; 
data.nongen_ind = nongen_ind;
data.Pgenbuses = d.PIndex;
data.Ybus = d.Ybus;
data.machine.H = d.H;
data.machine.Xdp = d.Xdp;
data.machine.Xqp = d.Xqp;
data.machine.TE = d.TE;
data.machine.Ae = d.Ae;
data.machine.Be = d.Be;
data.fr = d.fr;
data.to = d.to;
data.IndexI = d.IndexI;
data.IndexJ = d.IndexJ;
data.ImpAngle = d.ImpAngle;
data.Yabs = d.Yabs;
data.M = M;

var.theta = theta;
var.delta = delta;
var.Id = Id;
var.Iq = Iq;
var.V = V;

var.Efd = Efd;

var.eigen.D1 = D1;
var.eigen.D2 = D2;
var.eigen.D3 = D3;
var.eigen.D4 = D4;
var.eigen.D5 = D5;
var.eigen.D6 = D6;
var.eigen.D7 = D7;
var.eigen.C1 = C1;
var.eigen.C2 = C2;
var.eigen.B1 = B1;
var.eigen.B2 = B2;

var.eigen.A = A;
var.eigen.B = B;
var.eigen.C = C;
var.eigen.D = D;

var.eigen.zero1 = zero1;
var.eigen.zero2 = zero2;
var.eigen.zero3 = zero3;
var.eigen.zero4 = zero4;
%% Closed-form sensitivity
%-------------------------------Hsp,Hsq------------------------------------
Hsp_v = -real(f(GenNode,GenNode));
Hsp_theta = -real(g(GenNode,GenNode));
Hsq_v = -imag(f(GenNode,GenNode));
Hsq_theta = -imag(g(GenNode,GenNode));
Hsp_delta = -Hsp_theta;
Hsq_delta = -Hsq_theta;
Hsp_Id = -diag(Vsin);
Hsp_Iq = -diag(Vcos);
Hsq_Id = -diag(Vcos);
Hsq_Iq = diag(Vsin);
%-------------------------------Hed/Heq------------------------------------
Hed_v = diag(-sin_delta_theta);
Heq_v = diag(-cos_delta_theta);
Hed_theta = diag(cos_delta_theta);
Heq_theta = diag(-sin_delta_theta);
Hed_delta = diag(-cos_delta_theta);
Heq_delta = diag(sin_delta_theta);
Hed_Id = diag(-d.Rs);
Heq_Id = diag(-d.Xdp);
Hed_Iq = diag(d.Xqp);
Heq_Iq = diag(-d.Rs);
Hed_Edp = diag(ones(GenNumber,1));
Heq_Eqp = diag(ones(GenNumber,1));
%------------------------------Hiq/Hfd-------------------------------------
Hiq_Iq = -diag(d.Xq - d.Xqp);
Hfd_Id = -diag(d.Xd - d.Xdp);
Hiq_Edp = diag(ones(GenNumber,1));
Hfd_Eqp = -diag(ones(GenNumber,1));
Hfd_Efd = diag(ones(GenNumber,1));
%------------------------------A、B-------------------------------------------
data.a11 = Hsp_v;
data.a12 = Hsp_theta;
data.a13 = Hsp_delta;
data.a14 = Hsp_Id;
data.a15 = Hsp_Iq;

data.a21 = Hsq_v;
data.a22 = Hsq_theta;
data.a23 = Hsq_delta;
data.a24 = Hsq_Id;
data.a25 = Hsq_Iq;

data.a31 = Hed_v;
data.a32 = Hed_theta;
data.a33 = Hed_delta;
data.a34 = Hed_Id;
data.a35 = Hed_Iq;
data.a36 = Hed_Edp;

data.a41 = Heq_v;
data.a42 = Heq_theta;
data.a43 = Heq_delta;
data.a44 = Heq_Id;
data.a45 = Heq_Iq;
data.a47 = Heq_Eqp;

data.a55 = Hiq_Iq;
data.a56 = Hiq_Edp;

data.a64 = Hfd_Id;
data.a67 = Hfd_Eqp;
data.a68 = Hfd_Efd;

data.b1 = -ones(GenNumber,1);
data.b2 = sparse(GenNumber,1);
data.b3 = sparse(GenNumber,1);
data.b4 = sparse(GenNumber,1);
data.b5 = sparse(GenNumber,1);
data.b6 = sparse(GenNumber,1);

data.b01 = sparse(GenNumber,1);
data.b02 = -ones(GenNumber,1);
data.b03 = sparse(GenNumber,1);
data.b04 = sparse(GenNumber,1);
data.b05 = sparse(GenNumber,1);
data.b06 = sparse(GenNumber,1);

 end

