function Jacobian = matrixF(Node,BalanceNode,Y,U,y1,y2,PVNode)
dU = sparse(diag(U));                                                      
HIJ = -dU*conj(Y*dU);  
ImagHIJ = imag(HIJ);
RealHIJ = real(HIJ);
y4 = sparse(diag(y1));
y5 = sparse(diag(y2));
H = ImagHIJ + y5;
N = RealHIJ - y4;
j = -RealHIJ - y4;
L = ImagHIJ - y5;
%% PV node equation correction
L(:,PVNode) = 0;
L(PVNode,:) = 0;
j(PVNode,:) = 0;
N(:,PVNode) = 0;
L = L + sparse(PVNode,PVNode,1,Node,Node);
%% Slack node equation correction
H(:,BalanceNode) = 0;
H(BalanceNode,:) = 0;
N(:,BalanceNode) = 0;
N(BalanceNode,:) = 0;
j(:,BalanceNode) = 0;
j(BalanceNode,:) = 0;
L(:,BalanceNode) = 0;
L(BalanceNode,:) = 0;
L = L + sparse(BalanceNode,BalanceNode,1,Node,Node);
H = H + sparse(BalanceNode,BalanceNode,1,Node,Node);
Jacobian = [H N;j L];
end
