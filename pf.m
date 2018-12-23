function [V,theta,d,P,Q,U]= pf(Node,MaxK,IterationPrecision,BalanceNode,PVNode,V,P,Q,theta,d,LineB,LineI,LineJ,LineY,TransformerI,TransformerJ,TransformerK,TransformerY,GroundI,GroundB)
%% Power Flow
[Y,d] = matrixY(d,Node,LineB,LineI,LineJ,LineY,TransformerI,TransformerJ,TransformerK,TransformerY,GroundI,GroundB);
for k = 1:MaxK
     [dp dq U y1 y2] = constant(BalanceNode,Y,V,theta,P,Q,PVNode);
     pq = [dp;dq];
     max_pq = max(abs(pq));
     if max_pq>IterationPrecision
         Jacobian = matrixF(Node,BalanceNode,Y,U,y1,y2,PVNode);
         [V theta] = iteration(Node,V,theta,pq,Jacobian);
     else
         if k==1
             Jacobian = matrixF(Node,BalanceNode,Y,U,y1,y2,PVNode);
         end
         break
     end
end
U = V.*exp(1j.*theta); 
s = U.*conj(Y*U);   
P = real(s);
Q = imag(s);
d.genP(BalanceNode) = P(BalanceNode); 
end

