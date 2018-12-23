function [dp dq U y1 y2]=constant(BalanceNode,Y,u,an,p,q,PVNode)
U = u.*exp(1i.*an);                                                        %polar form
s = U.*conj(Y*U);                                                          
y1 = real(s);                                                                   
y2 = imag(s);                                                              
dp = p-y1;                                                                 %active power imbalance equation_dP
dq = q-y2;                                                                 %reactive power imbalance equation_dQ
dp(BalanceNode) = 0;
dq(BalanceNode) = 0;
dq(PVNode) = 0;
end

