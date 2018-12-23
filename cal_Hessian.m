function[Pvv,Qvv,Pva,Qva,Pav,Qav,Paa,Qaa] = cal_Hessian(var,data,ind,Yp,Yq)
V = var.V;
DiagYp = sparse(ind,ind,Yp,data.nbuses,data.nbuses); 
DiagYq = sparse(ind,ind,Yq,data.nbuses,data.nbuses);  
Yp = diag(DiagYp);
Yq = diag(DiagYq);
non_ind =1:data.nbuses;
non_ind(ind) =[];
Ybus = data.Ybus;
Ybus(non_ind,:) = sparse(length(non_ind),data.nbuses);
[IndexI,IndexJ,Yval] = find(Ybus);
ImpAngle = angle(Yval);
Yabs = abs(Yval);
EleAngle2 = var.theta(IndexI) - var.theta(IndexJ) - ImpAngle;
TempCos = Yabs.*cos(EleAngle2);      
TempSin = Yabs.*sin(EleAngle2);       
Ycos = sparse(IndexI,IndexJ,TempCos,data.nbuses,data.nbuses);              %Y*Cos
Ysin = sparse(IndexI,IndexJ,TempSin,data.nbuses,data.nbuses);              %Y*Sin 
%------------------Pvv,Qvv ----------------------------------------- 
Pvv = Ycos.'*DiagYp + DiagYp*Ycos;                                         %[YcosT]*diag(Yq) + diag(Yq)*[Ycos]
Qvv = Ysin.'*DiagYq + DiagYq*Ysin;                                         %[YsinT]*diag(Yp) + diag(Yp)*[Ysin]                        
%------------------Pva,Qva -----------------------------------------    
Pva = (-diag(Ysin*V) + diag(V)*Ysin.')*DiagYp - DiagYp*diag(V)*Ysin...
    + diag(Ysin.'*diag(V)*Yp);                                             %(-diag(Ysin*V) + diag(V)*[YsinT])*diag(Yp) - diag(Yp)*diag(V)*[Ysin] + diag([YsinT]*diag(V)*Yp)
Qva = (diag(Ycos*V) - diag(V)*Ycos.')*DiagYq + DiagYq*diag(V)*Ycos...
    - diag(Ycos.'*diag(V)*Yq);                                             %(diag(Ycos*V) - diag(V)*[YcosT])*diag(Yq) + diag(Yq)*[Ycos] - diag([YcosT]*diag(V)*Yq)
%------------------Pav,Qav -----------------------------------------  
Pav = Pva.';
Qav = Qva.';
%------------------Paa,Qaa -----------------------------------------  
Paa = (-diag(Ycos.'*diag(V)*Yp) + diag(V)*DiagYp*Ycos)*diag(V) - ...
    (diag(Ycos*V) - diag(V)*Ycos.')*diag(V)*DiagYp;                        %(-diag([YcosT]*diag(V)*Yp) + diag(V)*diag(Yp)*[Ycos])*diag(V) - (diag(Ycos*V) - diag(V)*[YcosT])*diag(V)*diag(Yp)
Qaa = (-diag(Ysin*V) + diag(V)*Ysin.')*diag(V)*DiagYq ...
    - (-diag(V)*DiagYq*Ysin + diag(Ysin.'*diag(V)*Yq))*diag(V);            %(-diag([Ysin]*V) + diag(V)*[Ysin.'])*diag(V)*diag(Yq) - (-diag(V)*diag(Yq)*[Ysin] + diag([Ysin.']*diag(V)*Yq))*diag(V)
end