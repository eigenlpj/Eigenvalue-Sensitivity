function [u an]=iteration(n,u,an,S,F)
an_u=S;
J=-F;
an_u = J\S;
u=u+an_u(1+n:2*n).*u;   %voltage amplitude
an=an+an_u(1:n);        %voltage phase variable
end

