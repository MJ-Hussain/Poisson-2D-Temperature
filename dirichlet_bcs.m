%dirichlet Boundary Function
function phi_b=dirichlet_bcs(x,y,L,H)
if x==0 && y~=0
    phi_b=1;
elseif x==L && y~=0
    phi_b=y;
elseif x~=0 && y==0
    phi_b=-2;
elseif x~=0 && y==H
    phi_b=x^2; 
else
    phi_b=0;
end
end