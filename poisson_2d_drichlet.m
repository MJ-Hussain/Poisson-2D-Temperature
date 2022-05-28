%poisson 2d with dirichlet boundary
clc; clear;

%                  Input                   %

L=1;        %Lenght in x
H=1;        %Height in y
num_pts=50; %number of mesh points in x and y

%              Source Functions            %

%src_func=@(x,y) 0;                      %Zeros Source
src_func=@(x,y) 1;                       %Uniform Source
%src_func=@(x,y) x./L;                   %Linear Source in x
%src_func=@(x,y) x./L+2*y./H;            %Linear Source in both x and y
%src_func=@(x,y) (x./L).^2+sin(10*y./H); %Non-Linear source in both x and y

%        Calling Fucntion for solution      %

[PHI,S,X,Y, soln_time,mem_A, boundary_pts]=poisson(src_func,L,H,num_pts);

%                Plots                  %

%Plot of solution 
figure(1)
contourf(X,Y,PHI,100,'LineStyle','none')
shading interp
colorbar
title('Plot of Solution')
xlabel('x'); ylabel('y');

%Plots of source
figure(2)
imagesc(X,Y,S');
set(gca,'YDir','normal')
colorbar
title('Plot of Source')
xlabel('x'); ylabel('y');

%Plot of Dirichlet Boundary
figure(3)
imagesc(X,Y,boundary_pts');
set(gca,'YDir','normal')
colorbar
title('Plot of Boundary')
xlabel('x'); ylabel('y');

%                      Functions                       %

%Dirichlet Boundary function 
% function phi_b=dirichlet_bcs(x,y,L,H)
% if x==0 && y~=0
%    phi_b=1;
%    
% elseif x==L && y~=0
%    phi_b=y;
%    
% elseif x~=0 && y==0
%    phi_b=-2;
%     
% elseif x~=0 && y==H
%     phi_b=x^2;
%    
% else
%     phi_b=0;
% end
% 
% end
function [phi,Z,x,y, sol_time, inf_A, BD_pts]=poisson(src_func,L,H,num_pts)
tic;      %Time start

n=num_pts+1; nn=n*n;
h=L/num_pts;h2=h^2;k=H/num_pts;k2=k^2;
x=(0:num_pts)*h; % mesh values in x
y=(0:num_pts)*k; % mesh values in y

%Coefficient matrix A and LHS matrix b for system Ax=b
A=zeros(nn,nn);b=zeros(nn,1);

%Boundary point matrix
BD_pts=zeros(num_pts);
%Loop in interior points
for i=2:n-1 
    for j=2:n-1
        A(i+(j-1)*n,i-1+(j-1)*n)=1/h2;
        A(i+(j-1)*n,i+1+(j-1)*n)=1/h2;
        A(i+(j-1)*n,i+(j-1)*n)=-2/h2-2/k2;
        A(i+(j-1)*n,i+(j-2)*n)=1/k2;
        A(i+(j-1)*n,i+j*n)=1/k2;
        b(i+(j-1)*n)=-src_func(x(i),y(j));
    end
end

% bottom and top boundary points
for i=1:n 
    j=1;
    A(i+(j-1)*n,i+(j-1)*n)=1;
    b(i+(j-1)*n)=dirichlet_bcs(x(i),y(j),L,H);
    BD_pts(i,j)=b(i+(j-1)*n);
    j=n;
    A(i+(j-1)*n,i+(j-1)*n)=1;
    b(i+(j-1)*n)=dirichlet_bcs(x(i),y(j),L,H);
    BD_pts(i,j)=b(i+(j-1)*n);
end

% left and right boundary points
for j=2:n-1 
    i=1;
    A(i+(j-1)*n,i+(j-1)*n)=1;
    b(i+(j-1)*n)=dirichlet_bcs(x(i),y(j),L,H);
    BD_pts(i,j)=b(i+(j-1)*n);
    i=n;
    A(i+(j-1)*n,i+(j-1)*n)=1;
    b(i+(j-1)*n)=dirichlet_bcs(x(i),y(j),L,H);
    BD_pts(i,j)=b(i+(j-1)*n);
end

% solve of the system Ax=b
v=A\b; 
%reshaping from labeling
phi=reshape(v(1:nn),n,n)';


%Time taken
sol_time=toc;

%Information of matrix A
inf_A=whos('A');

%Source matrix
Z=zeros(length(x),length(y));
for i=1:length(x)
    for j=1:length(y)
        Z(i,j)=src_func(x(i),y(i));
    end
end

end