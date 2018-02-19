# velocity-between-plates
Determine the flow between plates while changing nodes, pressure gradiant, and also determine reynolds number from the velocity.

%Cory Cresswell
%1/30/18
%ME3156-308

clear
clc

h = 1;
n = 100; %number of nodes
dy = h/n;
dpdx = [-10 -5 0 5 10];
mu=8.9e-3;
u=zeros(n+1,length(dpdx));

%BC
u(1)=0;
u(n+1,:)=0;

% pressure gradiant

a=ones(n,5);
b=ones(n,5).*-2;
c=ones(n,5);
d= ones(n,5).*dy.^2.*dpdx./mu;

%CFD online Thomas Algorithm
for k = 1:5
    d(n,k)=d(n,k)-u(n+1,k);
%forward sub
     for i = 2:n
         m = (a(i,k))./(b(i - 1,k));
         b(i,k) = b(i,k) - m.*c(i - 1,k);
         d(i,k) = d(i,k) - m.*d(i - 1,k);
     end 

%Backward substitution phase

u(n,k) = d(n,k)./ b(n,k);
for j = n-1:-1:2
    u(j,k) = (d(j,k) - c(j,k) .*u(j + 1,k))./(b(j,k)) ;
end 
end

plot(u,0:n)
xlabel('velocity')
ylabel('nodes')
legend('dpdx=-10','dpdx=-5','dpdx=0','dpdx=5','dpdx=10')

%renolds number
p=1000;
l=0.001;
re=p.*u*l/mu;
