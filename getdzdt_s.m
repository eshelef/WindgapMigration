function dzdt=getdzdt_s(z, A, U, K, m, n, D, divide_id, dx)

%Computes the rate of elevation change using linear diffusion and stream power equations (e.g., Tucker and Hancock, 2010) to produce model results shown in Shelef and Goren 2021.

%Inputs:
%z-vector of elevaion values in each node [L], x-vector of distance from boundary for each node [L],
%A-vector drainage area for each node, U-uplift rate [L/t], K-erodibility
%coefficient [L^1-2m)t^-1], m,n - exponents, D - diffusion coefficient [L^2/t],
%divide_id-index of divide location, dx-node spacing

%Output:
%dzdt-vector of rate of change in elevatio for each node

%Authors: Eitan Shelef and Liran Goren, 2021.

s1=[0;z(2:divide_id)-z(1:divide_id-1)]/dx;%left of divide
s2=[z(divide_id:end-1)-z(divide_id+1:end);0]/dx;%right of divide
sd=max(abs([s1(end), s2(1)]));%steepest descent at divide
s=abs([s1(1:end-1);sd(1);s2(2:end)]);
dzdtf=K.*A.^m.*s.^n;

dzdth=D/dx^2 * (z(3:end)-2*z(2:end-1)+z(1:end-2));
dzdth=[0;dzdth;0];

dzdt=U-(dzdtf-dzdth);
dzdt([1,end])=0;%maintain constant elevation at boundaries