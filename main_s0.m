%main code 
%sample code that aims to produce model results shown in Shelef and Goren 2021.


%% prepare model inputs
%%upload template with 1 (junction) and zero all other nodes
close all, clear all, clc

load 'model_data'%this loads initial conditions elevation vector (zv0), distance vector (xv), and local area vector (av)
dx=xv(2)-xv(1);%distance between nodes

%set model parameters
m=0.45;%area exponent
n=1;%slope exponent
K=1e-5;%erodibility coeff [m^(1-2m) yr^-1]
U=1e-3;% uplift [m/yr]
D=0.24;% diffusion coeff [m^2/yr]



%plot initial conditions
close all
plot(xv,zv0,'-k');
hold on
f=find(lav==max(lav));
plot(xv(f), zv0(f), 'ok', 'markerfacecolor', [1 1 1]*0.5)
xlabel('Distance [m]')
ylabel('Elevation [m]')
shg


%constrain time step
dtD= dx^2/D;
av=cumsum(lav(2:end));
dtF= dx/(K*max(av).^m);%this is for n=1
dt=100;
if dt>dtF/2 | dt>dtD/2
    error('dt too large')
end

%% run a fixed confluence model to steady state
[tvo1, divvo1, zo1, Ao1]=runFixedConflluenceLEM2SS_s(zv0, xv, lav, U, K, m, n, D, dx, dt);

%plot simulated topography
close all
plot(xv,zo1,'-k');
hold on
plot(xv(f), zo1(f), 'ok', 'markerfacecolor', [1 1 1]*0.5)
xlabel('Distance [m]')
ylabel('Elevation [m]')
shg

%% run avulsion model to quasi-steady state
tend=dt*round(max(zv0)/U*10/dt);%approximated time constraint to quasi steady state
Ndt=1;%number of time steps between avulsions
[tvo2, divvo2, zo2, Ao2]=runAvulsionLEM_s(zv0, xv, lav, U, K, m, n, D, dx, tend, dt, Ndt);

%% plot simulated topography 
close all
plot(xv,zo2,'--k');
hold on
plot(xv(f), zo2(f), 'ok', 'markerfacecolor', [1 1 1]*0.5)
xlabel('Distance [m]')
ylabel('Elevation [m]')
shg
%% plot divide location vs time for the two simulations above
close all
plot(tvo1, divvo1, '-k'),shg
hold on
plot(tvo2, divvo2, '--k'),shg
xlabel('Time [yr]')
ylabel('Divide Location [m]')
set(gca, 'xlim', [0 1.1e7])



