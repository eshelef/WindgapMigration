function [tvo, divvo, zo, Ao]=runFixedConflluenceLEM2SS_s(z, x, lav, U, K, m, n, D, dx, dt)

%Simulate landscape evolution to steady state using linear diffusion 
%and stream power equations (e.g., Tucker and Hancock, 2010). The code is 
%designed specifically to produce model results shown in Shelef and Goren 2021. 

%Inputs:
%z-vector of initial elevation  [L], x-vector of distance from boundary [L],
%lav - vector of local area at a node , U-uplift [L/t], K- erodibility
%coefficient [L^1-2m)t^-1], m,n - exponents, D - diffusion coefficient [L^2/t],
%dt - time step [t]

%Outputs:
%tvo-vector of time steps, %divv0-vector of divide location index at each time step,
%zo-vector of output elevations, Ao-vector of output drainage area at each
%node.

%Authors: Eitan Shelef and Liran Goren, 2021.



%prepare output vectors
tv=zeros(ceil(max(z)/U*50/dt), 1)*NaN;%coarse approximation for length of output vectors
divv=tv;
av=tv;
N=length(x);

%set parameters
counter=1;
e_diff_counter=0;
p=1.1;%exponent for area fractioning(Freeman, 1991)
t=0;

while true
    
    
    divide_id = find(z == max(z),1);%get divide location
    
    
    if divide_id>=length(x)%stop if divide at boundary
        warning('divide at boundary')
        break
    end
    
    
    %partition drainage area at divide following Freeman 1991
    S_div_p=(z(divide_id)-z(divide_id+1))/dx;
    S_div_m=(z(divide_id)-z(divide_id-1))/dx;
    frac_p=S_div_p^p/(S_div_p^p+S_div_m^p);
    frac_m=1-frac_p;
    
    A = zeros(N,1)*NaN;
    A(divide_id) = lav(divide_id);
    Ap=A(divide_id)*frac_p;
    Am=A(divide_id)*frac_m;
    A(divide_id+1) = lav(divide_id+1)+Ap;
    A(divide_id-1) = lav(divide_id-1)+Am;
    
    
    for i = divide_id-2:-1:1%set A left of divide
        A(i) = A(i+1) + lav(i);
    end
    for i = divide_id+2:1:N %set A right of divide
        A(i) = A(i-1) + lav(i);
    end
    
    %record for outoputs
    tv(counter)=t;%run time
    divv(counter)=x(divide_id);%divide location
    av(counter)=A(divide_id);%drainage area
    counter=counter+1;
    
    %solve using RK45
    dzdt=getdzdt_s(z,A,U,K,m,n,D,divide_id,dx);
    
    %step 1:
    ztemp=z + dzdt*(dt/3);
    k2=getdzdt_s(ztemp,A,U,K,m,n,D,divide_id,dx);
    
    %step 2:
    ztemp=z + (dzdt+k2)*dt/6;
    k3=getdzdt_s(ztemp,A,U,K,m,n,D,divide_id,dx);
    
    
    %step 3:
    ztemp=z + (dzdt+k3*3.0)*dt/8;
    k4=getdzdt_s(ztemp,A,U,K,m,n,D,divide_id,dx);
    
    %step 4:
    ztemp=z + ( dzdt-k3*3.0+k4*4.0)*dt/2;
    k5=getdzdt_s(ztemp,A,U,K,m,n,D,divide_id,dx);
    
    %step 5:
    z_new=z + (dzdt+k4*4.0+k5)*dt/6;
    
    
    %stop after N times of 0 erosion, should be he same as N=1 and added to eliminate numerical effects
    sdzdt=sum(abs(z_new-z));
    if(sdzdt==0)
        e_diff_counter=e_diff_counter+1;
        if e_diff_counter>=100
            break
        end
    else
        e_diff_counter=0;
    end
    
    
    %outout every 1000 dts to track progress
    if mod(t,1000*dt)== 0
        disp('divide index, sum dzdt:')
        disp(divide_id)
        disp(sdzdt)
    end
    
    t = t+dt;
    z = z_new;
    
    
end

%set outputs
zo=z;
Ao=A;
tvo=tv(1:counter-1);
divvo=divv(1:counter-1);
disp('done');
end