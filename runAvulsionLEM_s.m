function [tvo, divvo, zo, Ao]=runAvulsionLEM_s(z, x, lav, U, K, m, n, D, dx, tend, dt, Ndt)

%Simulates landscape evolution with avulsions using linear diffusion and 
%stream power equations (e.g., Tucker and Hancock, 2010). The code is designed specifically to 
%produce model results shown in Shelef and Goren 2021. 

%Inputs:
%z- vector of initial elevation  [L], x- vector of distance from boundary [L],
%lav- vector of local area at a node , U- uplift [L/t], K- erodibility
%coefficient [L^1-2m)t^-1], m,n- exponents, D- diffusion coefficient [L^2/t],
%tend- simulation duration, dt- time step [t], Ndt- number of time steps
%between avulsions

%Outputs:
%tvo- vector of time steps, divv0- vector of divide location index at each time step,
%zo- vector of output elevations, Ao- vector of output drainage area at each
%node.

%Authors: Eitan Shelef and Liran Goren, 2021.



%output vectors
tv=zeros(ceil(tend/dt),1)*NaN;
divv=tv;


%parameters
p=1.1;%exponent for area fractioning(Freeman, 1991)
N=length(x);
f0=find(lav==max(lav));%confluences index
avulsion_dist=mean(diff(f0));%node distance between confluences
lav0=lav;%local area
counter=1;
t=0;
tlast_avulse=t;

while t < tend
    
    
    divide_id = find(z == max(z),1);
    
    if t-tlast_avulse>=Ndt*dt%set avulsions every ~Ndt
        tlast_avulse=t;
        
        rand_av=round( (rand(length(f0),1)-0.5)*avulsion_dist);
        newidx=f0+rand_av;%avulsion along all channel
        
        lav=lav0;        
        newidx(newidx<=1 | newidx>=N)= f0(newidx<=1 | newidx>=N);%eliminate exceptions
        
        %change confluence location (avulsions)
        lav(newidx)=lav0(f0);
        lav(f0)=lav0(newidx);
        
    end
    
    
    
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
    
    %record for outputs
    tv(counter)=t;%run time
    divv(counter)=x(divide_id);%divide location
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
    
    
 %outout every 1000 dts to track progress 
    if mod(t,1000*dt)== 0
        disp('divide index:')
        disp(divide_id)
        disp(t)
    end
    
    t = t+dt;
    z = z_new;
    
    
end


zo=z;
Ao=A;
tvo=tv(1:counter-1);
divvo=divv(1:counter-1);
disp('done');
end