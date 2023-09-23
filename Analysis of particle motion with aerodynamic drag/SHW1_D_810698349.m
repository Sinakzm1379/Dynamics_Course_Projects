clc
clear all
% request Data
p= input('Please enter your problem ( a or b ) : ','s');
v0= input('Please enter the velocity at t=0 : ');
theta0= input('Please enter the theta0 in degree : ');
g=9.81;
theta0= (theta0*pi)/180;
t_Max=2.*((v0.*sin(theta0))./(g));

if p=='a'
    H_max=((v0.*sin(theta0)).^2)./(2*g);
    t_t =t_Max;
    t=0:0.01:t_t;
    m=length(t);
    v=zeros(m,2);
    x=zeros(m,1);
    y=zeros(m,1);
    
% Calculate the requested items
    for i=1:m
        vx(i)=v0.*cos(theta0);
        vy(i)=v0.*sin(theta0)-g.*t(1,i);
        x(i,1)=v0.*cos(theta0).*t(1,i);
        y(i,1)=(-1/2).*(g).*((t(1,i).^2))+v0.*sin(theta0).*t(1,i);
        R(i,1)=sqrt((x(i,1).^2)+(y(i,1).^2));
        theta(i,1)=atan(y(i,1)./x(i,1));
        thetaV(i,1)=atan(vy(i)./vx(i));
        v_t(i,1)=sqrt(vx(i).^2+vy(i).^2);
        v_r(i,1)=vx(i).*cos(theta(i,1))+vy(i).*sin(theta(i,1));
        v_theta(i,1)=vy(i).*cos(theta(i,1))-vx(i).*sin(theta(i,1));
        a_n(i,1)=(g).*(cos(thetaV(i,1)));
        rho(i,1)=abs((v_t(i,1).^2)./(a_n(i,1)));
    end
    
else if (p=='b')
    k = input('Please enter the K : ');
    setGlobalx(k);
    t=0:0.01:t_Max;
    v0=[v0.*cos(theta0);v0.*sin(theta0)];
% Calculate the velocity differential fuction
    [t,v]=ode45(@velocity,t,v0);
    m=length(t);
    x=zeros(m,1);
    y=zeros(m,1);
    vx(1)=v(1,1);
    vy(1)=v(1,2);
    
% Calculate the position integral fuction
    for j=2:m
       x(j)=trapz(t(1:j,1),v(1:j,1))+x(1,1);
       y(j)=trapz(t(1:j,1),v(1:j,2))+y(1,1);
       vx(j)=v(j,1);
       vy(j)=v(j,2);
    end
    w=1;
% Find the total time
    while y(w)>=0
       w=w+1;
    end
e=w;
% Remove extra items
    while e<=m
       t(w)=[];
       x(w)=[];
       y(w)=[];
       vx(w)=[];
       vy(w)=[];
       e=e+1;
    end
    m=length(t);
    t_t=t(m);
    
% Calculate the requested items
    for i=1:m
        a(i,1)=(-k).*(vx(i)).*(sqrt(vx(i).^2+vy(i).^2));
        a(i,2)=(-k).*(vy(i)).*(sqrt(vx(i).^2+vy(i).^2))-g;
        R(i,1)=sqrt((x(i,1).^2)+(y(i,1).^2));
        theta(i,1)=atan(y(i,1)./x(i,1));
        thetaV(i,1)=atan(vy(i)./vx(i));
        v_t(i,1)=sqrt(vx(i).^2+vy(i).^2);
        v_r(i,1)=vx(i).*cos(theta(i,1))+vy(i).*sin(theta(i,1));
        v_theta(i,1)=vy(i).*cos(theta(i,1))-vx(i).*sin(theta(i,1));
        a_n(i,1)=(g).*(cos(thetaV(i,1)));
        rho(i,1)=abs((v_t(i,1).^2)./(a_n(i,1)));
    end
% clear workspace
    clear e w v j
    end
end

% request the tp
A=['Please enter the time you want ( smaller than ',num2str(t_t),' seconds ) : '];
t_p= input(A);

%  plot
figure(1)
plot(x(:,1),y(:,1));
title('X-Y plot')
xlabel('X')
ylabel('Y')
figure(2)
plot(t,rho(:,1));
title('rho-t plot')
xlabel('Time')
ylabel('Rho')
figure(3)
plot(t,vx,t,vy);
title('V-t plot (x-y)')
xlabel('Time')
ylabel('Velocity')
legend('V-x','V-y')
figure(4)
plot(t,v_r(:,1),t,v_theta(:,1));
title('V-t plot (r-theta)')
xlabel('Time')
ylabel('Velocity')
legend('V-r','V-theta')

% Final requested values 
    n=(t_p./0.01)+1;
    H_max=max(y(:,1));
    X_max=x(m,1);
    V_p=v_t(n,1);
    thetav_p=(thetaV(n,1)).*180./pi;
    X_p=x(n,1);
    Y_p=y(n,1);
    Vr_p=v_r(n,1);
    Vtheta_p=v_theta(n,1);
    rho_p=rho(n,1);
    
% Display numbers
    S(1).name = 'Total Range : ';
    S(1).num = X_max;
    S(2).name = 'Maximum H : ';
    S(2).num = H_max;
    S(3).name = 'Velocity(number) at tp : ';
    S(3).num = V_p;
    S(4).name = 'Velocity(degrees) at tp : ';
    S(4).num = thetav_p;
    S(5).name = 'X at tp : ';
    S(5).num = X_p;
    S(6).name = 'Y at tp : ';
    S(6).num = Y_p;
    S(7).name = 'Vr at tp : ';
    S(7).num = Vr_p;
    S(8).name = 'Vtheta at tp : ';
    S(8).num = Vtheta_p;
    S(9).name = 'Radius of curvature at tp : ';
    S(9).num = rho_p;
    B = [{S.name};{S.num}];
    fprintf('\n')
    fprintf('%10s %.2f\n ',B{:})
    fprintf('\n')

% clear workspace
clear A B i m p t_Max 

% velocity differential fuction
function dt = velocity(t,v)
g=9.81;
k=getGlobalk;
dt=zeros(2,1);
dt(1)=(-k).*(v(1)).*(sqrt(v(1).^2+v(2).^2));
dt(2)=(-k).*(v(2)).*(sqrt(v(1).^2+v(2).^2))-g;
end
% make k global for using in velocity function
function setGlobalx(val)
global ck
ck = val;
end
function r = getGlobalk
global ck
r = ck;
end