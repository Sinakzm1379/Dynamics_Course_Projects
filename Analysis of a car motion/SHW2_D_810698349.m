clc
clear all
%Default
g=9.81;
syms F1 F2 N1 N2 Ax Ay Bx By alpha
% type of application
p1= input('Do you want the plots? ( Y/N ) : ','s');

if p1=='N' || p1=='n'
% request Data
p2= input('Do the wheels slip? ( Y/N ) : ','s');
m= input('Please enter the m (for car) : ');
mw= input('Please enter the m (for each wheel) : ');
a1= input('Please enter the a1 (in meter) : ');
a2= input('Please enter the a2 (in meter) : ');
h= input('Please enter the h (in meter) : ');
r= input('Please enter the radius of each wheel (in meter) : ');
T= input('Please enter the couple on the rear wheels (in N.m) : ');
phi= input('Please enter the \Phi (in degree) : ');
phi=(phi*pi)/180;
% Calculate the requested items
eq1= 2*Ax+2*F1-2*mw*g*sin(phi)-2*mw*r*alpha==0;
eq2= 2*Ay+2*N1-2*mw*g*cos(phi)==0;
eq3=2*T-2*F1*r-mw*(r^2)*alpha==0;
eq4=2*Bx-2*F2-2*mw*g*sin(phi)-2*mw*r*alpha==0;
eq5= 2*By+2*N2-2*mw*g*cos(phi)==0;
eq6=2*F2*r-mw*(r^2)*alpha==0;
eq7=-2*Ax-2*Bx-m*g*sin(phi)-m*r*alpha==0;
eq8=-2*Ay-2*By-m*g*cos(phi)==0;
eq9=-2*T+(2*Ax+2*Bx)*(h-r)+a1*2*By-a2*2*Ay==0;
[A,B] = equationsToMatrix([eq1,eq2,eq3,eq4,eq5,eq6,eq7,eq8,eq9],[F1,F2,N1,N2,Ax,Ay,Bx,By,alpha]);
Y = linsolve(A,B);

if p2=='N' || p2=='n'
    c=0;
    X=double(Y);
    F1=X(1);
    F2=X(2);
    N1=X(3);
    N2=X(4);
    Ax=X(5);
    Ay=X(6);
    Bx=X(7);
    By=X(8);
    alpha=X(9);
    a=r*alpha;
    FA=sqrt(Ax^2+Ay^2);
    FB=sqrt(Bx^2+By^2);

elseif p2=='Y' || p2=='y'
    mus = input('Please enter the \mu s : ');
    muk = input('Please enter the \mu k : ');
    nF1=double(Y(1));
    nN1=double(Y(3));
    Fmax=mus*nN1;
    if nF1<=Fmax
        c=0;
        X=double(Y);
        F1=X(1);
        F2=X(2);
        N1=X(3);
        N2=X(4);
        Ax=X(5);
        Ay=X(6);
        Bx=X(7);
        By=X(8);
        alpha=X(9);
        a=r*alpha;
        FA=sqrt(Ax^2+Ay^2);
        FB=sqrt(Bx^2+By^2);
        
    elseif nF1>Fmax
        c=1;
        syms a
        eq1= 2*Ax+2*(muk*N1)-2*mw*g*sin(phi)-2*mw*a==0;
        eq2= 2*Ay+2*N1-2*mw*g*cos(phi)==0;
        eq3=2*T-2*(muk*N1)*r-mw*(r^2)*alpha==0;
        eq4=2*Bx-2*F2-2*mw*g*sin(phi)-2*mw*a==0;
        eq5= 2*By+2*N2-2*mw*g*cos(phi)==0;
        eq6=2*F2*r-mw*(r^2)*(a/r)==0;
        eq7=-2*Ax-2*Bx-m*g*sin(phi)-m*a==0;
        eq8=-2*Ay-2*By-m*g*cos(phi)==0;
        eq9=-2*T+(2*Ax+2*Bx)*(h-r)+a1*2*By-a2*2*Ay==0;
        [A,B] = equationsToMatrix([eq1,eq2,eq3,eq4,eq5,eq6,eq7,eq8,eq9],[a,F2,N1,N2,Ax,Ay,Bx,By,alpha]);
        Y = linsolve(A,B);
        X=double(Y);
        a=X(1);
        F2=X(2);
        N1=X(3);
        N2=X(4);
        Ax=X(5);
        Ay=X(6);
        Bx=X(7);
        By=X(8);
        alpha=X(9);
        F1=muk*N1;
        FA=sqrt(Ax^2+Ay^2);
        FB=sqrt(Bx^2+By^2);
    end
end

% Display numbers
    S(1).name = 'Friction between the road and rear wheels : ';
    S(1).num = F1;
    S(2).name = 'Friction between the road and front wheels : ';
    S(2).num = F2;
    S(3).name = 'Normal reactions on rear wheels by the road : ';
    S(3).num = N1;
    S(4).name = 'Normal reactions on front wheels by the road : ';
    S(4).num = N2;
    S(5).name = 'Acceleration of the car : ';
    S(5).num = a;
    if c==1
    S(6).name = 'Angular acceleration of the rear wheels : ';
    S(6).num = alpha;
    S(7).name = 'Angular acceleration of the front wheels : ';
    S(7).num = (a/r);
    S(8).name = 'Forces exerted on rear axles : ';
    S(8).num = FA;
    S(9).name = 'Forces exerted on front axles : ';
    S(9).num = FB;
    elseif c==0
    S(6).name = 'Angular acceleration of the wheels : ';
    S(6).num = alpha;
    S(7).name = 'Forces exerted on rear axles : ';
    S(7).num = FA;
    S(8).name = 'Forces exerted on front axles : ';
    S(8).num = FB;
    end
    Ans = [{S.name};{S.num}];
    fprintf('\n')
    if c==1
        fprintf('The Car slips')
        fprintf('\n')
    elseif c==0
        fprintf('The car doesnt slip')
        fprintf('\n')
    end
    fprintf('%10s%.2f\n ',Ans{:})
    fprintf('\n')
    
% clear workspace
clear eq1 eq2 eq3 eq4 eq5 eq6 eq7 eq8 eq9 A B nF1 nN1 g p1 p2 X Y S Ans

elseif p1=='Y' || p1=='y'
fprintf('Please wait')
fprintf('\n')
%data
m=500;
mw=15;
a1=1.8;
a2=1.5;
h=0.8;
r=0.35;
T=200;
phi=2;
phi=(phi*pi)/180;

%Plot 1
fprintf('Plot 1 :')
mp=300:1:800;
cc=length(mp);
for i=1:cc
    m=mp(i);
    eq1= 2*Ax+2*F1-2*mw*g*sin(phi)-2*mw*r*alpha==0;
    eq2= 2*Ay+2*N1-2*mw*g*cos(phi)==0;
    eq3=2*T-2*F1*r-mw*(r^2)*alpha==0;
    eq4=2*Bx-2*F2-2*mw*g*sin(phi)-2*mw*r*alpha==0;
    eq5= 2*By+2*N2-2*mw*g*cos(phi)==0;
    eq6=2*F2*r-mw*(r^2)*alpha==0;
    eq7=-2*Ax-2*Bx-m*g*sin(phi)-m*r*alpha==0;
    eq8=-2*Ay-2*By-m*g*cos(phi)==0;
    eq9=-2*T+(2*Ax+2*Bx)*(h-r)+a1*2*By-a2*2*Ay==0;
    [A,B] = equationsToMatrix([eq1,eq2,eq3,eq4,eq5,eq6,eq7,eq8,eq9],[F1,F2,N1,N2,Ax,Ay,Bx,By,alpha]);
    Y = linsolve(A,B);
    aF1(i)=r*(double(Y(9)));
    fprintf('.')
end
fprintf('\n')
subplot(3,2,1);
plot(mp,aF1,'color','#A2142F');
title('a(acceleration) versus m(body mass)')
xlabel('m')
ylabel('a')

%Plot 2
fprintf('Plot 2 :')
m=500;
rp=0.25:0.001:0.5;
cc=length(rp);
for i=1:cc
    r=rp(i);
    eq1= 2*Ax+2*F1-2*mw*g*sin(phi)-2*mw*r*alpha==0;
    eq2= 2*Ay+2*N1-2*mw*g*cos(phi)==0;
    eq3=2*T-2*F1*r-mw*(r^2)*alpha==0;
    eq4=2*Bx-2*F2-2*mw*g*sin(phi)-2*mw*r*alpha==0;
    eq5= 2*By+2*N2-2*mw*g*cos(phi)==0;
    eq6=2*F2*r-mw*(r^2)*alpha==0;
    eq7=-2*Ax-2*Bx-m*g*sin(phi)-m*r*alpha==0;
    eq8=-2*Ay-2*By-m*g*cos(phi)==0;
    eq9=-2*T+(2*Ax+2*Bx)*(h-r)+a1*2*By-a2*2*Ay==0;
    [A,B] = equationsToMatrix([eq1,eq2,eq3,eq4,eq5,eq6,eq7,eq8,eq9],[F1,F2,N1,N2,Ax,Ay,Bx,By,alpha]);
    Y = linsolve(A,B);
    aF2(i)=r*(double(Y(9)));
    fprintf('.')
end
fprintf('\n')
subplot(3,2,2);
plot(rp,aF2,'color','#7E2F8E');
title('a(acceleration) versus r(wheel radius)')
xlabel('r')
ylabel('a')

%Plot 3
fprintf('Plot 3 :')
r=0.35;
np=0.5:0.01:2;
cc=length(np);
for i=1:cc
    a2=3./(np(i)+1);
    a1=a2*np(i);
    eq1= 2*Ax+2*F1-2*mw*g*sin(phi)-2*mw*r*alpha==0;
    eq2= 2*Ay+2*N1-2*mw*g*cos(phi)==0;
    eq3=2*T-2*F1*r-mw*(r^2)*alpha==0;
    eq4=2*Bx-2*F2-2*mw*g*sin(phi)-2*mw*r*alpha==0;
    eq5= 2*By+2*N2-2*mw*g*cos(phi)==0;
    eq6=2*F2*r-mw*(r^2)*alpha==0;
    eq7=-2*Ax-2*Bx-m*g*sin(phi)-m*r*alpha==0;
    eq8=-2*Ay-2*By-m*g*cos(phi)==0;
    eq9=-2*T+(2*Ax+2*Bx)*(h-r)+a1*2*By-a2*2*Ay==0;
    [A,B] = equationsToMatrix([eq1,eq2,eq3,eq4,eq5,eq6,eq7,eq8,eq9],[F1,F2,N1,N2,Ax,Ay,Bx,By,alpha]);
    Y = linsolve(A,B);
    X=double(Y);
    nAx=X(5);
    nAy=X(6);
    nBx=X(7);
    nBy=X(8);
    nFA=sqrt(nAx^2+nAy^2);
    nFB=sqrt(nBx^2+nBy^2);
    q(i)=nFB/nFA;
    fprintf('.')
end
fprintf('\n')
subplot(3,2,3);
plot(np,q,'color','#77AC30');
title('FA/FB versus a1/a2')
xlabel('n')
ylabel('q')

%Plot 4
fprintf('Plot 4 :')
a1=1.8;
a2=1.5;
nT=300:1:1000;
cc=length(nT);
j=0;
for i=1:cc
    T=nT(i);
    eq1= 2*Ax+2*F1-2*mw*g*sin(phi)-2*mw*r*alpha==0;
    eq2= 2*Ay+2*N1-2*mw*g*cos(phi)==0;
    eq3=2*T-2*F1*r-mw*(r^2)*alpha==0;
    eq4=2*Bx-2*F2-2*mw*g*sin(phi)-2*mw*r*alpha==0;
    eq5= 2*By+2*N2-2*mw*g*cos(phi)==0;
    eq6=2*F2*r-mw*(r^2)*alpha==0;
    eq7=-2*Ax-2*Bx-m*g*sin(phi)-m*r*alpha==0;
    eq8=-2*Ay-2*By-m*g*cos(phi)==0;
    eq9=-2*T+(2*Ax+2*Bx)*(h-r)+a1*2*By-a2*2*Ay==0;
    [A,B] = equationsToMatrix([eq1,eq2,eq3,eq4,eq5,eq6,eq7,eq8,eq9],[F1,F2,N1,N2,Ax,Ay,Bx,By,alpha]);
    Y = linsolve(A,B);
    nF2(i)=double(Y(1));
    nMuN(i)=(0.8)*double(Y(3));
    if nMuN(i)<nF2(i)
        Ts=nT(i);
        j=j+1;
    end
    fprintf('.')
end
fprintf('\n')
k=i-j;
Ts=nT(k);
subplot(3,2,4);
plot(nT,nF2,'color','#D95319');
hold on
plot(nT,nMuN,'color','#0072BD');
title('\mu s.N(rear) and f(rear wheel) versus T')
xlabel('T')
ylabel('f')
legend('f','\mu s.N')

%Plot 5
fprintf('Plot 5 :')
TF=linspace(300,Ts,31);
cc=length(TF);
v(1)=0;
j=2;
for i=2:cc
    T=TF(i);
    eq1= 2*Ax+2*F1-2*mw*g*sin(phi)-2*mw*r*alpha==0;
    eq2= 2*Ay+2*N1-2*mw*g*cos(phi)==0;
    eq3=2*T-2*F1*r-mw*(r^2)*alpha==0;
    eq4=2*Bx-2*F2-2*mw*g*sin(phi)-2*mw*r*alpha==0;
    eq5= 2*By+2*N2-2*mw*g*cos(phi)==0;
    eq6=2*F2*r-mw*(r^2)*alpha==0;
    eq7=-2*Ax-2*Bx-m*g*sin(phi)-m*r*alpha==0;
    eq8=-2*Ay-2*By-m*g*cos(phi)==0;
    eq9=-2*T+(2*Ax+2*Bx)*(h-r)+a1*2*By-a2*2*Ay==0;
    [A,B] = equationsToMatrix([eq1,eq2,eq3,eq4,eq5,eq6,eq7,eq8,eq9],[F1,F2,N1,N2,Ax,Ay,Bx,By,alpha]);
    Y = linsolve(A,B);
    aF3(i)=r*(double(Y(9)));
    j=i-1;
    v(i)=aF3(i)+v(j);
    fprintf('.')
end
fprintf('\n')
subplot(3,2,[5,6]);
time=0:1:30;
plot(time,v,'color','#4DBEEE');
title('v(velocity) versus t(time)')
xlabel('t')
ylabel('v')

% clear workspace
clear eq1 eq2 eq3 eq4 eq5 eq6 eq7 eq8 eq9 A B g p1 X Y cc time F1 F2 Ax Ay Bx By N1 N2 alpha i j k
end