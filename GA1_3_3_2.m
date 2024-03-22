function GA1_3_3_2
%% Minimum value of multivariate function
clc
clear all
close all
%%
figure(1);
lbx = -2;ubx = 2;
lby = -2;uby = 2;
ezmesh('y*sin(2*pi*x)+x*cos(2*pi*y)',[lbx,ubx,lby,uby],50);
hold on;
%%
NIND = 40;
MAXGEN = 50;
PRECI = 20;
GGAP = 0.95;
px = 0.7;
pm = 0.01;
trace = zeros(3,MAXGEN);
FieldD=[PRECI PRECI;lbx lby;ubx uby;1 1;0 0;1 1;1 1];
Chrom = crtbp(NIND,PRECI*2);
%%
gen = 0;
XY = bs2rv(Chrom,FieldD);
X = XY(:,1);Y = XY(:,2);
ObjV = Y.*sin(2*pi*X)+X.*cos(2*pi*Y);
while gen<MAXGEN
    FitnV = ranking(-ObjV);
    SelCh = select('sus',Chrom,FitnV,GGAP);
    SelCh = recombin('xovsp',SelCh,px);
    SelCh = mut(SelCh,pm);
    XY = bs2rv(SelCh,FieldD);
    X = XY(:,1);Y = XY(:,2);
    ObjVSel = Y.*sin(2*pi*X)+X.*cos(2*pi*Y);
    [Chrom,ObjV] = reins(Chrom,SelCh,1,1,ObjV,ObjVSel);
    XY = bs2rv(Chrom,FieldD);
    gen = gen+1;
    [Y,I] = max(ObjV);
    trace(1:2,gen) = XY(I,:);
    trace(3,gen) = Y;
end
plot3(trace(1,:),trace(2,:),trace(3,:),'bo');
grid on;
plot3(XY(:,1),XY(:,2),ObjV,'bo');
hold off
%%
figure(2)
plot(1:MAXGEN,trace(3,:));
grid on
xlabel('遗传代数')
ylabel('解的变化')
title('进化过程')
bestZ = trace(3,end);
bestX = trace(1,end);
bestY = trace(2,end);
fprintf(['最优解\nx=',num2str(bestX),'\ny=',num2str(bestY),'\nz=',num2str(bestZ)])