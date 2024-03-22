function GA1_3_3
%% Minimum value of a unary function use GA
%% find  min f(x) = sin(10*pi*x)/x x=[1,2] 
clc
clear all
close all
%% plot f(x) 
figure(1)
hold on;
lb = 1;
ub = 2;                                                         % f(x) independent variable x=[1,2]
ezplot('sin(10*pi*X) / X',[lb,ub]);                             % f(x) plot curve
xlabel('自变量/X');
ylabel('函数值/Y');
%% define GA algorithm parameter
NIND = 40;                                                      % population size
MAXGEN = 20;                                                    % max genetic algebra
PRECI = 20;                                                     % individual length
GGAP = 0.95;                                                    % generation gap
px = 0.7;                                                       % cross probability
pm = 0.01;                                                      % mutation probability
trace = zeros(2,MAXGEN);                                        % initial vaule optimal result
FieldD = [PRECI;lb;ub;1;0;1;1];                                 % field descriptor
Chrom = crtbp(NIND,PRECI);                                      % create any discrete random population
%% optimzation
gen = 0;                                                        % generation count
X = bs2rv(Chrom,FieldD);                                        % initial population BIN to DEC
ObjV = sin(10*pi*X)./X;                                         % calculate objective function value
while gen < MAXGEN
    FitnV = ranking(ObjV);                                      % assign fitness values
    SelCh = select('sus',Chrom,FitnV,GGAP);                     % select
    SelCh = recombin('xovsp',SelCh,px);                         % recombination
    SelCh = mut(SelCh,pm);                                      % mutation
    X = bs2rv(SelCh,FieldD);                                    % descendant individual DEC convert
    ObjVSel = sin(10 * pi * X)./X;                              % cal descendant objective fun value
    [Chrom,ObjV] = reins(Chrom,SelCh,1,1,ObjV,ObjVSel);         % reinsert child to parent, get new population
    X = bs2rv(Chrom,FieldD);
    gen = gen+1;                                                % add gen count
    % obtain optimal solution and sequence for each generation
    % Y:opt solution,I: sequence
    [Y,I] = min(ObjV);
    trace(1,gen) = X(I);                                        % every gen opt
    trace(2,gen) = Y;
end
plot(trace(1,:),trace(2,:),'bo');                               % plot every gen opt
grid on
plot(X,ObjV,'b*');                                              % plot last gen population hold off
%% plot evolutionary diagram
figure(2)
plot(1:MAXGEN,trace(2,:));
grid on
xlabel('遗传代数')
ylabel('解的变化')
title('进化过程')
bestY = trace(2,end);
bestX = trace(1,end);
fprintf(['最优解:\nX=',num2str(bestX),'\nY=',num2str(bestY),'\n'])

