% This script is partly based on open source code available at: 
%       github.com/saucermanlab/netflux
% for simulating logic-based signaling networks as described originally in: 
%       Kraeutler, M.J., Soltis, A.R., & Saucerman, J.J. (2010). 'Modeling 
%       cardiac B-adrenergic signaling with normalized-Hill differential
%       equations: comparison with a biochemical model.' BMC Systems 
%       Biology.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Comparison_Wu_Script.m
% Simulates Col1mRNA and Col3mRNA levels at 3 levels of Stress, with and without
% exogenous AngII, and compares model predictions to experimental data from 
% Wu et al (2014), as described in our accompanying publication: 
% Irons & Humphrey (2020): Cell signaling model for arterial mechanobiology,
% PLOS Computational Biology. (Reproduces Fig 6)
%-----------------------------------------------
% Created by Linda Irons: linda.irons@yale.edu
% Last modified by Linda Irons, July 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
clear global

Col1mRNA=45;
Col3mRNA=46;
AngIIin=3;
Stress=1;
p38=17;

load('reactions_final.mat');
ODEfilename='ODElist_final';

%Uniform Parameters
w=1; 
n=1.25;
EC50=0.55;
p0=[w;n;EC50];                  
        
%Initial conditions (those that are fixed for all simulations) 
y0=zeros(1,length(speciesNames));
ival=0.2; 
y0([2,4,5])=[0.5,ival,ival];  %Wss, SACs, Integrins

%%% Low Stress %%%
%Stress=0.2, AngII=0
y0(1)=0.2;
y0(3)=0;
ymax(p38)=1;

%Generate data
[~,y_ref]=ODE_master(p0,y0,speciesNames,tau,ymax,  reactionRules,ODEfilename);

BasalCol1=y_ref(end,Col1mRNA);
BasalCol3=y_ref(end,Col3mRNA);
Fold_Col1=y_ref(end,Col1mRNA)/BasalCol1;
Fold_Col3=y_ref(end,Col3mRNA)/BasalCol3;
 
A1=Fold_Col1;
B1=Fold_Col3;
 
%Stress=0.2, AngII=0.2
y0(1)=0.2;
y0(3)=0.2;
ymax(p38)=1;

%Generate data
[~,y_ref]=ODE_master(p0,y0,speciesNames,tau,ymax,  reactionRules,ODEfilename); 

Fold_Col1=y_ref(end,Col1mRNA)/BasalCol1;
Fold_Col3=y_ref(end,Col3mRNA)/BasalCol3;
 
A4=Fold_Col1;
B4=Fold_Col3;

%Stress=0.2, AngII=0, p38 knockdown
y0(1)=0.2;
y0(3)=0;
ymax(p38)=0.1;

%Generate data
[~,y_ref]=ODE_master(p0,y0,speciesNames,tau,ymax,  reactionRules,ODEfilename); 

Fold_Col1=y_ref(end,Col1mRNA)/BasalCol1;
Fold_Col3=y_ref(end,Col3mRNA)/BasalCol3;
 
A7=Fold_Col1;
B7=Fold_Col3;
 
%%% Mid Stress %%%
%Stress=0.3, AngII=0
y0(1)=0.3;
y0(3)=0;
ymax(p38)=1;

%Generate data
[~,y_ref]=ODE_master(p0,y0,speciesNames,tau,ymax,  reactionRules,ODEfilename);
   
Fold_Col1=y_ref(end,Col1mRNA)/BasalCol1;
Fold_Col3=y_ref(end,Col3mRNA)/BasalCol3;
 
A2=Fold_Col1;
B2=Fold_Col3;

%Stress=0.3, AngII=0.2
y0(1)=0.3;
y0(3)=0.2;
ymax(p38)=1;

%Generate data
[~,y_ref]=ODE_master(p0,y0,speciesNames,tau,ymax,  reactionRules,ODEfilename);
  
Fold_Col1=y_ref(end,Col1mRNA)/BasalCol1;
Fold_Col3=y_ref(end,Col3mRNA)/BasalCol3;
 
A5=Fold_Col1;
B5=Fold_Col3;
 
%Stress=0.3, AngII=0, p38 knockdown
y0(1)=0.3;
y0(3)=0;
ymax(p38)=0.1;

%Generate data
[~,y_ref]=ODE_master(p0,y0,speciesNames,tau,ymax,  reactionRules,ODEfilename);
  
Fold_Col1=y_ref(end,Col1mRNA)/BasalCol1;
Fold_Col3=y_ref(end,Col3mRNA)/BasalCol3;
 
A8=Fold_Col1;
B8=Fold_Col3;
 
%Stress=0.4, AngII=0
y0(1)=0.4;
y0(3)=0;
ymax(p38)=1;

%Generate data
[~,y_ref]=ODE_master(p0,y0,speciesNames,tau,ymax,  reactionRules,ODEfilename);

Fold_Col1=y_ref(end,Col1mRNA)/BasalCol1;
Fold_Col3=y_ref(end,Col3mRNA)/BasalCol3;
 
A3=Fold_Col1;
B3=Fold_Col3;

%Stress=0.4, AngII=0.2
y0(1)=0.4;
y0(3)=0.2;
ymax(p38)=1;

%Generate data
[~,y_ref]=ODE_master(p0,y0,speciesNames,tau,ymax,  reactionRules,ODEfilename);
 
Fold_Col1=y_ref(end,Col1mRNA)/BasalCol1;
Fold_Col3=y_ref(end,Col3mRNA)/BasalCol3;
 
A6=Fold_Col1;
B6=Fold_Col3;

%Stress=0.4, AngII=0, p38 knockdown
y0(1)=0.4;
y0(3)=0;
ymax(p38)=0.1;
 
%Generate data
[~,y_ref]=ODE_master(p0,y0,speciesNames,tau,ymax,  reactionRules,ODEfilename);
  
Fold_Col1=y_ref(end,Col1mRNA)/BasalCol1;
Fold_Col3=y_ref(end,Col3mRNA)/BasalCol3;
 
A9=Fold_Col1;
B9=Fold_Col3;
 
%% Plot
 
figure(); hold on;
subplot(2,1,1);
bar([1,2,3],[A1,A2,A3]); hold on; bar([5,6,7],[A4,A5,A6]); bar([9,10,11],[A7,A8,A9]);
plot(1,1,'ko','LineWidth',1,'MarkerSize',6,'MarkerFaceColor','k')
plot(2,1.33,'ko','LineWidth',1,'MarkerSize',6,'MarkerFaceColor','k');
plot(3,2.24,'ko','LineWidth',1,'MarkerSize',6,'MarkerFaceColor','k');
plot(5,1.03,'ko','LineWidth',1,'MarkerSize',6,'MarkerFaceColor','k');
plot(6,1.12,'ko','LineWidth',1,'MarkerSize',6,'MarkerFaceColor','k');
plot(7,2.8,'ko','LineWidth',1,'MarkerSize',6,'MarkerFaceColor','k');
plot(11,1.36,'ko','LineWidth',1,'MarkerSize',6,'MarkerFaceColor','k');
plot([-0.5,12.5],[1,1],'k--');
xlim([0,12]);
ylabel('Col1mRNA');

subplot(2,1,2);
bar([1,2,3],[B1,B2,B3]); hold on; bar([5,6,7],[B4,B5,B6]); bar([9,10,11],[B7,B8,B9]);
plot(1,1,'ko','LineWidth',1,'MarkerSize',5,'MarkerFaceColor','k');
plot(2,1.98,'ko','LineWidth',1,'MarkerSize',5,'MarkerFaceColor','k');
plot(3,3.4,'ko','LineWidth',1,'MarkerSize',5,'MarkerFaceColor','k');
plot(5,1.59,'ko','LineWidth',1,'MarkerSize',5,'MarkerFaceColor','k');
plot(6,4.28,'ko','LineWidth',1,'MarkerSize',5,'MarkerFaceColor','k');
plot(7,7.5,'ko','LineWidth',1,'MarkerSize',5,'MarkerFaceColor','k');
plot(11,1.56,'ko','LineWidth',1,'MarkerSize',6,'MarkerFaceColor','k');
plot([-0.5,12.5],[1,1],'k--');
xlim([0,12]);
ylabel('Col3mRNA');