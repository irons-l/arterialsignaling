%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DefineOrdering.m
% Used in Knockdown_Script and IO_Script to reorder speciesNames
%-----------------------------------------------
% Created by Linda Irons: linda.irons@yale.edu
% Last modified by Linda Irons, July 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Ordering,OrderingStr]=DefineOrdering()

Stress = 1;
Wss = 2;
AngIIin = 3; 
SACs = 4;
Integrins = 5;
PDGF = 6;
AngII = 7;
latentTGFb1 = 8;
TGFb1 = 9;
TGFbR2 = 10;
TGFbR1 = 11;
pSmad23 = 12;
Smad4 = 13;
Smad7 = 14;
TSP1 = 15;
TIMP = 16;
p38 = 17;
JNK = 18;
ERK = 19;
MMP1 = 20;
MMP2 = 21;
MMP9 = 22;
AT1R = 23;
AT2R = 24;
PDGFR = 25;
NO = 26;
ET1 = 27;
ETAR = 28;
ETBR = 29;
PI3K = 30;
Akt = 31;
mTOR = 32;
mTORC1 = 33;
mTORC2 = 34;
p70S6K = 35;
Ca = 36;
MLCK = 37;
Myosin = 38;
FAK = 39;
Cdc42 = 40;
Arp23 = 41;
RhoA = 42;
ROCK = 43;
Actin = 44;
Col1mRNA = 45;
Col3mRNA = 46;
Col1 = 47;
Col3 = 48;
ActomyosinActivity = 49;
SMCproliferation = 50;

% Desired order
OrderingStr={'Stress';'Wss';'AngIIin';'SACs';'Integrins';'AngII';'AT1R';'AT2R';'PDGF';'PDGFR';'latentTGFb1';'TGFb1';'TGFbR2';'TGFbR1';
    'pSmad23';'Smad4';'Smad7';'p38';'JNK';'ERK';'TSP1'; 'MMP1';'MMP2';'MMP9';'TIMP';'Col1mRNA';'Col3mRNA';'Col1';'Col3';
    'NO';'ET1';'ETAR';'ETBR';'PI3K';'Akt';'mTOR';'mTORC1';'mTORC2';'p70S6K';
    'Ca';'MLCK';'Myosin';'FAK';'Cdc42';'Arp23';'RhoA';'ROCK';'Actin';'ActomyosinActivity';'SMCproliferation'};

Ordering=zeros(1,length(OrderingStr));
for i=1:length(OrderingStr)
    Ordering(i)=eval(OrderingStr{i});
end

end