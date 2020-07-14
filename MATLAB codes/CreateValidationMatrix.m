%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CreateValidationMatrix.m
% Used in IO_Script to store qualitative experimental observations,
% as described in our accompanying publication (Fig 2): 
% Irons & Humphrey (2020): Cell signaling model for arterial mechanobiology,
% PLOS Computational Biology. 
%-----------------------------------------------
% Created by Linda Irons: linda.irons@yale.edu
% Last modified by Linda Irons, July 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Exp_validation=CreateValidationMatrix(speciesNames, Ordering)

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

Exp_validation=0.1*ones(3,length(speciesNames)); %a non-zero value so 'no change' can be recorded as 0.

Exp_validation(Stress,[TGFb1, TSP1, TIMP, MMP2, MMP9, Akt, p70S6K, Col1mRNA, Col3mRNA, ActomyosinActivity, SMCproliferation])=1; %Stress increases
Exp_validation(Stress,[])=0; %Stress doesn't change
Exp_validation(Stress,[])=-1; %Stress decreases
Exp_validation(Stress,[ET1])=-0.1; %mixed

Exp_validation(Wss,[NO, TGFb1])=1; 
Exp_validation(Wss,MMP9)=0; 
Exp_validation(Wss,[MMP2, ET1, Akt, SMCproliferation])=-1;

Exp_validation(AngIIin,[TSP1, Akt, MMP1, p70S6K, Col1mRNA, Col3mRNA, ActomyosinActivity, SMCproliferation])=1;
Exp_validation(AngIIin,[])=0;
Exp_validation(AngIIin,[])=-1;
Exp_validation(AngIIin,[TIMP,TGFb1, MMP2, MMP9])=-0.1; 

Exp_validation=Exp_validation(:,Ordering);

end