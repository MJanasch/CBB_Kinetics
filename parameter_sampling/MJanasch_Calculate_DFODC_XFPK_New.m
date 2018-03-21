%% Function to calculate dfodc for the XFPK model
% Markus Janasch, Ph.D. Student, KTH
% Created: 2017-10-26, last modified: 2017-10-26

% This function reassigns the right parameter values to the constants used
% for promiscous enzymes and calculates the numerical values of the partial
% derivatives of the rate equations

function [dfodc,ParameterSet] = MJanasch_Calculate_DFODC_XFPK_New(N,ParameterSet,V_K_Indeces,Fluxes,SFull)


%% Define Metabolite Concentrations

    BPG         = N.species(1).initialConcentration;
    P2G         = N.species(2).initialConcentration;
    P3G         = N.species(3).initialConcentration;
    ADP         = N.species(4).initialConcentration;
    DHAP        = N.species(5).initialConcentration;
    E4P         = N.species(6).initialConcentration;
    F6P         = N.species(7).initialConcentration;
    FBP         = N.species(8).initialConcentration;
    GAP         = N.species(9).initialConcentration;
    NADPH       = N.species(10).initialConcentration;
    PEP         = N.species(11).initialConcentration;
    PHI         = N.species(12).initialConcentration;
    R5P         = N.species(13).initialConcentration;
    Ru5P        = N.species(14).initialConcentration;
    RuBP        = N.species(15).initialConcentration;
    S7P         = N.species(16).initialConcentration;
    SBP         = N.species(17).initialConcentration;
    Xu5P        = N.species(18).initialConcentration;
    PYR         = N.species(19).initialConcentration;
    ACETP       = N.species(20).initialConcentration;
    ATP         = N.species(21).initialConcentration;
    NADP        = N.species(22).initialConcentration;
    O2          = N.species(23).initialConcentration;
    CO2_cax     = N.species(24).initialConcentration;
    CO2_cyt     = N.species(25).initialConcentration;
    ACCOA       = N.species(26).initialConcentration;
    NADH        = N.species(27).initialConcentration;
    NAD         = N.species(28).initialConcentration;
    COA         = N.species(29).initialConcentration;
    BioMass_R5P = N.species(30).initialConcentration;
    BioMass_F6P = N.species(31).initialConcentration;
    BioMass_E4P = N.species(32).initialConcentration;
    BioMass_PEP = N.species(33).initialConcentration;
    BioMass_PYR = N.species(34).initialConcentration;
    BioMass_P3G = N.species(35).initialConcentration;
    PPool       = N.species(36).initialConcentration;

%% Modify Parameters for Enzyme Promiscuity
% Set parameter value of later parameters working in the same enzyme to the
% same value as the parameters are they firstly appear in the parameter
% list

% For TKT2 (same values as in TKT1)
ParameterSet(34) = ParameterSet(42); % KmF6P
ParameterSet(35) = ParameterSet(43); % KmGAP
ParameterSet(36) = ParameterSet(44); % KmE4P
ParameterSet(37) = ParameterSet(45); % KmXuP
ParameterSet(38) = ParameterSet(46); % KmS7P
ParameterSet(39) = ParameterSet(47); % KmR5P

% For FBA (same values as in ALD)
ParameterSet(27) = ParameterSet(54); % KmFBP
ParameterSet(28) = ParameterSet(55); % KmDHAP
ParameterSet(29) = ParameterSet(56); % KmGAP
ParameterSet(30) = ParameterSet(57); % KmSBP
ParameterSet(31) = ParameterSet(58); % KmE4P

% For SBPase (same values as FBPase)
ParameterSet(50) = ParameterSet(61); % KmFBP
ParameterSet(51) = ParameterSet(62); % KmSBP

% For XFPK2 (same values as XFPK1)
ParameterSet(107) = ParameterSet(115); % KmF6P
ParameterSet(108) = ParameterSet(116); % KmPHI
ParameterSet(109) = ParameterSet(117); % KmE4P
ParameterSet(110) = ParameterSet(118); % KmACETP
ParameterSet(111) = ParameterSet(119); % KmXu5P
ParameterSet(112) = ParameterSet(120); % KmGAP


%% Define Parameters

% Set Vmax's and K to 1, if not already
    for j = 1:length(V_K_Indeces)
        ParameterSet(V_K_Indeces(j)) = 1;     % Set Vmax & K to 1
    end
    
    
    V_maxv1     = ParameterSet(1);
    KmCO2_caxv1 = ParameterSet(2);
    KmO2v1      = ParameterSet(3);
    KmRuBPv1    = ParameterSet(4);
    K_eqv1      = ParameterSet(5);
    KiPHIv1     = ParameterSet(6);
    KiNADPHv1   = ParameterSet(7);
    KaPHIv1     = ParameterSet(8);
    
    V_maxv2     = ParameterSet(9);
    KmP3Gv2     = ParameterSet(10);
    KmATPv2     = ParameterSet(11);
    KmBPGv2     = ParameterSet(12);
    KmADPv2     = ParameterSet(13);
    K_eqv2      = ParameterSet(14);
    
    V_maxv3     = ParameterSet(15);
    KmNADPHv3   = ParameterSet(16);
    KmBPGv3     = ParameterSet(17);
    KmNADPv3    = ParameterSet(18);
    KmGAPv3     = ParameterSet(19);
    KmPHIv3     = ParameterSet(20);
    K_eqv3      = ParameterSet(21);
    
    V_maxv4     = ParameterSet(22);
    KmDHAPv4    = ParameterSet(23);
    KmGAPv4     = ParameterSet(24);
    K_eqv4      = ParameterSet(25);
    
    V_maxv5     = ParameterSet(26);
    KmFBPv5     = ParameterSet(27);
    KmDHAPv5    = ParameterSet(28);
    KmGAPv5     = ParameterSet(29);
    KmSBPv5     = ParameterSet(30);
    KmE4Pv5     = ParameterSet(31);
    K_eqv5      = ParameterSet(32);
    
    V_maxv6     = ParameterSet(33);
    KmF6Pv6     = ParameterSet(34);
    KmGAPv6     = ParameterSet(35);
    KmE4Pv6     = ParameterSet(36);
    KmXu5Pv6    = ParameterSet(37);
    KmS7Pv6     = ParameterSet(38);
    KmR5Pv6     = ParameterSet(39);
    K_eqv6      = ParameterSet(40);
    
    V_maxv7     = ParameterSet(41);
    KmF6Pv7     = ParameterSet(42); % like KmF6Pv6
    KmGAPv7     = ParameterSet(43); % like KmGAPv6
    KmE4Pv7     = ParameterSet(44); % like KmE4Pv6
    KmXu5Pv7    = ParameterSet(45); % like KmXu5Pv6
    KmS7Pv7     = ParameterSet(46); % like KmS7Pv6
    KmR5Pv7     = ParameterSet(47); % like KmR5Pv6
    K_eqv7      = ParameterSet(48);
    
    V_maxv8     = ParameterSet(49);
    KmFBPv8     = ParameterSet(50);
    KmSBPv8     = ParameterSet(51);
    K_eqv8      = ParameterSet(52);
    
        
    V_maxv9    = ParameterSet(53);
    KmFBPv9    = ParameterSet(54); % like KmFBPv5
    KmDHAPv9   = ParameterSet(55); % like KmDHAPv5
    KmGAPv9    = ParameterSet(56); % like KmGAPv5
    KmSBPv9    = ParameterSet(57); % like KmSBPv5
    KmE4Pv9    = ParameterSet(58); % like KmE4Pv5
    K_eqv9     = ParameterSet(59);
    
    V_maxv10    = ParameterSet(60);
    KmFBPv10    = ParameterSet(61); % like KmFBPv8
    KmSBPv10    = ParameterSet(62); % like KmSBPv8
    K_eqv10     = ParameterSet(63);
    
    V_maxv11    = ParameterSet(64);
    KmR5Pv11    = ParameterSet(65);
    KmRu5Pv11   = ParameterSet(66);
    K_eqv11     = ParameterSet(67);
   
    V_maxv12    = ParameterSet(68);
    KmXu5Pv12   = ParameterSet(69);
    KmRu5Pv12   = ParameterSet(70);
    K_eqv12     = ParameterSet(71);
    
    V_maxv13    = ParameterSet(72);
    KmRu5Pv13   = ParameterSet(73);
    KmATPv13    = ParameterSet(74);
    KmRuBPv13   = ParameterSet(75);
    KmADPv13    = ParameterSet(76);
    K_eqv13     = ParameterSet(77);
    KiPEPv13    = ParameterSet(78);
    KiADPv13    = ParameterSet(79);
    
    V_maxv14    = ParameterSet(80);
    KmP3Gv14    = ParameterSet(81);
    KmP2Gv14    = ParameterSet(82);
    K_eqv14     = ParameterSet(83);
    
    V_maxv15    = ParameterSet(84);
    KmP2Gv15    = ParameterSet(85);
    KmPEPv15    = ParameterSet(86);
    K_eqv15     = ParameterSet(87);
    
    V_maxv16    = ParameterSet(88);
    KmPEPv16    = ParameterSet(89);
    KmPYRv16    = ParameterSet(90);
    KmADPv16    = ParameterSet(91);
    KmATPv16    = ParameterSet(92);
    K_eqv16     = ParameterSet(93);
    KaR5Pv16    = ParameterSet(94);
    KiPHIv16    = ParameterSet(95);
    KiFBPv16    = ParameterSet(96);
    KiATPv16    = ParameterSet(97);
   
    V_maxv17    = ParameterSet(98);
    KmPYRv17    = ParameterSet(99);
    KmNADv17    = ParameterSet(100);
    KmCOAv17    = ParameterSet(101);
    KmACCOAv17  = ParameterSet(102);
    KmNADHv17   = ParameterSet(103);
    KmCO2_cytv17= ParameterSet(104);
    K_eqv17     = ParameterSet(105);
    
    V_maxv18    = ParameterSet(106);
    KmF6Pv18    = ParameterSet(107);
    KmPHIv18    = ParameterSet(108);
    KmE4Pv18    = ParameterSet(109);
    KmACETPv18  = ParameterSet(110);
    KmXu5Pv18   = ParameterSet(111);
    KmGAPv18    = ParameterSet(112);
    K_eqv18     = ParameterSet(113);
    
    V_maxv19    = ParameterSet(114);
    KmF6Pv19    = ParameterSet(115); % like KmF6Pv18 
    KmPHIv19    = ParameterSet(116); % like KmPHIv18
    KmE4Pv19    = ParameterSet(117); % like KmE4Pv18
    KmACETPv19  = ParameterSet(118); % like KmACETPv18
    KmXu5Pv19   = ParameterSet(119); % like KmXu5Pv18
    KmGAPv19    = ParameterSet(120); % like KmGAPv18
    K_eqv19     = ParameterSet(121);
    
    V_maxv20    = ParameterSet(122);
    KmACETPv20  = ParameterSet(123);
    KmCOAv20    = ParameterSet(124);
    KmACCOAv20  = ParameterSet(125);
    KmPHIv20    = ParameterSet(126);
    K_eqv20     = ParameterSet(127);
    
    V_maxv21    = ParameterSet(128);
    KmADPv21    = ParameterSet(129);
    KmPHIv21    = ParameterSet(130);
    KmATPv21    = ParameterSet(131);
    K_eqv21     = ParameterSet(132);
    
    V_maxv22    = ParameterSet(133);
    KmNADPv22   = ParameterSet(134);
    KmNADPHv22  = ParameterSet(135);
    K_eqv22     = ParameterSet(136);
    
    V_maxv23    = ParameterSet(137);
    KmR5Pv23    = ParameterSet(138);
    
    V_maxv24    = ParameterSet(139);
    KmF6Pv24    = ParameterSet(140);
    
    V_maxv25    = ParameterSet(141);
    KmE4Pv25    = ParameterSet(142);
    
    V_maxv26    = ParameterSet(143);
    KmPEPv26    = ParameterSet(144);

    V_maxv27    = ParameterSet(145);
    KmPYRv27    = ParameterSet(146);

    V_maxv28    = ParameterSet(147);
    KmP3Gv28    = ParameterSet(148);
    
    K_PPoolv29  = ParameterSet(149);


%% Rate Equations    
    
    pippo(1,1) = (PHI/(PHI+KaPHIv1))*V_maxv1*(RuBP*CO2_cax/(CO2_cax+KmCO2_caxv1*(1+O2/KmO2v1)))/(RuBP+KmRuBPv1*(1+PHI/KiPHIv1+NADPH/KiNADPHv1));
    pippo(2,1) = V_maxv2*(P3G*ATP/(KmP3Gv2*KmATPv2))*(1-BPG*ADP/(P3G*ATP)/K_eqv2)/((1+P3G/KmP3Gv2+BPG/KmBPGv2)*(1+ATP/KmATPv2+ADP/KmADPv2));
    pippo(3,1) = V_maxv3*(BPG*NADPH/(KmBPGv3*KmNADPHv3))*(1-GAP*NADP*PHI/(BPG*NADPH)/K_eqv3)/((1+NADP/KmNADPv3+NADPH/KmNADPHv3)*(1+GAP/KmGAPv3+BPG/KmBPGv3+PHI/KmPHIv3));
    pippo(4,1) = V_maxv4*(GAP/KmGAPv4)*(1-DHAP/GAP/K_eqv4)/(1+GAP/KmGAPv4+DHAP/KmDHAPv4);
    pippo(5,1) = V_maxv5*((DHAP/KmDHAPv5) * (GAP/KmGAPv5))*(1-FBP/(DHAP*GAP)/K_eqv5)/((1+FBP/KmFBPv5)*(1+SBP/KmSBPv5)+(1+DHAP/KmDHAPv5)*(1+E4P/KmE4Pv5)*(1+GAP/KmGAPv5)-1);
    
    pippo(6,1) = V_maxv6*(F6P*GAP/(KmF6Pv6*KmGAPv6))*(1-E4P*Xu5P/(F6P*GAP)/K_eqv6)/((1+F6P/KmF6Pv6+E4P/KmE4Pv6)*(1+GAP/KmGAPv6+Xu5P/KmXu5Pv6)*(1+S7P/KmS7Pv6+R5P/KmR5Pv6));
    pippo(7,1) = V_maxv7*(S7P*GAP/(KmS7Pv7*KmGAPv7))*(1-R5P*Xu5P/(S7P*GAP)/K_eqv7)/((1+F6P/KmF6Pv7+E4P/KmE4Pv7)*(1+GAP/KmGAPv7+Xu5P/KmXu5Pv7)*(1+S7P/KmS7Pv7+R5P/KmR5Pv7));
    pippo(8,1) = V_maxv8 * FBP/(FBP + KmFBPv8*(1+(SBP/KmSBPv8)));
    pippo(9,1) = V_maxv9*((DHAP/KmDHAPv9) * (E4P/KmE4Pv9))*(1-SBP/(DHAP*E4P)/K_eqv9)/((1+FBP/KmFBPv9)*(1+SBP/KmSBPv9)+(1+DHAP/KmDHAPv9)*(1+E4P/KmE4Pv9)*(1+GAP/KmGAPv9)-1);
    pippo(10,1)= V_maxv10 * SBP/(SBP + KmSBPv10*(1+(FBP/KmFBPv10)));
    
    pippo(11,1)= V_maxv11*(R5P/KmR5Pv11)*(1-Ru5P/R5P/K_eqv11)/(1+Ru5P/KmRu5Pv11+R5P/KmR5Pv11);
    pippo(12,1)= V_maxv12*(Xu5P/KmXu5Pv12)*(1-Ru5P/Xu5P/K_eqv12)/(1+Ru5P/KmRu5Pv12+Xu5P/KmXu5Pv12);
    pippo(13,1)= (KiPEPv13/(PEP+KiPEPv13))*(KiADPv13/(ADP+KiADPv13)) * V_maxv13*(Ru5P*ATP/(KmRu5Pv13*KmATPv13))*(1-RuBP*ADP/(Ru5P*ATP)/K_eqv13)/((1+Ru5P/KmRu5Pv13+RuBP/KmRuBPv13)*(1+ATP/KmATPv13+ADP/KmADPv13));
    pippo(14,1)= V_maxv14*(P3G/KmP3Gv14)*(1-P2G/P3G/K_eqv14)/(1+P3G/KmP3Gv14+P2G/KmP2Gv14);
    pippo(15,1)= V_maxv15*(P2G/KmP2Gv15)*(1-PEP/P2G/K_eqv15)/(1+P2G/KmP2Gv15+PEP/KmPEPv15);
    
    pippo(16,1)= (R5P/(R5P+KaR5Pv16))*(KiPHIv16/(PHI+KiPHIv16))*(KiFBPv16/(FBP+KiFBPv16))*(KiATPv16/(ATP+KiATPv16))*V_maxv16*(PEP*ADP/(KmPEPv16*KmADPv16))*(1-PYR*ATP/(PEP*ADP)/K_eqv16)/((1+PEP/KmPEPv16+PYR/KmPYRv16)*(1+ADP/KmADPv16+ATP/KmATPv16));
    pippo(17,1)= V_maxv17 * (PYR * NAD * COA/(KmPYRv17 * KmNADv17 * KmCOAv17))*(1-(ACCOA * NADH * CO2_cyt)/(PYR*NAD*COA*K_eqv17))/((1+ NAD/KmNADv17 + NADH/KmNADHv17)*(1 + PYR/KmPYRv17 + CO2_cyt/KmCO2_cytv17) * (ACCOA/KmACCOAv17 + COA/KmCOAv17));
    
    pippo(18,1)= V_maxv18*(F6P*PHI/(KmF6Pv18*KmPHIv18))*(1-E4P*ACETP/(F6P*PHI)/K_eqv18)/((1+F6P/KmF6Pv18+E4P/KmE4Pv18)*(1+PHI/KmPHIv18+ACETP/KmACETPv18)*(1+Xu5P/KmXu5Pv18+GAP/KmGAPv18));
    pippo(19,1)= V_maxv19*(Xu5P*PHI/(KmXu5Pv19*KmPHIv19))*(1-GAP*ACETP/(Xu5P*PHI)/K_eqv19)/((1+F6P/KmF6Pv19+E4P/KmE4Pv19)*(1+PHI/KmPHIv19+ACETP/KmACETPv19)*(1+Xu5P/KmXu5Pv19+GAP/KmGAPv19));
    pippo(20,1)= V_maxv20*(ACETP*COA/(KmACETPv20*KmCOAv20))*(1-ACCOA*PHI/(ACETP*COA)/K_eqv20)/((1+ACETP/KmACETPv20+ACCOA/KmACCOAv20)*(1+COA/KmCOAv20+PHI/KmPHIv20));
    
    pippo(21,1)= V_maxv21*(ADP*PHI-ATP/K_eqv21)/(KmADPv21*KmPHIv21*(1+ADP/KmADPv21+PHI/KmPHIv21+ATP/KmATPv21+(ADP*PHI)/(KmADPv21*KmPHIv21)));
    pippo(22,1)= V_maxv22 * (NADP / KmNADPv22) * (1 - NADPH / NADP / K_eqv22) / (1 + NADP / KmNADPv22 + NADPH / KmNADPHv22);
    
    pippo(23,1)= V_maxv23 * R5P/(KmR5Pv23 + R5P);
    pippo(24,1)= V_maxv24 * F6P/(KmF6Pv24 + F6P);
    pippo(25,1)= V_maxv25 * E4P/(KmE4Pv25 + E4P);
    pippo(26,1)= V_maxv26 * PEP/(KmPEPv26 + PEP);
    pippo(27,1)= V_maxv27 * PYR/(KmPYRv27 + PYR);
    pippo(28,1)= V_maxv28 * P3G/(KmP3Gv28 + P3G);
    pippo(29,1)= K_PPoolv29*(PPool - PHI);
    
    
%% Calculate Vmax and K
     
    for j=1:length(N.reaction)
        Fluxes(j);                                             
        Vmax_K = Fluxes(j)/pippo(j,1);                   
        
        ParameterSet(V_K_Indeces(j)) = Vmax_K;
    end
    
    V_maxv1     = ParameterSet(V_K_Indeces(1));
    V_maxv2     = ParameterSet(V_K_Indeces(2));
    V_maxv3     = ParameterSet(V_K_Indeces(3));
    V_maxv4     = ParameterSet(V_K_Indeces(4));
    V_maxv5     = ParameterSet(V_K_Indeces(5));
    V_maxv6     = ParameterSet(V_K_Indeces(6));
    V_maxv7     = ParameterSet(V_K_Indeces(7));
    V_maxv8     = ParameterSet(V_K_Indeces(8));
    V_maxv9     = ParameterSet(V_K_Indeces(9));
    V_maxv10    = ParameterSet(V_K_Indeces(10));
    V_maxv11    = ParameterSet(V_K_Indeces(11));
    V_maxv12    = ParameterSet(V_K_Indeces(12));
    V_maxv13    = ParameterSet(V_K_Indeces(13));
    V_maxv14    = ParameterSet(V_K_Indeces(14));
    V_maxv15    = ParameterSet(V_K_Indeces(15));
    V_maxv16    = ParameterSet(V_K_Indeces(16));
    V_maxv17    = ParameterSet(V_K_Indeces(17));
    V_maxv18    = ParameterSet(V_K_Indeces(18));
    V_maxv19    = ParameterSet(V_K_Indeces(19));
    V_maxv20    = ParameterSet(V_K_Indeces(20));
    V_maxv21    = ParameterSet(V_K_Indeces(21));
    V_maxv22    = ParameterSet(V_K_Indeces(22));
    V_maxv23    = ParameterSet(V_K_Indeces(23));
    V_maxv24    = ParameterSet(V_K_Indeces(24));
    V_maxv25    = ParameterSet(V_K_Indeces(25));
    V_maxv26    = ParameterSet(V_K_Indeces(26));
    V_maxv27    = ParameterSet(V_K_Indeces(27));
    V_maxv28    = ParameterSet(V_K_Indeces(28));
    K_PPoolv29  = ParameterSet(V_K_Indeces(29));
    
    
%% Derivatives for DFODC matrix (Numerical values of the derivatives)
S = size(SFull,1);
R = size(SFull,2);

% Initialize dfodc matrix
dfodc = zeros(S,R);

% Row 1
dfodc(1,2)  = (ATP*P3G*V_maxv2*((ADP*BPG)/(ATP*K_eqv2*P3G) - 1))/(KmATPv2*KmBPGv2*KmP3Gv2*(BPG/KmBPGv2 + P3G/KmP3Gv2 + 1)^2*(ADP/KmADPv2 + ATP/KmATPv2 + 1)) - (ADP*V_maxv2)/(K_eqv2*KmATPv2*KmP3Gv2*(BPG/KmBPGv2 + P3G/KmP3Gv2 + 1)*(ADP/KmADPv2 + ATP/KmATPv2 + 1));
dfodc(1,3)  = (BPG*NADPH*V_maxv3*((GAP*NADP*PHI)/(BPG*K_eqv3*NADPH) - 1))/(KmBPGv3^2*KmNADPHv3*(NADP/KmNADPv3 + NADPH/KmNADPHv3 + 1)*(BPG/KmBPGv3 + GAP/KmGAPv3 + PHI/KmPHIv3 + 1)^2) - (NADPH*V_maxv3*((GAP*NADP*PHI)/(BPG*K_eqv3*NADPH) - 1))/(KmBPGv3*KmNADPHv3*(NADP/KmNADPv3 + NADPH/KmNADPHv3 + 1)*(BPG/KmBPGv3 + GAP/KmGAPv3 + PHI/KmPHIv3 + 1)) + (GAP*NADP*PHI*V_maxv3)/(BPG*K_eqv3*KmBPGv3*KmNADPHv3*(NADP/KmNADPv3 + NADPH/KmNADPHv3 + 1)*(BPG/KmBPGv3 + GAP/KmGAPv3 + PHI/KmPHIv3 + 1));

% Row 2
dfodc(2,14) = (P3G*V_maxv14*(P2G/(K_eqv14*P3G) - 1))/(KmP2Gv14*KmP3Gv14*(P2G/KmP2Gv14 + P3G/KmP3Gv14 + 1)^2) - V_maxv14/(K_eqv14*KmP3Gv14*(P2G/KmP2Gv14 + P3G/KmP3Gv14 + 1));
dfodc(2,15) = (P2G*V_maxv15*(PEP/(K_eqv15*P2G) - 1))/(KmP2Gv15^2*(P2G/KmP2Gv15 + PEP/KmPEPv15 + 1)^2) - (V_maxv15*(PEP/(K_eqv15*P2G) - 1))/(KmP2Gv15*(P2G/KmP2Gv15 + PEP/KmPEPv15 + 1)) + (PEP*V_maxv15)/(K_eqv15*KmP2Gv15*P2G*(P2G/KmP2Gv15 + PEP/KmPEPv15 + 1));

% Row 3
dfodc(3,1)  = 0;
dfodc(3,2)  = (ATP*P3G*V_maxv2*((ADP*BPG)/(ATP*K_eqv2*P3G) - 1))/(KmATPv2*KmP3Gv2^2*(BPG/KmBPGv2 + P3G/KmP3Gv2 + 1)^2*(ADP/KmADPv2 + ATP/KmATPv2 + 1)) - (ATP*V_maxv2*((ADP*BPG)/(ATP*K_eqv2*P3G) - 1))/(KmATPv2*KmP3Gv2*(BPG/KmBPGv2 + P3G/KmP3Gv2 + 1)*(ADP/KmADPv2 + ATP/KmATPv2 + 1)) + (ADP*BPG*V_maxv2)/(K_eqv2*KmATPv2*KmP3Gv2*P3G*(BPG/KmBPGv2 + P3G/KmP3Gv2 + 1)*(ADP/KmADPv2 + ATP/KmATPv2 + 1));
dfodc(3,14) = (P3G*V_maxv14*(P2G/(K_eqv14*P3G) - 1))/(KmP3Gv14^2*(P2G/KmP2Gv14 + P3G/KmP3Gv14 + 1)^2) - (V_maxv14*(P2G/(K_eqv14*P3G) - 1))/(KmP3Gv14*(P2G/KmP2Gv14 + P3G/KmP3Gv14 + 1)) + (P2G*V_maxv14)/(K_eqv14*KmP3Gv14*P3G*(P2G/KmP2Gv14 + P3G/KmP3Gv14 + 1));
dfodc(3,28) = V_maxv28/(KmP3Gv28 + P3G) - (P3G*V_maxv28)/(KmP3Gv28 + P3G)^2;

% Row 4
dfodc(4,2)  = (ATP*P3G*V_maxv2*((ADP*BPG)/(ATP*K_eqv2*P3G) - 1))/(KmADPv2*KmATPv2*KmP3Gv2*(BPG/KmBPGv2 + P3G/KmP3Gv2 + 1)*(ADP/KmADPv2 + ATP/KmATPv2 + 1)^2) - (BPG*V_maxv2)/(K_eqv2*KmATPv2*KmP3Gv2*(BPG/KmBPGv2 + P3G/KmP3Gv2 + 1)*(ADP/KmADPv2 + ATP/KmATPv2 + 1));
dfodc(4,13) = (ATP*KiADPv13*KiPEPv13*Ru5P*V_maxv13*((ADP*RuBP)/(ATP*K_eqv13*Ru5P) - 1))/(KmATPv13*KmRu5Pv13*(ADP + KiADPv13)^2*(KiPEPv13 + PEP)*(Ru5P/KmRu5Pv13 + RuBP/KmRuBPv13 + 1)*(ADP/KmADPv13 + ATP/KmATPv13 + 1)) - (KiADPv13*KiPEPv13*RuBP*V_maxv13)/(K_eqv13*KmATPv13*KmRu5Pv13*(ADP + KiADPv13)*(KiPEPv13 + PEP)*(Ru5P/KmRu5Pv13 + RuBP/KmRuBPv13 + 1)*(ADP/KmADPv13 + ATP/KmATPv13 + 1)) + (ATP*KiADPv13*KiPEPv13*Ru5P*V_maxv13*((ADP*RuBP)/(ATP*K_eqv13*Ru5P) - 1))/(KmADPv13*KmATPv13*KmRu5Pv13*(ADP + KiADPv13)*(KiPEPv13 + PEP)*(Ru5P/KmRu5Pv13 + RuBP/KmRuBPv13 + 1)*(ADP/KmADPv13 + ATP/KmATPv13 + 1)^2);
dfodc(4,16) = (ADP*KiATPv16*KiFBPv16*KiPHIv16*PEP*R5P*V_maxv16*((ATP*PYR)/(ADP*K_eqv16*PEP) - 1))/(KmADPv16^2*KmPEPv16*(ATP + KiATPv16)*(FBP + KiFBPv16)*(KiPHIv16 + PHI)*(KaR5Pv16 + R5P)*(PEP/KmPEPv16 + PYR/KmPYRv16 + 1)*(ADP/KmADPv16 + ATP/KmATPv16 + 1)^2) - (KiATPv16*KiFBPv16*KiPHIv16*PEP*R5P*V_maxv16*((ATP*PYR)/(ADP*K_eqv16*PEP) - 1))/(KmADPv16*KmPEPv16*(ATP + KiATPv16)*(FBP + KiFBPv16)*(KiPHIv16 + PHI)*(KaR5Pv16 + R5P)*(PEP/KmPEPv16 + PYR/KmPYRv16 + 1)*(ADP/KmADPv16 + ATP/KmATPv16 + 1)) + (ATP*KiATPv16*KiFBPv16*KiPHIv16*PYR*R5P*V_maxv16)/(ADP*K_eqv16*KmADPv16*KmPEPv16*(ATP + KiATPv16)*(FBP + KiFBPv16)*(KiPHIv16 + PHI)*(KaR5Pv16 + R5P)*(PEP/KmPEPv16 + PYR/KmPYRv16 + 1)*(ADP/KmADPv16 + ATP/KmATPv16 + 1));
dfodc(4,21) = (PHI*V_maxv21)/(KmADPv21*KmPHIv21*(ADP/KmADPv21 + ATP/KmATPv21 + PHI/KmPHIv21 + (ADP*PHI)/(KmADPv21*KmPHIv21) + 1)) + (V_maxv21*(1/KmADPv21 + PHI/(KmADPv21*KmPHIv21))*(ATP/K_eqv21 - ADP*PHI))/(KmADPv21*KmPHIv21*(ADP/KmADPv21 + ATP/KmATPv21 + PHI/KmPHIv21 + (ADP*PHI)/(KmADPv21*KmPHIv21) + 1)^2);

% Row 5
dfodc(5,4)  = (GAP*V_maxv4*(DHAP/(GAP*K_eqv4) - 1))/(KmDHAPv4*KmGAPv4*(DHAP/KmDHAPv4 + GAP/KmGAPv4 + 1)^2) - V_maxv4/(K_eqv4*KmGAPv4*(DHAP/KmDHAPv4 + GAP/KmGAPv4 + 1));
dfodc(5,5)  = (FBP*V_maxv5)/(DHAP*K_eqv5*KmDHAPv5*KmGAPv5*((FBP/KmFBPv5 + 1)*(SBP/KmSBPv5 + 1) + (DHAP/KmDHAPv5 + 1)*(E4P/KmE4Pv5 + 1)*(GAP/KmGAPv5 + 1) - 1)) - (GAP*V_maxv5*(FBP/(DHAP*GAP*K_eqv5) - 1))/(KmDHAPv5*KmGAPv5*((FBP/KmFBPv5 + 1)*(SBP/KmSBPv5 + 1) + (DHAP/KmDHAPv5 + 1)*(E4P/KmE4Pv5 + 1)*(GAP/KmGAPv5 + 1) - 1)) + (DHAP*GAP*V_maxv5*(E4P/KmE4Pv5 + 1)*(GAP/KmGAPv5 + 1)*(FBP/(DHAP*GAP*K_eqv5) - 1))/(KmDHAPv5^2*KmGAPv5*((FBP/KmFBPv5 + 1)*(SBP/KmSBPv5 + 1) + (DHAP/KmDHAPv5 + 1)*(E4P/KmE4Pv5 + 1)*(GAP/KmGAPv5 + 1) - 1)^2);
dfodc(5,9)  = (SBP*V_maxv9)/(DHAP*K_eqv9*KmDHAPv9*KmE4Pv9*((FBP/KmFBPv9 + 1)*(SBP/KmSBPv9 + 1) + (DHAP/KmDHAPv9 + 1)*(E4P/KmE4Pv9 + 1)*(GAP/KmGAPv9 + 1) - 1)) - (E4P*V_maxv9*(SBP/(DHAP*E4P*K_eqv9) - 1))/(KmDHAPv9*KmE4Pv9*((FBP/KmFBPv9 + 1)*(SBP/KmSBPv9 + 1) + (DHAP/KmDHAPv9 + 1)*(E4P/KmE4Pv9 + 1)*(GAP/KmGAPv9 + 1) - 1)) + (DHAP*E4P*V_maxv9*(SBP/(DHAP*E4P*K_eqv9) - 1)*(E4P/KmE4Pv9 + 1)*(GAP/KmGAPv9 + 1))/(KmDHAPv9^2*KmE4Pv9*((FBP/KmFBPv9 + 1)*(SBP/KmSBPv9 + 1) + (DHAP/KmDHAPv9 + 1)*(E4P/KmE4Pv9 + 1)*(GAP/KmGAPv9 + 1) - 1)^2);

% Row 6
dfodc(6,5)  = (DHAP*GAP*V_maxv5*(DHAP/KmDHAPv5 + 1)*(GAP/KmGAPv5 + 1)*(FBP/(DHAP*GAP*K_eqv5) - 1))/(KmDHAPv5*KmE4Pv5*KmGAPv5*((FBP/KmFBPv5 + 1)*(SBP/KmSBPv5 + 1) + (DHAP/KmDHAPv5 + 1)*(E4P/KmE4Pv5 + 1)*(GAP/KmGAPv5 + 1) - 1)^2);
dfodc(6,6)  = (F6P*GAP*V_maxv6*((E4P*Xu5P)/(F6P*GAP*K_eqv6) - 1))/(KmE4Pv6*KmF6Pv6*KmGAPv6*(GAP/KmGAPv6 + Xu5P/KmXu5Pv6 + 1)*(R5P/KmR5Pv6 + S7P/KmS7Pv6 + 1)*(E4P/KmE4Pv6 + F6P/KmF6Pv6 + 1)^2) - (V_maxv6*Xu5P)/(K_eqv6*KmF6Pv6*KmGAPv6*(GAP/KmGAPv6 + Xu5P/KmXu5Pv6 + 1)*(R5P/KmR5Pv6 + S7P/KmS7Pv6 + 1)*(E4P/KmE4Pv6 + F6P/KmF6Pv6 + 1));
dfodc(6,7)  = (GAP*S7P*V_maxv7*((R5P*Xu5P)/(GAP*K_eqv7*S7P) - 1))/(KmE4Pv7*KmGAPv7*KmS7Pv7*(GAP/KmGAPv7 + Xu5P/KmXu5Pv7 + 1)*(R5P/KmR5Pv7 + S7P/KmS7Pv7 + 1)*(E4P/KmE4Pv7 + F6P/KmF6Pv7 + 1)^2);
dfodc(6,9)  = (SBP*V_maxv9)/(E4P*K_eqv9*KmDHAPv9*KmE4Pv9*((FBP/KmFBPv9 + 1)*(SBP/KmSBPv9 + 1) + (DHAP/KmDHAPv9 + 1)*(E4P/KmE4Pv9 + 1)*(GAP/KmGAPv9 + 1) - 1)) - (DHAP*V_maxv9*(SBP/(DHAP*E4P*K_eqv9) - 1))/(KmDHAPv9*KmE4Pv9*((FBP/KmFBPv9 + 1)*(SBP/KmSBPv9 + 1) + (DHAP/KmDHAPv9 + 1)*(E4P/KmE4Pv9 + 1)*(GAP/KmGAPv9 + 1) - 1)) + (DHAP*E4P*V_maxv9*(SBP/(DHAP*E4P*K_eqv9) - 1)*(DHAP/KmDHAPv9 + 1)*(GAP/KmGAPv9 + 1))/(KmDHAPv9*KmE4Pv9^2*((FBP/KmFBPv9 + 1)*(SBP/KmSBPv9 + 1) + (DHAP/KmDHAPv9 + 1)*(E4P/KmE4Pv9 + 1)*(GAP/KmGAPv9 + 1) - 1)^2);
dfodc(6,18) = (F6P*PHI*V_maxv18*((ACETP*E4P)/(F6P*K_eqv18*PHI) - 1))/(KmE4Pv18*KmF6Pv18*KmPHIv18*(ACETP/KmACETPv18 + PHI/KmPHIv18 + 1)*(GAP/KmGAPv18 + Xu5P/KmXu5Pv18 + 1)*(E4P/KmE4Pv18 + F6P/KmF6Pv18 + 1)^2) - (ACETP*V_maxv18)/(K_eqv18*KmF6Pv18*KmPHIv18*(ACETP/KmACETPv18 + PHI/KmPHIv18 + 1)*(GAP/KmGAPv18 + Xu5P/KmXu5Pv18 + 1)*(E4P/KmE4Pv18 + F6P/KmF6Pv18 + 1));
dfodc(6,19) = (PHI*V_maxv19*Xu5P*((ACETP*GAP)/(K_eqv19*PHI*Xu5P) - 1))/(KmE4Pv19*KmPHIv19*KmXu5Pv19*(ACETP/KmACETPv19 + PHI/KmPHIv19 + 1)*(GAP/KmGAPv19 + Xu5P/KmXu5Pv19 + 1)*(E4P/KmE4Pv19 + F6P/KmF6Pv19 + 1)^2);
dfodc(6,25) = V_maxv25/(E4P + KmE4Pv25) - (E4P*V_maxv25)/(E4P + KmE4Pv25)^2;

% Row 7 
dfodc(7,6)  = (F6P*GAP*V_maxv6*((E4P*Xu5P)/(F6P*GAP*K_eqv6) - 1))/(KmF6Pv6^2*KmGAPv6*(GAP/KmGAPv6 + Xu5P/KmXu5Pv6 + 1)*(R5P/KmR5Pv6 + S7P/KmS7Pv6 + 1)*(E4P/KmE4Pv6 + F6P/KmF6Pv6 + 1)^2) - (GAP*V_maxv6*((E4P*Xu5P)/(F6P*GAP*K_eqv6) - 1))/(KmF6Pv6*KmGAPv6*(GAP/KmGAPv6 + Xu5P/KmXu5Pv6 + 1)*(R5P/KmR5Pv6 + S7P/KmS7Pv6 + 1)*(E4P/KmE4Pv6 + F6P/KmF6Pv6 + 1)) + (E4P*V_maxv6*Xu5P)/(F6P*K_eqv6*KmF6Pv6*KmGAPv6*(GAP/KmGAPv6 + Xu5P/KmXu5Pv6 + 1)*(R5P/KmR5Pv6 + S7P/KmS7Pv6 + 1)*(E4P/KmE4Pv6 + F6P/KmF6Pv6 + 1));
dfodc(7,7)  = (GAP*S7P*V_maxv7*((R5P*Xu5P)/(GAP*K_eqv7*S7P) - 1))/(KmF6Pv7*KmGAPv7*KmS7Pv7*(GAP/KmGAPv7 + Xu5P/KmXu5Pv7 + 1)*(R5P/KmR5Pv7 + S7P/KmS7Pv7 + 1)*(E4P/KmE4Pv7 + F6P/KmF6Pv7 + 1)^2);
dfodc(7,8)  = 0;
dfodc(7,10) = 0;
dfodc(7,18) = (F6P*PHI*V_maxv18*((ACETP*E4P)/(F6P*K_eqv18*PHI) - 1))/(KmF6Pv18^2*KmPHIv18*(ACETP/KmACETPv18 + PHI/KmPHIv18 + 1)*(GAP/KmGAPv18 + Xu5P/KmXu5Pv18 + 1)*(E4P/KmE4Pv18 + F6P/KmF6Pv18 + 1)^2) - (PHI*V_maxv18*((ACETP*E4P)/(F6P*K_eqv18*PHI) - 1))/(KmF6Pv18*KmPHIv18*(ACETP/KmACETPv18 + PHI/KmPHIv18 + 1)*(GAP/KmGAPv18 + Xu5P/KmXu5Pv18 + 1)*(E4P/KmE4Pv18 + F6P/KmF6Pv18 + 1)) + (ACETP*E4P*V_maxv18)/(F6P*K_eqv18*KmF6Pv18*KmPHIv18*(ACETP/KmACETPv18 + PHI/KmPHIv18 + 1)*(GAP/KmGAPv18 + Xu5P/KmXu5Pv18 + 1)*(E4P/KmE4Pv18 + F6P/KmF6Pv18 + 1));
dfodc(7,19) = (PHI*V_maxv19*Xu5P*((ACETP*GAP)/(K_eqv19*PHI*Xu5P) - 1))/(KmF6Pv19*KmPHIv19*KmXu5Pv19*(ACETP/KmACETPv19 + PHI/KmPHIv19 + 1)*(GAP/KmGAPv19 + Xu5P/KmXu5Pv19 + 1)*(E4P/KmE4Pv19 + F6P/KmF6Pv19 + 1)^2);
dfodc(7,24) = V_maxv24/(F6P + KmF6Pv24) - (F6P*V_maxv24)/(F6P + KmF6Pv24)^2;

% Row 8
dfodc(8,5)  = (DHAP*GAP*V_maxv5*(SBP/KmSBPv5 + 1)*(FBP/(DHAP*GAP*K_eqv5) - 1))/(KmDHAPv5*KmFBPv5*KmGAPv5*((FBP/KmFBPv5 + 1)*(SBP/KmSBPv5 + 1) + (DHAP/KmDHAPv5 + 1)*(E4P/KmE4Pv5 + 1)*(GAP/KmGAPv5 + 1) - 1)^2) - V_maxv5/(K_eqv5*KmDHAPv5*KmGAPv5*((FBP/KmFBPv5 + 1)*(SBP/KmSBPv5 + 1) + (DHAP/KmDHAPv5 + 1)*(E4P/KmE4Pv5 + 1)*(GAP/KmGAPv5 + 1) - 1));
dfodc(8,8)  = V_maxv8/(FBP + KmFBPv8*(SBP/KmSBPv8 + 1)) - (FBP*V_maxv8)/(FBP + KmFBPv8*(SBP/KmSBPv8 + 1))^2;
dfodc(8,9)  = (DHAP*E4P*V_maxv9*(SBP/(DHAP*E4P*K_eqv9) - 1)*(SBP/KmSBPv9 + 1))/(KmDHAPv9*KmE4Pv9*KmFBPv9*((FBP/KmFBPv9 + 1)*(SBP/KmSBPv9 + 1) + (DHAP/KmDHAPv9 + 1)*(E4P/KmE4Pv9 + 1)*(GAP/KmGAPv9 + 1) - 1)^2);
dfodc(8,10) = -(KmSBPv10*SBP*V_maxv10)/(KmFBPv10*(SBP + KmSBPv10*(FBP/KmFBPv10 + 1))^2);
dfodc(8,16) = (ADP*KiATPv16*KiFBPv16*KiPHIv16*PEP*R5P*V_maxv16*((ATP*PYR)/(ADP*K_eqv16*PEP) - 1))/(KmADPv16*KmPEPv16*(ATP + KiATPv16)*(FBP + KiFBPv16)^2*(KiPHIv16 + PHI)*(KaR5Pv16 + R5P)*(PEP/KmPEPv16 + PYR/KmPYRv16 + 1)*(ADP/KmADPv16 + ATP/KmATPv16 + 1));

% Row 9
dfodc(9,3)	= (BPG*NADPH*V_maxv3*((GAP*NADP*PHI)/(BPG*K_eqv3*NADPH) - 1))/(KmBPGv3*KmGAPv3*KmNADPHv3*(NADP/KmNADPv3 + NADPH/KmNADPHv3 + 1)*(BPG/KmBPGv3 + GAP/KmGAPv3 + PHI/KmPHIv3 + 1)^2) - (NADP*PHI*V_maxv3)/(K_eqv3*KmBPGv3*KmNADPHv3*(NADP/KmNADPv3 + NADPH/KmNADPHv3 + 1)*(BPG/KmBPGv3 + GAP/KmGAPv3 + PHI/KmPHIv3 + 1));
dfodc(9,4)  = (GAP*V_maxv4*(DHAP/(GAP*K_eqv4) - 1))/(KmGAPv4^2*(DHAP/KmDHAPv4 + GAP/KmGAPv4 + 1)^2) - (V_maxv4*(DHAP/(GAP*K_eqv4) - 1))/(KmGAPv4*(DHAP/KmDHAPv4 + GAP/KmGAPv4 + 1)) + (DHAP*V_maxv4)/(GAP*K_eqv4*KmGAPv4*(DHAP/KmDHAPv4 + GAP/KmGAPv4 + 1));
dfodc(9,5)  = (FBP*V_maxv5)/(GAP*K_eqv5*KmDHAPv5*KmGAPv5*((FBP/KmFBPv5 + 1)*(SBP/KmSBPv5 + 1) + (DHAP/KmDHAPv5 + 1)*(E4P/KmE4Pv5 + 1)*(GAP/KmGAPv5 + 1) - 1)) - (DHAP*V_maxv5*(FBP/(DHAP*GAP*K_eqv5) - 1))/(KmDHAPv5*KmGAPv5*((FBP/KmFBPv5 + 1)*(SBP/KmSBPv5 + 1) + (DHAP/KmDHAPv5 + 1)*(E4P/KmE4Pv5 + 1)*(GAP/KmGAPv5 + 1) - 1)) + (DHAP*GAP*V_maxv5*(DHAP/KmDHAPv5 + 1)*(E4P/KmE4Pv5 + 1)*(FBP/(DHAP*GAP*K_eqv5) - 1))/(KmDHAPv5*KmGAPv5^2*((FBP/KmFBPv5 + 1)*(SBP/KmSBPv5 + 1) + (DHAP/KmDHAPv5 + 1)*(E4P/KmE4Pv5 + 1)*(GAP/KmGAPv5 + 1) - 1)^2);
dfodc(9,6)  = (F6P*GAP*V_maxv6*((E4P*Xu5P)/(F6P*GAP*K_eqv6) - 1))/(KmF6Pv6*KmGAPv6^2*(GAP/KmGAPv6 + Xu5P/KmXu5Pv6 + 1)^2*(R5P/KmR5Pv6 + S7P/KmS7Pv6 + 1)*(E4P/KmE4Pv6 + F6P/KmF6Pv6 + 1)) - (F6P*V_maxv6*((E4P*Xu5P)/(F6P*GAP*K_eqv6) - 1))/(KmF6Pv6*KmGAPv6*(GAP/KmGAPv6 + Xu5P/KmXu5Pv6 + 1)*(R5P/KmR5Pv6 + S7P/KmS7Pv6 + 1)*(E4P/KmE4Pv6 + F6P/KmF6Pv6 + 1)) + (E4P*V_maxv6*Xu5P)/(GAP*K_eqv6*KmF6Pv6*KmGAPv6*(GAP/KmGAPv6 + Xu5P/KmXu5Pv6 + 1)*(R5P/KmR5Pv6 + S7P/KmS7Pv6 + 1)*(E4P/KmE4Pv6 + F6P/KmF6Pv6 + 1));
dfodc(9,7)  = (GAP*S7P*V_maxv7*((R5P*Xu5P)/(GAP*K_eqv7*S7P) - 1))/(KmGAPv7^2*KmS7Pv7*(GAP/KmGAPv7 + Xu5P/KmXu5Pv7 + 1)^2*(R5P/KmR5Pv7 + S7P/KmS7Pv7 + 1)*(E4P/KmE4Pv7 + F6P/KmF6Pv7 + 1)) - (S7P*V_maxv7*((R5P*Xu5P)/(GAP*K_eqv7*S7P) - 1))/(KmGAPv7*KmS7Pv7*(GAP/KmGAPv7 + Xu5P/KmXu5Pv7 + 1)*(R5P/KmR5Pv7 + S7P/KmS7Pv7 + 1)*(E4P/KmE4Pv7 + F6P/KmF6Pv7 + 1)) + (R5P*V_maxv7*Xu5P)/(GAP*K_eqv7*KmGAPv7*KmS7Pv7*(GAP/KmGAPv7 + Xu5P/KmXu5Pv7 + 1)*(R5P/KmR5Pv7 + S7P/KmS7Pv7 + 1)*(E4P/KmE4Pv7 + F6P/KmF6Pv7 + 1));
dfodc(9,9)  = (DHAP*E4P*V_maxv9*(SBP/(DHAP*E4P*K_eqv9) - 1)*(DHAP/KmDHAPv9 + 1)*(E4P/KmE4Pv9 + 1))/(KmDHAPv9*KmE4Pv9*KmGAPv9*((FBP/KmFBPv9 + 1)*(SBP/KmSBPv9 + 1) + (DHAP/KmDHAPv9 + 1)*(E4P/KmE4Pv9 + 1)*(GAP/KmGAPv9 + 1) - 1)^2);
dfodc(9,18) = (F6P*PHI*V_maxv18*((ACETP*E4P)/(F6P*K_eqv18*PHI) - 1))/(KmGAPv18*KmF6Pv18*KmPHIv18*(ACETP/KmACETPv18 + PHI/KmPHIv18 + 1)*(GAP/KmGAPv18 + Xu5P/KmXu5Pv18 + 1)^2*(E4P/KmE4Pv18 + F6P/KmF6Pv18 + 1));
dfodc(9,19) = (PHI*V_maxv19*Xu5P*((ACETP*GAP)/(K_eqv19*PHI*Xu5P) - 1))/(KmGAPv19*KmPHIv19*KmXu5Pv19*(ACETP/KmACETPv19 + PHI/KmPHIv19 + 1)*(GAP/KmGAPv19 + Xu5P/KmXu5Pv19 + 1)^2*(E4P/KmE4Pv19 + F6P/KmF6Pv19 + 1)) - (ACETP*V_maxv19)/(K_eqv19*KmPHIv19*KmXu5Pv19*(ACETP/KmACETPv19 + PHI/KmPHIv19 + 1)*(GAP/KmGAPv19 + Xu5P/KmXu5Pv19 + 1)*(E4P/KmE4Pv19 + F6P/KmF6Pv19 + 1));

% Row 10
dfodc(10,1) = -(CO2_cax*KmRuBPv1*PHI*RuBP*V_maxv1)/(KiNADPHv1*(CO2_cax + KmCO2_caxv1*(O2/KmO2v1 + 1))*(RuBP + KmRuBPv1*(NADPH/KiNADPHv1 + PHI/KiPHIv1 + 1))^2*(KaPHIv1 + PHI));
dfodc(10,3) = (BPG*NADPH*V_maxv3*((GAP*NADP*PHI)/(BPG*K_eqv3*NADPH) - 1))/(KmBPGv3*KmNADPHv3^2*(NADP/KmNADPv3 + NADPH/KmNADPHv3 + 1)^2*(BPG/KmBPGv3 + GAP/KmGAPv3 + PHI/KmPHIv3 + 1)) - (BPG*V_maxv3*((GAP*NADP*PHI)/(BPG*K_eqv3*NADPH) - 1))/(KmBPGv3*KmNADPHv3*(NADP/KmNADPv3 + NADPH/KmNADPHv3 + 1)*(BPG/KmBPGv3 + GAP/KmGAPv3 + PHI/KmPHIv3 + 1)) + (GAP*NADP*PHI*V_maxv3)/(K_eqv3*KmBPGv3*KmNADPHv3*NADPH*(NADP/KmNADPv3 + NADPH/KmNADPHv3 + 1)*(BPG/KmBPGv3 + GAP/KmGAPv3 + PHI/KmPHIv3 + 1));
dfodc(10,22)= (NADP*V_maxv22*(NADPH/(K_eqv22*NADP) - 1))/(KmNADPv22*KmNADPHv22*(NADP/KmNADPv22 + NADPH/KmNADPHv22 + 1)^2) - V_maxv22/(K_eqv22*KmNADPv22*(NADP/KmNADPv22 + NADPH/KmNADPHv22 + 1));

% Row 11
dfodc(11,13)= (ATP*KiADPv13*KiPEPv13*Ru5P*V_maxv13*((ADP*RuBP)/(ATP*K_eqv13*Ru5P) - 1))/(KmATPv13*KmRu5Pv13*(ADP + KiADPv13)*(KiPEPv13 + PEP)^2*(Ru5P/KmRu5Pv13 + RuBP/KmRuBPv13 + 1)*(ADP/KmADPv13 + ATP/KmATPv13 + 1));
dfodc(11,15)= (P2G*V_maxv15*(PEP/(K_eqv15*P2G) - 1))/(KmPEPv15*KmP2Gv15*(P2G/KmP2Gv15 + PEP/KmPEPv15 + 1)^2) - V_maxv15/(K_eqv15*KmP2Gv15*(P2G/KmP2Gv15 + PEP/KmPEPv15 + 1));
dfodc(11,16)= (ADP*KiATPv16*KiFBPv16*KiPHIv16*PEP*R5P*V_maxv16*((ATP*PYR)/(ADP*K_eqv16*PEP) - 1))/(KmADPv16*KmPEPv16^2*(ATP + KiATPv16)*(FBP + KiFBPv16)*(KiPHIv16 + PHI)*(KaR5Pv16 + R5P)*(PEP/KmPEPv16 + PYR/KmPYRv16 + 1)^2*(ADP/KmADPv16 + ATP/KmATPv16 + 1)) - (ADP*KiATPv16*KiFBPv16*KiPHIv16*R5P*V_maxv16*((ATP*PYR)/(ADP*K_eqv16*PEP) - 1))/(KmADPv16*KmPEPv16*(ATP + KiATPv16)*(FBP + KiFBPv16)*(KiPHIv16 + PHI)*(KaR5Pv16 + R5P)*(PEP/KmPEPv16 + PYR/KmPYRv16 + 1)*(ADP/KmADPv16 + ATP/KmATPv16 + 1)) + (ATP*KiATPv16*KiFBPv16*KiPHIv16*PYR*R5P*V_maxv16)/(K_eqv16*KmADPv16*KmPEPv16*PEP*(ATP + KiATPv16)*(FBP + KiFBPv16)*(KiPHIv16 + PHI)*(KaR5Pv16 + R5P)*(PEP/KmPEPv16 + PYR/KmPYRv16 + 1)*(ADP/KmADPv16 + ATP/KmATPv16 + 1));
dfodc(11,26)= V_maxv26/(KmPEPv26 + PEP) - (PEP*V_maxv26)/(KmPEPv26 + PEP)^2;

% Row 12
dfodc(12,1) = (CO2_cax*RuBP*V_maxv1)/((CO2_cax + KmCO2_caxv1*(O2/KmO2v1 + 1))*(RuBP + KmRuBPv1*(NADPH/KiNADPHv1 + PHI/KiPHIv1 + 1))*(KaPHIv1 + PHI)) - (CO2_cax*PHI*RuBP*V_maxv1)/((CO2_cax + KmCO2_caxv1*(O2/KmO2v1 + 1))*(RuBP + KmRuBPv1*(NADPH/KiNADPHv1 + PHI/KiPHIv1 + 1))*(KaPHIv1 + PHI)^2) - (CO2_cax*KmRuBPv1*PHI*RuBP*V_maxv1)/(KiPHIv1*(CO2_cax + KmCO2_caxv1*(O2/KmO2v1 + 1))*(RuBP + KmRuBPv1*(NADPH/KiNADPHv1 + PHI/KiPHIv1 + 1))^2*(KaPHIv1 + PHI));
dfodc(12,3) = (BPG*NADPH*V_maxv3*((GAP*NADP*PHI)/(BPG*K_eqv3*NADPH) - 1))/(KmBPGv3*KmNADPHv3*KmPHIv3*(NADP/KmNADPv3 + NADPH/KmNADPHv3 + 1)*(BPG/KmBPGv3 + GAP/KmGAPv3 + PHI/KmPHIv3 + 1)^2) - (GAP*NADP*V_maxv3)/(K_eqv3*KmBPGv3*KmNADPHv3*(NADP/KmNADPv3 + NADPH/KmNADPHv3 + 1)*(BPG/KmBPGv3 + GAP/KmGAPv3 + PHI/KmPHIv3 + 1));
dfodc(12,8) = 0;
dfodc(12,10)= 0;
dfodc(12,16)= (ADP*KiATPv16*KiFBPv16*KiPHIv16*PEP*R5P*V_maxv16*((ATP*PYR)/(ADP*K_eqv16*PEP) - 1))/(KmADPv16*KmPEPv16*(ATP + KiATPv16)*(FBP + KiFBPv16)*(KiPHIv16 + PHI)^2*(KaR5Pv16 + R5P)*(PEP/KmPEPv16 + PYR/KmPYRv16 + 1)*(ADP/KmADPv16 + ATP/KmATPv16 + 1));
dfodc(12,18)= (F6P*PHI*V_maxv18*((ACETP*E4P)/(F6P*K_eqv18*PHI) - 1))/(KmF6Pv18*KmPHIv18^2*(ACETP/KmACETPv18 + PHI/KmPHIv18 + 1)^2*(GAP/KmGAPv18 + Xu5P/KmXu5Pv18 + 1)*(E4P/KmE4Pv18 + F6P/KmF6Pv18 + 1)) - (F6P*V_maxv18*((ACETP*E4P)/(F6P*K_eqv18*PHI) - 1))/(KmF6Pv18*KmPHIv18*(ACETP/KmACETPv18 + PHI/KmPHIv18 + 1)*(GAP/KmGAPv18 + Xu5P/KmXu5Pv18 + 1)*(E4P/KmE4Pv18 + F6P/KmF6Pv18 + 1)) + (ACETP*E4P*V_maxv18)/(K_eqv18*KmF6Pv18*KmPHIv18*PHI*(ACETP/KmACETPv18 + PHI/KmPHIv18 + 1)*(GAP/KmGAPv18 + Xu5P/KmXu5Pv18 + 1)*(E4P/KmE4Pv18 + F6P/KmF6Pv18 + 1));
dfodc(12,19)= (PHI*V_maxv19*Xu5P*((ACETP*GAP)/(K_eqv19*PHI*Xu5P) - 1))/(KmPHIv19^2*KmXu5Pv19*(ACETP/KmACETPv19 + PHI/KmPHIv19 + 1)^2*(GAP/KmGAPv19 + Xu5P/KmXu5Pv19 + 1)*(E4P/KmE4Pv19 + F6P/KmF6Pv19 + 1)) - (V_maxv19*Xu5P*((ACETP*GAP)/(K_eqv19*PHI*Xu5P) - 1))/(KmPHIv19*KmXu5Pv19*(ACETP/KmACETPv19 + PHI/KmPHIv19 + 1)*(GAP/KmGAPv19 + Xu5P/KmXu5Pv19 + 1)*(E4P/KmE4Pv19 + F6P/KmF6Pv19 + 1)) + (ACETP*GAP*V_maxv19)/(K_eqv19*KmPHIv19*KmXu5Pv19*PHI*(ACETP/KmACETPv19 + PHI/KmPHIv19 + 1)*(GAP/KmGAPv19 + Xu5P/KmXu5Pv19 + 1)*(E4P/KmE4Pv19 + F6P/KmF6Pv19 + 1));
dfodc(12,20)= (ACETP*COA*V_maxv20*((ACCOA*PHI)/(ACETP*COA*K_eqv20) - 1))/(KmACETPv20*KmCOAv20*KmPHIv20*(COA/KmCOAv20 + PHI/KmPHIv20 + 1)^2*(ACCOA/KmACCOAv20 + ACETP/KmACETPv20 + 1)) - (ACCOA*V_maxv20)/(K_eqv20*KmACETPv20*KmCOAv20*(COA/KmCOAv20 + PHI/KmPHIv20 + 1)*(ACCOA/KmACCOAv20 + ACETP/KmACETPv20 + 1));
dfodc(12,21)= (ADP*V_maxv21)/(KmADPv21*KmPHIv21*(ADP/KmADPv21 + ATP/KmATPv21 + PHI/KmPHIv21 + (ADP*PHI)/(KmADPv21*KmPHIv21) + 1)) + (V_maxv21*(1/KmPHIv21 + ADP/(KmADPv21*KmPHIv21))*(ATP/K_eqv21 - ADP*PHI))/(KmADPv21*KmPHIv21*(ADP/KmADPv21 + ATP/KmATPv21 + PHI/KmPHIv21 + (ADP*PHI)/(KmADPv21*KmPHIv21) + 1)^2);
dfodc(12,29)= -K_PPoolv29;

% Row 13
dfodc(13,6) = (F6P*GAP*V_maxv6*((E4P*Xu5P)/(F6P*GAP*K_eqv6) - 1))/(KmF6Pv6*KmGAPv6*KmR5Pv6*(GAP/KmGAPv6 + Xu5P/KmXu5Pv6 + 1)*(R5P/KmR5Pv6 + S7P/KmS7Pv6 + 1)^2*(E4P/KmE4Pv6 + F6P/KmF6Pv6 + 1));
dfodc(13,7) = (GAP*S7P*V_maxv7*((R5P*Xu5P)/(GAP*K_eqv7*S7P) - 1))/(KmGAPv7*KmR5Pv7*KmS7Pv7*(GAP/KmGAPv7 + Xu5P/KmXu5Pv7 + 1)*(R5P/KmR5Pv7 + S7P/KmS7Pv7 + 1)^2*(E4P/KmE4Pv7 + F6P/KmF6Pv7 + 1)) - (V_maxv7*Xu5P)/(K_eqv7*KmGAPv7*KmS7Pv7*(GAP/KmGAPv7 + Xu5P/KmXu5Pv7 + 1)*(R5P/KmR5Pv7 + S7P/KmS7Pv7 + 1)*(E4P/KmE4Pv7 + F6P/KmF6Pv7 + 1));
dfodc(13,11)= (R5P*V_maxv11*(Ru5P/(K_eqv11*R5P) - 1))/(KmR5Pv11^2*(R5P/KmR5Pv11 + Ru5P/KmRu5Pv11 + 1)^2) - (V_maxv11*(Ru5P/(K_eqv11*R5P) - 1))/(KmR5Pv11*(R5P/KmR5Pv11 + Ru5P/KmRu5Pv11 + 1)) + (Ru5P*V_maxv11)/(K_eqv11*KmR5Pv11*R5P*(R5P/KmR5Pv11 + Ru5P/KmRu5Pv11 + 1));
dfodc(13,16)= (ADP*KiATPv16*KiFBPv16*KiPHIv16*PEP*R5P*V_maxv16*((ATP*PYR)/(ADP*K_eqv16*PEP) - 1))/(KmADPv16*KmPEPv16*(ATP + KiATPv16)*(FBP + KiFBPv16)*(KiPHIv16 + PHI)*(KaR5Pv16 + R5P)^2*(PEP/KmPEPv16 + PYR/KmPYRv16 + 1)*(ADP/KmADPv16 + ATP/KmATPv16 + 1)) - (ADP*KiATPv16*KiFBPv16*KiPHIv16*PEP*V_maxv16*((ATP*PYR)/(ADP*K_eqv16*PEP) - 1))/(KmADPv16*KmPEPv16*(ATP + KiATPv16)*(FBP + KiFBPv16)*(KiPHIv16 + PHI)*(KaR5Pv16 + R5P)*(PEP/KmPEPv16 + PYR/KmPYRv16 + 1)*(ADP/KmADPv16 + ATP/KmATPv16 + 1));
dfodc(13,23)= V_maxv23/(KmR5Pv23 + R5P) - (R5P*V_maxv23)/(KmR5Pv23 + R5P)^2;

% Row 14
dfodc(14,11)= (R5P*V_maxv11*(Ru5P/(K_eqv11*R5P) - 1))/(KmR5Pv11*KmRu5Pv11*(R5P/KmR5Pv11 + Ru5P/KmRu5Pv11 + 1)^2) - V_maxv11/(K_eqv11*KmR5Pv11*(R5P/KmR5Pv11 + Ru5P/KmRu5Pv11 + 1));
dfodc(14,12)= (V_maxv12*Xu5P*(Ru5P/(K_eqv12*Xu5P) - 1))/(KmRu5Pv12*KmXu5Pv12*(Ru5P/KmRu5Pv12 + Xu5P/KmXu5Pv12 + 1)^2) - V_maxv12/(K_eqv12*KmXu5Pv12*(Ru5P/KmRu5Pv12 + Xu5P/KmXu5Pv12 + 1));
dfodc(14,13)= (ATP*KiADPv13*KiPEPv13*Ru5P*V_maxv13*((ADP*RuBP)/(ATP*K_eqv13*Ru5P) - 1))/(KmATPv13*KmRu5Pv13^2*(ADP + KiADPv13)*(KiPEPv13 + PEP)*(Ru5P/KmRu5Pv13 + RuBP/KmRuBPv13 + 1)^2*(ADP/KmADPv13 + ATP/KmATPv13 + 1)) - (ATP*KiADPv13*KiPEPv13*V_maxv13*((ADP*RuBP)/(ATP*K_eqv13*Ru5P) - 1))/(KmATPv13*KmRu5Pv13*(ADP + KiADPv13)*(KiPEPv13 + PEP)*(Ru5P/KmRu5Pv13 + RuBP/KmRuBPv13 + 1)*(ADP/KmADPv13 + ATP/KmATPv13 + 1)) + (ADP*KiADPv13*KiPEPv13*RuBP*V_maxv13)/(K_eqv13*KmATPv13*KmRu5Pv13*Ru5P*(ADP + KiADPv13)*(KiPEPv13 + PEP)*(Ru5P/KmRu5Pv13 + RuBP/KmRuBPv13 + 1)*(ADP/KmADPv13 + ATP/KmATPv13 + 1));

% Row 15
dfodc(15,1) = (CO2_cax*PHI*V_maxv1)/((CO2_cax + KmCO2_caxv1*(O2/KmO2v1 + 1))*(RuBP + KmRuBPv1*(NADPH/KiNADPHv1 + PHI/KiPHIv1 + 1))*(KaPHIv1 + PHI)) - (CO2_cax*PHI*RuBP*V_maxv1)/((CO2_cax + KmCO2_caxv1*(O2/KmO2v1 + 1))*(RuBP + KmRuBPv1*(NADPH/KiNADPHv1 + PHI/KiPHIv1 + 1))^2*(KaPHIv1 + PHI));
dfodc(15,13)= (ATP*KiADPv13*KiPEPv13*Ru5P*V_maxv13*((ADP*RuBP)/(ATP*K_eqv13*Ru5P) - 1))/(KmATPv13*KmRuBPv13*KmRu5Pv13*(ADP + KiADPv13)*(KiPEPv13 + PEP)*(Ru5P/KmRu5Pv13 + RuBP/KmRuBPv13 + 1)^2*(ADP/KmADPv13 + ATP/KmATPv13 + 1)) - (ADP*KiADPv13*KiPEPv13*V_maxv13)/(K_eqv13*KmATPv13*KmRu5Pv13*(ADP + KiADPv13)*(KiPEPv13 + PEP)*(Ru5P/KmRu5Pv13 + RuBP/KmRuBPv13 + 1)*(ADP/KmADPv13 + ATP/KmATPv13 + 1));

% Row 16
dfodc(16,6) = (F6P*GAP*V_maxv6*((E4P*Xu5P)/(F6P*GAP*K_eqv6) - 1))/(KmF6Pv6*KmGAPv6*KmS7Pv6*(GAP/KmGAPv6 + Xu5P/KmXu5Pv6 + 1)*(R5P/KmR5Pv6 + S7P/KmS7Pv6 + 1)^2*(E4P/KmE4Pv6 + F6P/KmF6Pv6 + 1));
dfodc(16,7) = (GAP*S7P*V_maxv7*((R5P*Xu5P)/(GAP*K_eqv7*S7P) - 1))/(KmGAPv7*KmS7Pv7^2*(GAP/KmGAPv7 + Xu5P/KmXu5Pv7 + 1)*(R5P/KmR5Pv7 + S7P/KmS7Pv7 + 1)^2*(E4P/KmE4Pv7 + F6P/KmF6Pv7 + 1)) - (GAP*V_maxv7*((R5P*Xu5P)/(GAP*K_eqv7*S7P) - 1))/(KmGAPv7*KmS7Pv7*(GAP/KmGAPv7 + Xu5P/KmXu5Pv7 + 1)*(R5P/KmR5Pv7 + S7P/KmS7Pv7 + 1)*(E4P/KmE4Pv7 + F6P/KmF6Pv7 + 1)) + (R5P*V_maxv7*Xu5P)/(K_eqv7*KmGAPv7*KmS7Pv7*S7P*(GAP/KmGAPv7 + Xu5P/KmXu5Pv7 + 1)*(R5P/KmR5Pv7 + S7P/KmS7Pv7 + 1)*(E4P/KmE4Pv7 + F6P/KmF6Pv7 + 1));
dfodc(16,8) = 0;
dfodc(16,10)= 0;

% Row 17
dfodc(17,5) = (DHAP*GAP*V_maxv5*(FBP/KmFBPv5 + 1)*(FBP/(DHAP*GAP*K_eqv5) - 1))/(KmDHAPv5*KmGAPv5*KmSBPv5*((FBP/KmFBPv5 + 1)*(SBP/KmSBPv5 + 1) + (DHAP/KmDHAPv5 + 1)*(E4P/KmE4Pv5 + 1)*(GAP/KmGAPv5 + 1) - 1)^2);
dfodc(17,8) = -(FBP*KmFBPv8*V_maxv8)/(KmSBPv8*(FBP + KmFBPv8*(SBP/KmSBPv8 + 1))^2);
dfodc(17,9) = (DHAP*E4P*V_maxv9*(SBP/(DHAP*E4P*K_eqv9) - 1)*(FBP/KmFBPv9 + 1))/(KmDHAPv9*KmE4Pv9*KmSBPv9*((FBP/KmFBPv9 + 1)*(SBP/KmSBPv9 + 1) + (DHAP/KmDHAPv9 + 1)*(E4P/KmE4Pv9 + 1)*(GAP/KmGAPv9 + 1) - 1)^2) - V_maxv9/(K_eqv9*KmDHAPv9*KmE4Pv9*((FBP/KmFBPv9 + 1)*(SBP/KmSBPv9 + 1) + (DHAP/KmDHAPv9 + 1)*(E4P/KmE4Pv9 + 1)*(GAP/KmGAPv9 + 1) - 1));
dfodc(17,10)= V_maxv10/(SBP + KmSBPv10*(FBP/KmFBPv10 + 1)) - (SBP*V_maxv10)/(SBP + KmSBPv10*(FBP/KmFBPv10 + 1))^2;

% Row 18
dfodc(18,6) = (F6P*GAP*V_maxv6*((E4P*Xu5P)/(F6P*GAP*K_eqv6) - 1))/(KmF6Pv6*KmGAPv6*KmXu5Pv6*(GAP/KmGAPv6 + Xu5P/KmXu5Pv6 + 1)^2*(R5P/KmR5Pv6 + S7P/KmS7Pv6 + 1)*(E4P/KmE4Pv6 + F6P/KmF6Pv6 + 1)) - (E4P*V_maxv6)/(K_eqv6*KmF6Pv6*KmGAPv6*(GAP/KmGAPv6 + Xu5P/KmXu5Pv6 + 1)*(R5P/KmR5Pv6 + S7P/KmS7Pv6 + 1)*(E4P/KmE4Pv6 + F6P/KmF6Pv6 + 1));
dfodc(18,7) = (GAP*S7P*V_maxv7*((R5P*Xu5P)/(GAP*K_eqv7*S7P) - 1))/(KmGAPv7*KmS7Pv7*KmXu5Pv7*(GAP/KmGAPv7 + Xu5P/KmXu5Pv7 + 1)^2*(R5P/KmR5Pv7 + S7P/KmS7Pv7 + 1)*(E4P/KmE4Pv7 + F6P/KmF6Pv7 + 1)) - (R5P*V_maxv7)/(K_eqv7*KmGAPv7*KmS7Pv7*(GAP/KmGAPv7 + Xu5P/KmXu5Pv7 + 1)*(R5P/KmR5Pv7 + S7P/KmS7Pv7 + 1)*(E4P/KmE4Pv7 + F6P/KmF6Pv7 + 1));
dfodc(18,12)= (V_maxv12*Xu5P*(Ru5P/(K_eqv12*Xu5P) - 1))/(KmXu5Pv12^2*(Ru5P/KmRu5Pv12 + Xu5P/KmXu5Pv12 + 1)^2) - (V_maxv12*(Ru5P/(K_eqv12*Xu5P) - 1))/(KmXu5Pv12*(Ru5P/KmRu5Pv12 + Xu5P/KmXu5Pv12 + 1)) + (Ru5P*V_maxv12)/(K_eqv12*KmXu5Pv12*Xu5P*(Ru5P/KmRu5Pv12 + Xu5P/KmXu5Pv12 + 1));
dfodc(18,18)= (F6P*PHI*V_maxv18*((ACETP*E4P)/(F6P*K_eqv18*PHI) - 1))/(KmF6Pv18*KmPHIv18*KmXu5Pv18*(ACETP/KmACETPv18 + PHI/KmPHIv18 + 1)*(GAP/KmGAPv18 + Xu5P/KmXu5Pv18 + 1)^2*(E4P/KmE4Pv18 + F6P/KmF6Pv18 + 1));
dfodc(18,19)= (PHI*V_maxv19*Xu5P*((ACETP*GAP)/(K_eqv19*PHI*Xu5P) - 1))/(KmPHIv19*KmXu5Pv19^2*(ACETP/KmACETPv19 + PHI/KmPHIv19 + 1)*(GAP/KmGAPv19 + Xu5P/KmXu5Pv19 + 1)^2*(E4P/KmE4Pv19 + F6P/KmF6Pv19 + 1)) - (PHI*V_maxv19*((ACETP*GAP)/(K_eqv19*PHI*Xu5P) - 1))/(KmPHIv19*KmXu5Pv19*(ACETP/KmACETPv19 + PHI/KmPHIv19 + 1)*(GAP/KmGAPv19 + Xu5P/KmXu5Pv19 + 1)*(E4P/KmE4Pv19 + F6P/KmF6Pv19 + 1)) + (ACETP*GAP*V_maxv19)/(K_eqv19*KmPHIv19*KmXu5Pv19*Xu5P*(ACETP/KmACETPv19 + PHI/KmPHIv19 + 1)*(GAP/KmGAPv19 + Xu5P/KmXu5Pv19 + 1)*(E4P/KmE4Pv19 + F6P/KmF6Pv19 + 1));

% Row 19
dfodc(19,16)= (ADP*KiATPv16*KiFBPv16*KiPHIv16*PEP*R5P*V_maxv16*((ATP*PYR)/(ADP*K_eqv16*PEP) - 1))/(KmADPv16*KmPEPv16*KmPYRv16*(ATP + KiATPv16)*(FBP + KiFBPv16)*(KiPHIv16 + PHI)*(KaR5Pv16 + R5P)*(PEP/KmPEPv16 + PYR/KmPYRv16 + 1)^2*(ADP/KmADPv16 + ATP/KmATPv16 + 1)) - (ATP*KiATPv16*KiFBPv16*KiPHIv16*R5P*V_maxv16)/(K_eqv16*KmADPv16*KmPEPv16*(ATP + KiATPv16)*(FBP + KiFBPv16)*(KiPHIv16 + PHI)*(KaR5Pv16 + R5P)*(PEP/KmPEPv16 + PYR/KmPYRv16 + 1)*(ADP/KmADPv16 + ATP/KmATPv16 + 1));
dfodc(19,17)= (COA*NAD*PYR*V_maxv17*((ACCOA*CO2_cyt*NADH)/(COA*K_eqv17*NAD*PYR) - 1))/(KmCOAv17*KmNADv17*KmPYRv17^2*(ACCOA/KmACCOAv17 + COA/KmCOAv17)*(CO2_cyt/KmCO2_cytv17 + PYR/KmPYRv17 + 1)^2*(NAD/KmNADv17 + NADH/KmNADHv17 + 1)) - (COA*NAD*V_maxv17*((ACCOA*CO2_cyt*NADH)/(COA*K_eqv17*NAD*PYR) - 1))/(KmCOAv17*KmNADv17*KmPYRv17*(ACCOA/KmACCOAv17 + COA/KmCOAv17)*(CO2_cyt/KmCO2_cytv17 + PYR/KmPYRv17 + 1)*(NAD/KmNADv17 + NADH/KmNADHv17 + 1)) + (ACCOA*CO2_cyt*NADH*V_maxv17)/(K_eqv17*KmCOAv17*KmNADv17*KmPYRv17*PYR*(ACCOA/KmACCOAv17 + COA/KmCOAv17)*(CO2_cyt/KmCO2_cytv17 + PYR/KmPYRv17 + 1)*(NAD/KmNADv17 + NADH/KmNADHv17 + 1));
dfodc(19,27)= V_maxv27/(KmPYRv27 + PYR) - (PYR*V_maxv27)/(KmPYRv27 + PYR)^2;

% Row 20
dfodc(20,18)= (F6P*PHI*V_maxv18*((ACETP*E4P)/(F6P*K_eqv18*PHI) - 1))/(KmACETPv18*KmF6Pv18*KmPHIv18*(ACETP/KmACETPv18 + PHI/KmPHIv18 + 1)^2*(GAP/KmGAPv18 + Xu5P/KmXu5Pv18 + 1)*(E4P/KmE4Pv18 + F6P/KmF6Pv18 + 1)) - (E4P*V_maxv18)/(K_eqv18*KmF6Pv18*KmPHIv18*(ACETP/KmACETPv18 + PHI/KmPHIv18 + 1)*(GAP/KmGAPv18 + Xu5P/KmXu5Pv18 + 1)*(E4P/KmE4Pv18 + F6P/KmF6Pv18 + 1));
dfodc(20,19)= (PHI*V_maxv19*Xu5P*((ACETP*GAP)/(K_eqv19*PHI*Xu5P) - 1))/(KmACETPv19*KmPHIv19*KmXu5Pv19*(ACETP/KmACETPv19 + PHI/KmPHIv19 + 1)^2*(GAP/KmGAPv19 + Xu5P/KmXu5Pv19 + 1)*(E4P/KmE4Pv19 + F6P/KmF6Pv19 + 1)) - (GAP*V_maxv19)/(K_eqv19*KmPHIv19*KmXu5Pv19*(ACETP/KmACETPv19 + PHI/KmPHIv19 + 1)*(GAP/KmGAPv19 + Xu5P/KmXu5Pv19 + 1)*(E4P/KmE4Pv19 + F6P/KmF6Pv19 + 1));
dfodc(20,20)= (ACETP*COA*V_maxv20*((ACCOA*PHI)/(ACETP*COA*K_eqv20) - 1))/(KmACETPv20^2*KmCOAv20*(COA/KmCOAv20 + PHI/KmPHIv20 + 1)*(ACCOA/KmACCOAv20 + ACETP/KmACETPv20 + 1)^2) - (COA*V_maxv20*((ACCOA*PHI)/(ACETP*COA*K_eqv20) - 1))/(KmACETPv20*KmCOAv20*(COA/KmCOAv20 + PHI/KmPHIv20 + 1)*(ACCOA/KmACCOAv20 + ACETP/KmACETPv20 + 1)) + (ACCOA*PHI*V_maxv20)/(ACETP*K_eqv20*KmACETPv20*KmCOAv20*(COA/KmCOAv20 + PHI/KmPHIv20 + 1)*(ACCOA/KmACCOAv20 + ACETP/KmACETPv20 + 1));

% Row 21
dfodc(21,2) = (ATP*P3G*V_maxv2*((ADP*BPG)/(ATP*K_eqv2*P3G) - 1))/(KmATPv2^2*KmP3Gv2*(BPG/KmBPGv2 + P3G/KmP3Gv2 + 1)*(ADP/KmADPv2 + ATP/KmATPv2 + 1)^2) - (P3G*V_maxv2*((ADP*BPG)/(ATP*K_eqv2*P3G) - 1))/(KmATPv2*KmP3Gv2*(BPG/KmBPGv2 + P3G/KmP3Gv2 + 1)*(ADP/KmADPv2 + ATP/KmATPv2 + 1)) + (ADP*BPG*V_maxv2)/(ATP*K_eqv2*KmATPv2*KmP3Gv2*(BPG/KmBPGv2 + P3G/KmP3Gv2 + 1)*(ADP/KmADPv2 + ATP/KmATPv2 + 1));
dfodc(21,13)= (ATP*KiADPv13*KiPEPv13*Ru5P*V_maxv13*((ADP*RuBP)/(ATP*K_eqv13*Ru5P) - 1))/(KmATPv13^2*KmRu5Pv13*(ADP + KiADPv13)*(KiPEPv13 + PEP)*(Ru5P/KmRu5Pv13 + RuBP/KmRuBPv13 + 1)*(ADP/KmADPv13 + ATP/KmATPv13 + 1)^2) - (KiADPv13*KiPEPv13*Ru5P*V_maxv13*((ADP*RuBP)/(ATP*K_eqv13*Ru5P) - 1))/(KmATPv13*KmRu5Pv13*(ADP + KiADPv13)*(KiPEPv13 + PEP)*(Ru5P/KmRu5Pv13 + RuBP/KmRuBPv13 + 1)*(ADP/KmADPv13 + ATP/KmATPv13 + 1)) + (ADP*KiADPv13*KiPEPv13*RuBP*V_maxv13)/(ATP*K_eqv13*KmATPv13*KmRu5Pv13*(ADP + KiADPv13)*(KiPEPv13 + PEP)*(Ru5P/KmRu5Pv13 + RuBP/KmRuBPv13 + 1)*(ADP/KmADPv13 + ATP/KmATPv13 + 1));
dfodc(21,16)= (ADP*KiATPv16*KiFBPv16*KiPHIv16*PEP*R5P*V_maxv16*((ATP*PYR)/(ADP*K_eqv16*PEP) - 1))/(KmADPv16*KmPEPv16*(ATP + KiATPv16)^2*(FBP + KiFBPv16)*(KiPHIv16 + PHI)*(KaR5Pv16 + R5P)*(PEP/KmPEPv16 + PYR/KmPYRv16 + 1)*(ADP/KmADPv16 + ATP/KmATPv16 + 1)) - (KiATPv16*KiFBPv16*KiPHIv16*PYR*R5P*V_maxv16)/(K_eqv16*KmADPv16*KmPEPv16*(ATP + KiATPv16)*(FBP + KiFBPv16)*(KiPHIv16 + PHI)*(KaR5Pv16 + R5P)*(PEP/KmPEPv16 + PYR/KmPYRv16 + 1)*(ADP/KmADPv16 + ATP/KmATPv16 + 1)) + (ADP*KiATPv16*KiFBPv16*KiPHIv16*PEP*R5P*V_maxv16*((ATP*PYR)/(ADP*K_eqv16*PEP) - 1))/(KmADPv16*KmATPv16*KmPEPv16*(ATP + KiATPv16)*(FBP + KiFBPv16)*(KiPHIv16 + PHI)*(KaR5Pv16 + R5P)*(PEP/KmPEPv16 + PYR/KmPYRv16 + 1)*(ADP/KmADPv16 + ATP/KmATPv16 + 1)^2);
dfodc(21,21)= (V_maxv21*(ATP/K_eqv21 - ADP*PHI))/(KmADPv21*KmATPv21*KmPHIv21*(ADP/KmADPv21 + ATP/KmATPv21 + PHI/KmPHIv21 + (ADP*PHI)/(KmADPv21*KmPHIv21) + 1)^2) - V_maxv21/(K_eqv21*KmADPv21*KmPHIv21*(ADP/KmADPv21 + ATP/KmATPv21 + PHI/KmPHIv21 + (ADP*PHI)/(KmADPv21*KmPHIv21) + 1));

% Row 22
dfodc(22,3) = (BPG*NADPH*V_maxv3*((GAP*NADP*PHI)/(BPG*K_eqv3*NADPH) - 1))/(KmBPGv3*KmNADPv3*KmNADPHv3*(NADP/KmNADPv3 + NADPH/KmNADPHv3 + 1)^2*(BPG/KmBPGv3 + GAP/KmGAPv3 + PHI/KmPHIv3 + 1)) - (GAP*PHI*V_maxv3)/(K_eqv3*KmBPGv3*KmNADPHv3*(NADP/KmNADPv3 + NADPH/KmNADPHv3 + 1)*(BPG/KmBPGv3 + GAP/KmGAPv3 + PHI/KmPHIv3 + 1));
dfodc(22,22)= (NADP*V_maxv22*(NADPH/(K_eqv22*NADP) - 1))/(KmNADPv22^2*(NADP/KmNADPv22 + NADPH/KmNADPHv22 + 1)^2) - (V_maxv22*(NADPH/(K_eqv22*NADP) - 1))/(KmNADPv22*(NADP/KmNADPv22 + NADPH/KmNADPHv22 + 1)) + (NADPH*V_maxv22)/(K_eqv22*KmNADPv22*NADP*(NADP/KmNADPv22 + NADPH/KmNADPHv22 + 1));

dfodc = transpose(dfodc);
    end