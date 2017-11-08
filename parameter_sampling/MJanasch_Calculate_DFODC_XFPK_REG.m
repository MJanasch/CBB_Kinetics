%% Function to calculate dfodc for the XFPK model
% Markus Janasch, Ph.D. Student, KTH
% Created: 2017-10-26, last modified: 2017-10-26

% This function reassigns the right parameter values to the constants used
% for promiscous enzymes and calculates the numerical values of the partial
% derivatives of the rate equations

function [dfodc,ParameterSet] = MJanasch_Calculate_DFODC_XFPK_REG(N,ParameterSet,V_K_Indeces,Fluxes,SFull)


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
ParameterSet(42) = ParameterSet(34); % KmF6P
ParameterSet(43) = ParameterSet(35); % KmGAP
ParameterSet(44) = ParameterSet(36); % KmE4P
ParameterSet(45) = ParameterSet(37); % KmXuP
ParameterSet(46) = ParameterSet(38); % KmS7P
ParameterSet(47) = ParameterSet(39); % KmR5P

% For FBA (same values as in ALD)
ParameterSet(61) = ParameterSet(27); % KmFBP
ParameterSet(62) = ParameterSet(28); % KmDHAP
ParameterSet(63) = ParameterSet(29); % KmGAP
ParameterSet(64) = ParameterSet(30); % KmSBP
ParameterSet(65) = ParameterSet(31); % KmE4P

% For SBPase (same values as FBPrec)
ParameterSet(68) = ParameterSet(50); % KmFBP
ParameterSet(69) = ParameterSet(51); % KmSBP

% For XFPK2 (same values as XFPK1)
ParameterSet(124) = ParameterSet(114); % KmF6P
ParameterSet(125) = ParameterSet(115); % KmPHI
ParameterSet(126) = ParameterSet(116); % KmE4P
ParameterSet(127) = ParameterSet(117); % KmACETP
ParameterSet(128) = ParameterSet(118); % KmXu5P
ParameterSet(129) = ParameterSet(119); % KmGAP


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
    
    V_maxv9     = ParameterSet(53);
    KmF6Pv9     = ParameterSet(54);
    KmATPv9     = ParameterSet(55);
    KmFBPv9     = ParameterSet(56);
    KmADPv9     = ParameterSet(57);
    K_eqv9      = ParameterSet(58);
    KiATPv9     = ParameterSet(59);
    
    V_maxv10    = ParameterSet(60);
    KmFBPv10    = ParameterSet(61); % like KmFBPv5
    KmDHAPv10   = ParameterSet(62); % like KmDHAPv5
    KmGAPv10    = ParameterSet(63); % like KmGAPv5
    KmSBPv10    = ParameterSet(64); % like KmSBPv5
    KmE4Pv10    = ParameterSet(65); % like KmE4Pv5
    K_eqv10     = ParameterSet(66);
    
    V_maxv11    = ParameterSet(67);
    KmFBPv11    = ParameterSet(68); % like KmFBPv8
    KmSBPv11    = ParameterSet(69); % like KmSBPv8
    K_eqv11     = ParameterSet(70);
    
    V_maxv12    = ParameterSet(71);
    KmR5Pv12    = ParameterSet(72);
    KmRu5Pv12   = ParameterSet(73);
    K_eqv12     = ParameterSet(74);
   
    V_maxv13    = ParameterSet(75);
    KmXu5Pv13   = ParameterSet(76);
    KmRu5Pv13   = ParameterSet(77);
    K_eqv13     = ParameterSet(78);
    
    V_maxv14    = ParameterSet(79);
    KmRu5Pv14   = ParameterSet(80);
    KmATPv14    = ParameterSet(81);
    KmRuBPv14   = ParameterSet(82);
    KmADPv14    = ParameterSet(83);
    K_eqv14     = ParameterSet(84);
    KiPEPv14    = ParameterSet(85);
    KiADPv14    = ParameterSet(86);
    
    V_maxv15    = ParameterSet(87);
    KmP3Gv15    = ParameterSet(88);
    KmP2Gv15    = ParameterSet(89);
    K_eqv15     = ParameterSet(90);
    
    V_maxv16    = ParameterSet(91);
    KmP2Gv16    = ParameterSet(92);
    KmPEPv16    = ParameterSet(93);
    K_eqv16     = ParameterSet(94);
    
    V_maxv17    = ParameterSet(95);
    KmPEPv17    = ParameterSet(96);
    KmPYRv17    = ParameterSet(97);
    KmADPv17    = ParameterSet(98);
    KmATPv17    = ParameterSet(99);
    K_eqv17     = ParameterSet(100);
    KaR5Pv17    = ParameterSet(101);
    KiPHIv17    = ParameterSet(102);
    KiFBPv17    = ParameterSet(103);
    KiATPv17    = ParameterSet(104);
   
    V_maxv18    = ParameterSet(105);
    KmPYRv18    = ParameterSet(106);
    KmNADv18    = ParameterSet(107);
    KmCOAv18    = ParameterSet(108);
    KmACCOAv18  = ParameterSet(109);
    KmNADHv18   = ParameterSet(110);
    KmCO2_cytv18= ParameterSet(111);
    K_eqv18     = ParameterSet(112);
    
    V_maxv19    = ParameterSet(113);
    KmF6Pv19    = ParameterSet(114);
    KmPHIv19    = ParameterSet(115);
    KmE4Pv19    = ParameterSet(116);
    KmACETPv19  = ParameterSet(117);
    KmXu5Pv19   = ParameterSet(118);
    KmGAPv19    = ParameterSet(119);
    K_eqv19     = ParameterSet(120);
    KiATPv19    = ParameterSet(121);
    KaADPv19    = ParameterSet(122);
    
    
    V_maxv20    = ParameterSet(123);
    KmF6Pv20    = ParameterSet(124); % like KmF6Pv19 
    KmPHIv20    = ParameterSet(125); % like KmPHIv19
    KmE4Pv20    = ParameterSet(126); % like KmE4Pv19
    KmACETPv20  = ParameterSet(127); % like KmACETPv19
    KmXu5Pv20   = ParameterSet(128); % like KmXu5Pv19
    KmGAPv20    = ParameterSet(129); % like KmGAPv19
    K_eqv20     = ParameterSet(130);
    KiATPv20    = ParameterSet(131);
    KaADPv20    = ParameterSet(132);
    
    V_maxv21    = ParameterSet(133);
    KmACETPv21  = ParameterSet(134);
    KmCOAv21    = ParameterSet(135);
    KmACCOAv21  = ParameterSet(136);
    KmPHIv21    = ParameterSet(137);
    K_eqv21     = ParameterSet(138);
    
    V_maxv22    = ParameterSet(139);
    KmADPv22    = ParameterSet(140);
    KmPHIv22    = ParameterSet(141);
    KmATPv22    = ParameterSet(142);
    K_eqv22     = ParameterSet(143);
    
    V_maxv23    = ParameterSet(144);
    KmNADPv23   = ParameterSet(145);
    KmNADPHv23  = ParameterSet(146);
    K_eqv23     = ParameterSet(147);
    
    V_maxv24    = ParameterSet(148);
    KmR5Pv24    = ParameterSet(149);
    
    V_maxv25    = ParameterSet(150);
    KmF6Pv25    = ParameterSet(151);
    
    V_maxv26    = ParameterSet(152);
    KmE4Pv26    = ParameterSet(153);
    
    V_maxv27    = ParameterSet(154);
    KmPEPv27    = ParameterSet(155);

    V_maxv28    = ParameterSet(156);
    KmPYRv28    = ParameterSet(157);

    V_maxv29    = ParameterSet(158);
    KmP3Gv29    = ParameterSet(159);
    
    K_PPoolv30  = ParameterSet(160);


%% Rate Equations    
    
    pippo(1,1) = (PHI/(PHI+KaPHIv1))*V_maxv1*(RuBP*CO2_cax/(CO2_cax+KmCO2_caxv1*(1+O2/KmO2v1)))/(RuBP+KmRuBPv1*(1+PHI/KiPHIv1+NADPH/KiNADPHv1));
    pippo(2,1) = V_maxv2*(P3G*ATP/(KmP3Gv2*KmATPv2))*(1-BPG*ADP/(P3G*ATP)/K_eqv2)/((1+P3G/KmP3Gv2+BPG/KmBPGv2)*(1+ATP/KmATPv2+ADP/KmADPv2));
    pippo(3,1) = V_maxv3*(BPG*NADPH/(KmBPGv3*KmNADPHv3))*(1-GAP*NADP*PHI/(BPG*NADPH)/K_eqv3)/((1+NADP/KmNADPv3+NADPH/KmNADPHv3)*(1+GAP/KmGAPv3+BPG/KmBPGv3+PHI/KmPHIv3));
    pippo(4,1) = V_maxv4*(GAP/KmGAPv4)*(1-DHAP/GAP/K_eqv4)/(1+GAP/KmGAPv4+DHAP/KmDHAPv4);
    pippo(5,1) = V_maxv5*((DHAP/KmDHAPv5) * (GAP/KmGAPv5))*(1-FBP/(DHAP*GAP)/K_eqv5)/((1+FBP/KmFBPv5)*(1+SBP/KmSBPv5)+(1+DHAP/KmDHAPv5)*(1+E4P/KmE4Pv5)*(1+GAP/KmGAPv5)-1);
    
    pippo(6,1) = V_maxv6*(F6P*GAP/(KmF6Pv6*KmGAPv6))*(1-E4P*Xu5P/(F6P*GAP)/K_eqv6)/((1+F6P/KmF6Pv6+E4P/KmE4Pv6)*(1+GAP/KmGAPv6+Xu5P/KmXu5Pv6)*(1+S7P/KmS7Pv6+R5P/KmR5Pv6));
    pippo(7,1) = V_maxv7*(S7P*GAP/(KmS7Pv7*KmGAPv7))*(1-R5P*Xu5P/(S7P*GAP)/K_eqv7)/((1+F6P/KmF6Pv7+E4P/KmE4Pv7)*(1+GAP/KmGAPv7+Xu5P/KmXu5Pv7)*(1+S7P/KmS7Pv7+R5P/KmR5Pv7));
    pippo(8,1) = V_maxv8 * FBP/(FBP + KmFBPv8*(1+(SBP/KmSBPv8)));
    pippo(9,1) = (KiATPv9/(ATP+KiATPv9)) * V_maxv9*(F6P*ATP/(KmF6Pv9*KmATPv9))*(1-FBP*ADP/(F6P*ATP)/K_eqv9)/((1+F6P/KmF6Pv9+FBP/KmFBPv9)*(1+ATP/KmATPv9+ADP/KmADPv9));
    pippo(10,1)= V_maxv10*((DHAP/KmDHAPv10) * (E4P/KmE4Pv10))*(1-SBP/(DHAP*E4P)/K_eqv10)/((1+FBP/KmFBPv10)*(1+SBP/KmSBPv10)+(1+DHAP/KmDHAPv10)*(1+E4P/KmE4Pv10)*(1+GAP/KmGAPv10)-1);
    
    pippo(11,1)= V_maxv11 * SBP/(SBP + KmSBPv11*(1+(FBP/KmFBPv11)));
    pippo(12,1)= V_maxv12*(R5P/KmR5Pv12)*(1-Ru5P/R5P/K_eqv12)/(1+Ru5P/KmRu5Pv12+R5P/KmR5Pv12);
    pippo(13,1)= V_maxv13*(Xu5P/KmXu5Pv13)*(1-Ru5P/Xu5P/K_eqv13)/(1+Ru5P/KmRu5Pv13+Xu5P/KmXu5Pv13);
    pippo(14,1)= (KiPEPv14/(PEP+KiPEPv14))*(KiADPv14/(ADP+KiADPv14)) * V_maxv14*(Ru5P*ATP/(KmRu5Pv14*KmATPv14))*(1-RuBP*ADP/(Ru5P*ATP)/K_eqv14)/((1+Ru5P/KmRu5Pv14+RuBP/KmRuBPv14)*(1+ATP/KmATPv14+ADP/KmADPv14));
    pippo(15,1)= V_maxv15*(P3G/KmP3Gv15)*(1-P2G/P3G/K_eqv15)/(1+P3G/KmP3Gv15+P2G/KmP2Gv15);
    
    pippo(16,1)= V_maxv16*(P2G/KmP2Gv16)*(1-PEP/P2G/K_eqv16)/(1+P2G/KmP2Gv16+PEP/KmPEPv16);
    pippo(17,1)= (R5P/(R5P+KaR5Pv17))*(KiPHIv17/(PHI+KiPHIv17))*(KiFBPv17/(FBP+KiFBPv17))*(KiATPv17/(ATP+KiATPv17))*V_maxv17*(PEP*ADP/(KmPEPv17*KmADPv17))*(1-PYR*ATP/(PEP*ADP)/K_eqv17)/((1+PEP/KmPEPv17+PYR/KmPYRv17)*(1+ADP/KmADPv17+ATP/KmATPv17));
    pippo(18,1)= V_maxv18 * (PYR * NAD * COA/(KmPYRv18 * KmNADv18 * KmCOAv18))*(1-(ACCOA * NADH * CO2_cyt)/(PYR*NAD*COA*K_eqv18))/((1+ NAD/KmNADv18 + NADH/KmNADHv18)*(1 + PYR/KmPYRv18 + CO2_cyt/KmCO2_cytv18) * (ACCOA/KmACCOAv18 + COA/KmCOAv18));
    
    pippo(19,1)= (ADP/(ADP+KaADPv19))*(KiATPv19/(ATP+KiATPv19))*V_maxv19*(F6P*PHI/(KmF6Pv19*KmPHIv19))*(1-E4P*ACETP/(F6P*PHI)/K_eqv19)/((1+F6P/KmF6Pv19+E4P/KmE4Pv19)*(1+PHI/KmPHIv19+ACETP/KmACETPv19)*(1+Xu5P/KmXu5Pv19+GAP/KmGAPv19));
    pippo(20,1)= (ADP/(ADP+KaADPv20))*(KiATPv20/(ATP+KiATPv20))*V_maxv20*(Xu5P*PHI/(KmXu5Pv20*KmPHIv20))*(1-GAP*ACETP/(Xu5P*PHI)/K_eqv20)/((1+F6P/KmF6Pv20+E4P/KmE4Pv20)*(1+PHI/KmPHIv20+ACETP/KmACETPv20)*(1+Xu5P/KmXu5Pv20+GAP/KmGAPv20));
    pippo(21,1)= V_maxv21*(ACETP*COA/(KmACETPv21*KmCOAv21))*(1-ACCOA*PHI/(ACETP*COA)/K_eqv21)/((1+ACETP/KmACETPv21+ACCOA/KmACCOAv21)*(1+COA/KmCOAv21+PHI/KmPHIv21));
    
    pippo(22,1)= V_maxv22*(ADP*PHI-ATP/K_eqv22)/(KmADPv22*KmPHIv22*(1+ADP/KmADPv22+PHI/KmPHIv22+ATP/KmATPv22+(ADP*PHI)/(KmADPv22*KmPHIv22)));
    pippo(23,1)= V_maxv23 * (NADP / KmNADPv23) * (1 - NADPH / NADP / K_eqv23) / (1 + NADP / KmNADPv23 + NADPH / KmNADPHv23);
    
    pippo(24,1)= V_maxv24 * R5P/(KmR5Pv24 + R5P);
    pippo(25,1)= V_maxv25*F6P/(KmF6Pv25+F6P);
    pippo(26,1)= V_maxv26*E4P/(KmE4Pv26+E4P);
    pippo(27,1)= V_maxv27*PEP/(KmPEPv27+PEP);
    pippo(28,1)= V_maxv28*PYR/(KmPYRv28+PYR);
    pippo(29,1)= V_maxv29*P3G/(KmP3Gv29+P3G);
    pippo(30,1)= K_PPoolv30*(PPool - PHI);    

    
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
    V_maxv29    = ParameterSet(V_K_Indeces(29));
    K_PPoolv30  = ParameterSet(V_K_Indeces(30));
    
    
%% Derivatives for DFODC matrix (Numerical values of the derivatives)
S = size(SFull,1);
R = size(SFull,2);

% Initialize dfodc matrix
dfodc = zeros(S,R);

% Row 1
dfodc(1,2)  = (ATP*P3G*V_maxv2*((ADP*BPG)/(ATP*K_eqv2*P3G) - 1))/(KmATPv2*KmBPGv2*KmP3Gv2*(BPG/KmBPGv2 + P3G/KmP3Gv2 + 1)^2*(ADP/KmADPv2 + ATP/KmATPv2 + 1)) - (ADP*V_maxv2)/(K_eqv2*KmATPv2*KmP3Gv2*(BPG/KmBPGv2 + P3G/KmP3Gv2 + 1)*(ADP/KmADPv2 + ATP/KmATPv2 + 1));
dfodc(1,3)  = (BPG*NADPH*V_maxv3*((GAP*NADP*PHI)/(BPG*K_eqv3*NADPH) - 1))/(KmBPGv3^2*KmNADPHv3*(NADP/KmNADPv3 + NADPH/KmNADPHv3 + 1)*(BPG/KmBPGv3 + GAP/KmGAPv3 + PHI/KmPHIv3 + 1)^2) - (NADPH*V_maxv3*((GAP*NADP*PHI)/(BPG*K_eqv3*NADPH) - 1))/(KmBPGv3*KmNADPHv3*(NADP/KmNADPv3 + NADPH/KmNADPHv3 + 1)*(BPG/KmBPGv3 + GAP/KmGAPv3 + PHI/KmPHIv3 + 1)) + (GAP*NADP*PHI*V_maxv3)/(BPG*K_eqv3*KmBPGv3*KmNADPHv3*(NADP/KmNADPv3 + NADPH/KmNADPHv3 + 1)*(BPG/KmBPGv3 + GAP/KmGAPv3 + PHI/KmPHIv3 + 1));

% Row 2
dfodc(2,15) = (P3G*V_maxv15*(P2G/(K_eqv15*P3G) - 1))/(KmP2Gv15*KmP3Gv15*(P2G/KmP2Gv15 + P3G/KmP3Gv15 + 1)^2) - V_maxv15/(K_eqv15*KmP3Gv15*(P2G/KmP2Gv15 + P3G/KmP3Gv15 + 1));
dfodc(2,16) = (P2G*V_maxv16*(PEP/(K_eqv16*P2G) - 1))/(KmP2Gv16^2*(P2G/KmP2Gv16 + PEP/KmPEPv16 + 1)^2) - (V_maxv16*(PEP/(K_eqv16*P2G) - 1))/(KmP2Gv16*(P2G/KmP2Gv16 + PEP/KmPEPv16 + 1)) + (PEP*V_maxv16)/(K_eqv16*KmP2Gv16*P2G*(P2G/KmP2Gv16 + PEP/KmPEPv16 + 1));

% Row 3
dfodc(3,1)  = 0;
dfodc(3,2)  = (ATP*P3G*V_maxv2*((ADP*BPG)/(ATP*K_eqv2*P3G) - 1))/(KmATPv2*KmP3Gv2^2*(BPG/KmBPGv2 + P3G/KmP3Gv2 + 1)^2*(ADP/KmADPv2 + ATP/KmATPv2 + 1)) - (ATP*V_maxv2*((ADP*BPG)/(ATP*K_eqv2*P3G) - 1))/(KmATPv2*KmP3Gv2*(BPG/KmBPGv2 + P3G/KmP3Gv2 + 1)*(ADP/KmADPv2 + ATP/KmATPv2 + 1)) + (ADP*BPG*V_maxv2)/(K_eqv2*KmATPv2*KmP3Gv2*P3G*(BPG/KmBPGv2 + P3G/KmP3Gv2 + 1)*(ADP/KmADPv2 + ATP/KmATPv2 + 1));
dfodc(3,15) = (P3G*V_maxv15*(P2G/(K_eqv15*P3G) - 1))/(KmP3Gv15^2*(P2G/KmP2Gv15 + P3G/KmP3Gv15 + 1)^2) - (V_maxv15*(P2G/(K_eqv15*P3G) - 1))/(KmP3Gv15*(P2G/KmP2Gv15 + P3G/KmP3Gv15 + 1)) + (P2G*V_maxv15)/(K_eqv15*KmP3Gv15*P3G*(P2G/KmP2Gv15 + P3G/KmP3Gv15 + 1));
dfodc(3,29) = V_maxv29/(KmP3Gv29 + P3G) - (P3G*V_maxv29)/(KmP3Gv29 + P3G)^2;

% Row 4
dfodc(4,2)  = (ATP*P3G*V_maxv2*((ADP*BPG)/(ATP*K_eqv2*P3G) - 1))/(KmADPv2*KmATPv2*KmP3Gv2*(BPG/KmBPGv2 + P3G/KmP3Gv2 + 1)*(ADP/KmADPv2 + ATP/KmATPv2 + 1)^2) - (BPG*V_maxv2)/(K_eqv2*KmATPv2*KmP3Gv2*(BPG/KmBPGv2 + P3G/KmP3Gv2 + 1)*(ADP/KmADPv2 + ATP/KmATPv2 + 1));
dfodc(4,9)  = (ATP*F6P*KiATPv9*V_maxv9*((ADP*FBP)/(ATP*F6P*K_eqv9) - 1))/(KmADPv9*KmATPv9*KmF6Pv9*(ATP + KiATPv9)*(F6P/KmF6Pv9 + FBP/KmFBPv9 + 1)*(ADP/KmADPv9 + ATP/KmATPv9 + 1)^2) - (FBP*KiATPv9*V_maxv9)/(K_eqv9*KmATPv9*KmF6Pv9*(ATP + KiATPv9)*(F6P/KmF6Pv9 + FBP/KmFBPv9 + 1)*(ADP/KmADPv9 + ATP/KmATPv9 + 1));
dfodc(4,14) = (ATP*KiADPv14*KiPEPv14*Ru5P*V_maxv14*((ADP*RuBP)/(ATP*K_eqv14*Ru5P) - 1))/(KmATPv14*KmRu5Pv14*(ADP + KiADPv14)^2*(KiPEPv14 + PEP)*(Ru5P/KmRu5Pv14 + RuBP/KmRuBPv14 + 1)*(ADP/KmADPv14 + ATP/KmATPv14 + 1)) - (KiADPv14*KiPEPv14*RuBP*V_maxv14)/(K_eqv14*KmATPv14*KmRu5Pv14*(ADP + KiADPv14)*(KiPEPv14 + PEP)*(Ru5P/KmRu5Pv14 + RuBP/KmRuBPv14 + 1)*(ADP/KmADPv14 + ATP/KmATPv14 + 1)) + (ATP*KiADPv14*KiPEPv14*Ru5P*V_maxv14*((ADP*RuBP)/(ATP*K_eqv14*Ru5P) - 1))/(KmADPv14*KmATPv14*KmRu5Pv14*(ADP + KiADPv14)*(KiPEPv14 + PEP)*(Ru5P/KmRu5Pv14 + RuBP/KmRuBPv14 + 1)*(ADP/KmADPv14 + ATP/KmATPv14 + 1)^2);
dfodc(4,17) = (ADP*KiATPv17*KiFBPv17*KiPHIv17*PEP*R5P*V_maxv17*((ATP*PYR)/(ADP*K_eqv17*PEP) - 1))/(KmADPv17^2*KmPEPv17*(ATP + KiATPv17)*(FBP + KiFBPv17)*(KiPHIv17 + PHI)*(KaR5Pv17 + R5P)*(PEP/KmPEPv17 + PYR/KmPYRv17 + 1)*(ADP/KmADPv17 + ATP/KmATPv17 + 1)^2) - (KiATPv17*KiFBPv17*KiPHIv17*PEP*R5P*V_maxv17*((ATP*PYR)/(ADP*K_eqv17*PEP) - 1))/(KmADPv17*KmPEPv17*(ATP + KiATPv17)*(FBP + KiFBPv17)*(KiPHIv17 + PHI)*(KaR5Pv17 + R5P)*(PEP/KmPEPv17 + PYR/KmPYRv17 + 1)*(ADP/KmADPv17 + ATP/KmATPv17 + 1)) + (ATP*KiATPv17*KiFBPv17*KiPHIv17*PYR*R5P*V_maxv17)/(ADP*K_eqv17*KmADPv17*KmPEPv17*(ATP + KiATPv17)*(FBP + KiFBPv17)*(KiPHIv17 + PHI)*(KaR5Pv17 + R5P)*(PEP/KmPEPv17 + PYR/KmPYRv17 + 1)*(ADP/KmADPv17 + ATP/KmATPv17 + 1));
dfodc(4,19) = (ADP*F6P*KiATPv19*PHI*V_maxv19*((ACETP*E4P)/(F6P*K_eqv19*PHI) - 1))/(KmF6Pv19*KmPHIv19*(ADP + KaADPv19)^2*(ATP + KiATPv19)*(ACETP/KmACETPv19 + PHI/KmPHIv19 + 1)*(GAP/KmGAPv19 + Xu5P/KmXu5Pv19 + 1)*(E4P/KmE4Pv19 + F6P/KmF6Pv19 + 1)) - (F6P*KiATPv19*PHI*V_maxv19*((ACETP*E4P)/(F6P*K_eqv19*PHI) - 1))/(KmF6Pv19*KmPHIv19*(ADP + KaADPv19)*(ATP + KiATPv19)*(ACETP/KmACETPv19 + PHI/KmPHIv19 + 1)*(GAP/KmGAPv19 + Xu5P/KmXu5Pv19 + 1)*(E4P/KmE4Pv19 + F6P/KmF6Pv19 + 1));
dfodc(4,20) = (ADP*KiATPv20*PHI*V_maxv20*Xu5P*((ACETP*GAP)/(K_eqv20*PHI*Xu5P) - 1))/(KmPHIv20*KmXu5Pv20*(ADP + KaADPv20)^2*(ATP + KiATPv20)*(ACETP/KmACETPv20 + PHI/KmPHIv20 + 1)*(GAP/KmGAPv20 + Xu5P/KmXu5Pv20 + 1)*(E4P/KmE4Pv20 + F6P/KmF6Pv20 + 1)) - (KiATPv20*PHI*V_maxv20*Xu5P*((ACETP*GAP)/(K_eqv20*PHI*Xu5P) - 1))/(KmPHIv20*KmXu5Pv20*(ADP + KaADPv20)*(ATP + KiATPv20)*(ACETP/KmACETPv20 + PHI/KmPHIv20 + 1)*(GAP/KmGAPv20 + Xu5P/KmXu5Pv20 + 1)*(E4P/KmE4Pv20 + F6P/KmF6Pv20 + 1));
dfodc(4,22) = (PHI*V_maxv22)/(KmADPv22*KmPHIv22*(ADP/KmADPv22 + ATP/KmATPv22 + PHI/KmPHIv22 + (ADP*PHI)/(KmADPv22*KmPHIv22) + 1)) + (V_maxv22*(1/KmADPv22 + PHI/(KmADPv22*KmPHIv22))*(ATP/K_eqv22 - ADP*PHI))/(KmADPv22*KmPHIv22*(ADP/KmADPv22 + ATP/KmATPv22 + PHI/KmPHIv22 + (ADP*PHI)/(KmADPv22*KmPHIv22) + 1)^2);


% Row 5
dfodc(5,4)  = (GAP*V_maxv4*(DHAP/(GAP*K_eqv4) - 1))/(KmDHAPv4*KmGAPv4*(DHAP/KmDHAPv4 + GAP/KmGAPv4 + 1)^2) - V_maxv4/(K_eqv4*KmGAPv4*(DHAP/KmDHAPv4 + GAP/KmGAPv4 + 1));
dfodc(5,5)  = (FBP*V_maxv5)/(DHAP*K_eqv5*KmDHAPv5*KmGAPv5*((FBP/KmFBPv5 + 1)*(SBP/KmSBPv5 + 1) + (DHAP/KmDHAPv5 + 1)*(E4P/KmE4Pv5 + 1)*(GAP/KmGAPv5 + 1) - 1)) - (GAP*V_maxv5*(FBP/(DHAP*GAP*K_eqv5) - 1))/(KmDHAPv5*KmGAPv5*((FBP/KmFBPv5 + 1)*(SBP/KmSBPv5 + 1) + (DHAP/KmDHAPv5 + 1)*(E4P/KmE4Pv5 + 1)*(GAP/KmGAPv5 + 1) - 1)) + (DHAP*GAP*V_maxv5*(E4P/KmE4Pv5 + 1)*(GAP/KmGAPv5 + 1)*(FBP/(DHAP*GAP*K_eqv5) - 1))/(KmDHAPv5^2*KmGAPv5*((FBP/KmFBPv5 + 1)*(SBP/KmSBPv5 + 1) + (DHAP/KmDHAPv5 + 1)*(E4P/KmE4Pv5 + 1)*(GAP/KmGAPv5 + 1) - 1)^2);
dfodc(5,10) = (SBP*V_maxv10)/(DHAP*K_eqv10*KmDHAPv10*KmE4Pv10*((FBP/KmFBPv10 + 1)*(SBP/KmSBPv10 + 1) + (DHAP/KmDHAPv10 + 1)*(E4P/KmE4Pv10 + 1)*(GAP/KmGAPv10 + 1) - 1)) - (E4P*V_maxv10*(SBP/(DHAP*E4P*K_eqv10) - 1))/(KmDHAPv10*KmE4Pv10*((FBP/KmFBPv10 + 1)*(SBP/KmSBPv10 + 1) + (DHAP/KmDHAPv10 + 1)*(E4P/KmE4Pv10 + 1)*(GAP/KmGAPv10 + 1) - 1)) + (DHAP*E4P*V_maxv10*(SBP/(DHAP*E4P*K_eqv10) - 1)*(E4P/KmE4Pv10 + 1)*(GAP/KmGAPv10 + 1))/(KmDHAPv10^2*KmE4Pv10*((FBP/KmFBPv10 + 1)*(SBP/KmSBPv10 + 1) + (DHAP/KmDHAPv10 + 1)*(E4P/KmE4Pv10 + 1)*(GAP/KmGAPv10 + 1) - 1)^2);

% Row 6
dfodc(6,5)  = (DHAP*GAP*V_maxv5*(DHAP/KmDHAPv5 + 1)*(GAP/KmGAPv5 + 1)*(FBP/(DHAP*GAP*K_eqv5) - 1))/(KmDHAPv5*KmE4Pv5*KmGAPv5*((FBP/KmFBPv5 + 1)*(SBP/KmSBPv5 + 1) + (DHAP/KmDHAPv5 + 1)*(E4P/KmE4Pv5 + 1)*(GAP/KmGAPv5 + 1) - 1)^2);
dfodc(6,6)  = (F6P*GAP*V_maxv6*((E4P*Xu5P)/(F6P*GAP*K_eqv6) - 1))/(KmE4Pv6*KmF6Pv6*KmGAPv6*(GAP/KmGAPv6 + Xu5P/KmXu5Pv6 + 1)*(R5P/KmR5Pv6 + S7P/KmS7Pv6 + 1)*(E4P/KmE4Pv6 + F6P/KmF6Pv6 + 1)^2) - (V_maxv6*Xu5P)/(K_eqv6*KmF6Pv6*KmGAPv6*(GAP/KmGAPv6 + Xu5P/KmXu5Pv6 + 1)*(R5P/KmR5Pv6 + S7P/KmS7Pv6 + 1)*(E4P/KmE4Pv6 + F6P/KmF6Pv6 + 1));
dfodc(6,7)  = (GAP*S7P*V_maxv7*((R5P*Xu5P)/(GAP*K_eqv7*S7P) - 1))/(KmE4Pv7*KmGAPv7*KmS7Pv7*(GAP/KmGAPv7 + Xu5P/KmXu5Pv7 + 1)*(R5P/KmR5Pv7 + S7P/KmS7Pv7 + 1)*(E4P/KmE4Pv7 + F6P/KmF6Pv7 + 1)^2);
dfodc(6,10) = (SBP*V_maxv10)/(E4P*K_eqv10*KmDHAPv10*KmE4Pv10*((FBP/KmFBPv10 + 1)*(SBP/KmSBPv10 + 1) + (DHAP/KmDHAPv10 + 1)*(E4P/KmE4Pv10 + 1)*(GAP/KmGAPv10 + 1) - 1)) - (DHAP*V_maxv10*(SBP/(DHAP*E4P*K_eqv10) - 1))/(KmDHAPv10*KmE4Pv10*((FBP/KmFBPv10 + 1)*(SBP/KmSBPv10 + 1) + (DHAP/KmDHAPv10 + 1)*(E4P/KmE4Pv10 + 1)*(GAP/KmGAPv10 + 1) - 1)) + (DHAP*E4P*V_maxv10*(SBP/(DHAP*E4P*K_eqv10) - 1)*(DHAP/KmDHAPv10 + 1)*(GAP/KmGAPv10 + 1))/(KmDHAPv10*KmE4Pv10^2*((FBP/KmFBPv10 + 1)*(SBP/KmSBPv10 + 1) + (DHAP/KmDHAPv10 + 1)*(E4P/KmE4Pv10 + 1)*(GAP/KmGAPv10 + 1) - 1)^2);
dfodc(6,19) = (ADP*F6P*KiATPv19*PHI*V_maxv19*((ACETP*E4P)/(F6P*K_eqv19*PHI) - 1))/(KmE4Pv19*KmF6Pv19*KmPHIv19*(ADP + KaADPv19)*(ATP + KiATPv19)*(ACETP/KmACETPv19 + PHI/KmPHIv19 + 1)*(GAP/KmGAPv19 + Xu5P/KmXu5Pv19 + 1)*(E4P/KmE4Pv19 + F6P/KmF6Pv19 + 1)^2) - (ACETP*ADP*KiATPv19*V_maxv19)/(K_eqv19*KmF6Pv19*KmPHIv19*(ADP + KaADPv19)*(ATP + KiATPv19)*(ACETP/KmACETPv19 + PHI/KmPHIv19 + 1)*(GAP/KmGAPv19 + Xu5P/KmXu5Pv19 + 1)*(E4P/KmE4Pv19 + F6P/KmF6Pv19 + 1));
dfodc(6,20) = (ADP*KiATPv20*PHI*V_maxv20*Xu5P*((ACETP*GAP)/(K_eqv20*PHI*Xu5P) - 1))/(KmE4Pv20*KmPHIv20*KmXu5Pv20*(ADP + KaADPv20)*(ATP + KiATPv20)*(ACETP/KmACETPv20 + PHI/KmPHIv20 + 1)*(GAP/KmGAPv20 + Xu5P/KmXu5Pv20 + 1)*(E4P/KmE4Pv20 + F6P/KmF6Pv20 + 1)^2);
dfodc(6,26) = V_maxv26/(E4P + KmE4Pv26) - (E4P*V_maxv26)/(E4P + KmE4Pv26)^2;

% Row 7 
dfodc(7,6)  = (F6P*GAP*V_maxv6*((E4P*Xu5P)/(F6P*GAP*K_eqv6) - 1))/(KmF6Pv6^2*KmGAPv6*(GAP/KmGAPv6 + Xu5P/KmXu5Pv6 + 1)*(R5P/KmR5Pv6 + S7P/KmS7Pv6 + 1)*(E4P/KmE4Pv6 + F6P/KmF6Pv6 + 1)^2) - (GAP*V_maxv6*((E4P*Xu5P)/(F6P*GAP*K_eqv6) - 1))/(KmF6Pv6*KmGAPv6*(GAP/KmGAPv6 + Xu5P/KmXu5Pv6 + 1)*(R5P/KmR5Pv6 + S7P/KmS7Pv6 + 1)*(E4P/KmE4Pv6 + F6P/KmF6Pv6 + 1)) + (E4P*V_maxv6*Xu5P)/(F6P*K_eqv6*KmF6Pv6*KmGAPv6*(GAP/KmGAPv6 + Xu5P/KmXu5Pv6 + 1)*(R5P/KmR5Pv6 + S7P/KmS7Pv6 + 1)*(E4P/KmE4Pv6 + F6P/KmF6Pv6 + 1));
dfodc(7,7)  = (GAP*S7P*V_maxv7*((R5P*Xu5P)/(GAP*K_eqv7*S7P) - 1))/(KmF6Pv7*KmGAPv7*KmS7Pv7*(GAP/KmGAPv7 + Xu5P/KmXu5Pv7 + 1)*(R5P/KmR5Pv7 + S7P/KmS7Pv7 + 1)*(E4P/KmE4Pv7 + F6P/KmF6Pv7 + 1)^2);
dfodc(7,8)  = 0;
dfodc(7,9)  = (ATP*F6P*KiATPv9*V_maxv9*((ADP*FBP)/(ATP*F6P*K_eqv9) - 1))/(KmATPv9*KmF6Pv9^2*(ATP + KiATPv9)*(F6P/KmF6Pv9 + FBP/KmFBPv9 + 1)^2*(ADP/KmADPv9 + ATP/KmATPv9 + 1)) - (ATP*KiATPv9*V_maxv9*((ADP*FBP)/(ATP*F6P*K_eqv9) - 1))/(KmATPv9*KmF6Pv9*(ATP + KiATPv9)*(F6P/KmF6Pv9 + FBP/KmFBPv9 + 1)*(ADP/KmADPv9 + ATP/KmATPv9 + 1)) + (ADP*FBP*KiATPv9*V_maxv9)/(F6P*K_eqv9*KmATPv9*KmF6Pv9*(ATP + KiATPv9)*(F6P/KmF6Pv9 + FBP/KmFBPv9 + 1)*(ADP/KmADPv9 + ATP/KmATPv9 + 1));
dfodc(7,11) = 0;
dfodc(7,19) = (ADP*F6P*KiATPv19*PHI*V_maxv19*((ACETP*E4P)/(F6P*K_eqv19*PHI) - 1))/(KmF6Pv19^2*KmPHIv19*(ADP + KaADPv19)*(ATP + KiATPv19)*(ACETP/KmACETPv19 + PHI/KmPHIv19 + 1)*(GAP/KmGAPv19 + Xu5P/KmXu5Pv19 + 1)*(E4P/KmE4Pv19 + F6P/KmF6Pv19 + 1)^2) - (ADP*KiATPv19*PHI*V_maxv19*((ACETP*E4P)/(F6P*K_eqv19*PHI) - 1))/(KmF6Pv19*KmPHIv19*(ADP + KaADPv19)*(ATP + KiATPv19)*(ACETP/KmACETPv19 + PHI/KmPHIv19 + 1)*(GAP/KmGAPv19 + Xu5P/KmXu5Pv19 + 1)*(E4P/KmE4Pv19 + F6P/KmF6Pv19 + 1)) + (ACETP*ADP*E4P*KiATPv19*V_maxv19)/(F6P*K_eqv19*KmF6Pv19*KmPHIv19*(ADP + KaADPv19)*(ATP + KiATPv19)*(ACETP/KmACETPv19 + PHI/KmPHIv19 + 1)*(GAP/KmGAPv19 + Xu5P/KmXu5Pv19 + 1)*(E4P/KmE4Pv19 + F6P/KmF6Pv19 + 1));
dfodc(7,20) = (ADP*KiATPv20*PHI*V_maxv20*Xu5P*((ACETP*GAP)/(K_eqv20*PHI*Xu5P) - 1))/(KmF6Pv20*KmPHIv20*KmXu5Pv20*(ADP + KaADPv20)*(ATP + KiATPv20)*(ACETP/KmACETPv20 + PHI/KmPHIv20 + 1)*(GAP/KmGAPv20 + Xu5P/KmXu5Pv20 + 1)*(E4P/KmE4Pv20 + F6P/KmF6Pv20 + 1)^2);
dfodc(7,25) = V_maxv25/(F6P + KmF6Pv25) - (F6P*V_maxv25)/(F6P + KmF6Pv25)^2;

% Row 8
dfodc(8,5)  = (DHAP*GAP*V_maxv5*(SBP/KmSBPv5 + 1)*(FBP/(DHAP*GAP*K_eqv5) - 1))/(KmDHAPv5*KmFBPv5*KmGAPv5*((FBP/KmFBPv5 + 1)*(SBP/KmSBPv5 + 1) + (DHAP/KmDHAPv5 + 1)*(E4P/KmE4Pv5 + 1)*(GAP/KmGAPv5 + 1) - 1)^2) - V_maxv5/(K_eqv5*KmDHAPv5*KmGAPv5*((FBP/KmFBPv5 + 1)*(SBP/KmSBPv5 + 1) + (DHAP/KmDHAPv5 + 1)*(E4P/KmE4Pv5 + 1)*(GAP/KmGAPv5 + 1) - 1));
dfodc(8,8)  = V_maxv8/(FBP + KmFBPv8*(SBP/KmSBPv8 + 1)) - (FBP*V_maxv8)/(FBP + KmFBPv8*(SBP/KmSBPv8 + 1))^2;
dfodc(8,9)  = (ATP*F6P*KiATPv9*V_maxv9*((ADP*FBP)/(ATP*F6P*K_eqv9) - 1))/(KmATPv9*KmF6Pv9*KmFBPv9*(ATP + KiATPv9)*(F6P/KmF6Pv9 + FBP/KmFBPv9 + 1)^2*(ADP/KmADPv9 + ATP/KmATPv9 + 1)) - (ADP*KiATPv9*V_maxv9)/(K_eqv9*KmATPv9*KmF6Pv9*(ATP + KiATPv9)*(F6P/KmF6Pv9 + FBP/KmFBPv9 + 1)*(ADP/KmADPv9 + ATP/KmATPv9 + 1));
dfodc(8,10) = (DHAP*E4P*V_maxv10*(SBP/(DHAP*E4P*K_eqv10) - 1)*(SBP/KmSBPv10 + 1))/(KmDHAPv10*KmFBPv10*KmE4Pv10*((FBP/KmFBPv10 + 1)*(SBP/KmSBPv10 + 1) + (DHAP/KmDHAPv10 + 1)*(E4P/KmE4Pv10 + 1)*(GAP/KmGAPv10 + 1) - 1)^2);
dfodc(8,11) = -(KmSBPv11*SBP*V_maxv11)/(KmFBPv11*(SBP + KmSBPv11*(FBP/KmFBPv11 + 1))^2);
dfodc(8,17) = (ADP*KiATPv17*KiFBPv17*KiPHIv17*PEP*R5P*V_maxv17*((ATP*PYR)/(ADP*K_eqv17*PEP) - 1))/(KmADPv17*KmPEPv17*(ATP + KiATPv17)*(FBP + KiFBPv17)^2*(KiPHIv17 + PHI)*(KaR5Pv17 + R5P)*(PEP/KmPEPv17 + PYR/KmPYRv17 + 1)*(ADP/KmADPv17 + ATP/KmATPv17 + 1));

% Row 9
dfodc(9,3)	= (BPG*NADPH*V_maxv3*((GAP*NADP*PHI)/(BPG*K_eqv3*NADPH) - 1))/(KmBPGv3*KmGAPv3*KmNADPHv3*(NADP/KmNADPv3 + NADPH/KmNADPHv3 + 1)*(BPG/KmBPGv3 + GAP/KmGAPv3 + PHI/KmPHIv3 + 1)^2) - (NADP*PHI*V_maxv3)/(K_eqv3*KmBPGv3*KmNADPHv3*(NADP/KmNADPv3 + NADPH/KmNADPHv3 + 1)*(BPG/KmBPGv3 + GAP/KmGAPv3 + PHI/KmPHIv3 + 1));
dfodc(9,4)  = (GAP*V_maxv4*(DHAP/(GAP*K_eqv4) - 1))/(KmGAPv4^2*(DHAP/KmDHAPv4 + GAP/KmGAPv4 + 1)^2) - (V_maxv4*(DHAP/(GAP*K_eqv4) - 1))/(KmGAPv4*(DHAP/KmDHAPv4 + GAP/KmGAPv4 + 1)) + (DHAP*V_maxv4)/(GAP*K_eqv4*KmGAPv4*(DHAP/KmDHAPv4 + GAP/KmGAPv4 + 1));
dfodc(9,5)  = (FBP*V_maxv5)/(GAP*K_eqv5*KmDHAPv5*KmGAPv5*((FBP/KmFBPv5 + 1)*(SBP/KmSBPv5 + 1) + (DHAP/KmDHAPv5 + 1)*(E4P/KmE4Pv5 + 1)*(GAP/KmGAPv5 + 1) - 1)) - (DHAP*V_maxv5*(FBP/(DHAP*GAP*K_eqv5) - 1))/(KmDHAPv5*KmGAPv5*((FBP/KmFBPv5 + 1)*(SBP/KmSBPv5 + 1) + (DHAP/KmDHAPv5 + 1)*(E4P/KmE4Pv5 + 1)*(GAP/KmGAPv5 + 1) - 1)) + (DHAP*GAP*V_maxv5*(DHAP/KmDHAPv5 + 1)*(E4P/KmE4Pv5 + 1)*(FBP/(DHAP*GAP*K_eqv5) - 1))/(KmDHAPv5*KmGAPv5^2*((FBP/KmFBPv5 + 1)*(SBP/KmSBPv5 + 1) + (DHAP/KmDHAPv5 + 1)*(E4P/KmE4Pv5 + 1)*(GAP/KmGAPv5 + 1) - 1)^2);
dfodc(9,6)  = (F6P*GAP*V_maxv6*((E4P*Xu5P)/(F6P*GAP*K_eqv6) - 1))/(KmF6Pv6*KmGAPv6^2*(GAP/KmGAPv6 + Xu5P/KmXu5Pv6 + 1)^2*(R5P/KmR5Pv6 + S7P/KmS7Pv6 + 1)*(E4P/KmE4Pv6 + F6P/KmF6Pv6 + 1)) - (F6P*V_maxv6*((E4P*Xu5P)/(F6P*GAP*K_eqv6) - 1))/(KmF6Pv6*KmGAPv6*(GAP/KmGAPv6 + Xu5P/KmXu5Pv6 + 1)*(R5P/KmR5Pv6 + S7P/KmS7Pv6 + 1)*(E4P/KmE4Pv6 + F6P/KmF6Pv6 + 1)) + (E4P*V_maxv6*Xu5P)/(GAP*K_eqv6*KmF6Pv6*KmGAPv6*(GAP/KmGAPv6 + Xu5P/KmXu5Pv6 + 1)*(R5P/KmR5Pv6 + S7P/KmS7Pv6 + 1)*(E4P/KmE4Pv6 + F6P/KmF6Pv6 + 1));
dfodc(9,7)  = (GAP*S7P*V_maxv7*((R5P*Xu5P)/(GAP*K_eqv7*S7P) - 1))/(KmGAPv7^2*KmS7Pv7*(GAP/KmGAPv7 + Xu5P/KmXu5Pv7 + 1)^2*(R5P/KmR5Pv7 + S7P/KmS7Pv7 + 1)*(E4P/KmE4Pv7 + F6P/KmF6Pv7 + 1)) - (S7P*V_maxv7*((R5P*Xu5P)/(GAP*K_eqv7*S7P) - 1))/(KmGAPv7*KmS7Pv7*(GAP/KmGAPv7 + Xu5P/KmXu5Pv7 + 1)*(R5P/KmR5Pv7 + S7P/KmS7Pv7 + 1)*(E4P/KmE4Pv7 + F6P/KmF6Pv7 + 1)) + (R5P*V_maxv7*Xu5P)/(GAP*K_eqv7*KmGAPv7*KmS7Pv7*(GAP/KmGAPv7 + Xu5P/KmXu5Pv7 + 1)*(R5P/KmR5Pv7 + S7P/KmS7Pv7 + 1)*(E4P/KmE4Pv7 + F6P/KmF6Pv7 + 1));
dfodc(9,10) = (DHAP*E4P*V_maxv10*(SBP/(DHAP*E4P*K_eqv10) - 1)*(DHAP/KmDHAPv10 + 1)*(E4P/KmE4Pv10 + 1))/(KmDHAPv10*KmE4Pv10*KmGAPv10*((FBP/KmFBPv10 + 1)*(SBP/KmSBPv10 + 1) + (DHAP/KmDHAPv10 + 1)*(E4P/KmE4Pv10 + 1)*(GAP/KmGAPv10 + 1) - 1)^2);
dfodc(9,19) = (ADP*F6P*KiATPv19*PHI*V_maxv19*((ACETP*E4P)/(F6P*K_eqv19*PHI) - 1))/(KmGAPv19*KmF6Pv19*KmPHIv19*(ADP + KaADPv19)*(ATP + KiATPv19)*(ACETP/KmACETPv19 + PHI/KmPHIv19 + 1)*(GAP/KmGAPv19 + Xu5P/KmXu5Pv19 + 1)^2*(E4P/KmE4Pv19 + F6P/KmF6Pv19 + 1));
dfodc(9,20) = (ADP*KiATPv20*PHI*V_maxv20*Xu5P*((ACETP*GAP)/(K_eqv20*PHI*Xu5P) - 1))/(KmGAPv20*KmPHIv20*KmXu5Pv20*(ADP + KaADPv20)*(ATP + KiATPv20)*(ACETP/KmACETPv20 + PHI/KmPHIv20 + 1)*(GAP/KmGAPv20 + Xu5P/KmXu5Pv20 + 1)^2*(E4P/KmE4Pv20 + F6P/KmF6Pv20 + 1)) - (ACETP*ADP*KiATPv20*V_maxv20)/(K_eqv20*KmPHIv20*KmXu5Pv20*(ADP + KaADPv20)*(ATP + KiATPv20)*(ACETP/KmACETPv20 + PHI/KmPHIv20 + 1)*(GAP/KmGAPv20 + Xu5P/KmXu5Pv20 + 1)*(E4P/KmE4Pv20 + F6P/KmF6Pv20 + 1));

% Row 10
dfodc(10,1) = -(CO2_cax*KmRuBPv1*PHI*RuBP*V_maxv1)/(KiNADPHv1*(CO2_cax + KmCO2_caxv1*(O2/KmO2v1 + 1))*(RuBP + KmRuBPv1*(NADPH/KiNADPHv1 + PHI/KiPHIv1 + 1))^2*(KaPHIv1 + PHI));
dfodc(10,3) = (BPG*NADPH*V_maxv3*((GAP*NADP*PHI)/(BPG*K_eqv3*NADPH) - 1))/(KmBPGv3*KmNADPHv3^2*(NADP/KmNADPv3 + NADPH/KmNADPHv3 + 1)^2*(BPG/KmBPGv3 + GAP/KmGAPv3 + PHI/KmPHIv3 + 1)) - (BPG*V_maxv3*((GAP*NADP*PHI)/(BPG*K_eqv3*NADPH) - 1))/(KmBPGv3*KmNADPHv3*(NADP/KmNADPv3 + NADPH/KmNADPHv3 + 1)*(BPG/KmBPGv3 + GAP/KmGAPv3 + PHI/KmPHIv3 + 1)) + (GAP*NADP*PHI*V_maxv3)/(K_eqv3*KmBPGv3*KmNADPHv3*NADPH*(NADP/KmNADPv3 + NADPH/KmNADPHv3 + 1)*(BPG/KmBPGv3 + GAP/KmGAPv3 + PHI/KmPHIv3 + 1));
dfodc(10,23)= (NADP*V_maxv23*(NADPH/(K_eqv23*NADP) - 1))/(KmNADPv23*KmNADPHv23*(NADP/KmNADPv23 + NADPH/KmNADPHv23 + 1)^2) - V_maxv23/(K_eqv23*KmNADPv23*(NADP/KmNADPv23 + NADPH/KmNADPHv23 + 1));

% Row 11
dfodc(11,14)= (ATP*KiADPv14*KiPEPv14*Ru5P*V_maxv14*((ADP*RuBP)/(ATP*K_eqv14*Ru5P) - 1))/(KmATPv14*KmRu5Pv14*(ADP + KiADPv14)*(KiPEPv14 + PEP)^2*(Ru5P/KmRu5Pv14 + RuBP/KmRuBPv14 + 1)*(ADP/KmADPv14 + ATP/KmATPv14 + 1));
dfodc(11,16)= (P2G*V_maxv16*(PEP/(K_eqv16*P2G) - 1))/(KmPEPv16*KmP2Gv16*(P2G/KmP2Gv16 + PEP/KmPEPv16 + 1)^2) - V_maxv16/(K_eqv16*KmP2Gv16*(P2G/KmP2Gv16 + PEP/KmPEPv16 + 1));
dfodc(11,17)= (ADP*KiATPv17*KiFBPv17*KiPHIv17*PEP*R5P*V_maxv17*((ATP*PYR)/(ADP*K_eqv17*PEP) - 1))/(KmADPv17*KmPEPv17^2*(ATP + KiATPv17)*(FBP + KiFBPv17)*(KiPHIv17 + PHI)*(KaR5Pv17 + R5P)*(PEP/KmPEPv17 + PYR/KmPYRv17 + 1)^2*(ADP/KmADPv17 + ATP/KmATPv17 + 1)) - (ADP*KiATPv17*KiFBPv17*KiPHIv17*R5P*V_maxv17*((ATP*PYR)/(ADP*K_eqv17*PEP) - 1))/(KmADPv17*KmPEPv17*(ATP + KiATPv17)*(FBP + KiFBPv17)*(KiPHIv17 + PHI)*(KaR5Pv17 + R5P)*(PEP/KmPEPv17 + PYR/KmPYRv17 + 1)*(ADP/KmADPv17 + ATP/KmATPv17 + 1)) + (ATP*KiATPv17*KiFBPv17*KiPHIv17*PYR*R5P*V_maxv17)/(K_eqv17*KmADPv17*KmPEPv17*PEP*(ATP + KiATPv17)*(FBP + KiFBPv17)*(KiPHIv17 + PHI)*(KaR5Pv17 + R5P)*(PEP/KmPEPv17 + PYR/KmPYRv17 + 1)*(ADP/KmADPv17 + ATP/KmATPv17 + 1));
dfodc(11,27)= V_maxv27/(KmPEPv27 + PEP) - (PEP*V_maxv27)/(KmPEPv27 + PEP)^2;

% Row 12
dfodc(12,1) = (CO2_cax*RuBP*V_maxv1)/((CO2_cax + KmCO2_caxv1*(O2/KmO2v1 + 1))*(RuBP + KmRuBPv1*(NADPH/KiNADPHv1 + PHI/KiPHIv1 + 1))*(KaPHIv1 + PHI)) - (CO2_cax*PHI*RuBP*V_maxv1)/((CO2_cax + KmCO2_caxv1*(O2/KmO2v1 + 1))*(RuBP + KmRuBPv1*(NADPH/KiNADPHv1 + PHI/KiPHIv1 + 1))*(KaPHIv1 + PHI)^2) - (CO2_cax*KmRuBPv1*PHI*RuBP*V_maxv1)/(KiPHIv1*(CO2_cax + KmCO2_caxv1*(O2/KmO2v1 + 1))*(RuBP + KmRuBPv1*(NADPH/KiNADPHv1 + PHI/KiPHIv1 + 1))^2*(KaPHIv1 + PHI));
dfodc(12,3) = (BPG*NADPH*V_maxv3*((GAP*NADP*PHI)/(BPG*K_eqv3*NADPH) - 1))/(KmBPGv3*KmNADPHv3*KmPHIv3*(NADP/KmNADPv3 + NADPH/KmNADPHv3 + 1)*(BPG/KmBPGv3 + GAP/KmGAPv3 + PHI/KmPHIv3 + 1)^2) - (GAP*NADP*V_maxv3)/(K_eqv3*KmBPGv3*KmNADPHv3*(NADP/KmNADPv3 + NADPH/KmNADPHv3 + 1)*(BPG/KmBPGv3 + GAP/KmGAPv3 + PHI/KmPHIv3 + 1));
dfodc(12,8) = 0;
dfodc(12,11)= 0;
dfodc(12,17)= (ADP*KiATPv17*KiFBPv17*KiPHIv17*PEP*R5P*V_maxv17*((ATP*PYR)/(ADP*K_eqv17*PEP) - 1))/(KmADPv17*KmPEPv17*(ATP + KiATPv17)*(FBP + KiFBPv17)*(KiPHIv17 + PHI)^2*(KaR5Pv17 + R5P)*(PEP/KmPEPv17 + PYR/KmPYRv17 + 1)*(ADP/KmADPv17 + ATP/KmATPv17 + 1));
dfodc(12,19)= (ADP*F6P*KiATPv19*PHI*V_maxv19*((ACETP*E4P)/(F6P*K_eqv19*PHI) - 1))/(KmF6Pv19*KmPHIv19^2*(ADP + KaADPv19)*(ATP + KiATPv19)*(ACETP/KmACETPv19 + PHI/KmPHIv19 + 1)^2*(GAP/KmGAPv19 + Xu5P/KmXu5Pv19 + 1)*(E4P/KmE4Pv19 + F6P/KmF6Pv19 + 1)) - (ADP*F6P*KiATPv19*V_maxv19*((ACETP*E4P)/(F6P*K_eqv19*PHI) - 1))/(KmF6Pv19*KmPHIv19*(ADP + KaADPv19)*(ATP + KiATPv19)*(ACETP/KmACETPv19 + PHI/KmPHIv19 + 1)*(GAP/KmGAPv19 + Xu5P/KmXu5Pv19 + 1)*(E4P/KmE4Pv19 + F6P/KmF6Pv19 + 1)) + (ACETP*ADP*E4P*KiATPv19*V_maxv19)/(K_eqv19*KmF6Pv19*KmPHIv19*PHI*(ADP + KaADPv19)*(ATP + KiATPv19)*(ACETP/KmACETPv19 + PHI/KmPHIv19 + 1)*(GAP/KmGAPv19 + Xu5P/KmXu5Pv19 + 1)*(E4P/KmE4Pv19 + F6P/KmF6Pv19 + 1));
dfodc(12,20)= (ADP*KiATPv20*PHI*V_maxv20*Xu5P*((ACETP*GAP)/(K_eqv20*PHI*Xu5P) - 1))/(KmPHIv20^2*KmXu5Pv20*(ADP + KaADPv20)*(ATP + KiATPv20)*(ACETP/KmACETPv20 + PHI/KmPHIv20 + 1)^2*(GAP/KmGAPv20 + Xu5P/KmXu5Pv20 + 1)*(E4P/KmE4Pv20 + F6P/KmF6Pv20 + 1)) - (ADP*KiATPv20*V_maxv20*Xu5P*((ACETP*GAP)/(K_eqv20*PHI*Xu5P) - 1))/(KmPHIv20*KmXu5Pv20*(ADP + KaADPv20)*(ATP + KiATPv20)*(ACETP/KmACETPv20 + PHI/KmPHIv20 + 1)*(GAP/KmGAPv20 + Xu5P/KmXu5Pv20 + 1)*(E4P/KmE4Pv20 + F6P/KmF6Pv20 + 1)) + (ACETP*ADP*GAP*KiATPv20*V_maxv20)/(K_eqv20*KmPHIv20*KmXu5Pv20*PHI*(ADP + KaADPv20)*(ATP + KiATPv20)*(ACETP/KmACETPv20 + PHI/KmPHIv20 + 1)*(GAP/KmGAPv20 + Xu5P/KmXu5Pv20 + 1)*(E4P/KmE4Pv20 + F6P/KmF6Pv20 + 1));
dfodc(12,21)= (ACETP*COA*V_maxv21*((ACCOA*PHI)/(ACETP*COA*K_eqv21) - 1))/(KmACETPv21*KmCOAv21*KmPHIv21*(COA/KmCOAv21 + PHI/KmPHIv21 + 1)^2*(ACCOA/KmACCOAv21 + ACETP/KmACETPv21 + 1)) - (ACCOA*V_maxv21)/(K_eqv21*KmACETPv21*KmCOAv21*(COA/KmCOAv21 + PHI/KmPHIv21 + 1)*(ACCOA/KmACCOAv21 + ACETP/KmACETPv21 + 1));
dfodc(12,22)= (ADP*V_maxv22)/(KmADPv22*KmPHIv22*(ADP/KmADPv22 + ATP/KmATPv22 + PHI/KmPHIv22 + (ADP*PHI)/(KmADPv22*KmPHIv22) + 1)) + (V_maxv22*(1/KmPHIv22 + ADP/(KmADPv22*KmPHIv22))*(ATP/K_eqv22 - ADP*PHI))/(KmADPv22*KmPHIv22*(ADP/KmADPv22 + ATP/KmATPv22 + PHI/KmPHIv22 + (ADP*PHI)/(KmADPv22*KmPHIv22) + 1)^2);
dfodc(12,30)= -K_PPoolv30;

% Row 13
dfodc(13,6) = (F6P*GAP*V_maxv6*((E4P*Xu5P)/(F6P*GAP*K_eqv6) - 1))/(KmF6Pv6*KmGAPv6*KmR5Pv6*(GAP/KmGAPv6 + Xu5P/KmXu5Pv6 + 1)*(R5P/KmR5Pv6 + S7P/KmS7Pv6 + 1)^2*(E4P/KmE4Pv6 + F6P/KmF6Pv6 + 1));
dfodc(13,7) = (GAP*S7P*V_maxv7*((R5P*Xu5P)/(GAP*K_eqv7*S7P) - 1))/(KmGAPv7*KmR5Pv7*KmS7Pv7*(GAP/KmGAPv7 + Xu5P/KmXu5Pv7 + 1)*(R5P/KmR5Pv7 + S7P/KmS7Pv7 + 1)^2*(E4P/KmE4Pv7 + F6P/KmF6Pv7 + 1)) - (V_maxv7*Xu5P)/(K_eqv7*KmGAPv7*KmS7Pv7*(GAP/KmGAPv7 + Xu5P/KmXu5Pv7 + 1)*(R5P/KmR5Pv7 + S7P/KmS7Pv7 + 1)*(E4P/KmE4Pv7 + F6P/KmF6Pv7 + 1));
dfodc(13,12)= (R5P*V_maxv12*(Ru5P/(K_eqv12*R5P) - 1))/(KmR5Pv12^2*(R5P/KmR5Pv12 + Ru5P/KmRu5Pv12 + 1)^2) - (V_maxv12*(Ru5P/(K_eqv12*R5P) - 1))/(KmR5Pv12*(R5P/KmR5Pv12 + Ru5P/KmRu5Pv12 + 1)) + (Ru5P*V_maxv12)/(K_eqv12*KmR5Pv12*R5P*(R5P/KmR5Pv12 + Ru5P/KmRu5Pv12 + 1));
dfodc(13,17)= (ADP*KiATPv17*KiFBPv17*KiPHIv17*PEP*R5P*V_maxv17*((ATP*PYR)/(ADP*K_eqv17*PEP) - 1))/(KmADPv17*KmPEPv17*(ATP + KiATPv17)*(FBP + KiFBPv17)*(KiPHIv17 + PHI)*(KaR5Pv17 + R5P)^2*(PEP/KmPEPv17 + PYR/KmPYRv17 + 1)*(ADP/KmADPv17 + ATP/KmATPv17 + 1)) - (ADP*KiATPv17*KiFBPv17*KiPHIv17*PEP*V_maxv17*((ATP*PYR)/(ADP*K_eqv17*PEP) - 1))/(KmADPv17*KmPEPv17*(ATP + KiATPv17)*(FBP + KiFBPv17)*(KiPHIv17 + PHI)*(KaR5Pv17 + R5P)*(PEP/KmPEPv17 + PYR/KmPYRv17 + 1)*(ADP/KmADPv17 + ATP/KmATPv17 + 1));
dfodc(13,24)= V_maxv24/(KmR5Pv24 + R5P) - (R5P*V_maxv24)/(KmR5Pv24 + R5P)^2;

% Row 14
dfodc(14,12)= (R5P*V_maxv12*(Ru5P/(K_eqv12*R5P) - 1))/(KmR5Pv12*KmRu5Pv12*(R5P/KmR5Pv12 + Ru5P/KmRu5Pv12 + 1)^2) - V_maxv12/(K_eqv12*KmR5Pv12*(R5P/KmR5Pv12 + Ru5P/KmRu5Pv12 + 1));
dfodc(14,13)= (V_maxv13*Xu5P*(Ru5P/(K_eqv13*Xu5P) - 1))/(KmRu5Pv13*KmXu5Pv13*(Ru5P/KmRu5Pv13 + Xu5P/KmXu5Pv13 + 1)^2) - V_maxv13/(K_eqv13*KmXu5Pv13*(Ru5P/KmRu5Pv13 + Xu5P/KmXu5Pv13 + 1));
dfodc(14,14)= (ATP*KiADPv14*KiPEPv14*Ru5P*V_maxv14*((ADP*RuBP)/(ATP*K_eqv14*Ru5P) - 1))/(KmATPv14*KmRu5Pv14^2*(ADP + KiADPv14)*(KiPEPv14 + PEP)*(Ru5P/KmRu5Pv14 + RuBP/KmRuBPv14 + 1)^2*(ADP/KmADPv14 + ATP/KmATPv14 + 1)) - (ATP*KiADPv14*KiPEPv14*V_maxv14*((ADP*RuBP)/(ATP*K_eqv14*Ru5P) - 1))/(KmATPv14*KmRu5Pv14*(ADP + KiADPv14)*(KiPEPv14 + PEP)*(Ru5P/KmRu5Pv14 + RuBP/KmRuBPv14 + 1)*(ADP/KmADPv14 + ATP/KmATPv14 + 1)) + (ADP*KiADPv14*KiPEPv14*RuBP*V_maxv14)/(K_eqv14*KmATPv14*KmRu5Pv14*Ru5P*(ADP + KiADPv14)*(KiPEPv14 + PEP)*(Ru5P/KmRu5Pv14 + RuBP/KmRuBPv14 + 1)*(ADP/KmADPv14 + ATP/KmATPv14 + 1));

% Row 15
dfodc(15,1) = (CO2_cax*PHI*V_maxv1)/((CO2_cax + KmCO2_caxv1*(O2/KmO2v1 + 1))*(RuBP + KmRuBPv1*(NADPH/KiNADPHv1 + PHI/KiPHIv1 + 1))*(KaPHIv1 + PHI)) - (CO2_cax*PHI*RuBP*V_maxv1)/((CO2_cax + KmCO2_caxv1*(O2/KmO2v1 + 1))*(RuBP + KmRuBPv1*(NADPH/KiNADPHv1 + PHI/KiPHIv1 + 1))^2*(KaPHIv1 + PHI));
dfodc(15,14)= (ATP*KiADPv14*KiPEPv14*Ru5P*V_maxv14*((ADP*RuBP)/(ATP*K_eqv14*Ru5P) - 1))/(KmATPv14*KmRuBPv14*KmRu5Pv14*(ADP + KiADPv14)*(KiPEPv14 + PEP)*(Ru5P/KmRu5Pv14 + RuBP/KmRuBPv14 + 1)^2*(ADP/KmADPv14 + ATP/KmATPv14 + 1)) - (ADP*KiADPv14*KiPEPv14*V_maxv14)/(K_eqv14*KmATPv14*KmRu5Pv14*(ADP + KiADPv14)*(KiPEPv14 + PEP)*(Ru5P/KmRu5Pv14 + RuBP/KmRuBPv14 + 1)*(ADP/KmADPv14 + ATP/KmATPv14 + 1));

% Row 16
dfodc(16,6) = (F6P*GAP*V_maxv6*((E4P*Xu5P)/(F6P*GAP*K_eqv6) - 1))/(KmF6Pv6*KmGAPv6*KmS7Pv6*(GAP/KmGAPv6 + Xu5P/KmXu5Pv6 + 1)*(R5P/KmR5Pv6 + S7P/KmS7Pv6 + 1)^2*(E4P/KmE4Pv6 + F6P/KmF6Pv6 + 1));
dfodc(16,7) = (GAP*S7P*V_maxv7*((R5P*Xu5P)/(GAP*K_eqv7*S7P) - 1))/(KmGAPv7*KmS7Pv7^2*(GAP/KmGAPv7 + Xu5P/KmXu5Pv7 + 1)*(R5P/KmR5Pv7 + S7P/KmS7Pv7 + 1)^2*(E4P/KmE4Pv7 + F6P/KmF6Pv7 + 1)) - (GAP*V_maxv7*((R5P*Xu5P)/(GAP*K_eqv7*S7P) - 1))/(KmGAPv7*KmS7Pv7*(GAP/KmGAPv7 + Xu5P/KmXu5Pv7 + 1)*(R5P/KmR5Pv7 + S7P/KmS7Pv7 + 1)*(E4P/KmE4Pv7 + F6P/KmF6Pv7 + 1)) + (R5P*V_maxv7*Xu5P)/(K_eqv7*KmGAPv7*KmS7Pv7*S7P*(GAP/KmGAPv7 + Xu5P/KmXu5Pv7 + 1)*(R5P/KmR5Pv7 + S7P/KmS7Pv7 + 1)*(E4P/KmE4Pv7 + F6P/KmF6Pv7 + 1));
dfodc(16,8) = 0;
dfodc(16,11)= 0;

% Row 17
dfodc(17,5) = (DHAP*GAP*V_maxv5*(FBP/KmFBPv5 + 1)*(FBP/(DHAP*GAP*K_eqv5) - 1))/(KmDHAPv5*KmGAPv5*KmSBPv5*((FBP/KmFBPv5 + 1)*(SBP/KmSBPv5 + 1) + (DHAP/KmDHAPv5 + 1)*(E4P/KmE4Pv5 + 1)*(GAP/KmGAPv5 + 1) - 1)^2);
dfodc(17,8) = -(FBP*KmFBPv8*V_maxv8)/(KmSBPv8*(FBP + KmFBPv8*(SBP/KmSBPv8 + 1))^2);
dfodc(17,10)= (DHAP*E4P*V_maxv10*(SBP/(DHAP*E4P*K_eqv10) - 1)*(FBP/KmFBPv10 + 1))/(KmDHAPv10*KmE4Pv10*KmSBPv10*((FBP/KmFBPv10 + 1)*(SBP/KmSBPv10 + 1) + (DHAP/KmDHAPv10 + 1)*(E4P/KmE4Pv10 + 1)*(GAP/KmGAPv10 + 1) - 1)^2) - V_maxv10/(K_eqv10*KmDHAPv10*KmE4Pv10*((FBP/KmFBPv10 + 1)*(SBP/KmSBPv10 + 1) + (DHAP/KmDHAPv10 + 1)*(E4P/KmE4Pv10 + 1)*(GAP/KmGAPv10 + 1) - 1));
dfodc(17,11)= V_maxv11/(SBP + KmSBPv11*(FBP/KmFBPv11 + 1)) - (SBP*V_maxv11)/(SBP + KmSBPv11*(FBP/KmFBPv11 + 1))^2;

% Row 18
dfodc(18,6) = (F6P*GAP*V_maxv6*((E4P*Xu5P)/(F6P*GAP*K_eqv6) - 1))/(KmF6Pv6*KmGAPv6*KmXu5Pv6*(GAP/KmGAPv6 + Xu5P/KmXu5Pv6 + 1)^2*(R5P/KmR5Pv6 + S7P/KmS7Pv6 + 1)*(E4P/KmE4Pv6 + F6P/KmF6Pv6 + 1)) - (E4P*V_maxv6)/(K_eqv6*KmF6Pv6*KmGAPv6*(GAP/KmGAPv6 + Xu5P/KmXu5Pv6 + 1)*(R5P/KmR5Pv6 + S7P/KmS7Pv6 + 1)*(E4P/KmE4Pv6 + F6P/KmF6Pv6 + 1));
dfodc(18,7) = (GAP*S7P*V_maxv7*((R5P*Xu5P)/(GAP*K_eqv7*S7P) - 1))/(KmGAPv7*KmS7Pv7*KmXu5Pv7*(GAP/KmGAPv7 + Xu5P/KmXu5Pv7 + 1)^2*(R5P/KmR5Pv7 + S7P/KmS7Pv7 + 1)*(E4P/KmE4Pv7 + F6P/KmF6Pv7 + 1)) - (R5P*V_maxv7)/(K_eqv7*KmGAPv7*KmS7Pv7*(GAP/KmGAPv7 + Xu5P/KmXu5Pv7 + 1)*(R5P/KmR5Pv7 + S7P/KmS7Pv7 + 1)*(E4P/KmE4Pv7 + F6P/KmF6Pv7 + 1));
dfodc(18,13)= (V_maxv13*Xu5P*(Ru5P/(K_eqv13*Xu5P) - 1))/(KmXu5Pv13^2*(Ru5P/KmRu5Pv13 + Xu5P/KmXu5Pv13 + 1)^2) - (V_maxv13*(Ru5P/(K_eqv13*Xu5P) - 1))/(KmXu5Pv13*(Ru5P/KmRu5Pv13 + Xu5P/KmXu5Pv13 + 1)) + (Ru5P*V_maxv13)/(K_eqv13*KmXu5Pv13*Xu5P*(Ru5P/KmRu5Pv13 + Xu5P/KmXu5Pv13 + 1));
dfodc(18,19)= (ADP*F6P*KiATPv19*PHI*V_maxv19*((ACETP*E4P)/(F6P*K_eqv19*PHI) - 1))/(KmF6Pv19*KmPHIv19*KmXu5Pv19*(ADP + KaADPv19)*(ATP + KiATPv19)*(ACETP/KmACETPv19 + PHI/KmPHIv19 + 1)*(GAP/KmGAPv19 + Xu5P/KmXu5Pv19 + 1)^2*(E4P/KmE4Pv19 + F6P/KmF6Pv19 + 1));
dfodc(18,20)= (ADP*KiATPv20*PHI*V_maxv20*Xu5P*((ACETP*GAP)/(K_eqv20*PHI*Xu5P) - 1))/(KmPHIv20*KmXu5Pv20^2*(ADP + KaADPv20)*(ATP + KiATPv20)*(ACETP/KmACETPv20 + PHI/KmPHIv20 + 1)*(GAP/KmGAPv20 + Xu5P/KmXu5Pv20 + 1)^2*(E4P/KmE4Pv20 + F6P/KmF6Pv20 + 1)) - (ADP*KiATPv20*PHI*V_maxv20*((ACETP*GAP)/(K_eqv20*PHI*Xu5P) - 1))/(KmPHIv20*KmXu5Pv20*(ADP + KaADPv20)*(ATP + KiATPv20)*(ACETP/KmACETPv20 + PHI/KmPHIv20 + 1)*(GAP/KmGAPv20 + Xu5P/KmXu5Pv20 + 1)*(E4P/KmE4Pv20 + F6P/KmF6Pv20 + 1)) + (ACETP*ADP*GAP*KiATPv20*V_maxv20)/(K_eqv20*KmPHIv20*KmXu5Pv20*Xu5P*(ADP + KaADPv20)*(ATP + KiATPv20)*(ACETP/KmACETPv20 + PHI/KmPHIv20 + 1)*(GAP/KmGAPv20 + Xu5P/KmXu5Pv20 + 1)*(E4P/KmE4Pv20 + F6P/KmF6Pv20 + 1));

% Row 19
dfodc(19,17)= (ADP*KiATPv17*KiFBPv17*KiPHIv17*PEP*R5P*V_maxv17*((ATP*PYR)/(ADP*K_eqv17*PEP) - 1))/(KmADPv17*KmPEPv17*KmPYRv17*(ATP + KiATPv17)*(FBP + KiFBPv17)*(KiPHIv17 + PHI)*(KaR5Pv17 + R5P)*(PEP/KmPEPv17 + PYR/KmPYRv17 + 1)^2*(ADP/KmADPv17 + ATP/KmATPv17 + 1)) - (ATP*KiATPv17*KiFBPv17*KiPHIv17*R5P*V_maxv17)/(K_eqv17*KmADPv17*KmPEPv17*(ATP + KiATPv17)*(FBP + KiFBPv17)*(KiPHIv17 + PHI)*(KaR5Pv17 + R5P)*(PEP/KmPEPv17 + PYR/KmPYRv17 + 1)*(ADP/KmADPv17 + ATP/KmATPv17 + 1));
dfodc(19,18)= (COA*NAD*PYR*V_maxv18*((ACCOA*CO2_cyt*NADH)/(COA*K_eqv18*NAD*PYR) - 1))/(KmCOAv18*KmNADv18*KmPYRv18^2*(ACCOA/KmACCOAv18 + COA/KmCOAv18)*(CO2_cyt/KmCO2_cytv18 + PYR/KmPYRv18 + 1)^2*(NAD/KmNADv18 + NADH/KmNADHv18 + 1)) - (COA*NAD*V_maxv18*((ACCOA*CO2_cyt*NADH)/(COA*K_eqv18*NAD*PYR) - 1))/(KmCOAv18*KmNADv18*KmPYRv18*(ACCOA/KmACCOAv18 + COA/KmCOAv18)*(CO2_cyt/KmCO2_cytv18 + PYR/KmPYRv18 + 1)*(NAD/KmNADv18 + NADH/KmNADHv18 + 1)) + (ACCOA*CO2_cyt*NADH*V_maxv18)/(K_eqv18*KmCOAv18*KmNADv18*KmPYRv18*PYR*(ACCOA/KmACCOAv18 + COA/KmCOAv18)*(CO2_cyt/KmCO2_cytv18 + PYR/KmPYRv18 + 1)*(NAD/KmNADv18 + NADH/KmNADHv18 + 1));
dfodc(19,28)= V_maxv28/(KmPYRv28 + PYR) - (PYR*V_maxv28)/(KmPYRv28 + PYR)^2;

% Row 20
dfodc(20,19)= (ADP*F6P*KiATPv19*PHI*V_maxv19*((ACETP*E4P)/(F6P*K_eqv19*PHI) - 1))/(KmACETPv19*KmF6Pv19*KmPHIv19*(ADP + KaADPv19)*(ATP + KiATPv19)*(ACETP/KmACETPv19 + PHI/KmPHIv19 + 1)^2*(GAP/KmGAPv19 + Xu5P/KmXu5Pv19 + 1)*(E4P/KmE4Pv19 + F6P/KmF6Pv19 + 1)) - (ADP*E4P*KiATPv19*V_maxv19)/(K_eqv19*KmF6Pv19*KmPHIv19*(ADP + KaADPv19)*(ATP + KiATPv19)*(ACETP/KmACETPv19 + PHI/KmPHIv19 + 1)*(GAP/KmGAPv19 + Xu5P/KmXu5Pv19 + 1)*(E4P/KmE4Pv19 + F6P/KmF6Pv19 + 1));
dfodc(20,20)= (ADP*KiATPv20*PHI*V_maxv20*Xu5P*((ACETP*GAP)/(K_eqv20*PHI*Xu5P) - 1))/(KmACETPv20*KmPHIv20*KmXu5Pv20*(ADP + KaADPv20)*(ATP + KiATPv20)*(ACETP/KmACETPv20 + PHI/KmPHIv20 + 1)^2*(GAP/KmGAPv20 + Xu5P/KmXu5Pv20 + 1)*(E4P/KmE4Pv20 + F6P/KmF6Pv20 + 1)) - (ADP*GAP*KiATPv20*V_maxv20)/(K_eqv20*KmPHIv20*KmXu5Pv20*(ADP + KaADPv20)*(ATP + KiATPv20)*(ACETP/KmACETPv20 + PHI/KmPHIv20 + 1)*(GAP/KmGAPv20 + Xu5P/KmXu5Pv20 + 1)*(E4P/KmE4Pv20 + F6P/KmF6Pv20 + 1));
dfodc(20,21)= (ACETP*COA*V_maxv21*((ACCOA*PHI)/(ACETP*COA*K_eqv21) - 1))/(KmACETPv21^2*KmCOAv21*(COA/KmCOAv21 + PHI/KmPHIv21 + 1)*(ACCOA/KmACCOAv21 + ACETP/KmACETPv21 + 1)^2) - (COA*V_maxv21*((ACCOA*PHI)/(ACETP*COA*K_eqv21) - 1))/(KmACETPv21*KmCOAv21*(COA/KmCOAv21 + PHI/KmPHIv21 + 1)*(ACCOA/KmACCOAv21 + ACETP/KmACETPv21 + 1)) + (ACCOA*PHI*V_maxv21)/(ACETP*K_eqv21*KmACETPv21*KmCOAv21*(COA/KmCOAv21 + PHI/KmPHIv21 + 1)*(ACCOA/KmACCOAv21 + ACETP/KmACETPv21 + 1));

% Row 21
dfodc(21,2) = (ATP*P3G*V_maxv2*((ADP*BPG)/(ATP*K_eqv2*P3G) - 1))/(KmATPv2^2*KmP3Gv2*(BPG/KmBPGv2 + P3G/KmP3Gv2 + 1)*(ADP/KmADPv2 + ATP/KmATPv2 + 1)^2) - (P3G*V_maxv2*((ADP*BPG)/(ATP*K_eqv2*P3G) - 1))/(KmATPv2*KmP3Gv2*(BPG/KmBPGv2 + P3G/KmP3Gv2 + 1)*(ADP/KmADPv2 + ATP/KmATPv2 + 1)) + (ADP*BPG*V_maxv2)/(ATP*K_eqv2*KmATPv2*KmP3Gv2*(BPG/KmBPGv2 + P3G/KmP3Gv2 + 1)*(ADP/KmADPv2 + ATP/KmATPv2 + 1));
dfodc(21,9) = (ATP*F6P*KiATPv9*V_maxv9*((ADP*FBP)/(ATP*F6P*K_eqv9) - 1))/(KmATPv9*KmF6Pv9*(ATP + KiATPv9)^2*(F6P/KmF6Pv9 + FBP/KmFBPv9 + 1)*(ADP/KmADPv9 + ATP/KmATPv9 + 1)) - (F6P*KiATPv9*V_maxv9*((ADP*FBP)/(ATP*F6P*K_eqv9) - 1))/(KmATPv9*KmF6Pv9*(ATP + KiATPv9)*(F6P/KmF6Pv9 + FBP/KmFBPv9 + 1)*(ADP/KmADPv9 + ATP/KmATPv9 + 1)) + (ATP*F6P*KiATPv9*V_maxv9*((ADP*FBP)/(ATP*F6P*K_eqv9) - 1))/(KmATPv9^2*KmF6Pv9*(ATP + KiATPv9)*(F6P/KmF6Pv9 + FBP/KmFBPv9 + 1)*(ADP/KmADPv9 + ATP/KmATPv9 + 1)^2) + (ADP*FBP*KiATPv9*V_maxv9)/(ATP*K_eqv9*KmATPv9*KmF6Pv9*(ATP + KiATPv9)*(F6P/KmF6Pv9 + FBP/KmFBPv9 + 1)*(ADP/KmADPv9 + ATP/KmATPv9 + 1));
dfodc(21,14)= (ATP*KiADPv14*KiPEPv14*Ru5P*V_maxv14*((ADP*RuBP)/(ATP*K_eqv14*Ru5P) - 1))/(KmATPv14^2*KmRu5Pv14*(ADP + KiADPv14)*(KiPEPv14 + PEP)*(Ru5P/KmRu5Pv14 + RuBP/KmRuBPv14 + 1)*(ADP/KmADPv14 + ATP/KmATPv14 + 1)^2) - (KiADPv14*KiPEPv14*Ru5P*V_maxv14*((ADP*RuBP)/(ATP*K_eqv14*Ru5P) - 1))/(KmATPv14*KmRu5Pv14*(ADP + KiADPv14)*(KiPEPv14 + PEP)*(Ru5P/KmRu5Pv14 + RuBP/KmRuBPv14 + 1)*(ADP/KmADPv14 + ATP/KmATPv14 + 1)) + (ADP*KiADPv14*KiPEPv14*RuBP*V_maxv14)/(ATP*K_eqv14*KmATPv14*KmRu5Pv14*(ADP + KiADPv14)*(KiPEPv14 + PEP)*(Ru5P/KmRu5Pv14 + RuBP/KmRuBPv14 + 1)*(ADP/KmADPv14 + ATP/KmATPv14 + 1));
dfodc(21,17)= (ADP*KiATPv17*KiFBPv17*KiPHIv17*PEP*R5P*V_maxv17*((ATP*PYR)/(ADP*K_eqv17*PEP) - 1))/(KmADPv17*KmPEPv17*(ATP + KiATPv17)^2*(FBP + KiFBPv17)*(KiPHIv17 + PHI)*(KaR5Pv17 + R5P)*(PEP/KmPEPv17 + PYR/KmPYRv17 + 1)*(ADP/KmADPv17 + ATP/KmATPv17 + 1)) - (KiATPv17*KiFBPv17*KiPHIv17*PYR*R5P*V_maxv17)/(K_eqv17*KmADPv17*KmPEPv17*(ATP + KiATPv17)*(FBP + KiFBPv17)*(KiPHIv17 + PHI)*(KaR5Pv17 + R5P)*(PEP/KmPEPv17 + PYR/KmPYRv17 + 1)*(ADP/KmADPv17 + ATP/KmATPv17 + 1)) + (ADP*KiATPv17*KiFBPv17*KiPHIv17*PEP*R5P*V_maxv17*((ATP*PYR)/(ADP*K_eqv17*PEP) - 1))/(KmADPv17*KmATPv17*KmPEPv17*(ATP + KiATPv17)*(FBP + KiFBPv17)*(KiPHIv17 + PHI)*(KaR5Pv17 + R5P)*(PEP/KmPEPv17 + PYR/KmPYRv17 + 1)*(ADP/KmADPv17 + ATP/KmATPv17 + 1)^2);
dfodc(21,19)= (ADP*F6P*KiATPv19*PHI*V_maxv19*((ACETP*E4P)/(F6P*K_eqv19*PHI) - 1))/(KmF6Pv19*KmPHIv19*(ADP + KaADPv19)*(ATP + KiATPv19)^2*(ACETP/KmACETPv19 + PHI/KmPHIv19 + 1)*(GAP/KmGAPv19 + Xu5P/KmXu5Pv19 + 1)*(E4P/KmE4Pv19 + F6P/KmF6Pv19 + 1));
dfodc(21,20)= (ADP*KiATPv20*PHI*V_maxv20*Xu5P*((ACETP*GAP)/(K_eqv20*PHI*Xu5P) - 1))/(KmPHIv20*KmXu5Pv20*(ADP + KaADPv20)*(ATP + KiATPv20)^2*(ACETP/KmACETPv20 + PHI/KmPHIv20 + 1)*(GAP/KmGAPv20 + Xu5P/KmXu5Pv20 + 1)*(E4P/KmE4Pv20 + F6P/KmF6Pv20 + 1));
dfodc(21,22)= (V_maxv22*(ATP/K_eqv22 - ADP*PHI))/(KmADPv22*KmATPv22*KmPHIv22*(ADP/KmADPv22 + ATP/KmATPv22 + PHI/KmPHIv22 + (ADP*PHI)/(KmADPv22*KmPHIv22) + 1)^2) - V_maxv22/(K_eqv22*KmADPv22*KmPHIv22*(ADP/KmADPv22 + ATP/KmATPv22 + PHI/KmPHIv22 + (ADP*PHI)/(KmADPv22*KmPHIv22) + 1));

% Row 22
dfodc(21,3) = (BPG*NADPH*V_maxv3*((GAP*NADP*PHI)/(BPG*K_eqv3*NADPH) - 1))/(KmBPGv3*KmNADPv3*KmNADPHv3*(NADP/KmNADPv3 + NADPH/KmNADPHv3 + 1)^2*(BPG/KmBPGv3 + GAP/KmGAPv3 + PHI/KmPHIv3 + 1)) - (GAP*PHI*V_maxv3)/(K_eqv3*KmBPGv3*KmNADPHv3*(NADP/KmNADPv3 + NADPH/KmNADPHv3 + 1)*(BPG/KmBPGv3 + GAP/KmGAPv3 + PHI/KmPHIv3 + 1));
dfodc(21,23)= (NADP*V_maxv23*(NADPH/(K_eqv23*NADP) - 1))/(KmNADPv23^2*(NADP/KmNADPv23 + NADPH/KmNADPHv23 + 1)^2) - (V_maxv23*(NADPH/(K_eqv23*NADP) - 1))/(KmNADPv23*(NADP/KmNADPv23 + NADPH/KmNADPHv23 + 1)) + (NADPH*V_maxv23)/(K_eqv23*KmNADPv23*NADP*(NADP/KmNADPv23 + NADPH/KmNADPHv23 + 1));

dfodc = transpose(dfodc);
    end