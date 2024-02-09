function [freq,rho0,lambda0,mu0,rho1,lambda1,mu1] = Choix_Materiaux(Mat_n, Inc_n,freq, varargin)
% Data concernant divers matériaux: solides, liquides, gas
% Extraits des tables de Onda Corporation ou autres
%
% choix des matériaux :
%
% 'PDMS'
% 'silicone'
% 'epoxy'
% 'AQUALENE'
% 'zz_pure'
% 'eau_mer'
% 'gel_agar'
% '47V20'
% 'gel'
% 'gel_test'
% 'PMMA'
% 'plomb'
% 'verre'
% 'FC40'
% 'air_sec'
% 'SF6'
% 'FC40+SF6'
% 'Sonovue'
% 'mat_fortedensité'
% 'aerogel'

%% -------------- vevTeur frequence ---------------------------------------
% nFREQ = 1000; %nFREQ = 1000;
% fMIN = fMAX/nFREQ;
% freq = linspace(fMIN,fMAX,nFREQ).';


%% ------------ propriétés des  materiaux ---------------------------------
[rho0, vL0, vT0, alphaL0, alphaT0] = material_No(Mat_n, freq);
[rho1, vL1, vT1, alphaL1, alphaT1] = material_No(Inc_n, freq);



%% ---- Passage en vL, mu et lambda complexes -----------------------------
% !!! attention mu n'est pas complexe car alphaT n'est pas pris en compte

kL0     = freq*2*pi./vL0 + 1i*alphaL0;
kT0     = freq*2*pi./vT0 + 1i*alphaT0;
mu0     = (rho0.*(freq*2*pi).^2) ./ kT0.^2 ;
lambda0 = (rho0.*(freq*2*pi).^2) ./ kL0.^2 - 2.*mu0 ;

kL1     = freq*2*pi./vL1 + 1i*alphaL1;
kT1     = freq*2*pi./vT1 + 1i*alphaT1;
mu1     = (rho1.*(freq*2*pi).^2) ./ kT1.^2 ;
lambda1 = (rho1.*(freq*2*pi).^2) ./ kL1.^2 - 2.*mu1 ;

%---------------------------------------
%       DATA materiaux
%---------------------------------------

    function [rho, vL, vT, alphaL, alphaT ] = material_No(Mat_n,freq)
        %  Data concernant divers matériaux: solides, liquides, gas
        %  Extraits des tables de Onda Corporation ou autres
        %  Les vitesses vL et vT sont données en mm.µs^(-1)
        %  Les modules (lambda et mu) sont donnés en MPa
        %  Les atténuation (alpha) sont données en Np.MHz^(-n).mm^(-1)
        
        tempArray = ones(1,length(freq));
        
        switch Mat_n
            
            case 'PDMS'
                rho = 1.00 * tempArray;
                vL = 1.000 * tempArray;
                alphaL = 0.1 * freq.^(1.5);
                vT = 0.030 * tempArray;
                alphaT = 100 * freq.^(1.5);
                
            case 'Silicone'
                %%% Silicone Bleu --- Mesures du 7 mars 2011 syr 5MHz_0_empty
                rho = 1.11 * tempArray;                   %%% Mesuré au labo
                vL = 1.00 * tempArray;                    %%% Mesures du 11/04/2011
                alphaL = 0.0158 * freq.^(1.7);
                mu = 1;                                   %%% (en MPa)
                vT = sqrt(mu *1000)/1000 * tempArray;     %%% Valeur hypothetique
                alphaT = 1000 * alphaL;                   %%% Valeur hypothetique
                
            case 'PVC'
                %%% PVC de l'envLume couleur ivoire pour manip G' et G''
                rho = 1.4 * tempArray;                    %%% Mesuré au labo
                vL = 2.34 * tempArray;                    %%% Mesures du 11/04/2011
                alphaL = 0.0158 * freq.^(1.7);
                vT = 1.04 * tempArray;                    %%% Valeur hypothetique
                alphaT = 100 * alphaL;
                
            case 'Epoxy'
                %%% cf Sayers, J. Phys. D (1983)
                rho = 1.20 * tempArray;
                vL = 2.64 * tempArray;
                alphaL = 0.0 * freq;
                vT = 1.20 * tempArray;
                alphaT = 0 * alphaL;
                
            case 'Aqualene'
                %%% cf Olympus
                rho = 0.92 * tempArray;
                vL = 1.59 * tempArray;
                alphaL = 0.003 * freq.^(1.5);
                vT = 0.80 * tempArray;
                alphaT = 2 * alphaL;                      %%% Valeur hypothetique
                
            case 'Cemcat'
                %%% Polyuréthane --- Mesures du 26 juin 2012 1MHz
                rho = 1.02 * tempArray;
                vL = 1.41 * tempArray;
                alphaL = freq.^(2);
                vT = 0.110 * tempArray;                   %%% Cf These G. Lepert
                alphaT = 200 * alphaL;                    %%% Valeur hypothetique
                %vT = freq.^(2);
                %mu = rho.*(vT.^2).*(1-i*freq.^(1/6)/3);
                %alphaT = imag(2*pi*freq.*sqrt(rho./mu)) ;
                
            case 'Polystyrene'
                rho = 1.05 * tempArray;
                vL = 2.35 * tempArray;
                alphaL = 0.0 * freq;
                vT = 1.12 * tempArray;
                alphaT = 0.0 * freq;                      %%% Pas d'onde transverse
                
            case 'Cork'
                %%% avec E = 18.6 MPa & nu = 0
                rho = 0.23 * tempArray;
                vL = 0.284 * tempArray;
                alphaL = freq.^(2);                       %%% Valeur hypothetique
                vT = 0.201 * tempArray;
                alphaT = 200 * alphaL;                    %%% Valeur hypothetique
                
                
            case 'Water'
                rho = 1.00 * tempArray;
                vL = 1.5 * tempArray;
                alphaL = 0.000003 * freq.^(2);
                vT = 0.0 * tempArray;                     %%% Pas d'onde transverse
                alphaT = 0.0 * tempArray;                 %%% Pas d'onde transverse
                
                
            case 'Agar'
                %%% cf Bouazzak, J. Phys. IV (1994)
                rho = 1.005 * tempArray;
                vL = 1.490 * tempArray;
                alphaL = 0.00000633 * freq.^(2);
                vT = 0.0 * tempArray;                     %%% Pas d'onde transverse
                alphaT = 0.0 * freq;                      %%% Pas d'onde transverse
                
            case '47V20'
                %%% cf Samet, CFA 2010
                rho = 0.95 * tempArray;
                vL = 0.976 * tempArray;
                alphaL = 0.0007 * freq.^(1.87);
                vT = 0.0 * tempArray;                     %%% Pas d'onde transverse
                alphaT = 0.0 * freq;                      %%% Pas d'onde transverse
                
            case 'HMDSO'                                  %%% (Hexamethyldisiloxane)
                rho = 0.760 * tempArray;
                vL = 0.894 * tempArray;
                alphaL = 0.0007 * freq.^(1.87);
                vT = 0.0 * tempArray;                     %%% Pas d'onde transverse
                alphaT = 0.0 * freq;                      %%% Pas d'onde transverse
                
            case 'Hair gel'
                %%% Gel cheveux Auchan dilué
                rho = 1.005 * tempArray;
                vL = 1.48 * tempArray;%1.492 * tempArray;
                alphaL = 0.000083 * freq.^(2);
                vT = 0.0 * tempArray;                     %%% Pas d'onde transverse
                alphaT = 0.0 * tempArray;                 %%% Pas d'onde transverse
                
                
            case 'Xanthan'
                %%% Gel cheveux Auchan dilué
                rho = 1.005 * tempArray;
                vL = 1.52 * tempArray;%1.492 * tempArray;
                alphaL = 0.000083 * freq.^(2);
                vT = 0.0 * tempArray;                     %%% Pas d'onde transverse
                alphaT = 0.0 * tempArray;                 %%% Pas d'onde transverse
                
            case 'PMMA'
                %%% cf He, J. Acoust. Soc. Am. (2000)
                rho = 1.176 * tempArray;
                vL = 2.75 * tempArray;
                alphaL = 10.^(0.94 * freq/20)/10;
                vT = 1.34 * tempArray;                    %%% Leymarie, Thèse 2002
                alphaT = 2 * alphaL;                      %%% Valeur hypothetique
                
            case 'Glass'
                %%% cf K. Hildebrand's thesis
                rho = 2.50 * tempArray;
                vL = 5.60 * tempArray;
                alphaL = 7e-7 * freq.^(0);                %%% Valeur donnée à 2MHz
                vT = 3.4 * tempArray;
                alphaT = 0.0 * freq;                      %%% Attenuation negligee
                
            case 'Steel'
                %%% cf Onda_Solids ("Steel - mild")
                rho = 7.8 * tempArray;
                vL = 5.9 * tempArray;
                alphaL = 0.0 * freq;                      %%% Attenuation negligee
                vT = 3.2 * tempArray;
                alphaT = 0.0 * freq;                      %%% Attenuation negligee
                
            case 'Aluminum'
                %%% cf Onda_Solids ("Steel - mild")
                rho = 2.71 * tempArray;
                vL = 6.40 * tempArray;
                alphaL = 0.0 * freq;                      %%% Attenuation negligee
                vT = 3.12 * tempArray;
                alphaT = 0.0 * freq;                      %%% Attenuation negligee
                
            case 'Tungstene'
                %%% cf Onda_Solids ("Tungsten")
                rho = 19.4 * tempArray;
                vL = 5.2 * tempArray;
                alphaL = 0.0 * freq;                      %%% Attenuation negligee
                vT = 2.9 * tempArray;
                alphaT = 0.0 * freq;                      %%% Attenuation negligee
                
            case 'Lead'
                %%% cf Onda_Solids ("Lead")
                rho = 11.2 * tempArray;
                vL = 2.2 * tempArray;
                alphaL = 0.0 * freq;                      %%% Attenuation negligee
                vT = 0.7 * tempArray;
                alphaT = 0.0 * freq;                      %%% Attenuation negligee
                
            case 'Bismuth'
                %%% cf Onda_Solids ("Bismuth")
                rho = 9.8 * tempArray;
                vL = 2.2 * tempArray;
                alphaL = 0.0 * freq;                      %%% Attenuation negligee
                vT = 1.1 * tempArray;
                alphaT = 0.0 * freq;                      %%% Attenuation negligee
                
            case 'Baryte'
                rho = 4.3 * tempArray;
                vL = 7.0 * tempArray;
                alphaL = 0.0 * freq;                      %%% Attenuation negligee
                vT = 3.7 * tempArray;
                alphaT = 0.0 * freq;                      %%% Attenuation negligee
                
            case 'Tin'
                %%% cf Onda_Solids ("tin")
                rho = 7.3 * tempArray;
                vL = 3.3 * tempArray;
                alphaL = 0.0 * freq;                      %%% Attenuation negligee
                vT = 1.7 * tempArray;
                alphaT = 0.0 * freq;                      %%% Attenuation negligee
                
            case 'BiPbInSnCd'
                %%% cf Indium Corporation ("Indalloy 177")
                rho = 9.16 * tempArray;
                vL = 2.37 * tempArray;                    %%% mesure au labo
                alphaL = 0.0 * freq;                      %%% Attenuation negligee
                vT = 1.05 * tempArray;                    %%% mesure au labo
                alphaT = 0.0 * freq;                      %%% Attenuation negligee
                
            case 'Mercury'
                %%% cf Onda_Liquids (Mercury at 25.0°C)
                rho = 13.5 * tempArray;
                vL = 1.43 * tempArray;
                alphaL = 0.0 * freq;                      %%% Attenuation negligee
                vT = 0.0 * tempArray;
                alphaT = 0.0 * freq;                      %%% Attenuation negligee
                
            case 'GaIn'
                %%% cf Morley, Rev. Sci. Instrum. (2008)
                rho = 6.28 * tempArray;
                vL = 2.74 * tempArray;
                alphaL = 0.0 * freq;                      %%% Attenuation negligee
                vT = 0.0 * tempArray;
                alphaT = 0.0 * freq;                      %%% Attenuation negligee
                
            case 'FC40'
                %%% cf Fluorinert & Onda_Liquids
                rho = 1.85 * tempArray;
                vL = 0.64 * tempArray;
                alphaL = 0.001 * freq.^(2);               %%% cf These Tallon
                vT = 0.0 * tempArray;                     %%% Pas d'onde transverse
                alphaT = 0.0 * freq;                      %%% Pas d'onde transverse
                
            case 'FC75'
                %%% cf Fluorinert & Onda_Liquids
                rho = 1.76 * tempArray;  % Mesuré
                vL = 0.585 * tempArray;
                alphaL = 0.0 * freq.^(0);
                vT = 0.0 * tempArray;                     %%% Pas d'onde transverse
                alphaT = 0.0 * freq;                      %%% Pas d'onde transverse
                
            case 'Ferrofluid'
                %%% cf Brunet, Phys. Rev. Lett. (2013)
                rho = 1.9 * tempArray;
                vL = 0.670 * tempArray;
                alphaL = 0.0057 * freq.^(-2);
                vT = 0.0 * tempArray;                     %%% Pas d'onde transverse
                alphaT = 0.0 * freq;                      %%% Pas d'onde transverse
                
            case 'Air'
                %%% cf Onda_Gases ("Air - at 20°C")
                rho = 0.001293 * tempArray;
                vL = 0.344 * tempArray;
                alphaL = 0.0 * freq;
                vT = 0.0 * tempArray;                     %%% Pas d'onde transverse
                alphaT = 0.0 * freq;                      %%% Pas d'onde transverse
                
            case 'SF6'
                %%% cf http://www.sf6.fr
                rho = 0.00627 * tempArray;
                vL = 0.136 * tempArray;
                vT = 0.0 * tempArray;                     %%% Pas d'onde transverse
                alphaT = 0.0 * freq;                      %%% Pas d'onde transverse
                
            case 'Aerogel'
                %%% cf Iino, Acoust. Sci. & Tech. (2011)
                rho = 0.25 * tempArray;
                vL = 0.15 * tempArray;
                alphaL = 0.4949 * (2.5119) * freq.^(1.9);
                vT = 0.058 * tempArray;                   %%% d apres coeff Poisson
                alphaT = 10 * alphaL;                     %%% Valeur hypothetique
                
                
            case 'Xerogel'
                %%% Donneés extraites de fit de BM (réunion ANR 2014-02-13)
                rho = 0.6 * tempArray;                    %%% thèse Simon (2014)
                vL = 0.330 * tempArray;
                alphaL = 0.12 * freq.^(1.1);              %%% cf Gomez, APL (2002)
                vT = 0.183 * tempArray;                   %%% nu = 0.2780
                alphaT = 0.226 * freq.^(0.5);             %%% cf Gomez, APL (2002)
                
            case 'HIPE'
                %%% cf Brunet, Nature Materials 2015
                rho = 0.60 * tempArray;                   %%% porosite = 40%
                vL = 0.08 * tempArray;
                alphaL = 60 * freq.^(1.5);
                vT = 0.040 * tempArray;                   %%% nu = 0.3333
                alphaT = 200 * freq.^(1.5);
                
            case 'Rubber'
                %%% cf Liu, Science (2000)
                rho = 1.3 * tempArray;
                vL = 0.02 * tempArray;
                alphaL = 0 * freq.^(2);
                vT = 0.0 * tempArray;                     %%% nu = 0.5
                alphaT = 0.0 * freq;                      %%% Pas d'onde transverse
                
            case 'Nylon'
                %%% cf Croenne, AIP Advances (2011)
                rho = 1.15 * tempArray;
                vL = 2.4 * tempArray;
                alphaL = 0 * freq.^(2);
                vT = 0.95 * tempArray;                     %%% nu = 0.5
                alphaT = 0.0 * freq;                      %%% Pas d'onde transverse
                
        end
        
    end

end