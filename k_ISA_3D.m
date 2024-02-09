function [f,k0a,vl0,kISA,vph,vgr,ve,ltr] = k_ISA_3D(M0,M1,M2,N,f,Rm,P,phi,Nr,Rwidth)

% Mat = materiau pour la matrice
% Inc = materiau pour les inclusions
% N = nombre de modes
% fMAX = freq maximum pour le calcul
% amoy = rayon moyen sphere (mm)
% P = polydispersite
% phi = fraction volumique rho2,vl2,vt2,alphal2,alphat2

f(f==0) = 1e-5;

%% Entrees
phi = phi/100;            % fraction volumique en %
P = P/100;                % polydispersite en %

rho0 = M0(1); vl0 = M0(2); alphal0 = M0(3);
rho1 = M1(1); vl1 = M1(2); vt1 = M1(3); alphal1 = M1(4); alphat1 = M1(5);
rho2 = M2(1); vl2 = M2(2); vt2 = M2(3); alphal2 = M2(4); alphat2 = M1(5);

[kISA,ltr,vl0,vph,vgr,ve] = fct_Calc_k_ISA(Rm,f,N,phi,P,rho0,vl0,alphal0,rho1,vl1,vt1,alphal1,alphat1,rho2,vl2,vt2,alphal2,alphat2,Nr,Rwidth);

%%

k0a = real(2*pi*f./real(vl0)*Rm(1));
vl0 = real(vl0(1));

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Sub-routines

function [kISA,ltr,vl0,vph,vgr,ve] = fct_Calc_k_ISA(Rm,f,N,phi,P,rho0,vl0,alphal0,rho1,vl1,vt1,alphal1,alphat1,rho2,vl2,vt2,alphal2,alphat2,Nr,Rwidth)
%%

if P == 0
    
    sig=0;
    r=1;
    eta=phi/(4/3.*pi*(Rm(1)*r)^3);
    
else
    
    sig=P;
    r=linspace(1 - Rwidth*sig,1 + Rwidth*sig,Nr);
    
    %Distribution gaussienne
    
    gauss = exp(-0.5*((r-1)/sig).^2);
    gaussn =( 1/sum(gauss)) * gauss;
    
    eta = gaussn.* phi / sum( gaussn .* ( (4/3)*pi*((Rm(1)*r).^3) ) );
    
end

omega=2*pi*f;
Nr = length(r);
kL0 = omega./vl0 + 1i*alphal0.*f.^(2);

Nf=length(f);
A = zeros(N,Nf,Nr);
F0 = zeros(Nr,Nf);

for ir=1:Nr
    
    [S] = Calc_Scatt_Amp_Layers(vl0,vl1,vl2,vt1,vt2,rho0,rho1, ...
        rho2,alphal0,alphal1,alphat1,alphal2,alphat2,omega,Rm(1)*r(ir),Rm(2)*r(ir),Rm(3)*r(ir),N);
    
    A(:,:,ir) = S;
    
    F0(ir,:)=foncDiff(A(:,:,ir),0,kL0);                                    % entrees : (An,theta,kL0)
    
end

%% K effectif (ISA)

kISA=sqrt(kL0.^2+4*pi*eta*F0);

kpoly = sqrt(repmat(kL0,[Nr,1]).^2+4.*pi.*repmat(eta.',[1,Nf]).*F0);      % <!> premiere version kpoly etait definie a partir de kISA

vgr = gradient(omega,real(kISA)) ;
vgrpoly = diff(repmat(omega.',[1,Nr])).'./diff(real(kpoly).').';

%% Calcul du <cos>

%%%%%%%%%%%%%%%%%%%%%%%% Moyenne de <cos> sur les rayons %%%%%%%%%%%%%%%%%%

cosMoy = zeros(Nf,Nr);

for ir = 1:Nr
    
    cosMoy(:,ir) = calcul_cosMoy(A(:,:,ir),kL0);
    
end

ltr = 1./real(4*pi./kL0.*(eta*(imag(F0).*(1-cosMoy.'))));                  % l* Avec le theoreme optique

%% Vitesse de l'energie

dtheta = 1e-1;
theta = 0 : dtheta : pi;
Fdiff = zeros(length(f),Nr,length(theta));

for ir = 1 : Nr
    
    for it = 1 : length(theta)
        
        Fdiff(:,ir,it)= foncDiff(A(:,:,ir),theta(it),kL0);
        
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

vph = omega.'./real(kISA.');

eta = repmat(eta,[Nf-1,1]);
F0 = F0.';
vgrpoly = vgrpoly.';
omega = repmat(omega,[Nr,1]).';
kpoly = kpoly.';

Ind = vl0.^2./vph(1 : Nf - 1);

d1 = (diff(real(F0))./diff(omega))./real(kpoly(1 : Nf  - 1,:)).*vgrpoly;

theta = reshape(theta,1,1,length(theta));
theta = repmat(theta,[Nf - 1,Nr,1]);

d2 = (Fdiff(1:length(f)-1,:,:).*conj(Fdiff(1:length(f)-1,:,:)))...
    .*(diff(unwrap(angle(-Fdiff)))./diff(repmat(omega,[1,1,size(theta,3)])))...
    .*repmat(vgrpoly,[1,1,size(theta,3)]);

d2 = squeeze(sum(d2.*dtheta.*sin(theta),3));

dve = sum(2.*pi.*eta.*(d1 + d2),2);

ve = Ind./(1 + dve);
ve = cat(1,real(ve),real(ve(end)));

end


function En = Calc_Scatt_Amp_Layers(v0l,v1l,v2l,v1t,v2t,rho0,rho1,rho2,alphal0,alphal1,alphat1,alphal2,alphat2,om,R1,R2,R3,N)

% resolution du systeme [M] . [A] = [C]

%% Amplitude onde diffusee core-shell dans liquide.

Nf = length(om);

rho3 = 0; v3l = 0; alphal3 = 0;                                            % to be done

k0l = om./v0l + alphal0.*(om./(2*pi)).^2.*1i;
k1l = om./v1l + alphal1.*(om./(2*pi)).^2.*1i;
k1t = om./v1t + alphat1.*(om./(2*pi)).^2.*1i;
k2l = om./v2l + alphal2.*(om./(2*pi)).^2.*1i;
k2t = om./v2t + alphat2.*(om./(2*pi)).^2.*1i;
k3l = om./v3l + alphal3.*(om./(2*pi)).^2.*1i;

mu1 = om.^2.*rho1./k1t.^2;                                                 % shear modulus
mu2 = om.^2.*rho2./k2t.^2;                                                 % shear modulus

%% %%%%%%%%%%%%%%%%%%%%%Calcul scattering coefficients %%%%%%%%%%%%%%%%%%%%

En = zeros(N,Nf);                                                          % Coeff de scattering de l'onde diffusee

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Plain solid %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (R2 == 0) && (R3 == 0) && (v1t~=0)
    
    for n = 0:N-1
        
        %% Vecteur S
        
        S1 = k0l.*besselSp(n,k0l.*R1,1);
        S2 = -om.^2.*rho0.*besselS(n,k0l.*R1,1);
        
        %% Matrice A
        
        m11 = k1l.*besselSp(n,k1l.*R1,1);
        m12 = -n.*(n+1)/R1.*besselS(n,k1t.*R1,1);
        m13 = -k0l.*besselSp(n,k0l*R1,2);
        
        m21 = besselS(n,k1l.*R1,1).*(2.*mu1.*k1l.^2-rho1.*om.^2) + 2.*mu1.*k1l.^2.*besselSpp(n,k1l.*R1,1);
        m22 = 2.*mu1.*n.*(n+1)./R1.*(besselS(n,k1t.*R1,1)./R1 - k1t.*besselSp(n,k1t.*R1,1));
        m23 = om.^2.*rho0.*besselS(n,k0l.*R1,2);
        
        m31 = 2.*(k1l.*besselSp(n,k1l.*R1,1)./R1 - besselS(n,k1l.*R1,1)./(R1.^2));
        m32 = -(k1t.^2.*besselSpp(n,k1t.*R1,1) - besselS(n,k1t.*R1,1).*(2-n.*(n+1))./(R1.^2));
        
        %% Methode de Cramer
        
        for jj = 1 : length(om)
            
            S = [S1(jj),S2(jj),0];
            
            MM = [m11(jj),m12(jj),m13(jj); ...
                m21(jj),m22(jj),m23(jj); ...
                m31(jj),m32(jj),0];
            
            M9 = MM;
            M9(:,3) = S;
            
            En(n+1,jj) = det(M9)./det(MM);                                 % Cramer
            
        end
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%% Plain fluid %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
elseif (R2 == 0) && (R3 == 0) && (v1t==0)
    
    for n = 0:N-1
        
        %% Vecteur S
        
        S1 = k0l.*besselSp(n,k0l.*R1,1);
        S2 = rho0.*besselS(n,k0l.*R1,1);
        
        %% Matrice A
        
        m11 = k1l.*besselSp(n,k1l.*R1,1);
        m12 = -k0l.*besselSp(n,k0l*R1,2);
        
        m21 = besselS(n,k1l.*R1,1).*rho1;
        m22 = -rho0.*besselS(n,k0l.*R1,2);
        
        %% Methode de Cramer
        
        En(n+1,:) = (m11.*S2-m21.*S1)./(m11.*m22-m21.*m12);
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%% Core-shell %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
elseif (R2 ~= 0) && (R3 == 0)
    
    if (v1t~=0) && (v2t~=0)                                                % Solid/Solid
        
        IndSol1 = ones(1,Nf);
        IndSol2 = ones(1,Nf);
        
    elseif (v1t~=0) && (v2t==0)                                            % Solid/Fluid
        
        IndSol1 = ones(1,Nf);
        mu2 = zeros(1,Nf); k2t = ones(1,Nf); IndSol2 = zeros(1,Nf);
        
    elseif (v1t==0) && (v2t~=0)                                            % Fluid/Solid
        
        mu1 = zeros(1,Nf); k1t = ones(1,Nf); IndSol1 = zeros(1,Nf);
        IndSol2 = ones(1,Nf);
        
    elseif (v1t==0) && (v2t==0)                                            % Fluid/Fluid
        
        mu1 = zeros(1,Nf); k1t = ones(1,Nf); IndSol1 = zeros(1,Nf);
        mu2 = zeros(1,Nf); k2t = ones(1,Nf); IndSol2 = zeros(1,Nf);
        
        
    end
    
    for n = 0:N-1
        
        %% Vecteur S
        
        S1 = k0l.*besselSp(n,k0l.*R1,1);
        S2 = -om.^2.*rho0.*besselS(n,k0l.*R1,1);
        
        %% Matrice A
        
        m11 = k1l.*besselSp(n,k1l.*R1,1);
        m12 = k1l.*besselSp(n,k1l.*R1,2);
        m13 = -n.*(n+1)/R1.*besselS(n,k1t.*R1,1).*IndSol1;
        m14 = -n.*(n+1)/R1.*besselS(n,k1t.*R1,2).*IndSol1;
        m19 = -k0l.*besselSp(n,k0l*R1,2);
        
        m21 = besselS(n,k1l.*R1,1).*(2.*mu1.*k1l.^2-rho1.*om.^2) + 2.*mu1.*k1l.^2.*besselSpp(n,k1l.*R1,1);
        m22 = besselS(n,k1l.*R1,2).*(2.*mu1.*k1l.^2-rho1.*om.^2) + 2.*mu1.*k1l.^2.*besselSpp(n,k1l.*R1,2);
        m23 = 2.*mu1.*n.*(n+1)./R1.*(besselS(n,k1t.*R1,1)./R1 - k1t.*besselSp(n,k1t.*R1,1));
        m24 = 2.*mu1.*n.*(n+1)./R1.*(besselS(n,k1t.*R1,2)./R1 - k1t.*besselSp(n,k1t.*R1,2));
        m29 = om.^2.*rho0.*besselS(n,k0l.*R1,2);
        
        m31 = 2.*(k1l.*besselSp(n,k1l.*R1,1)./R1 - besselS(n,k1l.*R1,1)./(R1.^2));
        m32 = 2.*(k1l.*besselSp(n,k1l.*R1,2)./R1 - besselS(n,k1l.*R1,2)./(R1.^2));
        m33 = -(k1t.^2.*besselSpp(n,k1t.*R1,1) - besselS(n,k1t.*R1,1).*(2-n.*(n+1))./(R1.^2)).*IndSol1;
        m34 = -(k1t.^2.*besselSpp(n,k1t.*R1,2) - besselS(n,k1t.*R1,2).*(2-n.*(n+1))./(R1.^2)).*IndSol1;
        
        m41 = k1l.*besselSp(n,k1l.*R2,1);
        m42 = k1l.*besselSp(n,k1l.*R2,2);
        m43 = -n.*(n+1)/R2.*besselS(n,k1t.*R2,1).*IndSol1;
        m44 = -n.*(n+1)/R2.*besselS(n,k1t.*R2,2).*IndSol1;
        m45 = -k2l.*besselSp(n,k2l.*R2,1);
        m47 = n.*(n+1)/R2.*besselS(n,k2t.*R2,1).*IndSol2;
        
        m51 = besselS(n,k1l.*R2,1)./R2;
        m52 = besselS(n,k1l.*R2,2)./R2;
        m53 = -(k1t.*besselSp(n,k1t.*R2,1) + besselS(n,k1t.*R2,1)./R2).*IndSol1;
        m54 = -(k1t.*besselSp(n,k1t.*R2,2) + besselS(n,k1t.*R2,2)./R2).*IndSol1;
        m55 = -besselS(n,k2l.*R2,1)./R2;
        m57 = (k2t.*besselSp(n,k2t.*R2,1) + besselS(n,k2t.*R2,1)./R2).*IndSol2;
        
        m61 = besselS(n,k1l.*R2,1).*(2.*mu1.*k1l.^2-rho1.*om.^2) + 2.*mu1.*k1l.^2.*besselSpp(n,k1l.*R2,1);
        m62 = besselS(n,k1l.*R2,2).*(2.*mu1.*k1l.^2-rho1.*om.^2) + 2.*mu1.*k1l.^2.*besselSpp(n,k1l.*R2,2);
        m63 = 2.*mu1.*n.*(n+1)./R2.*(besselS(n,k1t.*R2,1)./R2 - k1t.*besselSp(n,k1t.*R2,1));
        m64 = 2.*mu1.*n.*(n+1)./R2.*(besselS(n,k1t.*R2,2)./R2 - k1t.*besselSp(n,k1t.*R2,2));
        m65 = -(besselS(n,k2l.*R2,1).*(2.*mu2.*k2l.^2-rho2.*om.^2) + 2.*mu2.*k2l.^2.*besselSpp(n,k2l.*R2,1));
        m67 = -2.*mu2.*n.*(n+1)./R2.*(besselS(n,k2t.*R2,1)./R2 - k2t.*besselSp(n,k2t.*R2,1)).*IndSol2;
        
        m71 = 2.*(k1l.*besselSp(n,k1l.*R2,1)./R2 - besselS(n,k1l.*R2,1)./(R2.^2)).*mu1;
        m72 = 2.*(k1l.*besselSp(n,k1l.*R2,2)./R2 - besselS(n,k1l.*R2,2)./(R2.^2)).*mu1;
        m73 = -(k1t.^2.*besselSpp(n,k1t.*R2,1) - besselS(n,k1t.*R2,1).*(2-n.*(n+1))./(R2.^2)).*mu1.*IndSol1;
        m74 = -(k1t.^2.*besselSpp(n,k1t.*R2,2) - besselS(n,k1t.*R2,2).*(2-n.*(n+1))./(R2.^2)).*mu1.*IndSol1;
        m75 = -2.*(k2l.*besselSp(n,k2l.*R2,1)./R2 - besselS(n,k2l.*R2,1)./(R2.^2)).*mu2;
        m77 = (k2t.^2.*besselSpp(n,k2t.*R2,1) - besselS(n,k2t.*R2,1).*(2-n.*(n+1))./(R2.^2)).*mu2.*IndSol2;
        
        %% Methode de Cramer
      
        for jj = 1 : length(om)
            
            if (v1t~=0) && (v2t~=0)                                        % Solid/Solid
                
                S = [S1(jj),S2(jj),0,0,0,0,0];
                
                MM = [m11(jj),m12(jj),m13(jj),m14(jj),0,0,m19(jj); ...
                    m21(jj),m22(jj),m23(jj),m24(jj),0,0,m29(jj); ...
                    m31(jj),m32(jj),m33(jj),m34(jj),0,0,0; ...
                    m41(jj),m42(jj),m43(jj),m44(jj),m45(jj),m47(jj),0; ...
                    m51(jj),m52(jj),m53(jj),m54(jj),m55(jj),m57(jj),0; ...
                    m61(jj),m62(jj),m63(jj),m64(jj),m65(jj),m67(jj),0; ...
                    m71(jj),m72(jj),m73(jj),m74(jj),m75(jj),m77(jj),0];
                
                M9 = MM;
                M9(:,7) = S;
                
                En(n+1,jj) = det(M9)./det(MM);                             % Cramer
                
            elseif (v1t~=0) && (v2t==0)                                    % Solid/Fluid
                
                S = [S1(jj),S2(jj),0,0,0,0];
                
                MM = [m11(jj),m12(jj),m13(jj),m14(jj),0,m19(jj); ...
                    m21(jj),m22(jj),m23(jj),m24(jj),0,m29(jj); ...
                    m31(jj),m32(jj),m33(jj),m34(jj),0,0; ...
                    m41(jj),m42(jj),m43(jj),m44(jj),m45(jj),0; ...
                    m61(jj),m62(jj),m63(jj),m64(jj),m65(jj),0; ...
                    m71(jj),m72(jj),m73(jj),m74(jj),0,0];
                
                M9 = MM;
                M9(:,6) = S;
                
                En(n+1,jj) = det(M9)./det(MM);                             % Cramer
                
            elseif (v1t==0) && (v2t~=0)                                    % Fluid/Solid
                
                S = [S1(jj),S2(jj),0,0,0];
                
                MM = [m11(jj),m12(jj),0,0,m19(jj); ...
                    m21(jj),m22(jj),0,0,m29(jj); ...
                    m41(jj),m42(jj),m45(jj),m47(jj),0; ...
                    m61(jj),m62(jj),m65(jj),m67(jj),0; ...
                    0,0,m75(jj),m77(jj),0];
                
                M9 = MM;
                M9(:,5) = S;
                
                En(n+1,jj) = det(M9)./det(MM);                             % Cramer
                
            elseif (v1t==0) && (v2t==0)                                    % Fluid/Fluid
                
                S = [S1(jj),S2(jj),0,0];
                
                MM = [m11(jj),m12(jj),0,m19(jj); ...
                    m21(jj),m22(jj),0,m29(jj); ...
                    m41(jj),m42(jj),m45(jj),0; ...
                    m61(jj),m62(jj),m65(jj),0];
                
                M9 = MM;
                M9(:,4) = S;
                
                En(n+1,jj) = det(M9)./det(MM);                             % Cramer
                
            end            
            
        end
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%% solide/fluide/vide %%%%%%%%%%%%%%%%%%%%%%%%
    
elseif (R1 ~= R2) && (R2 ~= R3) && (v1t~=0) && (v2t==0)
    
    for n = 0:N-1
        
        %% Vecteur S
        
        S1 = k0l.*besselSp(n,k0l.*R1,1);
        S2 = -om.^2.*rho0.*besselS(n,k0l.*R1,1);
        
        %% Matrice A
        
        m11 = k1l.*besselSp(n,k1l.*R1,1);
        m12 = k1l.*besselSp(n,k1l.*R1,2);
        m13 = -n.*(n+1)/R1.*besselS(n,k1t.*R1,1);
        m14 = -n.*(n+1)/R1.*besselS(n,k1t.*R1,2);
        m19 = -k0l.*besselSp(n,k0l*R1,2);
        
        m21 = besselS(n,k1l.*R1,1).*(2.*mu1.*k1l.^2-rho1.*om.^2) + 2.*mu1.*k1l.^2.*besselSpp(n,k1l.*R1,1);
        m22 = besselS(n,k1l.*R1,2).*(2.*mu1.*k1l.^2-rho1.*om.^2) + 2.*mu1.*k1l.^2.*besselSpp(n,k1l.*R1,2);
        m23 = 2.*mu1.*n.*(n+1)./R1.*(besselS(n,k1t.*R1,1)./R1 - k1t.*besselSp(n,k1t.*R1,1));
        m24 = 2.*mu1.*n.*(n+1)./R1.*(besselS(n,k1t.*R1,2)./R1 - k1t.*besselSp(n,k1t.*R1,2));
        m29 = om.^2.*rho0.*besselS(n,k0l.*R1,2);
        
        m31 = 2.*(k1l.*besselSp(n,k1l.*R1,1)./R1 - besselS(n,k1l.*R1,1)./(R1.^2));
        m32 = 2.*(k1l.*besselSp(n,k1l.*R1,2)./R1 - besselS(n,k1l.*R1,2)./(R1.^2));
        m33 = -(k1t.^2.*besselSpp(n,k1t.*R1,1) - besselS(n,k1t.*R1,1).*(2-n.*(n+1))./(R1.^2));
        m34 = -(k1t.^2.*besselSpp(n,k1t.*R1,2) - besselS(n,k1t.*R1,2).*(2-n.*(n+1))./(R1.^2));
        
        m41 = k1l.*besselSp(n,k1l.*R2,1);
        m42 = k1l.*besselSp(n,k1l.*R2,2);
        m43 = -n.*(n+1)/R2.*besselS(n,k1t.*R2,1);
        m44 = -n.*(n+1)/R2.*besselS(n,k1t.*R2,2);
        m45 = -k2l.*besselSp(n,k2l.*R2,1);
        m46 = -k2l.*besselSp(n,k2l.*R2,2);
        
        m61 = besselS(n,k1l.*R2,1).*(2.*mu1.*k1l.^2-rho1.*om.^2) + 2.*mu1.*k1l.^2.*besselSpp(n,k1l.*R2,1);
        m62 = besselS(n,k1l.*R2,2).*(2.*mu1.*k1l.^2-rho1.*om.^2) + 2.*mu1.*k1l.^2.*besselSpp(n,k1l.*R2,2);
        m63 = 2.*mu1.*n.*(n+1)./R2.*(besselS(n,k1t.*R2,1)./R2 - k1t.*besselSp(n,k1t.*R2,1));
        m64 = 2.*mu1.*n.*(n+1)./R2.*(besselS(n,k1t.*R2,2)./R2 - k1t.*besselSp(n,k1t.*R2,2));
        m65 = -(besselS(n,k2l.*R2,1).*(-rho2.*om.^2));
        m66 = -(besselS(n,k2l.*R2,2).*(-rho2.*om.^2));
        
        m85 = besselS(n,k2l.*R3,1).*(-rho2.*om.^2);
        m86 = besselS(n,k2l.*R3,2).*(-rho2.*om.^2);
        
        m95 = 2.*(k2l.*besselSp(n,k2l.*R3,1)./R3 - besselS(n,k2l.*R3,1)./(R3.^2));
        m96 = 2.*(k2l.*besselSp(n,k2l.*R3,2)./R3 - besselS(n,k2l.*R3,2)./(R3.^2));
        
        %% Methode de Cramer
        
        zz = zeros(1,Nf);
        
        hcat1 = cat(1,m11,m12,m13,m14,zz,zz,m19);
        hcat2 = cat(1,m21,m22,m23,m24,zz,zz,m29);
        hcat3 = cat(1,m31,m32,m33,m34,zz,zz,zz);
        hcat4 = cat(1,m41,m42,m43,m44,m45,m46,zz);
        hcat6 = cat(1,m61,m62,m63,m64,m65,m66,zz);
        hcat8 = cat(1,zz,zz,zz,zz,m85,m86,zz);
        hcat9 = cat(1,zz,zz,zz,zz,m95,m96,zz);
        
        MM = cat(3,hcat1,hcat2,hcat3,hcat4,hcat6,hcat8,hcat9);
        MM = permute(MM,[3,1,2]);
        
        M9 = MM;
        
        S = cat(1,S1,S2,zz,zz,zz,zz,zz);
        M9(:,7,:) = S;
        
        for jj = 1 : length(om)
            
            En(n+1,jj) = det(M9(:,:,jj))./det(MM(:,:,jj));                 % Cramer
            
        end
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%% solide/solide/fluide %%%%%%%%%%%%%%%%%%%%%%%%%
    
else
    
    for n = 0:N-1
        
        %% Vecteur S
        
        S1 = k0l.*besselSp(n,k0l.*R1,1);
        S2 = -om.^2.*rho0.*besselS(n,k0l.*R1,1);
        
        %% Matrice A
        
        m11 = k1l.*besselSp(n,k1l.*R1,1);
        m12 = k1l.*besselSp(n,k1l.*R1,2);
        m13 = -n.*(n+1)/R1.*besselS(n,k1t.*R1,1);
        m14 = -n.*(n+1)/R1.*besselS(n,k1t.*R1,2);
        m19 = -k0l.*besselSp(n,k0l*R1,2);
        
        m21 = besselS(n,k1l.*R1,1).*(2.*mu1.*k1l.^2-rho1.*om.^2) + 2.*mu1.*k1l.^2.*besselSpp(n,k1l.*R1,1);
        m22 = besselS(n,k1l.*R1,2).*(2.*mu1.*k1l.^2-rho1.*om.^2) + 2.*mu1.*k1l.^2.*besselSpp(n,k1l.*R1,2);
        m23 = 2.*mu1.*n.*(n+1)./R1.*(besselS(n,k1t.*R1,1)./R1 - k1t.*besselSp(n,k1t.*R1,1));
        m24 = 2.*mu1.*n.*(n+1)./R1.*(besselS(n,k1t.*R1,2)./R1 - k1t.*besselSp(n,k1t.*R1,2));
        m29 = om.^2.*rho0.*besselS(n,k0l.*R1,2);
        
        m31 = 2.*(k1l.*besselSp(n,k1l.*R1,1)./R1 - besselS(n,k1l.*R1,1)./(R1.^2));
        m32 = 2.*(k1l.*besselSp(n,k1l.*R1,2)./R1 - besselS(n,k1l.*R1,2)./(R1.^2));
        m33 = -(k1t.^2.*besselSpp(n,k1t.*R1,1) - besselS(n,k1t.*R1,1).*(2-n.*(n+1))./(R1.^2));
        m34 = -(k1t.^2.*besselSpp(n,k1t.*R1,2) - besselS(n,k1t.*R1,2).*(2-n.*(n+1))./(R1.^2));
        
        m41 = k1l.*besselSp(n,k1l.*R2,1);
        m42 = k1l.*besselSp(n,k1l.*R2,2);
        m43 = -n.*(n+1)/R2.*besselS(n,k1t.*R2,1);
        m44 = -n.*(n+1)/R2.*besselS(n,k1t.*R2,2);
        m45 = -k2l.*besselSp(n,k2l.*R2,1);
        m46 = -k2l.*besselSp(n,k2l.*R2,2);
        m47 = n.*(n+1)/R2.*besselS(n,k2t.*R2,1);
        m48 = n.*(n+1)/R2.*besselS(n,k2t.*R2,2);
        
        m51 = besselSp(n,k1l.*R2,1)./R2;
        m52 = besselSp(n,k1l.*R2,2)./R2;
        m53 = -(k1t.*besselSp(n,k1t.*R2,1) + besselS(n,k1t.*R2,1)./R2);
        m54 = -(k1t.*besselSp(n,k1t.*R2,2) + besselS(n,k1t.*R2,2)./R2);
        m55 = -besselSp(n,k2l.*R2,1)./R2;
        m56 = -besselSp(n,k2l.*R2,2)./R2;
        m57 = (k2t.*besselSp(n,k2t.*R2,1) + besselS(n,k2t.*R2,1)./R2);
        m58 = (k2t.*besselSp(n,k2t.*R2,2) + besselS(n,k2t.*R2,2)./R2);
        
        m61 = besselS(n,k1l.*R2,1).*(2.*mu1.*k1l.^2-rho1.*om.^2) + 2.*mu1.*k1l.^2.*besselSpp(n,k1l.*R2,1);
        m62 = besselS(n,k1l.*R2,2).*(2.*mu1.*k1l.^2-rho1.*om.^2) + 2.*mu1.*k1l.^2.*besselSpp(n,k1l.*R2,2);
        m63 = 2.*mu1.*n.*(n+1)./R2.*(besselS(n,k1t.*R2,1)./R2 - k1t.*besselSp(n,k1t.*R2,1));
        m64 = 2.*mu1.*n.*(n+1)./R2.*(besselS(n,k1t.*R2,2)./R2 - k1t.*besselSp(n,k1t.*R2,2));
        m65 = -(besselS(n,k2l.*R2,1).*(2.*mu2.*k2l.^2-rho2.*om.^2) + 2.*mu2.*k2l.^2.*besselSpp(n,k2l.*R2,1));
        m66 = -(besselS(n,k2l.*R2,2).*(2.*mu2.*k2l.^2-rho2.*om.^2) + 2.*mu2.*k2l.^2.*besselSpp(n,k2l.*R2,2));
        m67 = -2.*mu2.*n.*(n+1)./R2.*(besselS(n,k2t.*R2,1)./R2 - k2t.*besselSp(n,k2t.*R2,1));
        m68 = -2.*mu2.*n.*(n+1)./R2.*(besselS(n,k2t.*R2,2)./R2 - k2t.*besselSp(n,k2t.*R2,2));
        
        m71 = 2.*(k1l.*besselSp(n,k1l.*R2,1)./R2 - besselS(n,k1l.*R2,1)./(R2.^2));
        m72 = 2.*(k1l.*besselSp(n,k1l.*R2,2)./R2 - besselS(n,k1l.*R2,2)./(R2.^2));
        m73 = -(k1t.^2.*besselSpp(n,k1t.*R2,1) - besselS(n,k1t.*R2,1).*(2-n.*(n+1))./(R2.^2));
        m74 = -(k1t.^2.*besselSpp(n,k1t.*R2,2) - besselS(n,k1t.*R2,2).*(2-n.*(n+1))./(R2.^2));
        m75 = -2.*(k2l.*besselSp(n,k2l.*R2,1)./R2 - besselS(n,k2l.*R2,1)./(R2.^2));
        m76 = -2.*(k2l.*besselSp(n,k2l.*R2,2)./R2 - besselS(n,k2l.*R2,2)./(R2.^2));
        m77 = k2t.^2.*besselSpp(n,k2t.*R2,1) - besselS(n,k2t.*R2,1).*(2-n.*(n+1))./(R2.^2);
        m78 = k2t.^2.*besselSpp(n,k2t.*R2,2) - besselS(n,k2t.*R2,2).*(2-n.*(n+1))./(R2.^2);
        
        m85 = besselS(n,k2l.*R3,1).*(2.*mu2.*k2l.^2-rho2.*om.^2) + 2.*mu2.*k2l.^2.*besselSpp(n,k2l.*R3,1);
        m86 = besselS(n,k2l.*R3,2).*(2.*mu2.*k2l.^2-rho2.*om.^2) + 2.*mu2.*k2l.^2.*besselSpp(n,k2l.*R3,2);
        m87 = 2.*mu2.*n.*(n+1)./R3.*(besselS(n,k2t.*R3,1)./R3 - k2t.*besselSp(n,k2t.*R3,1));
        m88 = 2.*mu2.*n.*(n+1)./R3.*(besselS(n,k2t.*R3,2)./R3 - k2t.*besselSp(n,k2t.*R3,2));
        m810 = besselS(n,k3l.*R3,1).*(rho3.*om.^2);
        
        m95 = 2.*(k2l.*besselSp(n,k2l.*R3,1)./R3 - besselS(n,k2l.*R3,1)./(R3.^2));
        m96 = 2.*(k2l.*besselSp(n,k2l.*R3,2)./R3 - besselS(n,k2l.*R3,2)./(R3.^2));
        m97 = -(k2t.^2.*besselSpp(n,k2t.*R3,1) - besselS(n,k2t.*R3,1).*(2-n.*(n+1))./(R3.^2));
        m98 = -(k2t.^2.*besselSpp(n,k2t.*R3,2) - besselS(n,k2t.*R3,2).*(2-n.*(n+1))./(R3.^2));
        
        m105 = -k2l.*besselSp(n,k2l.*R3,1);
        m106 = -k2l.*besselSp(n,k2l.*R3,2);
        m107 = n.*(n+1)/R2.*besselS(n,k2t.*R3,1);
        m108 = n.*(n+1)/R2.*besselS(n,k2t.*R3,2);
        m1010 = -k3l.*besselSp(n,k3l.*R3,1);
        
        %% Methode de Cramer
        
        zz = zeros(1,Nf);
        D = zz;
        
        hcat1 = cat(1,m11,m12,m13,m14,zz,zz,zz,zz,m19,zz);
        hcat2 = cat(1,m21,m22,m23,m24,zz,zz,zz,zz,m29,zz);
        hcat3 = cat(1,m31,m32,m33,m34,zz,zz,zz,zz,zz,zz);
        hcat4 = cat(1,m41,m42,m43,m44,m45,m46,m47,m48,zz,zz);
        hcat5 = cat(1,m51,m52,m53,m54,m55,m56,m57,m58,zz,zz);
        hcat6 = cat(1,m61,m62,m63,m64,m65,m66,m67,m68,zz,zz);
        hcat7 = cat(1,m71,m72,m73,m74,m75,m76,m77,m78,zz,zz);
        hcat8 = cat(1,zz,zz,zz,zz,m85,m86,m87,m88,zz,m810);
        hcat9 = cat(1,zz,zz,zz,zz,m95,m96,m97,m98,zz,zz);
        hcat10 = cat(1,zz,zz,zz,zz,m105,m106,m107,m108,zz,m1010);
        
        MM = cat(3,hcat1,hcat2,hcat3,hcat4,hcat5,hcat6,hcat7,hcat8,hcat9,hcat10);
        MM = permute(MM,[3,1,2]);
        
        M9 = MM;
        
        S = cat(1,S1,S2,zz,zz,zz,zz,zz,zz,zz,zz);
        M9(:,9,:) = S;
        
        for jj = 1 : length(om)
            
            En(n+1,jj) = det(M9(:,:,jj))./det(MM(:,:,jj));                 % Cramer
            
        end
        
    end
    
end

end

%% Fonctions locales

%% Calcul fonction de diffusion
function fTheta=foncDiff(an,theta,k0)

N=size(an,1);
fTheta=0;
P = ones(length(theta),3);

for n=0:N-1
    
    if n == 1
        
        P(:,2) = ones(size(theta));
        P(:,3) = cos(theta).';
        
    elseif n > 1
        
        %         P(:,3) = LegendreB(cos(theta),n,squeeze(P(:,2)),squeeze(P(:,1)));
        N  = n-1;
        A1 = (2 * N + 1) ./ (N + 1);
        A2 = N ./ (N + 1);
        P(:,3) = A1.* cos(theta).'.* squeeze(P(:,2)) - A2.* squeeze(P(:,1));
        
    end
    
    %     fTheta=(2.*n+1).*an(n+1,:).*legN(n,cos(theta))+fTheta;
    fTheta=(2.*n+1).*an(n+1,:).*P(:,3)+fTheta;
    
    P(:,1) = P(:,2);
    P(:,2) = P(:,3);
    
end

fTheta=fTheta./(1i.*k0);

end

function w=besselS(n,z,u)

if u==1
    w=sqrt(pi./(2.*z)).*besselj(n+0.5,z);
elseif u==2
    w=sqrt(pi./(2.*z)).*besselh(n+0.5,1,z);
else
    w=sqrt(pi./(2.*z)).*besselh(n+0.5,2,z);
end

end

function w=besselSp(n,x,u) %Derivees

w = besselS(n-1,x,u) - (n+1)./x.*besselS(n,x,u);

end

function w=besselSpp(n,x,u) %Derivees sec. (v. HDR V. Leroy)

if u==1
    w=sqrt(pi./(2.*x)).*(2.*besselj(n+1.5,x)./x + besselj(n+0.5,x).*((n.*(n-1))./x.^2 - 1));
elseif u==2
    w=sqrt(pi./(2.*x)).*(2.*besselh(n+1.5,1,x)./x + besselh(n+0.5,1,x).*((n.*(n-1))./x.^2 - 1));
else
    w=sqrt(pi./(2.*x)).*(2.*besselh(n+1.5,2,x)./x + besselh(n+0.5,2,x).*((n.*(n-1))./x.^2 - 1));
end

end

%% Calcul cos Moy

function cosMoy = calcul_cosMoy (A,k0)


theta = 0 : 0.1 : pi ;
Ftheta = zeros(size(A,2),length(theta));
Nf = size(A,2);

for ii = 1 : length(theta)
    
    Ftheta(:,ii) = foncDiff(A,theta(ii),k0);
    
end

cosMoy = sum(abs(Ftheta).^2.*cos(repmat(theta,[Nf,1])).*sin(repmat(theta,[Nf,1])),2) ...
    ./(sum(abs(Ftheta).^2.*sin(repmat(theta,[Nf,1])),2));

end

