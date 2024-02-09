function [f,k0a,vl0,kISA,vph,vgr,ve,ltr] = k_ISA_2D_test(M0,M1,M2,N,f,Rm,P,phi,Nr,Rwidth)

% Mat = materiau pour la matrice
% Inc = materiau pour les inclusions
% N = nombre de modes
% fMAX = freq maximum pour le calcul
% amoy = rayon moyen sphere (mm)
% P = polydispersite
% phi = fraction volumique rho2,vl2,vt2,alphal2,alphat2

f(f==0) = eps;

%% Entrees

global rho0 rho1 rho2 vl0 vl1 vt1 vl2 vt2 alphal0 alphal1 alphat1 alphal2 alphat2

phi = phi/100;            % fraction volumique en %
P = P/100;                % polydispersite en %

rho0 = M0(1); vl0 = M0(2); alphal0 = M0(3);
rho1 = M1(1); vl1 = M1(2); vt1 = M1(3); alphal1 = M1(4); alphat1 = M1(5);
rho2 = M2(1); vl2 = M2(2); vt2 = M2(3); alphal2 = M2(4); alphat2 = M1(5);

[kISA,ltr,vph,vgr,ve] = fct_Calc_k_ISA(Rm,f,N,phi,P,Nr,Rwidth);

%%

k0a = real(2.*pi.*f./vl0)*Rm(1);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Sub-routines

function [kISA,ltr,vph,vgr,ve] = fct_Calc_k_ISA(Rm,f,N,phi,P,Nr,Rwidth)
%%

if P == 0

    sig=0;
    r=1;
    eta=phi/(pi*(Rm(1)*r)^2);

else

    sig=P;
    r=linspace(1 - Rwidth*sig,1 + Rwidth*sig,Nr);

    %Distribution gaussienne

    gauss = exp(-0.5*((r-1)/sig).^2);
    gaussn =( 1/sum(gauss)) * gauss;

    eta = gaussn.* phi / sum( gaussn .* ( pi*((Rm(1)*r).^2) ) );

end

global rho0 rho1 rho2 vl0 vl1 vt1 vl2 vt2 alphal0 alphal1 alphat1 alphal2 alphat2

omega=2*pi*f;
k0l = omega./vl0 + alphal0.*(omega./(2*pi)).^2.*1i;


Nr = length(r);
Nf=length(f);
A = zeros(N,Nf,Nr);
F0 = zeros(Nr,Nf);

ee = ones(N,1).*2;
ee(1) = 1;
ee = repmat(ee,[1,Nf]);
nn = repmat((1:N).',[1,Nf]);

for ir=1:Nr

    [S] = Calc_Scatt_Amp_Layers(omega,Rm(1)*r(ir),Rm(2)*r(ir),N);

    A(:,:,ir) = S;

    F0(ir,:) = -1i.*sum(ee.*squeeze(A(:,:,ir)));                                    % entrees : (An,theta,kL0)

end

%% K effectif (ISA)

kISA=sqrt(k0l.^2+4.*pi.*eta*F0);

kpoly = sqrt(repmat(k0l,[Nr,1]).^2+4.*pi.*repmat(eta.',[1,Nf]).*F0);       % <!> premiere version kpoly etait definie a partir de kISA

vgr = gradient(omega,real(kISA)) ;
vgrpoly = diff(repmat(omega.',[1,Nr])).'./diff(real(kpoly).').';

%% fonction de diffusion (D. Royer & T. Valier-brasier Vol. 2)

dtheta = 1e-1;
theta = 0 : dtheta : 2*pi;
Ntheta = length(theta);
Fdiff = zeros(length(f),Nr,Ntheta);

ee = repmat(ee,[1,1,Ntheta]);
nn = repmat(nn,[1,1,Ntheta]);
theta = repmat(theta.',[1,N,Nf]);
theta = permute(theta,[2,3,1]);

fTheta = zeros(Nf,Nr,Ntheta);

for ir = 1 : Nr

    fTheta(:,ir,:) = -1i.*squeeze(sum(ee.*repmat(squeeze(A(:,:,ir)),[1,1,Ntheta]).*cos(nn.*theta),1));

end

clear ee nn

%% Calcul du <cos>

theta = theta(1,:,:);
theta = repmat(theta,[Nr,1,1]);
theta = permute(theta,[2,1,3]);

cosMoy = squeeze(sum(abs(fTheta).^2.*cos(theta),3) ...
    ./(sum(abs(fTheta).^2,3)));

ltr = 1./real(4*pi./k0l.*(eta*(imag(F0).*(1-cosMoy.'))));                  % l* Avec le theoreme optique

%% Vitesse de l'énergie

vph = omega.'./real(kISA.');

eta = repmat(eta,[Nf-1,1]);
F0 = F0.';
vgrpoly = vgrpoly.';
omega = repmat(omega,[Nr,1]).';
kpoly = kpoly.';

Ind = vl0.^2./vph(1 : Nf - 1);

d1 = (diff(real(F0))./diff(omega))./real(kpoly(1 : Nf  - 1,:)).*vgrpoly;

d2 = (Fdiff(1:length(f)-1,:,:).*conj(Fdiff(1:length(f)-1,:,:)))...
    .*(diff(unwrap(angle(-Fdiff)))./diff(repmat(omega,[1,1,size(theta,3)])))...
    .*repmat(vgrpoly,[1,1,size(theta,3)]);

d2 = squeeze(sum(d2.*dtheta,3));

dve = sum(2.*pi.*eta.*(d1 + d2),2);

ve = Ind./(1 + dve);
ve = cat(1,real(ve),real(ve(end)));

end


function Sn = Calc_Scatt_Amp_Layers(om,R1,R2,N)

% resolution du systeme [M] . [A] = [C]

%%

global rho0 rho1 rho2 vl0 vl1 vt1 vl2 vt2 alphal0 alphal1 alphat1 alphal2 alphat2

issolid1 = max(vt1) ~= 0;
issolid2 = max(vt2) ~= 0;

k0l = om./vl0 + alphal0.*(om./(2*pi)).^2.*1i;
k1l = om./vl1 + alphal1.*(om./(2*pi)).^2.*1i;
k2l = om./vl2 + alphal2.*(om./(2*pi)).^2.*1i;

k1t = zeros(size(k1l));
k2t = zeros(size(k2l));

mu1 = zeros(size(k1l));
mu2 = zeros(size(k2l));

if issolid1

    k1t = om./vt1 + alphat1.*(om./(2*pi)).^2.*1i;
    mu1 = om.^2.*rho1./k1t.^2;                                             % shear modulus

end

if issolid2

    k2t = om./vt2 + alphat2.*(om./(2*pi)).^2.*1i;
    mu2 = om.^2.*rho2./k2t.^2;                                             % shear modulus

end


lambda0 = om.^2.*rho0./k1l.^2;                                             % lambda modulus
lambda1 = om.^2.*rho1./k1l.^2 - 2.*mu1;                                    % lambda modulus
lambda2 = om.^2.*rho2./k2l.^2 - 2.*mu2;                                    % lambda modulus

%% Amplitude onde diffusee core-shell dans liquide.

Nf = length(om);

%% %%%%%%%%%%%%%%%%%%%%%Calcul scattering coefficients %%%%%%%%%%%%%%%%%%%%

Sn = zeros(N,Nf);                                                          % Coeff de scattering de l'onde diffusee

if R2 == 0

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Plain %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    for n = 0:N-1

        m11 = Unl(k1l,n,R1,1); m13 = Unt(k1t,n,R1,1); m17 = -Unl(k0l,n,R1,2);
        m21 = Snl(lambda1,mu1,k1l,n,R1,1); m23 = Snt(mu1,k1t,n,R1,1); m27 = -Snl(lambda0,0,k0l,n,R1,2);
        m31 = Tnl(mu1,k1l,n,R1,1); m33 = Tnt(mu1,k1t,n,R1,1);

        S1 =  Unl(k0l,n,R1,1); S2 =  Snl(lambda0,0,k0l,n,R1,1);

        %% Methode de Cramer

        for jj = 1 : length(om)

            if issolid1

                S = [S1(jj),S2(jj),0];

                MM = [m11(jj),m13(jj),m17(jj); ...
                    m21(jj),m23(jj),m27(jj); ...
                    m31(jj),m33(jj),0];

            else

                S = [S1(jj),S2(jj)];

                MM = [m11(jj),m17(jj); ...
                    m21(jj),m27(jj)];


            end

            M = MM;
            M(:,end) = S;

            Sn(n+1,jj) = det(M)./det(MM);                                  % Cramer

        end

    end

elseif R2 ~= 0

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Core-shell %%%%%%%%%%%%%%%%%%%%%%%%%%%%

    for n = 0:N-1

        m11 = Unl(k1l,n,R1,1); m12 = Unl(k1l,n,R1,2); m13 = Unt(k1t,n,R1,1); m14 = Unt(k1t,n,R1,2); m17 = -Unl(k0l,n,R1,2);
        m21 = Snl(lambda1,mu1,k1l,n,R1,1); m22 = Snl(lambda1,mu1,k1l,n,R1,2); m23 = Snt(mu1,k1t,n,R1,1); m24 = Snt(mu1,k1t,n,R1,2); m27 = -Snl(lambda0,0,k0l,n,R1,2);
        m31 = Tnl(mu1,k1l,n,R1,1); m32 = Tnl(mu1,k1l,n,R1,2); m33 = Tnt(mu1,k1t,n,R1,1); m34 = Tnt(mu1,k1t,n,R1,2);
        m41 = Unl(k1l,n,R2,1); m42 = Unl(k1l,n,R2,2); m43 = Unt(k1t,n,R2,1); m44 = Unt(k1t,n,R2,2); m45 = -Unl(k2l,n,R2,1); m46 = -Unt(k2t,n,R2,1);
        m51 = Vnl(k1l,n,R2,1); m52 = Vnl(k1l,n,R2,2); m53 = Vnt(k1t,n,R2,1); m54 = Vnt(k1t,n,R2,2); m55 = -Vnl(k2l,n,R2,1); m56 = -Vnt(k2t,n,R2,1);
        m61 = Snl(lambda1,mu1,k1l,n,R2,1); m62 = Snl(lambda1,mu1,k1l,n,R2,2); m63 = Snt(mu1,k1t,n,R2,1); m64 = Snt(mu1,k1t,n,R2,2); m65 = -Snl(lambda2,mu2,k2l,n,R2,1); m66 = -Snt(mu2,k2t,n,R2,1);
        m71 = Tnl(mu1,k1l,n,R2,1); m72 = Tnl(mu1,k1l,n,R2,2); m73 = Tnt(mu1,k1t,n,R2,1); m74 = Tnt(mu1,k1t,n,R2,2); m75 = -Tnl(mu2,k2l,n,R2,1); m76 = -Tnt(mu2,k2t,n,R2,1);

        S1 =  Unl(k0l,n,R1,1); S2 =  Snl(lambda0,0,k0l,n,R1,1);

        %% Methode de Cramer

        for jj = 1 : length(om)

            if issolid1 && issolid2

                S = [S1(jj),S2(jj),0,0,0,0,0];

                MM = [m11(jj),m12(jj),m13(jj),m14(jj),0,0,m17(jj); ...
                    m21(jj),m22(jj),m23(jj),m24(jj),0,0,m27(jj); ...
                    m31(jj),m32(jj),m33(jj),m34(jj),0,0,0; ...
                    m41(jj),m42(jj),m43(jj),m44(jj),m45(jj),m46(jj),0; ...
                    m51(jj),m52(jj),m53(jj),m54(jj),m55(jj),m56(jj),0; ...
                    m61(jj),m62(jj),m63(jj),m64(jj),m65(jj),m66(jj),0; ...
                    m71(jj),m72(jj),m73(jj),m74(jj),m75(jj),m76(jj),0];

            elseif issolid1 && (issolid2 == 0)

                S = [S1(jj),S2(jj),0,0,0,0];

                MM = [m11(jj),m12(jj),m13(jj),m14(jj),0,m17(jj); ...
                    m21(jj),m22(jj),m23(jj),m24(jj),0,m27(jj); ...
                    m31(jj),m32(jj),m33(jj),m34(jj),0,0; ...
                    m41(jj),m42(jj),m43(jj),m44(jj),m45(jj),0; ...
                    m61(jj),m62(jj),m63(jj),m64(jj),m65(jj),0; ...
                    m71(jj),m72(jj),m73(jj),m74(jj),0,0];

            else

                S = [S1(jj),S2(jj),0,0];

                MM = [m11(jj),m12(jj),0,m17(jj); ...
                    m21(jj),m22(jj),0,m27(jj); ...
                    m41(jj),m42(jj),m45(jj),0; ...
                    m61(jj),m62(jj),m65(jj),0];

            end

            M = MM;
            M(:,end) = S;

            Sn(n+1,jj) = det(M)./det(MM);                                  % Cramer

        end

    end

end

end

%% Fonctions locales

function M = Unl(k,n,R,kind)

if kind == 1

    M = k.*bjp(k.*R,n,1);

elseif kind == 2

    M = k.*bjp(k.*R,n,2);

end

end


function M = Unt(k,n,R,kind)

if kind == 1

    M = n.*besselj(n,k.*R)./R;

elseif kind == 2

    M = n.*besselh(n,k.*R)./R;

end

end


function M = Vnl(k,n,R,kind)

if kind == 1

    M = -n.*besselj(n,k.*R)./R;

elseif kind == 2

    M = -n.*besselh(n,k.*R)./R;

end

end

function M = Vnt(k,n,R,kind)

if kind == 1

    M = -k.*bjp(k.*R,n,1);

elseif kind == 2

    M = -k.*bjp(k.*R,n,2);

end

end

function M = Snl(lam,mu,k,n,R,kind)

if kind == 1

    M = (lam + 2.*mu).*k.^2.*bjpp(k.*R,n,1) + lam./R.*(k.*bjp(k.*R,n,1) - ...
        n.^2./R.*besselj(n,k.*R));

elseif kind == 2

    M = (lam + 2.*mu).*k.^2.*bjpp(k.*R,n,2) + lam./R.*(k.*bjp(k.*R,n,2) - ...
        n.^2./R.*besselh(n,k.*R));

end

end

function M = Snt(mu,k,n,R,kind)

if kind == 1

    M = 2.*mu.*n./R.*(k.*bjp(k.*R,n,1) - besselj(n,k.*R)./R);

elseif kind == 2

    M = 2.*mu.*n./R.*(k.*bjp(k.*R,n,2) - besselh(n,k.*R)./R);

end

end

function M = Tnl(mu,k,n,R,kind)

if kind == 1

    M = 2.*mu.*n./R.*(k.*bjp(k.*R,n,1) - besselj(n,k.*R)./R);

elseif kind == 2

    M = 2.*mu.*n./R.*(k.*bjp(k.*R,n,2) - besselh(n,k.*R)./R);

end

end

function M = Tnt(mu,k,n,R,kind)

if kind == 1

    M = mu./R.*(k.*bjp(k.*R,n,1) - k.^2.*R.*bjpp(k.*R,n,1) - n.^2./R.*besselj(n,k.*R));

elseif kind == 2

    M = mu./R.*(k.*bjp(k.*R,n,2) - k.^2.*R.*bjpp(k.*R,n,2) - n.^2./R.*besselh(n,k.*R));

end

end

function M = bjp(x,n,kind)

if kind == 1

    M = 0.5.*(besselj(n-1,x) - besselj(n+1,x));

elseif kind == 2

    M = 0.5.*(besselh(n-1,x) - besselh(n+1,x));

end

end

function M = bjpp(x,n,kind)

if kind == 1

    M = 0.25.*(besselj(n-2,x) - 2.*besselj(n,x) + besselj(n+2,x));

elseif kind == 2

    M = 0.25.*(besselj(n-2,x) - 2.*besselj(n,x) + besselj(n+2,x));

end

end