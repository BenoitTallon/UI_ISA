function [x,y,P0it,P0st,P1t,r,theta] = fct_Simu_deplacements(a,f,N,c0,c1,rho0,rho1,alpha0,alpha1,lim,pas,t,M)

%% Calcul des amplitudes An, Bn

om0=2*pi*f;                                                                % Pulsation dans la matrice
k0=om0/c0 + 1i.*alpha0;                                                    % Nombre d'onde dans la matrice
k1=om0/c1 + 1i.*alpha1;                                                    % Nombre d'onde dans l'inclusion

A=zeros(N,length(f));
B=A;
q=rho1*k0/(rho0*k1);

for ii=1:length(f)
    
    A(:,ii)=-(q.*besselSp(0:N-1,k0(ii)*a,1).*besselS(0:N-1,k1(ii).*a,1) - ...
        besselS(0:N-1,k0(ii).*a,1).*besselSp(0:N-1,k1(ii)*a,1))./ ...
        (q.*besselSp(0:N-1,k0(ii)*a,2).*besselS(0:N-1,k1(ii).*a,1)- ...
        besselS(0:N-1,k0(ii).*a,2).*besselSp(0:N-1,k1(ii)*a,1));
    
    B(:,ii)=A(:,ii).'.*besselS(0:N-1,k0(ii)*a,2)./besselS(0:N-1,k1(ii)*a,1) + ...
        besselS(0:N-1,k0(ii)*a,1)./besselS(0:N-1,k1(ii)*a,1);
    
end

%% Calcul du champ de deplacements

xx=-lim/2:pas:lim/2;                                                       % Discretisation en x (mm)
yy=-lim/2:pas:lim/2;                                                       % Discretisation en y (mm)
[x y]=meshgrid(xx,yy);

Nx=length(x);
Ny=length(y);

r=sqrt(x.^2+y.^2);                                                         % Passage en coord. polaires
r1=r;
theta1=atan(y./x).*(x>=0).*(y>=0);
theta2=(atan(y./(-x))+pi/2)'.*(x<0).*(y>=0);
theta3=(atan(y./(x))+pi).*(x<0).*(y<0);
theta4=(atan(-y./(x))+3*pi/2)'.*(x>=0).*(y<0);

theta=theta1+theta2+theta3+theta4;

theta1=theta;

R=zeros(Nx,Ny,length(t));
T=R;

[Uri0,Urs0,Uts0,Uti0,P0i,P0s,Ur1,Ut1,P1] = CalculDep(Nx,Ny,k0,k1,r,theta,a,A,B,rho0,rho1,om0,M);  % sub-routine calcul des champs de deplacement

Ur = Uri0 + Urs0 + Ur1;
Utheta = Uts0 + Uti0 + Ut1;

for it=1:length(t)
    
    r=r1+real(Ur.*exp(-1i.*om0.*t(it)));                                   % Creation de la grille deformee (en coord. polaires)
    r(isnan(r)) = 0;
    theta=theta1+real(Utheta.*exp(-1i.*om0.*t(it)));
    
    R(:,:,it)=r;
    T(:,:,it)=theta;
    
    clc
    
end

%%

%% Champ de déplacements a la frontière

pas=2*pi/500;                                                              % Pas discretisation spatiale du contour de l'inclusion
theta=0:pas:2*pi;

Nt=length(theta);

P0it = zeros(Nx,Ny,length(t));  P0st = P0it;    P1t = P0it;

theta1=theta;

R(isnan(R)) = 0;
[xxx,yyy] = (find(r == min(min(r))));
T(xxx,yyy+1:end,:) = 2.*pi.*ones(size(T(xxx,yyy+1:end,:)));
T(xxx,1 : yyy-1,:) = pi.*ones(size(T(xxx,1 : yyy-1,:)));
T(xxx,yyy,:) = zeros(size(T(xxx,yyy,:)));

for it = 1 : length(t)
    
    P0it(:,:,it) = real(P0i.*exp(-1i.*om0.*t(it)));
    P0st(:,:,it) = real(P0s.*exp(-1i.*om0.*t(it)));
    P1t(:,:,it) = real(P1.*exp(-1i.*om0.*t(it)));
    
end

[Ura,Uta] = CalculDepFront(Nt,k1,theta,a,B,rho1,om0);

Ura = M.*Ura;
Uta = M.*Uta;

r = zeros(length(theta1),length(t));
theta = r;

for it=1:length(t)
    
    r(:,it)=a+real(Ura.*exp(-1i.*om0.*t(it)));
    
    theta(:,it)=theta1+real(Uta.*exp(-1i.*om0.*t(it)));
    
    x=R(:,:,it).*cos(T(:,:,it));
    y=R(:,:,it).*sin(T(:,:,it));
    
end

end

%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Sub-routines %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Calculation of pressure and displacement fields

function [Uri0,Urs0,Uts0,Uti0,P0i,P0s,Ur1,Ut1,P1] = CalculDep(Nx,Ny,k0,k1,r,theta,a,A,B,rho0,rho1,om0,M)

Ursi0 = zeros(Nx,Ny);   Urss0 = Ursi0;  Utsi0 = Ursi0;  Utss0 = Ursi0;
Psi0 = Ursi0;           Pss0 = Ursi0;   Urs1 = Ursi0;   Uts1 = Ursi0;
Ps1 = Ursi0;

N=length(A);

for n=0:N-1
    
    P=legN(n,cos(theta),Nx,Ny);
    
    if n == 0
        
        Pp = zeros(size(P));
        
    else
        
        P11=legN(n - 1,cos(theta),Nx,Ny);
        Pp=(n./(cos(theta).^2 - 1)).*(cos(theta).*P - P11);                     % Calcul polynômes de Legendre et derivees
        
    end
    
    %% Dans la matrice
    
    Uri0 = 1i.^n.*(2*n+1).*P.*k0.*besselSp(n,k0.*r,1)./ ...
        (rho0.*om0^2).*(r>a)+Ursi0;                                        % Deplacement Onde incidente (vaut 0 dans l'inclusion)
    
    Urs0 = 1i.^n.*(2*n+1).*P.*k0.*A(n+1).*besselSp(n,k0.*r,2)./ ...        % Deplacement Onde diffusee (vaut 0 dans l'inclusion)
        (rho0.*om0^2).*(r>a)+Urss0;
    
    Uti0=-sin(theta).*1i.^n.*(2*n+1).*Pp.*besselS(n,k0.*r,1)./...
        (r.*rho0.*om0^2).*(r>a)+Utsi0;                                         % idem
    
    Uts0=-sin(theta).*1i.^n.*(2*n+1).*Pp.*A(n+1).*besselS(n,k0.*r,2)./...
        (r.*rho0.*om0^2).*(r>a)+Utss0;                                     % idem
    
    P0i=1i.^n.*(2*n+1).*P.*besselS(n,k0.*r,1).*(r>a)+Psi0;                 % Pression incidente dans la matrice
    
    P0s=1i.^n.*(2*n+1).*P.*A(n+1).* besselS(n,k0.*r,2).*(r>a)+Pss0;        % Pression incidente dans la matrice
    
    %% Dans l'inclusion
    
    Ur1=(1i.^n.*(2*n+1).*P.*k1.*B(n+1).*besselSp(n,k1.*r,1))./...
        (rho1.*om0^2).*(r<=a).*M+Urs1;                                     % vaut 0 dans la matrice
    
    Ut1=-sin(theta).*(1i.^n.*(2*n+1).*Pp.*B(n+1).*besselS(n,k1.*r,1))./...
        (r.*rho1.*om0^2).*(r<=a).*M+Uts1;                                  %idem
    
    P1=1i.^n.*(2*n+1).*P.*B(n+1).*besselS(n,k1.*r,1).*(r<=a)+Ps1;          % Pression dans l'inclusion
    
    %% Somme
    
    Ursi0 = Uri0;   Utsi0 = Uti0;   Utss0 = Uts0;   Urss0 = Urs0;
    Psi0 = P0i;     Pss0 = P0s;     Urs1 = Ur1;     Uts1 = Ut1;
    Ps1 = P1;
    
end

end

%% Displacements of the droplet surface

function [Ura,Uta] = CalculDepFront(Nt,k0,theta,a,A,rho0,om0)

% Calcul des déplacements à la frontière à partir des pressions int et ext

Tra=zeros(1,Nt); Tta=zeros(1,Nt);

N=length(A);


for n=0:N-1
    
    P=legN(n,cos(theta),1,Nt);
    
    if n == 0
        
        Pp = zeros(size(P));
        
    else
        
        P11=legN(n - 1,cos(theta),1,Nt);
        Pp=(n./(cos(theta).^2 - 1)).*(cos(theta).*P - P11);                     % Calcul polynômes de Legendre et derivees
        
    end
    
    %% Interface
    
    Ura=(1i.^n.*(2*n+1).*P.*k0.*A(n+1).*besselSp(n,k0.*a,1))./...
        (rho0.*om0^2)+Tra;                                         % vaut 0 dans la matrice
    
    Uta=-sin(theta).*(1i.^n.*(2*n+1).*Pp.*A(n+1).*besselS(n,k0.*a,1))./...
        (a.*rho0.*om0^2)+Tta;
    
    %% Somme
    
    Tra=Ura;
    Tta=Uta;
    
end

end

%% Spherical Bessel functions

function w=besselS(n,z,u)

if u==1
    w=sqrt(pi./(2.*z)).*besselj(n+0.5,z);
else
    w=sqrt(pi./(2.*z)).*besselh(n+0.5,z);
end

end

%% Derivatives of spherical Bessel functions

function w=besselSp(n,x,u)

if u==1
    w=sqrt(pi./(2.*x)).*((n+0.5).*besselj(n+0.5,x)./...
        x-besselj(n+1.5,x)-besselj(n+0.5,x)./(2.*x));
else
    w=sqrt(pi./(2.*x)).*((n+0.5).*besselh(n+0.5,x)./...
        x-besselh(n+1.5,x)-besselh(n+0.5,x)./(2.*x));
end

end

%% Legendre polynomials

function P=legN(n,x,Nx,Ny)

P=legendre(n,x);

if n>0
    
    P=reshape(P(1,:,:),Nx,Ny);                                             % Suppression des tableaux m>0 (sauf pour n=0)
    
end

end
