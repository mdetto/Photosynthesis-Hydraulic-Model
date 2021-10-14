function [An,Ac,Aj,gs0,Rd,psiL,psix,psi0] = An_StomatOpt_v3(Vcmax25,k,p50,as,c0,t,dat)
% Net Photosyntesis with stomata optimization coupled with a 3-element hydraulic system
%
%
% Inputs
% Vcmax25:  maximum carboxillation velocity at reference temperature [mumol/m2/s]
% k:   vector [kr kx kL] max hydraulic conductance per unit of leaf area
%      for root, xylem and leaf [mmol m-2 s-1 Mpa-1]
% p50: vector [0 p50x p50L]  water potential at 50% loss of conductivity
%      for xylem and leaf (not defined for root) [Mpa]
% as:  vector [a asx asl] inflection of vulnerability curve for
%      xylem and leaf (not defined root) [Mpa-1]
% c0:  cost parameter associated with psiL [mumol CO2 / Mpa-t / m2]
% t:   nonlinear cost exponent of cost function
%
% dat: meteorological inputs structure
% dat.TL: leaf temperature (C)
% dat.D: VPD  (Kpa)
% dat.I: incident PAR (mumol m-2) (it can be a vector)
% dat.Ca: ambient CO2 concentration (PPM)
% dat.H: height of the leaf (m)
%
% Example:
%
% dat = struct('I',1000,'Ca',400,'TL',25,'D',1.5,'psiS',-0.5,'H',30);
% Vcmax25 = 40;
% Kmax = [5 10 10];
% p50 = [0 -1.5 -1];
% as = [0 6 3];
% c0 = 1;
% t = 2;
%   [An,Ac,Aj,gs0,Rd,psiL,psix,psi0]  = An_StomatOpt_v3(Vcmax25,Kmax,p50,as,c0,t,dat);
%
% Algorithm:
% The optimization is based on maximization of An - cost
% where cost = co x psiL^t
% psiL is computed from a three-element hydrualic model root-xylem-leaves
% soil water potential is the weighted average by root distribution
% For example, using Jackson 1996, psiS can be computed as
% beta = 0.961;f=zeros(4,1);
% f(1) = 1-exp(log(beta)*z1);
% f(2) = exp(log(beta)*z2)-exp(log(beta)*z3);
% f(3) = exp(log(beta)*z3)-exp(log(beta)*z4);
% f(4) = exp(log(beta)*z4);
% psiS = psi*f
% where psi is a vector of water potential at differnt depths z1, z2, z3 and z4 
% Photosyntesis is computed from the Farqhuar model with temperature
% dependent kynetics from Medelyn (2002)

%% parametrization
Tg = 25; %plant growing temperature
R = 8.414;            % J mol-1 K-1  
Kc25 = 404.9;         % mubar
E.Kc = 79430;         % J mol-1   
Ko25 = 278.4;         % mbar
E.Ko = 36380;         % J mol-1 
Gstar25 = 42.75;      % mumol/mol
E.Gstar = 37830;      % mubar
Rd25 = 0.015*Vcmax25; % mumol m-2 s-1
E.Rd = 46390;         % J mol-1 

% Vcmax (from Medelyn et al., 2002)
Ha.Vcmax = 66560;     % J mol-1  
Hd.Vcmax = 200000;    % J mol-1  
S.Vcmax = 638-1.07*Tg;% J mol-1 K-1 Kattge and Knorr (2007), table 3

% Jmax (from Medelyn et al., 2002)
Jmax25 = exp(1.01+0.89*log(Vcmax25)); % mumol/m2/s (Walker at al., 2014, table 4)
Ha.Jmax = 43965;      % J mol-1  
Hd.Jmax = 200000;     % J mol-1  
S.Jmax = 659-0.75*Tg;% J mol-1  Kattge and Knorr (2007), table 3
theta = .7;           % unitless
alpha = 0.36;         % unitless

% environmental conditions
Ca =  dat.Ca;         % mubar
O  =  210;            % mbar
D  =  dat.D*10;       % mbar
TK = dat.TL+273.15;   % Kelvin
psiS = dat.psiS;      % Mpa
H = dat.H;            % m
% temperature functions
Kc  =  Kc25*exp((TK-298)./(R*TK*298)*E.Kc);
Ko  =  Ko25*exp((TK-298)./(R*TK*298)*E.Ko);
Gstar  =  Gstar25*exp((TK-298)./(R*TK*298)*E.Gstar);
Rd  =  Rd25*exp((TK-298)./(R*TK*298)*E.Rd);

Vcmax  =  Vcmax25*exp((TK-298)./(R*TK*298)*Ha.Vcmax).*(1+exp((298*S.Vcmax-Hd.Vcmax)/(298*R)))...
       ./(1+exp((S.Vcmax*TK-Hd.Vcmax)./(R*TK)));
   
Jmax   =   Jmax25*exp((TK-298)./(R*TK*298)*Ha.Jmax).*(1+exp((298*S.Jmax-Hd.Jmax)/(298*R)))...
          ./(1+exp((S.Jmax*TK-Hd.Jmax)./(R*TK)));
      
% Vcmax  =  Vcmax_opt.*Hd.Vcmax.*exp(Ha.Vcmax.*(TK-Topt)/(TK*R*Topt))...
%        ./(Hd.Vcmax+exp((S.Vcmax*TK-Hd.Vcmax)./(R*TK)));
%    
% Jmax   =   Jmax25*exp((TK-298)./(R*TK*298)*Ha.Jmax).*(1+exp((298*S.Jmax-Hd.Jmax)/(298*R)))...
%           ./(1+exp((S.Jmax*TK-Hd.Jmax)./(R*TK)));      

Km = Kc.*(1+O./Ko); 

%% hydraulics

opts = optimoptions(@fmincon,'Algorithm','active-set','Display','off',...
    'SpecifyConstraintGradient',false);
% opts = optimset('TolX',1e-9,'TolFun',1e-9,'display','off');

An = zeros(length(dat.I),1);
gs0 = zeros(length(dat.I),1);
gs = zeros(2,1);
rho = 1e-3; %density of water (ton m-3)
g = 9.81; %gravity (m s-2)
pz = rho*g*H; %gravitational term

p0 = psiS;
p1 = p0-pz;
p2 = p1;

s20 = exp(as(2)*(p0+0-p50(2)));
s21 = exp(as(2)*(p1+pz-p50(2)));
s31 = exp(as(3)*(p1+pz-p50(3)));
s32 = exp(as(3)*(p2+pz-p50(3)));

k10 = k(1);
k20 = k(2).*s20./(s20+1);
k21 = k(2).*s21./(s21+1);
k31 = k(3).*s31./(s31+1);
k32 = k(3).*s32./(s32+1);

thetap = t*c0*abs(p2).^(t-1); 
LD = thetap./k32.*(1+k31./k21.*(1+k20./k10)).*D;

% lambda = -c0*psiS./Kmax./exp(as*(psiS-p50)).*(exp(as*(psiS-p50))+1);
%% carbon limited
a  =  (Vcmax-Rd)/2;
b  =  (Ca+Km)/2;
c  =  Rd./2.*(Ca+Km) + Vcmax./2.*(2*Gstar-Ca+Km);

gs(1) = fmincon(@hydraulics,0.001,-1,0,[],[],0,[],@nonlcon,opts);

% LD = c0*abs(psiS - pz).^t/Kmax*A2./(A2-1)*D;
% gmax  =  sqrt((2*a.*b-c).*(2*a.*b+c))./(2*b.^2).*(b-LD)./sqrt(2*b.*LD-LD.^2) - c./(2*b.^2);
% 

Ac(:,1) = a + b.*gs(1) - sqrt(b.^2.*gs(1).^2+c.*gs(1)+a.^2);

%% light limited
Aj = zeros(length(dat.I),1);
for i=1:length(dat.I)
    
J  =  (alpha.*dat.I(i) + Jmax - sqrt((alpha.*dat.I(i)+Jmax).^2 - 4*alpha*Jmax.*theta.*dat.I(i)))./(2*theta);
a  =  J/8-Rd/2;
b  =  Ca/2+Gstar;
c  =  Rd./2.*(Ca  + 2*Gstar) + J./2.*(Gstar - Ca/4);

A0  =  (2*a.*b-c)./(2*b) - sqrt(a.^2-(c./(2*b)).^2).*sqrt(LD./(2*b-LD));
if  real(A0)<0
    gmax = 0;
else
    gmax  =  sqrt((2*a.*b-c).*(2*a.*b+c))./(2*b.^2).*(b-LD)./sqrt(2*b.*LD-LD.^2) - c./(2*b.^2);
end

if gmax>0
    
  gs(2) = fmincon(@hydraulics,0.001,-1,0,[],[],0,[],@nonlcon,opts);

end

Aj(i) = a + b.*gs(2) - sqrt(b.^2.*gs(2).^2+c.*gs(2)+a.^2);
[An(i),index] = min([Ac Aj(i)]);
gs0(i) = 1.6*gs(index);
end

E = gs0*D;
psi0 = psiS - E/k(1);
psix = p50(2) - pz + 1/as(2)*log((exp(as(2)*(p0-   p50(2))) + 1).*exp(-as(2)*E/k(2)) - 1);
psiL = p50(3) - pz + 1/as(3)*log((exp(as(3)*(p1+pz-p50(3))) + 1).*exp(-as(3)*E/k(3)) - 1);

%% functions for solver
function [y,p0,p1,p2] = hydraulics(x)

E = 1.6*x*D;

p0 = psiS - E/k(1);
p1 = p50(2) - pz + 1/as(2)*log((exp(as(2)*(p0   -p50(2))) + 1).*exp(-as(2)*E/k(2)) - 1);
p2 = p50(3) - pz + 1/as(3)*log((exp(as(3)*(p1+pz-p50(3))) + 1).*exp(-as(3)*E/k(3)) - 1);

y = c0*abs(p2).^t - a - b.*x + sqrt(b.^2.*x.^2+c.*x+a.^2);

end


function [c,ceq] = nonlcon(x)

E = 1.6*x*D;

p0 = psiS - E/k(1);
p1 = p50(2) - pz + 1/as(2)*log((exp(as(2)*(p0-p50(2))) + 1).*exp(-as(2)*E/k(2)) - 1);

c = 1 - (exp(as(3)*(p1+pz-p50(3))) + 1).*exp(-as(3)*E/k(3));
ceq = [];

% W = exp(as(2)*(psiS-p50(2)-E/k(1)));
% dp1_dE = ((W+1)/k(2) + W/k(1))./(1+W-exp(as(2)*E/k(2)));
%      
%     
% dc_dE =   as(3)*exp(-as(3)*E/k(3)).*(exp(as(3)*(pz - p50(3) + p1)) + 1)/k(3) ...
%         - as(3)*exp(-as(3)*E/k(3)).* exp(as(3)*(pz - p50(3) + p1)).*dp1_dE;
% 
% gradc = dc_dE*1.6*D;
% gradceq = [];
end

end
%% co-limitation Vico et al., (2013)
% elseif strcmp(dat.model,'colimit')
%     
% J  =  (alpha.*dat.I + Jmax - sqrt((alpha.*dat.I+Jmax).^2 - 4*alpha*Jmax.*theta.*dat.I))./(2*theta);
% a  =  J/8-Rd/2;
% b  =  Ca+Km.*J./Vcmax/4;
% c  =  Rd./2.*(Ca+Km.*J./Vcmax/4) + J./8.*(2*Gstar-Ca + Km.*J./Vcmax/4);
% 
% 
% An  =  (2*a.*b-c)./(2*b) - sqrt(a.^2-(c./(2*b)).^2).*sqrt(LD./(2*b-LD));
% if  real(An)<0
%     gmax = 0;
% else
%     gmax  =  sqrt((2*a.*b-c).*(2*a.*b+c))./(2*b.^2).*(b-LD)./sqrt(2*b.*LD-LD.^2) - c./(2*b.^2);
% end
% 
% g0 = min(gmax,glim);
% if gmax>0
%     
%     A3 = exp(-as*D/Kmax);
%     A4 = (lambda0*D).^(as*psi0);
%     gs = fzero(@(x) ...
%         ((b - (b^2*x + c/2)/sqrt(b^2*x^2+c*x+a^2))).^(as*psi0) - ...
%         (A1*A3.^(1.6*x)-A2).*A4,...
%         [0 g0],opts);
% end
% 
% An = a + b.*gs - sqrt(b.^2.*gs.^2+c.*gs+a.^2);
% gs0 = 1.6*gs;
% psiL = p50 + 1/as*log((exp(as*(psiS-p50)) + 1)*exp(-as*gs0*D/Kmax) - 1);
% end


% check analitycal solutions
% gc  =  sqrt((a.*b-c).*(a.*b+c))./b.^2.*(b-2*LD)./sqrt(b.*LD-LD.^2) - 2*c./b.^2;
% cc(:,1)  =  Ca - Ac(:,1)./gc;
% Ac(:,2)  =  Vcmax./(Km+cc).*(cc-Gstar)-Rd;

% gj  =  sqrt((a.*b-c).*(a.*b+c))./b.^2.*(b-2*LD)./sqrt(b.*LD-LD.^2) - 2*c./b.^2;
% cj  =  Ca - Aj(:,1)./gj;
% Aj(:,2)  =  J./(4*cj+8*Gstar).*(cj-Gstar)-Rd;



%% Medelyn et al, 2002

% dat  =  readtable('C:\Users\mdetto\Dropbox (Smithsonian)\paper\Rubisco\Medelyn et al, 2002.xlsx');
% R = 8.414*1e-3;
% Ha.Jmax    =  median(dat.Ha);
% Hd.Jmax    =  median(dat.Hd);
% Topt.Jmax  =  median(dat.Topt);
% S.Jmax     =  median(R*log(dat.Ha./(dat.Hd-dat.Ha)) + dat.Hd./(273+dat.Topt));
% Ha.Vcmax   =  median(dat.Ha_1);
% Hd.Vcmax   =  median(dat.Hd_1);
% S.Vcmax    =  median(R*log(dat.Ha_1./(dat.Hd_1-dat.Ha_1)) + dat.Hd_1./(273+dat.Topt_1));

