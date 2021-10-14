function out = RT2S_v3(LAI,ze,f,N)
% Radiative transfer model, two-stream approximation
% Inputs:
%      LAI: leaf area index
%      I0:  total incoming radiation
%      ze: solar zenit angle (rad)
%      f: fraction of diffuse radiation
%      N: number of layers

rho = 0.1;   % leaf reflectance
tau = 0.05;  % leaf transmittance
rg = 0.1;    % soil reflectance
omega = rho+tau;
delta = rho-tau;
alfa = 1 - omega; % Absorption coefficient for incoming diffuse radiation


m = 10;    % number of |cos(r.rL)| classes of sun leaves (+1 for shaded leaves)

%% numerical solver

x0 = linspace(0,LAI,N);
[k0,J0] = LeafAngle(x0);

options = odeset('RelTol',1e-6,'AbsTol',1e-6,'Jacobian',@jac1,'Vectorized','off');%,...'InitialSlope',-k0(1));

temp1 = ode15s(@Direct,[0 LAI],1-f,options);


options = bvpset('RelTol',1e-6,'AbsTol',1e-6,'FJacobian',@fjac,...
    'BCJacobian',@fjacbc,'Vectorized','on','NMax', 1e+5, 'Stats', 'On');
% tic

solinit = bvpinit([0 LAI],[f,10]);
temp = bvp5c(@TwoStream,@bcfun,solinit,options);

% CPUTime=toc;

% display(['Solver converged, CPU Time: ' num2str(CPUTime,'%2.2f')])

sol = deval(temp,x0);
Idn(:,1) = sol(1,:);
Iup(:,1) = sol(2,:);
R0(:,1) = deval(temp1,x0);


%% absorbed radiation
% A =  (alfa + gamma)*Iup - gamma*Idn - sigma_bw*R0 + ...
%      (alfa + gamma)*Idn - gamma*Iup - sigma_fw*R0 + ...
%       k*R0;

gamma = 1/2*(omega+J0*delta);            % Backward scattering coefficient for incoming diffuse radiation
sigma_bw = 1/2*omega*k0 + 1/2*J0*delta;   % backward scattering for direct radiation
sigma_fw = omega*k0-sigma_bw;            % farward scattering for direct radiation

out.A0 =  alfa.*(Iup + Idn + k0.'.*R0);

% [out.A,out.p] = Absorbance(x0);

out.R0=R0;
out.Idn=Idn;
out.Iup=Iup;
out.x=x0;
out.alfa=alfa;
%% plottings
if nargout==0
subplot(121)
plot([Idn Iup R0],x)
set(gca,'ydir','reverse')
subplot(122)
% plot(cumsum(A)*x(2),x)
plot(A,x)
set(gca,'ydir','reverse')
end


%% two-stream equations
function dY = Direct(x,Y)

    
k  = spline(x0,k0,x);

dY(1,:) =  -k.*Y;


end

function J = jac1(x,~)
    k  = spline(x0,k0,x);
    J = -k;
end

%% two-stream equations
function dY = TwoStream(x,Y)

    
k  = spline(x0,k0,x);
J  = spline(x0,J0,x);
gamma = 1/2*(omega+J*delta);            % Backward scattering coefficient for incoming diffuse radiation
sigma_bw = 1/2*omega*k + 1/2*J*delta;   % backward scattering for direct radiation
sigma_fw = omega*k-sigma_bw;            % farward scattering for direct radiation

R = deval(temp1,x);

dY(1,:) =  -(alfa + gamma).*Y(1,:) + gamma.*Y(2,:) + sigma_fw.*R;
dY(2,:) =   (alfa + gamma).*Y(2,:) - gamma.*Y(1,:) - sigma_bw.*R;

end

    
%% boundaty conditions
function res = bcfun(ya,yb)

    R = deval(temp1,LAI);
    res = [ya(1)-f;yb(2) - rg*yb(1) - rg*R;];
end   

%% Jacobians
function [dbcdya,dbcdyb] = fjacbc(~,~)
    
dbcdya=[1 0
        0 0];

dbcdyb=[0   0
        -rg 1];

end
    
    
function Jac = fjac(x,~)

J  = spline(x0,J0,x);
gamma = 1/2*(omega+J*delta);            % Backward scattering coefficient for incoming diffuse radiation

Jac(1,1) =  -alfa-gamma;
Jac(1,2) =  +gamma;
Jac(2,1) =  -gamma;
Jac(2,2) =  +alfa+gamma;

end

%% Detto et al 2015 Leaf Angle Distribution
function [k,J] = LeafAngle(x)

n = length(x);
z = 55*(1-x/LAI);
ME=60-(60-22)*exp(-(z/25).^2);
SD=20-(20-13)*exp(-(z/25).^2);

tbar = ME./90;
st=(SD/90).^2;
s0=tbar.*(1-tbar);
nu = tbar.*(s0./st-1);
mu=(1-tbar).*(s0./st-1);
    J = zeros(1,n);
    k = zeros(1,n);
    for i=1:length(z)

    F = @(th) 2./pi./beta(mu(i),nu(i))*(1-2*th/pi).^(mu(i)-1).*(2*th/pi).^(nu(i)-1);
    
   if ze==0
        W1 = @(th) F(th).*cos(th);
        
        k(i) = integral(@(x) W1(x),0,pi/2);
    else
        W1 = @(th) F(th).*cos(th);
        W2 = @(th) F(th).*cos(th).*2./pi.*(sqrt(tan(th).^2.*tan(ze).^2-1)-asec(tan(th).*tan(ze)));

        k(i) = integral(@(x) W1(x),0,pi/2) + integral(@(x) W2(x),pi/2-ze,pi/2);

   end
    
    W1 = @(th) F(th).*cos(th).^2;
    J(i) = integral(@(x) W1(x),0,pi/2);

    end
end

% function [A,p] = Absorbance(x)
%     
%   
% n = length(x);
% z = 55*(1-x/LAI);
% ME=60-(60-22)*exp(-(z/25).^2);
% SD=20-(20-13)*exp(-(z/25).^2);
% 
% tbar = ME./90;
% st=(SD/90).^2;
% s0=tbar.*(1-tbar);
% nu = tbar.*(s0./st-1);
% mu=(1-tbar).*(s0./st-1);
% 
% 
% w=linspace(0,1,m+1);
% fhi = linspace(0,pi,2^10);
% u = linspace(0,1,2^10);
% p = zeros(n,m+1);
% wm = zeros(1,m);
% 
% A = zeros(n,m+1);
% A(:,1) = alfa*(Idn+Iup);
% gap = R0./(1-f);
% p(:,1) = 1-gap;
% for i=1:n
%     th = betaincinv(u,nu(i),mu(i))*pi/2;
%     
%     [y,x]=meshgrid(fhi,th);
%     
%     z = abs(sin(x).*sin(ze).*cos(y)+cos(x).*cos(ze));
%     
%     for j=1:m
%         use=z(:)>=w(j) & abs(z(:))<w(j+1);
%         p(i,j+1) = mean(use)*gap(i);
%         if p(i,j+1)>0
%             wm(i,j) = mean(z(use));
%         else
%             wm(i,j) = w(j)+w(2)/2;
%         end
%     end
%     
% end
% 
%     A(:,2:m+1)=(alfa*(1-f)/cos(ze).*wm + A(:,1));
%     
% end

end