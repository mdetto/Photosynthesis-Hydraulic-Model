%% simulations
clear
fname = 'Sim_v9.0.mat';

dat=load('inputs.mat');

% 
% fields = fieldnames(dat);
% 
% for i=1:length(fields)
%     if ~strcmp(fields{i},'x')
%         x = getfield(dat,fields{i});
%         [m,n,k]=size(x);
%         if n==1 && k==1
%             eval(['dat.' fields{i} ' = x(1:6:end);']);
%         elseif n>1 && k==1
%             eval(['dat.' fields{i} ' = x(:,1:6:end);']);
%         elseif n>1 && k>1
%             eval(['dat.' fields{i} ' = x(:,:,1:6:end);']);
%         end
%     end
% end

N = 10; %number of layers
n = length(dat.gpp0);
Kmax = [5 10 10]*2;
p50 = [0 -1.5 -1];
gamma = [0 6 3];
Ca = 400;
c0 = 1;
t = 2;
Vcmax25Top = 30;

gpp = zeros(n,1);
T = zeros(n,1);
gc = zeros(n,1);
psiL = zeros(n,1);
g1 = zeros(n,1);
P1 = zeros(n,1);
P2 = zeros(n,1);
parfor i=1:n
        n-i
    if dat.wet(i)
        use = dat.x<=6;
        Vcmax25 = Vcmax25Top*3.^(-dat.x(use)./6);
    else
        use = dat.x<=5;
        Vcmax25 = Vcmax25Top*3.^(-dat.x(use)./5);
    end
    
        

    An = zeros(N,11);
    gs = zeros(N,11);
    Rd = zeros(N,11);
    xxx = 0;

    for j=1:N
      Z = struct('I',dat.a(j,:,i)*dat.I0(i),'Ca',Ca,'TL',dat.tc(i),'D',dat.D(i),'psiS',dat.psiS(i),'H',(1-dat.x(j)/6)*30);       
       
        try
            [An(j,:),~,~,gs(j,:),Rd(j,:),xxx] = An_StomatOpt_v3(Vcmax25(j),Kmax,p50,gamma,c0,t,Z);

        catch ME
            if (strcmp(ME.identifier,'MATLAB:fzero:ValuesAtEndPtsSameSign'))
                An(j,:) = nan;
                gs(j,:) = nan;
                Rd(j,:) = nan;
                xxx = nan;
                
            end
        end
        if j==1
            psiL(i)=mean(xxx);
            g1(i)=mean(gs(1,:));
        end
    end
        

        
        gpp(i) = trapz(dat.x(use),sum((An(use,:)+Rd(use,:)).*dat.p(use,:,i),2));
        T(i) = trapz(dat.x(use),sum(gs(use,:).*dat.p(use,:,i),2))*dat.D(i)*10;
        gc(i) = T(i)./dat.D(i)/10;
        P1(i) = gpp(i)./sqrt(dat.D(i))./Ca;  
        P2(i) = gpp(i)./sqrt(dat.Da(i))./Ca;  
%     end
end

    disp(['error: ' num2str(sum(isnan(gpp)))])
    
save(fname,...
    'gpp','T','gc','P1','P2','Vcmax25Top','Kmax','p50','gamma','c0','Ca','dat','psiL')


%% plottings
addpath('C:\Users\mdetto\Dropbox (Smithsonian)\EC_processing');

figure(1);clf
for j=1:2
    subplot(2,2,j)
qt = linspace(0,1,9);
use0 = dat.wet==2-j;
sm = quantile(dat.psiS(use0),qt)';
col=flipud(parula(length(sm)-1))';
for i=length(sm)-1:-1:1
use2 = dat.psiS>=sm(i) & dat.psiS<sm(i+1) & dat.wet==2-j & ~isnan(gpp);

x = P1(use2);
y = gc(use2);

[xi,yi] = KernelRegression(x,y);
plot(xi,yi,'o-','markersize',4,'color',col(:,i));hold all

end

xlabel('{\boldmath$\frac{GPP}{C_a \sqrt{D}}$}','Interpreter','latex','fontsize',18)
ylabel('g_c  (mol m^{-2} s^{-1})')

if j==1
txt{1}  = 'wet';
title('wet season')
for i=2:length(sm)-2
txt{i} = ' ';
end
txt{length(sm)-1} = 'dry';
lgd = legend(txt,'location','northwest');
title(lgd,{'soil moisture','quantile'})
legend('boxoff')
annotation(gcf,'arrow',[0.23 0.23],...
    [0.87 0.73]);
elseif j==2
    title('dry season')
end
set(gca,'xtick',0.02:0.02:0.12)
% axis([0.02 0.12 0.1 1.1])
axis square
% axis([0.02 0.12 0.05 1.5])
end

% data from Wright and Cornejo (1990)
psimin = [-1.18971575 -1.793760878]; %wet&dry
SE = [0.063137398 0.11982874];


figure(2);clf
plot(dat.psiS(dat.wet),psiL(dat.wet),'.');hold all
plot(dat.psiS(~dat.wet),psiL(~dat.wet),'.')
clear h
for i=1:2
h(1) = refline(0,psimin(i)-SE(i));
h(2) = refline(0,psimin(i)+SE(i));
if i==1;set(h,'color','blue');else;set(h,'color','red');end
end
axis square
xlabel('\psi_{soil} (Mpa)')
ylabel('\psi_{leaf} (Mpa)')
clear h
h=get(gca,'children');
legend([h(6) h(5) h(3) h(1)],'wet','dry','Wright and Cornejo, 1990','')
legend('boxoff')
pause(.1)

%% LUE and WUE

figure(1)
qt = .75;
col =[0     0.8500
    0.4470  0.3250
    0.7410  0.0980];
subplot(223)
densityPlot(dat.I0(dat.wet),gpp(dat.wet),qt,col(:,1),0)
[xi,yi] = KernelRegression(dat.I0(dat.wet),gpp(dat.wet));
plot(xi,yi,'o-','markersize',4,'color',col(:,1))

densityPlot(dat.I0(~dat.wet),gpp(~dat.wet),qt,col(:,2),0)
[xi,yi] = KernelRegression(dat.I0(~dat.wet),gpp(~dat.wet));
plot(xi,yi,'o-','markersize',4,'color',col(:,2))

axis([0 2200 0 33])
xlabel('PAR (\mumol m^{-2} s^{-1})')
ylabel('GPP (\mumol m^{-2} s^{-1})')


h=get(gca,'children');
legend(h([1 4]),'wet','dry')
legend('boxoff')
title('Light use efficiency')
axis square

subplot(224)
densityPlot(T(dat.wet),gpp(dat.wet),qt,col(:,1),0)
[xi,yi] = KernelRegression(T(dat.wet),gpp(dat.wet));
plot(xi,yi,'o-','markersize',4,'color',col(:,1))

densityPlot(T(~dat.wet),gpp(~dat.wet),qt,col(:,2),0)
[xi,yi] = KernelRegression(T(~dat.wet),gpp(~dat.wet));
plot(xi,yi,'o-','markersize',4,'color',col(:,2))

xlabel('ET (\mumol m^{-2} s^{-1})')
ylabel('GPP (\mumol m^{-2} s^{-1})')

axis([-.3 11 0 33])
title('Water use efficiency')
axis square
sublabel

% quit