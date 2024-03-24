%- Script pipers_ngpt.m calculates the Paleothermometer parameters for the 
%- PIPERS data parameters are fsw,fgmw,fcdw,A,S

clear; close;

%- Freezing point temp.  This should be computed using in-situ salinity.
tfrz = -1.91;

%- LOAD PIPERS data
df=readtable('PIPERS_NG+Hydro.csv');
dF = table2array(df);

%- Column indices for He, Ne, Ar, Kr, Xe, T, and S.
idx = [5:9,16,15];

% - He, Ne,Ar,Kr,Xe,T,S
c = dF(:,idx); c(:,1:5) = c(:,1:5)*1e-3; 


%- Set the intial values for the lsqnonlin() iteration. 
x0(1) = 1;  %-f_surface water 
x0(2) = 1;  %-f_gmw
x0(3) = 1;  %-f_cdw
x0(4) = 0;  %-Air content, A.
x0(5) = 40; %- Salinity, psu
x0(6) = 1;  %- f_si.


%- these are precision estimates in %.  (std(c)/avg(c)*100).
cstd = std(c(:,1:end-2))';


for i = 1:size(c,1),

    cw = c(i,:)';
    cw(end+1,1) = 1;

    %- Adjust weight matrix to give same order weights to tracers
    %- and extra weight to continuity constraint.
    w = diag([1./cstd;-60;3;800]);
    w(1,1) = w(1,1)*2; 

    
    %- frw,fgmw,fcdw,Air,Sal,fsi
    xlb = [0 0 0 0 0 0];
    xub = [1 1 1 41428 40 1];

    %- Use the freezing point temp as the temp of dense shelf water.
    t = tfrz;

    opts = optimoptions(@lsqnonlin,'Display','off','Jacobian','off',...
        'tolfun',1e-10,'tolx',1e-10,'MaxFunEvals',1e5,'MaxIter',5e3);
    
    [x(i,:),RESNORM,RESIDUAL,EXITFLAG(i),OUTPUT,LAMBDA] = lsqnonlin(@(x) ngpaleofun(x,cw,w,t),x0,xlb,xub,opts);
end

%%
%- PLOT ALL MISFITS
rr = .892;
Airp = 0.97;
air = 0.11;
S = x(:,5); T = tfrz;

%- Equil. solubility concentration (umol/kg). 
Ceq = Airp*[LJSolMol(1,T,S),LJSolMol(2,T,S),LJSolMol(3,T,S),LJSolMol(4,T,S),LJSolMol(5,T,S)]'.*1e6;

%- Atmospheric partial pressures.
chi = Airp*[5.24e-6,18.18e-6, .00934,1.14e-6,87e-9]';

%- Air conntent in pure GMW.
Cpure = air/22.414*1e6*chi;

%- values from Posthelwaite, 2002 unpublished thesis and Top et al., (1988)
kiw = [1.33,.83,0.49,0.4,0.5];


%- This gives cdw end member
cdw = [0.0019,0.0080,16.2477,0.0039,0.0006,34.6290,1.6210];
cdw = cdw(:);

figure(3); clf;
Tmod = x(:,1).*T + x(:,2).*(-92) + x(:,3).*cdw(end)+ rr.*x(:,end).*3.35e5/4180*.1675 - rr.*T.*x(:,end);
Smod = x(:,1).*x(:,end-1) + x(:,3).*cdw(end-1)-x(:,end-1).*x(:,end).*.33;
mod = x(:,1)+x(:,2)+x(:,3)-rr.*x(:,end);

%- Noble gas reconstruction from NGPT model.
Xemod = (Ceq(5,:)'+chi(5).*x(:,4)).*x(:,1)-rr.*x(:,end).*kiw(5).*Ceq(5,:)' + x(:,2).*Cpure(5) + x(:,3).*cdw(5);
Krmod = (Ceq(4,:)'+chi(4).*x(:,4)).*x(:,1)-rr.*x(:,end).*kiw(4).*Ceq(4,:)' + x(:,2).*Cpure(4) + x(:,3).*cdw(4);
Armod = (Ceq(3,:)'+chi(3).*x(:,4)).*x(:,1)-rr.*x(:,end).*kiw(3).*Ceq(3,:)' + x(:,2).*Cpure(3) + x(:,3).*cdw(3);
Nemod = (Ceq(2,:)'+chi(2).*x(:,4)).*x(:,1)-rr.*x(:,end).*kiw(2).*Ceq(2,:)' + x(:,2).*Cpure(2) + x(:,3).*cdw(2);
Hemod = (Ceq(1,:)'+chi(1).*x(:,4)).*x(:,1)-rr.*x(:,end).*kiw(1).*Ceq(1,:)' + x(:,2).*Cpure(1) + x(:,3).*cdw(1);


scatter((Tmod-dF(:,16))./dF(:,15)*100,dF(:,10),'d','filled'); hold on;
scatter((Smod-dF(:,15))./dF(:,15)*100,dF(:,10),'s','filled'); hold on;
scatter((mod - 1)*100,dF(:,10),'^','filled');
set(gca,'ydir','rev'); set(get(gcf,'children'),'fontsize',14);
ylabel('Depth (m)');
xlabel('Model-data misfit in %');
leg = legend('tmod','smod','Water mass misfit','location','southwest');
grid on

figure(4); clf;
scatter((Hemod*1e3-dF(:,idx(1)))./dF(:,idx(1))*100,dF(:,10),'o','filled'); hold on;
scatter((Nemod*1e3-dF(:,idx(2)))./dF(:,idx(2))*100,dF(:,10),'s','filled'); hold on;
scatter((Armod*1e3-dF(:,idx(3)))./dF(:,idx(3))*100,dF(:,10),'d','filled');
scatter((Krmod*1e3-dF(:,idx(4)))./dF(:,idx(4))*100,dF(:,10),'^','filled');
scatter((Xemod*1e3-dF(:,idx(5)))./dF(:,idx(5))*100,dF(:,10),'p','filled');
set(gca,'ydir','rev'); set(get(gcf,'children'),'fontsize',14);
ylabel('Depth (m)');
xlabel('Model-data misfit in %');
leg = legend('He misfit','Ne misfit','Ar misfit','Kr misfit','Xe misfit','location','southwest');
grid on


%%
tnb = dF(:,2)<=35;
ris = dF(:,2) >35;

figure;  subplot(1,3,1);
scatter(x(tnb,5),dF(tnb,10),'d','filled'); set(gca,'ydir','rev');
hold on; scatter(x(ris,5),dF(ris,10),'d','filled','red'); 
xlabel('salinity');

subplot(1,3,2);
scatter(x(tnb,4),dF(tnb,10),'d','filled'); set(gca,'ydir','rev');
hold on; scatter(x(ris,4),dF(ris,10),'d','filled','red'); 
xlabel('air (mol/kg)');

figure; subplot(1,2,1);
scatter(x(tnb,6)*100,dF(tnb,10),'d','filled'); set(gca,'ydir','rev');
hold on; scatter(x(ris,6)*100,dF(ris,10),'d','filled','red'); 
xlabel('sea ice (%)');

subplot(1,2,2);
scatter(x(tnb,2)*100,dF(tnb,10),'d','filled'); set(gca,'ydir','rev');
hold on; scatter(x(ris,2)*100,dF(ris,10),'d','filled','red'); 
xlabel('GMW (%)');


figure; subplot(1,2,1);
scatter(x(tnb,1),dF(tnb,10),'d','filled'); set(gca,'ydir','rev');
hold on; scatter(x(ris,1),dF(ris,10),'d','filled','red'); 
xlabel('fsw');

subplot(1,2,2);
scatter(x(tnb,3),dF(tnb,10),'d','filled'); set(gca,'ydir','rev');
hold on; scatter(x(ris,3),dF(ris,10),'d','filled','red'); 
xlabel('cdw');





