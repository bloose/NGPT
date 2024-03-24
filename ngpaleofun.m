function [G] = ngpaleofun(x,cw,w,t),

%- function executed for computing the non-linear least squares solution
%- to the air injection + glacial meltwater.  Equation is:

air = 0.11;


fsw = x(1);
fgmw = x(2);
fcdw = x(3);
A = x(4);
S = x(5);
fsi = x(6);
T = t; %- Set temperature parameter to freeze point temp, based on S,p.

Airp = .97;

%- Noble gas partial pressures in air[He, Ne, Ar, Kr, Xe];
chi = Airp*[5.24e-6,18.18e-6, .00934,1.14e-6,87e-9]';

%- This gives umol of gas/kg of GMW
Cpure = air/22.414*1e6*chi;

%- This gives cdw end member
cdw = [0.0019,0.0080,16.2477,0.0039,0.0006,34.6290,1.6210];
cdw = cdw(:);


%- Use Lott-Jenkins solubility functions
Ceq = Airp.*[LJSolMol(1,T,S),LJSolMol(2,T,S),LJSolMol(3,T,S),LJSolMol(4,T,S),LJSolMol(5,T,S)]'.*1e6;

kiw = [1.38,0.83,0.49,0.4,0.5]'; %values from Top et al., (1988).


rsw = 1.028;

%- Noble gases in Obj. fun.
F = (Ceq+ A.*chi).*fsw - .917/rsw.*Ceq.*fsi.*kiw  + fgmw.*Cpure + fcdw.*cdw(1:5); 

%- Heat budget in Obj. fun
F(end+1,1) = T.*(fsw - .917/rsw.*fsi) + .917/rsw.*3.35e5/4180*.1675.*fsi + fgmw.*(-92) + fcdw.*cdw(end); %- Temp

%- Salt budget in Obj. fun
F(end+1,1) = fsw.*S + fgmw.*0 + fcdw.*cdw(end-1) - .917/rsw.*fsi*S.*.33;  %- Salinity 

%- Water conservation constraint.
F(end+1,1) = fsw + fgmw + fcdw - .917/rsw.*fsi;%- continuity


G = w*F - w*cw;

