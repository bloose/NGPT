%       function C = LJSolMol(nGas,T,S)
%
%       returns solubility equilibrium gas concentration
%       for any of the 5 stable noble gases in mol/kg
%       in equilibrium with air at one standard atmosphere
%       and 100% relative humidity at temperature T in deg C and
%       salinity S in PSS78 
%       Assuming the following atmospheric abundances:
%       1   He:  5.24 ppmV
%       2   Ne: 18.18 ppmV
%       3   Ar:  0.934%
%       4   Kr:  1.14 ppmV
%       5   Xe: 0.087 ppmV
%       6   He3: 5.24 x 1.384 pptV
%
%       S & T must be scalars, or vectors/arrays with same
%           dimensions. Otherwise an error will result
%       C is returned with same dimensions as S & T
%		Based on Jenkins, Lott, and Cahill (2019)A Determination 
%			of Atmospheric Helium, Neon, Argon, Krypton, and Xenon 
%			Solubility Concentrations in Water and Seawater. 
%			Marine Chemistry, Revised and resubmitted.
%	Requires: LJSolMol.mat
%
%   Example Usage: 
%        C = LJSolMol(2,10,35)
%			returns Ne solubility at 10C and 35 PSS78
%
function C = LJSolMol(nGas,T,S)
load LJSolMol;
cT=(T+273.15)/100;
if isnumeric(nGas)
    if (nGas > 0 ) & (nGas < 6)
        C=exp(Amol(1,nGas)+Amol(2,nGas)./cT+Amol(3,nGas)*log(cT)+ Amol(4,nGas)*cT ...
            + S.*(Amol(5,nGas)+Amol(6,nGas)*cT+Amol(7,nGas)*cT.^2)+Amol(8,nGas)*S.^2);
    else
                C=1.384e-6 * healph(T) .* exp(Amol(1,1)+Amol(2,1)./cT+Amol(3,1)*log(cT)+ Amol(4,1)*cT ...
            + S.*(Amol(5,1)+Amol(6,1)*cT+Amol(7,1)*cT.^2)+Amol(8,1)*S.^2);
    end
end