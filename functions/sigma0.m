function sigma = sigma0(Z,Gamma,p,d,w,Nim,Npx,NA,M)
d = abs(d);
if nargin==7
    sigma = p.*Gamma./d .* sqrt(2./(Z.*Nim.*w)) .* sqrt(log10(Npx.*Npx))/(8*sqrt(2));
else
    NA0 = sin(atan(1.22*Gamma*M/(4*d)));
    fracNA = 3.5379e-9*NA./(NA0-NA);
    sigma = (Z.*Gamma./d*p/(8*sqrt(2)) + fracNA) .* sqrt(2./(Nim.*w)) .* sqrt(log10(Npx.*Npx));
end
