function F = fermi_dirac(beta, mu, x)
% Fermi-dirac distribution
% F(x) = 1/(1+exp(beta(x-mu)))

F = 1 ./ (1+exp(beta*(x-mu)));

end