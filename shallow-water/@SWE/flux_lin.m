% Thu Apr 28 04:49:29 MSD 2011
% Karl Kästner
%
%% linearised st-venant equation
% SWE::flux_lin
function f = flux_lin(t, q, J0)
        f = J0*q;
end % SWE::flux_lin

