function Q = inv_hydrograph(F,Qmin,Qmax)
        Q = Qmin + 0.5*(Qmax-Qmin)*(1-cos(pi*F));
end
