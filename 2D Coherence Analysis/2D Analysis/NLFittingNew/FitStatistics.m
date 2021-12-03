function [ Stat ] = FitStatistics( hess, RSS, DoF, beta, alpha, N)
    Stat.RSS = RSS;
    Stat.DoF = DoF;
    Stat.SRSS = RSS/DoF;
    htol = eps*norm(hess,'fro');
    [hQ, hR, hE] = qr(hess,0);
    [~,hZ] = sort(hE);
    hk = sum(diag(abs(hR))>htol);
    invhess = rankKGeneralInv(hQ,hR,hk,hZ);
    Stat.cov = Stat.SRSS*invhess;
    Stat.stderr = sqrt(diag(Stat.cov));
    Stat.relerr = Stat.stderr./beta;
    n = length(beta) + numel(alpha); %total # params in model (including intercept)
    K = n+1; %you get 1 more for the variance of the residual
    Stat.AICc = N*log(RSS/N) + 2*K + 2*K*(K+1)/(N-K-1);
end

