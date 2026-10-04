function [lo, hi] = wilsonCI(p, n, z)
%WILSONCI  Wilson score confidence interval of a proportion.
%
%   [lo, hi] = wilsonCI(p, n)      95 % interval for an observed share p of n runs
%   [lo, hi] = wilsonCI(p, n, z)   other confidence level (z = 1.96 for 95 %)
%
% Unlike p +- z*sqrt(p(1-p)/n), the interval stays inside [0, 1] and is not
% degenerate for p = 0 or p = 1 (e.g. 0 of 20 runs: [0, 16 %]).

    if nargin < 3, z = 1.96; end
    den = 1 + z.^2./n;
    mid = (p + z.^2./(2*n))./den;
    half = z.*sqrt(p.*(1 - p)./n + z.^2./(4*n.^2))./den;
    lo = max(0, mid - half);
    hi = min(1, mid + half);
end
