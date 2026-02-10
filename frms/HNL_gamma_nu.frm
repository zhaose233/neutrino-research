Symbols M, d;
Vectors k, p1, p2, s;
Indices mu, nu;

Off Statistics;


Local AmpSq =
    - g_(1, p1) *
    
    (1 + g5_(1))/2 * g_(1, p2) * g_(1, mu) *
    
    ( g_(1, k) + M ) * (1 + g5_(1) * g_(1, s))/2 *
    
    g_(1, mu) * g_(1, p2) * (1 - g5_(1))/2;

Trace4, 1;

id p1 = k - p2;
.sort

id p2.p2 = 0;
id k.s = 0;
id s.s = -1;
id k.p2 = M^2 / 2;

.sort
Print;
.end