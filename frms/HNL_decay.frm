Vector k, p1, p2, p3;
Indices nu, mu;
Symbol sw, me, mN, p2m, p3m;

Local Ml2part1 = (g_(1, p2) + me) * g_(1, mu) * (2 * sw^2 - (1 - g5_(1)) / 2) * (g_(1, p3) - me) * g_(1, nu) * (2 * sw^2 - (1 - g5_(1)) / 2);
Local Ml2part2 = g_(2, p1) * g_(2, mu) * (1 - g5_(2)) / 2 * (g_(2, k) + mN) * g_(2, nu) * (1 - g5_(2)) / 2;

Local Mn2part1 = g_(1, p2) * g_(1, mu) * (1 - g5_(1)) / 2 * g_(1, p3)* g_(1, nu) * (1 - g5_(1)) / 2;
Local Mn2part2 = g_(2, p1) * g_(2, mu) * (1 - g5_(2)) / 2 * (g_(2, k) + mN) * g_(2, nu) * (1 - g5_(2)) / 2;

Local Ml2part = Ml2part1 * Ml2part2;
Local Mn2part = Mn2part1 * Mn2part2;

trace4, 1;
trace4, 2;

contract;

id k.p3 = mN * p3m;
id p1.p2 = mN * (mN - 2 * p3m);
id k.p2 = mN * p2m;
id p1.p3 = mN * (mN - 2 * p2m);
id me = 0;

Print;
.end