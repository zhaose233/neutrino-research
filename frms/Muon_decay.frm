Vector k, p1, p2, p3;
Indices nu, mu;
Symbol sw, me, mN, mmu, p2m, p3m;

Local M2part1 = (g_(1, p2) + me) * g_(1, mu) * (1 - g5_(1)) / 2 * g_(1, p3)* g_(1, nu) * (1 - g5_(1)) / 2;
Local M2part2 = (g_(2, p1) + mN) * g_(2, mu) * (1 - g5_(2)) / 2 * (g_(2, k) + mmu) * g_(2, nu) * (1 - g5_(2)) / 2;

Local Mpart = M2part1 * M2part2;

trace4, 1;
trace4, 2;

contract;

Print;
.end