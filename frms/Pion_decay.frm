Vector k, p1, p2, s;
Indices nu, mu;
Symbol ml, mN, mpi, p2m;

Local M2N = k(mu) * k(nu) * (g_(1, p1) - ml) * g_(1, mu) * (1 - g5_(1)) * (g_(1, p2) + mN) * g_(1, nu) * (1 - g5_(1));

Local M2NS = k(mu) * k(nu) * (g_(1, p1) - ml) * g_(1, mu) * (1 - g5_(1)) * (g_(1, p2) + mN) * (1 + g_(1, s)) / 2 * g_(1, nu) * (1 - g5_(1));

trace4, 1;


contract;

Print;
.end