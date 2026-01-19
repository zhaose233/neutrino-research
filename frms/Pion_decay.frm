Vector k, p1, p2;
Indices nu, mu;
Symbol ml, mN, mpi, p2m;

Local M2N = k(mu) * k(nu) * (g_(1, p2) + ml) * g_(1, mu) * (1 - g5_(1)) * (g_(1, p1) - mN) * g_(1, nu) * (1 - g5_(1));

trace4, 1;

id k.k  = mpi^2;
id k.p2 = mpi * (p2m^2 + mN^2)^(1/2);
id k.p1 = mpi * (p2m^2 + ml^2)^(1/2);
id p1.p2 = ((p2m^2 + mN^2)^(1/2) * (p2m^2 + ml^2)^(1/2) + p2m^2);

contract;

Print;
.end