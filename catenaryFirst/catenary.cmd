option solver ipopt;

solve;

printf {j in 0..N}: "%10.5f %10.5f \n", x[j], y[j];
