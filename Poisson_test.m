computed = load('pressure');
true = load('true');
error = abs(computed - true);
dx = 1/length(error);
norm_of_2 = norm(error, 2);
norm_of_2 = dx*norm_of_2;
