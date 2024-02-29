%x1 = randi([0, 3], 5)

%x = rand(7,9);
%n = 20;
%ndx = randperm(numel(x), n);
%[row,col] = ind2sub(size(x), ndx)

lambda = 2;
poissrnd(lambda, 1, 10)