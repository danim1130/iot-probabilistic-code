%redundancy_ratio_fix_probability(0.901^(1/10), 12, 10)

k = 5;
M = 0:1:30;
P = 0.05:0.05:0.95;
results = zeros(length(M), length(P));
for pi=1:length(P)
    for mi=1:length(M)
        results(mi, pi) = redundancy_ratio_fix_probability(P(pi) ^ (1/k), k + M(mi), k);
    end
end
fig = surf(P, M, results(:,:));
xlabel('Message probability');
ylabel('Extra packets');
zlabel('Energy ratio');
title(['Needed message parts = ' num2str(k)]);
%P_message_normal = 0.8;
%M = 10;
%N = 4;
%
%P_succ = P_message_normal^N;
%testfunc = @(x) binocdf(N - 1, M, x, "upper") - P_succ;
%P_solved = fzero(testfunc, 0.5);
%ratio = log(P_message_normal) / log(P_solved);


%function y = testfunc(P)
%    y = binocdf(N - 1, M, P_message_coded, "upper") - P_succ
%end