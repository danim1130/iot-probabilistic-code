function ratio = redundancy_ratio_fix_probability(P_message_normal, M, N)
    P_succ = P_message_normal^N;
    testfunc = @(x) special_binocdf(N - 1, M, x) - P_succ;
    P_solved = fzero(testfunc, 0.5);
    ratio = N / M * log(P_solved) / log(P_message_normal);
end

function y = special_binocdf(x, N, P)
    if P < 0
        y = 0;
    elseif P > 1
        y = 1;
    else
        y = binocdf(x, N, P, "upper");
    end
end