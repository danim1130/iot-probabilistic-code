function ratio = redundancy_ratio_fix_energy(P_message_normal, M, N)
    P_original = P_message_normal^N;
    P_message_new = P_message_normal^(M/N);
    P_new = binocdf(N - 1, M, P_message_new, "upper");
    ratio = P_new / P_original;
end