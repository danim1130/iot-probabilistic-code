error_epsilon = 0.1;
N = 2;
maxM = 6;

p_succ = 1-error_epsilon;
p_succ_message = zeros(1, maxM);
for m_num_calc=N:maxM
    p_search_low = 0;
    p_search_high = 1;
    while 1
        p_search_middle = (p_search_low + p_search_high) / 2;
        p_end_result = 0;
        for ii=N:m_num_calc
            p_end_result = p_end_result + nchoosek(m_num_calc, ii) * p_search_middle^ii * (1 - p_search_middle)^(m_num_calc - ii);
        end
        if p_end_result <= p_succ
           p_search_low = p_search_middle; 
        end
        if p_end_result >= p_succ
           p_search_high = p_search_middle;
        end
        if (abs(p_end_result - p_succ) / p_succ) < 0.0001
            break;
        end
    end
    p_succ_message(m_num_calc) = p_search_high;
end

result = zeros(1, maxM - N + 1);
for i=N:maxM
    result(i - N + 1) = (-1/log(p_succ_message(i))) * i;
end
plot([N:maxM], result)