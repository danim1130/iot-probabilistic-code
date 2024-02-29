function result = iot_main_grid_multipath(grid_width, grid_height, w_grid, energy, N, M, p_succ, sample_num)
    result = zeros(sample_num, 1);
    energy_messages_grid = zeros(grid_width, grid_height);
    for node_x = 1:grid_width
        for node_y = 1:grid_height
            if (node_x == grid_width && node_y == grid_height)
                continue;
            end
            if (node_x > node_y)
                energy_messages_grid(node_x, node_y) = energy_messages_grid(node_y, node_x);
                continue;
            end
            syms p;
            total_hops = (grid_width - node_x) + (grid_height - node_y);
            diagonal_hops = min(grid_width - node_x, grid_height - node_y);
            P_overall = (1 - (1 - p^2)^3)^diagonal_hops * p^(total_hops - diagonal_hops);
            
            p_end_result = 0;
            for ii=N:M
                p_end_result = p_end_result + nchoosek(M, ii) * P_overall^ii * (1 - P_overall)^(M - ii);
            end
            %p_search_low = 0;
            %p_search_high = 1;
            %while 1
            %    p_search_middle = (p_search_low + p_search_high) / 2;
            %    
            %    p = p_search_middle;
            %    p_message = subs(P_overall);
            %    
            %    p_end_result = 0;
            %    for ii=N:M
            %        p_end_result = p_end_result + nchoosek(M, ii) * p_message^ii * (1 - p_message)^(M - ii);
            %    end
            %    if p_end_result <= p_succ
            %       p_search_low = p_search_middle; 
            %    end
            %    if p_end_result >= p_succ
            %       p_search_high = p_search_middle;
            %    end
            %    if (abs(p_end_result - p_succ) / p_succ) < 0.0001
            %        break;
            %    end
            %end
            p_hop = vpasolve([p_end_result == p_succ], [0 1]);
            %p_hop = p_search_high;
            energy_messages_grid(node_x, node_y) = (w_grid / N^2) / log(p_hop);
        end
    end
    for i=1:sample_num
        result(i) = runGridTestOnlyClosest(grid_width, grid_height, energy, energy_messages_grid, N, M, p_succ);
    end
end

function message_count = runGridTestOnlyClosest(grid_width, grid_height, energy, energy_messages_grid, N, M, p_succ)
    node_energy_vector = ones(grid_height, grid_width) * energy;
    message_count = 0;
    
    while 1
        node_x = randi(grid_width);
        node_y = randi(grid_height);
        message_energy_near = energy_messages_grid(node_x, node_y);
        
        while(node_x ~= grid_width || node_y ~= grid_height)
            node_energy_vector(node_x, node_y) = node_energy_vector(node_x, node_y) - message_energy_near * M;
            if (node_x ~= grid_width && node_y ~= grid_height)
                node_energy_vector(node_x+1, node_y) = node_energy_vector(node_x+1, node_y) - message_energy_near * M;
                node_energy_vector(node_x, node_y+1) = node_energy_vector(node_x, node_y+1) - message_energy_near * M;
            end
            if (min(node_energy_vector,[],'all') < 0)
                break;
            end
            if (node_x ~= grid_width)
                node_x = node_x + 1;
            end
            if (node_y ~= grid_height)
                node_y = node_y + 1;
            end
        end
        if (min(node_energy_vector,[],'all') < 0)
            break;
        end
        message_count = message_count + 1;
    end
end