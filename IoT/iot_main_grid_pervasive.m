function result = iot_main_grid_pervasive(grid_width, grid_height, w_grid, energy, p_succ, sample_num)
    result = zeros(sample_num, 1);
    for i=1:sample_num
        result(i) = runGridTestOnlyClosest(grid_width, grid_height, w_grid, energy, p_succ);
    end
end

function message_count = runGridTestOnlyClosest(grid_width, grid_height, w_grid, energy, p_succ)
    node_energy_vector = ones(grid_height, grid_width) * energy;
    message_count = 0;
    
    while 1
        node_x = randi(grid_width);
        node_y = randi(grid_height);
        if node_x == grid_width && node_y == grid_height
            continue
        end
        
        syms var_g;
        energy_eq = 1;
        x_delta = grid_width - node_x;
        y_delta = grid_height - node_y;
        for i=0:x_delta
            for j=0:y_delta
                if i == 0 && j == 0
                    continue
                end
                energy_eq = energy_eq * (1 - exp(w_grid * (i^2 + j^2 + (x_delta - i)^2 + (y_delta - j)^2) / var_g));
            end
        end
        min_energy = 0;
        max_energy = min(min(node_energy_vector(node_y:grid_height, node_x:grid_width)));
        var_g = max_energy;
        p_calc = 1 - subs(energy_eq);
        if (p_calc < p_succ)
            break
        end
        while((max_energy - min_energy) / max_energy > 0.001 )
            var_g = (max_energy + min_energy) / 2;
            p_calc = 1 - subs(energy_eq);
            if p_calc < p_succ
                min_energy = (max_energy + min_energy) / 2;
            else
                max_energy = (max_energy + min_energy) / 2;
            end
        end
        
        node_energy_vector(node_y:grid_height, node_x:grid_width) = node_energy_vector(node_y:grid_height, node_x:grid_width) - max_energy;
        if (min(node_energy_vector,[],'all') < 0)
            break;
        end
        message_count = message_count + 1;
    end
end