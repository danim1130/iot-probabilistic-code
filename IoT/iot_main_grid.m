function result = iot_main_grid(grid_width, grid_height, w_grid, p_succ, energy, sample_num, use_diagonal)
    result = zeros(sample_num, 1);
    for i=1:sample_num
        result(i) = runGridTestOnlyClosest(grid_width, grid_height, w_grid, p_succ, energy, use_diagonal);
    end
end

function message_count = runGridTestOnlyClosest(grid_width, grid_height, w_grid, p_succ, energy, use_diagonal)
    node_energy_vector = ones(grid_height, grid_width) * energy;
    message_count = 0;
    while 1
        node_x = randi(grid_width);
        node_y = randi(grid_height);
        total_hops = (grid_width - node_x) + (grid_height - node_y);
        p_hop = p_succ ^ (1/total_hops);
        message_energy_near = w_grid / log(p_hop);
        if (node_x == grid_width && node_y == grid_height)
            continue;
        end
        while(node_x ~= grid_width || node_y ~= grid_height)
            node_energy_vector(node_x, node_y) = node_energy_vector(node_x, node_y) - message_energy_near;
            if (min(node_energy_vector,[],'all') < 0)
                break;
            end
            if (node_x == grid_width)
                node_y = node_y + 1;
            elseif (node_y == grid_height) 
                node_x = node_x + 1;
            else
                if use_diagonal
                    node_energies = [node_energy_vector(node_y + 1, node_x) node_energy_vector(node_y, node_x + 1) node_energy_vector(node_y + 1, node_x + 1)];
                else
                    node_energies = [node_energy_vector(node_y + 1, node_x) node_energy_vector(node_y, node_x + 1)];
                end
                [~, maxI] = max(node_energies);
                if maxI ~= 1
                    node_x = node_x + 1;
                end
                if maxI ~= 2
                    node_y = node_y + 1;
                end
            end
        end
        if (min(node_energy_vector,[],'all') < 0)
            break;
        end
        message_count = message_count + 1;
    end
end