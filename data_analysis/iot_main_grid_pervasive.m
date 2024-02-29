function result = iot_main_grid_pervasive(grid_width, grid_height, w_grid, d_grid, energy, p_succ, sample_num)
    %result = zeros(sample_num, 3);
    precomputed_energies = zeros(grid_height, grid_width);
    for node_y=1:grid_height
        for node_x=1:grid_width
            if (node_x == grid_width && node_y == grid_height)
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
            max_energy = energy;
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
            precomputed_energies(node_y, node_x) = max_energy;
        end
    end
    
    %for i=1:sample_num
        result = runGridTestOnlyClosest(grid_width, grid_height, w_grid, d_grid, energy, precomputed_energies);
    %end
end

function result_data = runGridTestOnlyClosest(grid_width, grid_height, w_grid, d_grid, energy, precomputed_energies)
    node_energy_vector = ones(grid_height, grid_width) * energy;
    result_data = ["Success", "Sent messages", "Starting distance", "One used energy"];
    for i=1:(grid_width * grid_height)
        result_data = [result_data, "Node energy " + int2str(i)];
    end
    while 1
        current_data_vector = zeros(1, length(result_data(1,:)));
        node_x = randi(grid_width);
        node_y = randi(grid_height);
        if node_x == grid_width && node_y == grid_height
            continue
        end
        
        max_energy = precomputed_energies(node_y, node_x);
        
        if (min(node_energy_vector(node_y:grid_height, node_x:grid_width) - max_energy,[],'all') < 0)
            break;
        end
        probability_matrix_1 = ones(grid_height, grid_width) * inf;
        probability_matrix_2 = ones(grid_height, grid_width) * inf;
        for i=node_x:grid_width
            for j=node_y:grid_height
                probability_matrix_1(j, i) = (i - node_x)^2 + (j - node_y)^2;
                probability_matrix_2(j, i) = (grid_width - i)^2 + (grid_height - j)^2;
            end
        end
        current_data_vector(4) = precomputed_energies(node_y, node_x);
        current_data_vector(3) = d_grid * sqrt((grid_width - node_x)^2 + (grid_height - node_y)^2);
        probability_matrix_1 = exp(w_grid * probability_matrix_1 / max_energy);
        probability_matrix_2 = exp(w_grid * probability_matrix_2 / max_energy);
        message_received_nodes = rand(grid_height, grid_width) < probability_matrix_1;
        current_data_vector(2) = current_data_vector(2) + 1;
        message_received_nodes(grid_height, grid_width) = 0; %Checked in second round
        current_data_vector(2) = current_data_vector(2) + sum(message_received_nodes,'all');
        node_energy_vector(message_received_nodes) = node_energy_vector(message_received_nodes) - max_energy;
        base_station_received_from = rand(grid_height, grid_width) < probability_matrix_2;
        if any(message_received_nodes .* base_station_received_from, 'all')
            current_data_vector(1) = 1;
        end
        node_energy_vector(grid_height, grid_width) = energy;
        current_data_vector(5:end) = reshape(node_energy_vector', 1, []);
        result_data = [result_data; current_data_vector];
    end
end