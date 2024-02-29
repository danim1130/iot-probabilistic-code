function [message_paths, new_energy_vector] = leach_step(energy_constants, base_station, node_energy_vector, p_succ, p_node_election, compression_ratio)
    node_count = length(energy_constants);
    p_sending_succ = sqrt(p_succ);
    node_sending_energy = energy_constants / log(p_sending_succ);
    
    selected_node_heads = rand(1, node_count) < p_node_election;
    message_counter(selected_node_heads) = 1;
    message_counter(base_station) = 0;
    selected_node_heads(base_station) = 1;
    message_path_pairs = [];
    for i = 1:node_count
        if (selected_node_heads(i) == 0)
            possible_energies = node_sending_energy(i, :);
            filtered_energies = possible_energies(selected_node_heads);
            sending_energy = min(filtered_energies);
            target_head = find(possible_energies == sending_energy, 1, 'first');
            node_energy_vector(i) = node_energy_vector(i) - sending_energy;
            if (node_energy_vector(i) >= 0)
                message_path_pairs(end + 1, :) = [i target_head];
            else
                node_energy_vector(i) = 0;
            end
        end
    end
    for i = 1:node_count
        if (selected_node_heads(i) == 1)
            spent_energy = node_sending_energy(i, base_station) * message_counter(i) * compression_ratio;
            node_energy_vector(i) = node_energy_vector(i) - spent_energy;
            if node_energy_vector(i) >= 0
                message_path_pairs(end + 1, :) = [i base_station];
            end
        end
    end
    message_paths = message_path_pairs;
    new_energy_vector = node_energy_vector;
end