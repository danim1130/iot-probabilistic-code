function result = leach_main(energy_constants, base_station, energy, error_epsilon, sample_num, p_node_election, compression_ratio)
    p_succ = 1-error_epsilon;
    result_message_count = [];
    for i=1:sample_num
        result_message_count(i) = runLEACHTest(energy_constants, base_station, energy, p_succ, p_node_election, compression_ratio);
    end
    result = result_message_count;
end

function message_count = runLEACHTest(energy_constants, base_station, energy, p_succ, p_node_election, compression_ratio)
    node_count = length(energy_constants);
    node_energy_vector = ones(1, node_count) * energy;
    p_sending_succ = sqrt(p_succ);
    node_sending_energy = energy_constants / log(p_sending_succ);
    has_dead_node = false;
    received_message = 0;
    
    while has_dead_node == false
        selected_node_heads = zeros(1, node_count);
        while selected_node_heads == zeros(1, node_count)
            selected_node_heads = rand(1, node_count) < p_node_election;    
        end
        message_counter = zeros(1, node_count);
        message_counter(selected_node_heads) = 1;
        message_counter(base_station) = 0;
        selected_node_heads(base_station) = 1;
        for i = 1:node_count
            if (selected_node_heads(i) == 0)
                possible_energies = node_sending_energy(i, :);
                filtered_energies = possible_energies(selected_node_heads);
                sending_energy = min(filtered_energies);
                target_head = find(possible_energies == sending_energy, 1, 'first');
                message_counter(target_head) = message_counter(target_head) + 1;
                node_energy_vector(i) = node_energy_vector(i) - sending_energy;
                if (node_energy_vector(i) < 0)
                    has_dead_node = true;
                    break
                end
            end
        end
        if (has_dead_node)
            break;
        end
        for i = 1:node_count
            if (selected_node_heads(i) == 1)
                spent_energy = node_sending_energy(i, base_station) * message_counter(i) * compression_ratio;
                node_energy_vector(i) = node_energy_vector(i) - spent_energy;
                if node_energy_vector(i) < 0
                    has_dead_node = true;
                    break
                end
                received_message = received_message + message_counter(i); 
            end
        end
    end
    message_count = received_message;
end