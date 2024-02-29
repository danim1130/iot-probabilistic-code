function result = pegasis_main(energy_constants, base_station, energy, error_epsilon, sample_num)
    p_succ = 1-error_epsilon;
    result_message_count = [];
    for i=1:sample_num
        result_message_count(i) = runPEGASISTest(energy_constants, base_station, energy, p_succ);
    end
    result = [min(result_message_count) mean(result_message_count) max(result_message_count)];
end

function message_count = runPEGASISTest(energy_constants, base_station, energy, p_succ)
    node_count = length(energy_constants);
    node_energy_vector = ones(1, node_count) * energy;
    p_succ_per_node = p_succ ^ (1/(node_count - 1));
    node_sending_energy = energy_constants / log(p_succ_per_node);
    message_count = 0;
    
    node_sending_energys_copy = -energy_constants;
    [~, max_index] = max(node_sending_energys_copy(:, base_station));
    node_sending_energys_copy(:, [base_station, max_index]) = Inf;
    chain = [max_index];
    
    while length(chain) ~= node_count - 1
        next_node = chain(end);
        [~, min_index] = min(node_sending_energys_copy(next_node, :));
        node_sending_energys_copy(:, min_index) = Inf;
        chain = [chain min_index];
    end
    
    base_head_index = 1;
    chain_end_index = length(chain);
    has_dead_node = false;
    
    while has_dead_node == false
        for i = 1:(base_head_index - 1)
            sending_message_count = i;
            current_node = chain(i);
            next_node = chain(i+1);
            message_energy = node_sending_energy(current_node, next_node) * sending_message_count;
            node_energy_vector(current_node) = node_energy_vector(current_node) - message_energy;
            if (node_energy_vector(current_node) < 0)
                has_dead_node = true;
                break;
            end
        end
        if has_dead_node
            break;
        end
        
        for i = chain_end_index:-1:(base_head_index + 1)
            sending_message_count = chain_end_index - i + 1;
            current_node = chain(i);
            next_node = chain(i-1);
            message_energy = node_sending_energy(current_node, next_node) * sending_message_count;
            node_energy_vector(current_node) = node_energy_vector(current_node) - message_energy;
            if (node_energy_vector(current_node) < 0)
                has_dead_node = true;
                break;
            end
        end
        if has_dead_node
            break;
        end
        head_sending_energy = node_sending_energy(base_head_index, base_station) * length(chain);
        node_energy_vector(base_head_index) = node_energy_vector(base_head_index) - head_sending_energy;
        if node_energy_vector(base_head_index) < 0
            has_dead_node = true;
            break;
        end
        message_count = message_count + length(chain);
        base_head_index = mod(base_head_index, length(chain)) + 1;
    end
end