function [message_paths, new_energy_vector] = k_hop_step(energy_constants, base_station, node_energy_vector, p_succ, hop_number)
    energy_diff_rel = 0.001;
    energy_diff_abs = max(node_energy_vector) * energy_diff_rel;
    
    node_count = length(energy_constants);
    possible_starting_nodes = 1:1:node_count;
    possible_starting_nodes(base_station) = [];
    energyLevelsMultiApprx = node_energy_vector;
    e = zeros(1, node_count);
    common_energy_apprx_low = 0;
    
    while ~isempty(possible_starting_nodes)
        source = possible_starting_nodes(randi(length(possible_starting_nodes)));
        common_energy_apprx_high = node_energy_vector(source);

        if (calcLargestPossibilityLimited(energy_constants, source, base_station, node_energy_vector, hop_number) < log(p_succ))
            possible_starting_nodes(possible_starting_nodes == source) = [];
        else 
            break;
        end
    end
    
    if isempty(possible_starting_nodes)
        message_paths = [];
        new_energy_vector = 0 * node_energy_vector;
        return;
    end

    while (common_energy_apprx_high - common_energy_apprx_low) > energy_diff_abs
        common_energy_test = (common_energy_apprx_low + common_energy_apprx_high) / 2;
        delta_energy_test = node_energy_vector - ones(1, node_count) * common_energy_test;
        delta_energy_test(delta_energy_test < 0) = 0;
        if calcLargestPossibilityLimited(energy_constants, source, base_station, delta_energy_test, hop_number) < log(p_succ)
            common_energy_apprx_high = common_energy_test;
        else
            common_energy_apprx_low = common_energy_test;    
        end
    end

    delta_energy = node_energy_vector - ones(1, node_count) * common_energy_apprx_low;
    delta_energy(delta_energy < 0) = 0;
    [logp, graph_path] = calcLargestPossibilityLimited(energy_constants, source, base_station, delta_energy, hop_number);
    assert(logp >= log(p_succ))
    node_energy_vector(graph_path) = common_energy_apprx_low;

    new_energy_vector = node_energy_vector;
    message_paths = [];
    for i=2:length(graph_path)
        message_paths(end + 1, :) = graph_path((i-1):i);
    end 
    message_paths(end + 1, :) = [graph_path(end) base_station];
end


function [logp_max, graph_path] = calcLargestPossibilityLimited(energy_constants, source, base_station, energy, hop_number)
    logp_matrix = -(energy_constants ./ (energy' * ones(1, length(energy))));
    logp_matrix(or(isnan(logp_matrix), isinf(logp_matrix))) = Inf;
    logp_matrix(logp_matrix < 0) = Inf;
    
    current_shortest_vector = logp_matrix(source, :);
    current_shortest_path = cell(1, length(current_shortest_vector));
    
    for i = 1:length(current_shortest_vector)
        current_shortest_path{i} = [];
    end
    
    for k = 2:hop_number
        temp_shortest_vector = current_shortest_vector;
        temp_shortest_path = current_shortest_path;
        for i = 1:length(current_shortest_vector)
            if i ~= source
                minIndex = -1;
                minValue = current_shortest_vector(i);
                for j = 1:length(current_shortest_vector)
                    if j ~= i
                        if (minValue > current_shortest_vector(j) + logp_matrix(j, i))
                            minIndex = j;
                            minValue = current_shortest_vector(j) + logp_matrix(j, i);
                        end
                    end
                end
                if (minIndex ~= -1)
                    temp_shortest_vector(i) = minValue;
                    temp_shortest_path{i} = [current_shortest_path{minIndex} minIndex];
                end
            end
        end
        current_shortest_vector = temp_shortest_vector;
        current_shortest_path = temp_shortest_path;
    end
    
    logp_max = -current_shortest_vector(base_station);
    graph_path = [source current_shortest_path{base_station}];
end