function result = iot_main(energy_constants, distance_matrix, base_station, energy, error_epsilon, sample_num, k_num)
    p_succ = 1-error_epsilon;
    result = zeros(sample_num, 3);
    
    %for i=1:sample_num
        result = runMultiHopApprxTestLimited(energy_constants, distance_matrix, base_station, energy, p_succ, k_num);
    %end
end

function energyLevelsMultiApprx = runMultiHopApprxTestLimited(energy_constants, distance_matrix, base_station, energy, p_succ, hop_number)
    energy_diff_rel = 0.0001;
    energy_diff_abs = energy * energy_diff_rel;
    
    node_count = length(energy_constants);
    node_energy_vector = ones(1, node_count) * energy;
    possible_starting_nodes = 1:1:node_count;
    possible_starting_nodes(base_station) = [];
    starting_array = ["Success", "Hop-number", "Starting distance", "Travel sum", "Min used energy", "Max used energy", "Sum used energy", "Optimal sum energy", "Path optimal sum energy"];
    energyLevelsMultiApprx = starting_array;
    for i=1:length(node_energy_vector)
        energyLevelsMultiApprx = [energyLevelsMultiApprx, "Node energy " + int2str(i)];
    end
    while 1
        current_data_vector = zeros(1, length(energyLevelsMultiApprx(1,:)));
        
        source = 1;%possible_starting_nodes(randi(length(possible_starting_nodes)));
        
        %Optimal energy calculation
        common_energy_apprx_low = 0;
        common_energy_apprx_high = energy;
        test_energy_vector = ones(1, node_count) * energy;
        
        while (common_energy_apprx_high - common_energy_apprx_low) > energy_diff_abs
            common_energy_test = (common_energy_apprx_low + common_energy_apprx_high) / 2;
            delta_energy_test = test_energy_vector - ones(1, node_count) * common_energy_test;
            delta_energy_test(delta_energy_test < 0) = 0;
            if calcLargestPossibilityLimited(energy_constants, source, base_station, delta_energy_test, hop_number) < log(p_succ)
                common_energy_apprx_high = common_energy_test;
            else
                common_energy_apprx_low = common_energy_test;    
            end
        end
        delta_energy = test_energy_vector - ones(1, node_count) * common_energy_apprx_low;
        delta_energy(delta_energy < 0) = 0;
        [logp, graph_path] = calcLargestPossibilityLimited(energy_constants, source, base_station, delta_energy, hop_number);
        sum_optimal_energy = length(graph_path) * (energy - common_energy_apprx_low);
        current_data_vector(8) = sum_optimal_energy;
        
        %Real energy calculation
        common_energy_apprx_low = 0;
        common_energy_apprx_high = node_energy_vector(source);
        current_data_vector(3) = distance_matrix(source, base_station);
        
        if (calcLargestPossibilityLimited(energy_constants, source, base_station, node_energy_vector, hop_number) < log(p_succ))
            break
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
        common_energy = common_energy_apprx_low;
        %[logp_2, graph_path_2] = calcLargestPossibility(energy_constants, source, base_station, delta_energy);
        [logp, graph_path] = calcLargestPossibilityLimited(energy_constants, source, base_station, delta_energy, hop_number);
        assert(logp >= log(p_succ))
        transfer_failed = false;
        
        %Path optimal energy calc
        common_energy_apprx_low = 0;
        common_energy_apprx_high = energy;
        test_energy_vector = zeros(1, node_count);
        test_energy_vector(graph_path) = energy;
        
        while (common_energy_apprx_high - common_energy_apprx_low) > energy_diff_abs
            common_energy_test = (common_energy_apprx_low + common_energy_apprx_high) / 2;
            delta_energy_test = test_energy_vector - ones(1, node_count) * common_energy_test;
            delta_energy_test(delta_energy_test < 0) = 0;
            if calcLargestPossibilityLimited(energy_constants, source, base_station, delta_energy_test, hop_number) < log(p_succ)
                common_energy_apprx_high = common_energy_test;
            else
                common_energy_apprx_low = common_energy_test;    
            end
        end
        test_delta_energy = test_energy_vector - ones(1, node_count) * common_energy_apprx_low;
        test_delta_energy(test_delta_energy < 0) = 0;
        [test_logp, test_graph_path] = calcLargestPossibilityLimited(energy_constants, source, base_station, test_delta_energy, hop_number);
        sum_optimal_energy_test = length(test_graph_path) * (energy - common_energy_apprx_low);
        current_data_vector(9) = sum_optimal_energy_test;
        
        for i = 1:length(graph_path)
            node = graph_path(i);
            energy_delta = node_energy_vector(node) - common_energy;
            if (current_data_vector(5) == 0)
                current_data_vector(5) = energy_delta;
            end
            current_data_vector(5) = min(current_data_vector(5), energy_delta);
            current_data_vector(6) = max(current_data_vector(6), energy_delta);
            current_data_vector(7) = current_data_vector(7) + energy_delta;
            
            node_energy_vector(node) = common_energy;
            distance = 0;
            if i ~= length(graph_path)
                node_path_constant = energy_constants(node, graph_path(i + 1));
                distance = distance_matrix(node, graph_path(i + 1));
            else
                node_path_constant = energy_constants(node, base_station);
                distance = distance_matrix(node, base_station);
            end
            message_probability = exp(node_path_constant/energy_delta);
            if (rand() > message_probability)
                transfer_failed = true;
                break;
            end
            current_data_vector(2) = current_data_vector(2) + 1;
            current_data_vector(4) = current_data_vector(4) + distance;
        end
        current_data_vector(1) = ~transfer_failed;
        current_data_vector((length(starting_array)+1):end) = node_energy_vector;
        energyLevelsMultiApprx = [energyLevelsMultiApprx; current_data_vector];
    end
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