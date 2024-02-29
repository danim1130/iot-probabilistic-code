function result = average_energy_test(energy_constants, base_station, error_epsilon, k_num)
    p_succ = 1-error_epsilon;
    
    result = runMultiHopApprxTestLimited(energy_constants, base_station, p_succ, k_num);
end

function energyLevelsMultiApprx = runMultiHopAverageEnergy(energy_constants, base_station, p_succ, hop_number)
    energy_diff_abs = 0.00001;
    %energy_diff_abs = energy * energy_diff_rel;
    
    node_count = length(energy_constants);
    possible_starting_nodes = 1:1:node_count;
    possible_starting_nodes(base_station) = [];
    for source = possible_starting_nodes
        if hop_number > 2
            
            common_energy_apprx_low = 0;
            common_energy_apprx_high = node_energy_vector(source);
            
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
            %[logp_2, graph_path_2] = calcLargestPossibility(energy_constants, source, base_station, delta_energy);
            [logp, graph_path] = calcLargestPossibilityLimited(energy_constants, source, base_station, delta_energy, hop_number);
            assert(logp >= log(p_succ))
        elseif hop_number == 2
            common_energy_apprx_low = -Inf;
            for middle = possible_starting_nodes
                c1 = energy_constants(source, middle);
                c2 = energy_constants(middle, base_station);
                x = node_energy_vector(source);
                y = node_energy_vector(middle);
                a = 1;
                b = -(x + y - (c1 + c2)/log(p_succ));
                c = x*y - (c1 * y + c2 * x) / log(p_succ);
                common_energy = (-b - sqrt(b^2 - 4 * a * c)) / (2 * a);
                %syms common_energy_var
                %eqn = c1/(x - common_energy_var) + c2/(y - common_energy_var) == log(P_Succ_Receive);
                %common_energy = solve([eqn, common_energy_var < x, common_energy_var < y]);
                
                if common_energy > common_energy_apprx_low
                    common_energy_apprx_low = common_energy;
                    if middle == source
                        graph_path = [source];
                    else
                        graph_path = [source, middle];
                    end
                end
            end
            if common_energy_apprx_low < 0
                break;
            end
        else
            common_energy_apprx_low = node_energy_vector(source) - energy_constants(source, base_station) / log(p_succ);
            if common_energy_apprx_low < 0
                break;
            end
            graph_path = [source];
        end
        energyLevelsMultiApprx(1) = energyLevelsMultiApprx(1) + 1;
        transfer_failed = false;
        for i = 1:length(graph_path)
            node = graph_path(i);
            energy_delta = node_energy_vector(node) - common_energy_apprx_low;
            node_energy_vector(node) = common_energy_apprx_low;
            if i ~= length(graph_path)
                node_path_constant = energy_constants(node, graph_path(i + 1));
            else
                node_path_constant = energy_constants(node, base_station);
            end
            message_probability = exp(node_path_constant/energy_delta);
            if (rand() > message_probability)
                transfer_failed = true;
                break;
            end
        end
        if ~transfer_failed
            energyLevelsMultiApprx(2) = energyLevelsMultiApprx(2) + 1; 
        end
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