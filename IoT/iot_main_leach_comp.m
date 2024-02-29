function result = iot_main_leach_comp(energy_constants, base_station, energy, error_epsilon, sample_num, k_nums)
    
    p_succ = 1-error_epsilon;
    result = zeros(sample_num, length(k_nums) + 1);
    
    for i=1:sample_num
        rng_saved = rng;
        for j=1:length(k_nums) + 1
            if j == 1
                message_count = runLEACHTest(energy_constants, base_station, energy, p_succ, 0.1, 0.4);
            else
                k_num = k_nums(j - 1);
                rng(rng_saved);
                if (k_num == 1)
                    message_count = runOneHopTest(energy_constants, base_station, energy, p_succ);
                elseif (k_num == 2)
                    message_count = runTwoHopTest(energy_constants, base_station, energy, p_succ);
                else
                    message_count = runMultiHopApprxTestLimited(energy_constants, base_station, energy, p_succ, k_num); 
                end
            end
            
            result(i, j) = message_count;
        end
    end
end

function message_count = runLEACHTest(energy_constants, base_station, energy, p_succ, p_node_election, compression_ratio)
    node_count = length(energy_constants);
    node_energy_vector = ones(1, node_count) * energy;
    p_sending_succ = sqrt(p_succ);
    node_sending_energy = energy_constants / log(p_sending_succ);
    has_dead_node = false;
    received_message = 0;
    
    while has_dead_node == false
        selected_node_heads = rand(1, node_count) < p_node_election;
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

function energyLevelsMultiApprx = runMultiHopApprxTestLimited(energy_constants, base_station, energy, p_succ, hop_number)
    energy_diff_rel = 0.0001;
    energy_diff_abs = energy * energy_diff_rel;
    
    node_count = length(energy_constants);
    node_energy_vector = ones(1, node_count) * energy;
    possible_starting_nodes = 1:1:node_count;
    possible_starting_nodes(base_station) = [];
    energyLevelsMultiApprx = 0;
    
    messageIndex = length(possible_starting_nodes) + 1;
    while 1
        if messageIndex > length(possible_starting_nodes)
            messageIndex = 1;
            message_perm = possible_starting_nodes(randperm(length(possible_starting_nodes)));
        end
        source = message_perm(messageIndex);%possible_starting_nodes(randi(length(possible_starting_nodes)));
        messageIndex = messageIndex + 1;
        
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
        node_energy_vector(graph_path) = common_energy_apprx_low;
        
        energyLevelsMultiApprx = energyLevelsMultiApprx + 1;
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

function energyLevelsOneHop = runOneHopTest(energy_constants, base_station, energy, p_succ)
    node_count = length(energy_constants);
    node_energy_vector = ones(1, node_count) * energy;
    possible_starting_nodes = 1:1:node_count;
    possible_starting_nodes(base_station) = [];
    energyLevelsOneHop = 0;
    messageIndex = length(possible_starting_nodes) + 1;
    while 1
        if messageIndex > length(possible_starting_nodes)
            messageIndex = 1;
            message_perm = possible_starting_nodes(randperm(length(possible_starting_nodes)));
        end
        source = message_perm(messageIndex);%possible_starting_nodes(randi(length(possible_starting_nodes)));
        messageIndex = messageIndex + 1;
        e = zeros(1, node_count);
        e(source) = energy_constants(source, base_station) / log(p_succ);
        node_energy_vector = node_energy_vector - e;
        if min(node_energy_vector) < 0
            break
        end
        energyLevelsOneHop = energyLevelsOneHop + 1;
    end
end

function energyLevelsTwoHop = runTwoHopTest(energy_constants, base_station, energy, p_succ)
    node_count = length(energy_constants);
    node_energy_vector = ones(1, node_count) * energy;
    possible_starting_nodes = 1:1:node_count;
    possible_starting_nodes(base_station) = [];
    possible_middle_nodes = possible_starting_nodes;
    energyLevelsTwoHop = 0;
    
    messageIndex = length(possible_starting_nodes) + 1;
    while 1
        if messageIndex > length(possible_starting_nodes)
            messageIndex = 1;
            message_perm = possible_starting_nodes(randperm(length(possible_starting_nodes)));
        end
        source = message_perm(messageIndex);%possible_starting_nodes(randi(length(possible_starting_nodes)));
        messageIndex = messageIndex + 1;
        best_delta_e = ones(1, node_count) * Inf;
        for middle = possible_middle_nodes
            
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

            e = zeros(1, node_count);
            e(source) = x - common_energy;
            e(middle) = y - common_energy;
            
            if best_delta_e(source) > e(source)
                best_delta_e = e;
            end
        end
        node_energy_vector = node_energy_vector - best_delta_e;
        if min(node_energy_vector) < 0
            break
        end
        energyLevelsTwoHop = energyLevelsTwoHop + 1;
    end
end

function e = calcTwoHopEnergy(START, MIDDLE, END)
    global energy_constants
    global P_Succ_Receive
    global NODE_COUNT
    global NODE_ENERGY_VECTOR
    c1 = energy_constants(START, MIDDLE);
    c2 = energy_constants(MIDDLE, END);
    x = NODE_ENERGY_VECTOR(START);
    y = NODE_ENERGY_VECTOR(MIDDLE);
    a = 1;
    b = -(x + y - (c1 + c2)/log(P_Succ_Receive));
    c = x*y - (c1 * y + c2 * x) / log(P_Succ_Receive);
    common_energy = (b + sqrt(b^2 - 4 * a * c)) / 2;
    
    %syms common_energy_var
    %eqn = c1/(x - common_energy_var) + c2/(y - common_energy_var) == log(P_Succ_Receive);
    %common_energy = solve([eqn, common_energy_var < x, common_energy_var < y]);
    
    e = zeros(1, NODE_COUNT);
    e(START) = x - common_energy;
    e(MIDDLE) = y - common_energy;
end

function energyLevelsMultiApprx = runMultiHopApprxTest(energy_constants, base_station, energy, p_succ)
    energy_diff_rel = 0.0001;
    energy_diff_abs = energy * energy_diff_rel;
    
    node_count = length(energy_constants);
    node_energy_vector = ones(1, node_count) * energy;
    possible_starting_nodes = 1:1:node_count;
    possible_starting_nodes(base_station) = [];
    energyLevelsMultiApprx = 0;
    while 1
        source = possible_starting_nodes(randi(length(possible_starting_nodes)));
        e = zeros(1, node_count);
        
        common_energy_apprx_low = 0;
        common_energy_apprx_high = node_energy_vector(source);
        
        if (calcLargestPossibility(energy_constants, source, base_station, node_energy_vector) < log(p_succ))
            break
        end
        
        while (common_energy_apprx_high - common_energy_apprx_low) > energy_diff_abs
            common_energy_test = (common_energy_apprx_low + common_energy_apprx_high) / 2;
            delta_energy_test = node_energy_vector - ones(1, node_count) * common_energy_test;
            delta_energy_test(delta_energy_test < 0) = 0;
            if calcLargestPossibility(energy_constants, source, base_station, delta_energy_test) < log(p_succ)
                common_energy_apprx_high = common_energy_test;
            else
                common_energy_apprx_low = common_energy_test;    
            end
        end
        
        delta_energy = node_energy_vector - ones(1, node_count) * common_energy_apprx_low;
        delta_energy(delta_energy < 0) = 0;
        [logp, graph_path] = calcLargestPossibility(energy_constants, source, base_station, delta_energy);
        assert(logp >= log(p_succ))
        graph_path(length(graph_path)) = [];
        node_energy_vector(graph_path) = common_energy_apprx_low;
        
        energyLevelsMultiApprx = energyLevelsMultiApprx + 1;
    end
end

function [logp_max, graph_path] = calcLargestPossibility(energy_constants, source, base_station, energy)
    logp_matrix = -(energy_constants ./ (energy' * ones(1, length(energy))));
    logp_matrix(or(isnan(logp_matrix), isinf(logp_matrix))) = 0;
    
    logp_graph = digraph(logp_matrix);
    [graph_path, logp] = shortestpath(logp_graph, source, base_station);
    logp_max = -logp;
end