energy_sum=200000;
error_epsilon=0.1;
sample_num=100;
distances_ratio = [1, 2, 5, 10, 20]; 

for distance=distances_ratio
    for num=1:1
        columns = 0:(1/distance):1;
        rows = columns;
        node_coordinates = combvec(rows, columns)';
        BS = length(node_coordinates);
        
        energy = energy_sum / BS;
        
        alpha = 2;
        teta = 1;
        sigma_power_z = 100;
        D=squareform(pdist(node_coordinates));
        node_energy_constants = -1 * (D.^alpha)*teta*sigma_power_z;
        for i=1:length(node_coordinates)
            for j=1:length(node_coordinates)
                node_a = node_coordinates(i,:);
                node_b = node_coordinates(j,:);
                %if node_a(1) ~= node_b(1) && node_a(2) ~= node_b(2)
                %    node_energy_constants(i,j) = node_energy_constants(i,j) * 1000;
                %end
            end
        end
        
        
        
        %BS = length(node_coordinates);
        twohop_res = iot_main(node_energy_constants, BS, energy, error_epsilon, sample_num, 2);
        
        %grid_res_false = iot_main_grid(length(columns), length(rows), node_energy_constants(1, 2), 1 - error_epsilon, energy, sample_num, false);
        
        %grid_res_true = iot_main_grid(length(columns), length(rows), node_energy_constants(1, 2), 1 - error_epsilon, energy, sample_num, true);
        
        %grid_res_multi = iot_main_grid_multipath(length(columns), length(rows), node_energy_constants(1, 2), energy, 1, 1, 1 - error_epsilon, sample_num);
        
        pervasive_res = iot_main_grid_pervasive(length(columns), length(rows), node_energy_constants(1, 2), energy, 1 - error_epsilon, sample_num);
        result_max = [twohop_res pervasive_res];
        mean(pervasive_res(:,2) ./ pervasive_res(:,1))
        writematrix(result_max, "comparison" + distance + "_" + num + ".txt");
    end
end