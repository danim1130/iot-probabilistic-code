rng('default');

node_counts = [5, 10, 20];
%sum_energy = 3000; for 10
sum_energy = 400;
error_epsilon = 0.2;
sample_num = 5;
k_nums = [2];
random_configuration_num = 5;

alpha = 4;
teta = 1.15;
sigma_power_z = 1.86;

for node_count=node_counts
    energy = sum_energy / node_count;
    for i=1:random_configuration_num
        state = rng;
        %save("rng_state_"+node_count+"_"+i+".mat", "state");
        
        BS = node_count;
        node_coordinates = rand(node_count,2);
        D=squareform(pdist(node_coordinates));
        node_energy_constants = -1 * (D.^alpha)*teta*sigma_power_z;
        
        %for k_num=k_nums
        %    result_part = iot_main(node_energy_constants, BS, energy, error_epsilon, sample_num, k_num);
        %    writematrix(result_part, "run_" + node_count + "_" + i + '_' + k_num + ".txt");
        %end

        for sample = 1:sample_num
            result_matrix = zeros(1, length(k_nums));
            for k_num=k_nums
                tic;
                result_part = iot_main(node_energy_constants, BS, energy, error_epsilon, 1, k_num);
                toc
                result_matrix(:, k_nums == k_num) = result_part(:, 2);
            end
            writematrix(result_matrix, "run_" + node_count + "_" + i + "_" + sample + ".txt");
        end
        
        %result_part = leach_main(node_energy_constants, BS, energy, error_epsilon, sample_num, 0.05, 0.5);
        %writematrix(result_part, "run_leach_hop_test_" + node_count + "_" + i + ".txt");
        %result_part = iot_main_nc(node_energy_constants, BS, energy, error_epsilon, 1, k_nums, 5, 5);
        %writematrix(result_part, "run_nc_" + node_count + "_" + i + ".txt");
    end
    node_count
end