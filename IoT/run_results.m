rng('default');

node_counts = [10];
energy = 20000;
error_epsilon = 0.05;
sample_num = 2;
k_nums = [1 2];
random_configuration_num = 5;

alpha = 2;
teta = 1;
sigma_power_z = 5;

for node_count=node_counts
    for i=1:random_configuration_num
        state = rng;
        %save("rng_state_"+node_count+"_"+i+".mat", "state");
        
        BS = node_count;
        node_coordinates = rand(node_count,2);
        D=squareform(pdist(node_coordinates));
        node_energy_constants = -1 * (D.^alpha)*teta*sigma_power_z;
        
        result_part = iot_main(node_energy_constants, BS, energy, error_epsilon, sample_num, k_nums);
        writematrix(result_part, "run_" + node_count + "_" + i + ".txt");
        
        %result_part = iot_main_leach_comp(node_energy_constants, BS, energy, error_epsilon, sample_num, k_nums);
        %writematrix(result_part, "run_leach_hop_" + node_count + "_" + i + ".txt");
        result_part = iot_main_nc(node_energy_constants, BS, energy, error_epsilon, 1, k_nums, 5, 5);
        writematrix(result_part, "run_nc_" + node_count + "_" + i + ".txt");
    end
end