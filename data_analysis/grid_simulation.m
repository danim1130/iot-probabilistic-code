my_seed = 0;
rand(my_seed);
energy_sum=10000;
error_epsilon=0.05;
sample_num=5;
distance = 2;

columns = 0:(1/distance):1;
rows = columns;
node_coordinates = [];
for i=1:length(rows)
    for j=1:length(columns)
        node_coordinates = [node_coordinates; rows(i) columns(j)];
    end
end
BS = length(node_coordinates);

energy = energy_sum;

alpha = 2;
teta = 1;
sigma_power_z = 5;
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
twohop_res = iot_main(node_energy_constants, D, BS, energy, error_epsilon, sample_num, 2);
writematrix(twohop_res, "test_twohop_result.csv");
%pervasive_res = iot_main_grid_pervasive(length(columns), length(rows), node_energy_constants(1, 2), 1/distance, energy, 1 - error_epsilon, sample_num);
%writematrix(pervasive_res, "pervasive_result.csv");