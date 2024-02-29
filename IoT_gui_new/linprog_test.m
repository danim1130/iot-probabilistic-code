node_count = 10;
BS = node_count;
node_coordinates = rand(node_count,2);

alpha = 2;
teta = 1;
sigma_power_z = 100;
D=squareform(pdist(node_coordinates));
node_energy_constants = -1 * (D.^alpha)*teta*sigma_power_z;
node_energy_vector = rand(1, node_count) * 10000 + 1500;

[path, new_energy] = k_hop_step(node_energy_constants, BS, node_energy_vector, 0.9, 4)