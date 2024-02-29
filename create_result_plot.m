node_count = 10;
prefix = 'new_network_coding/run_';
%files = dir(strcat('results/', 'run_topography.txt'));%prefix, int2str(node_count), '_', '*.txt'));
files = dir(strcat(prefix, int2str(node_count), '_', '*.txt'));
result_matrix = [];
for i=1:length(files)
    result_matrix = [result_matrix; readmatrix(strcat('new_network_coding/', files(i).name))];
end
%boxplot(result_matrix(:,(1:9)), "Labels", ["Leach" "Direct" "2-hop" "3-hop" "4-hop" "5-hop" "6-hop" "7-hop" "8-hop"]);
hold on;
bar(mean(result_matrix))
ylabel("Sent message");
xticks([1 2 ]);
xlim([0.5 2.5]);
xticklabels(["No NC" "15/10 NC"]);
title(strcat("Node count: ", int2str(node_count)));
hold off
%title("Average message count for 100 nodes", "FontWeight","normal", "FontSize", 10)