lambda = 2;
basenodeenergy = ones(1, 10)*2;
nodes = poissrnd(lambda, 1, 10) + basenodeenergy
%wsngen

for i=1:size(nodes,2)
    nodes = nodes-ones(1,10)
end

nodes;

plot(1:10, nodes)