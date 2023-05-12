function infected_neighborsget = infected_neighborsSEIR(neighborsAll, node, infected_list)
%Calculates the number of infected neighbors for the node
% return sum([infected_list[node_i] for node_i in G.neighbors(node)]) 

% neighbors = unique([G(find(G(:,1)==node),2); G(find(G(:,2)==node),1)]);
% infected_neighborsget=sum(ismember(infected_list,neighbors));
neighbors = neighborsAll{node};
infected_neighborsget = sum(infected_list(neighbors));