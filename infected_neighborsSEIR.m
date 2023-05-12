function infected_neighborsget = infected_neighborsSEIR(neighborsAll, node, infected_list)
neighbors = neighborsAll{node};
infected_neighborsget = sum(infected_list(neighbors));