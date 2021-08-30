public class RoutingGamesNode extends Node{

    public RoutingGamesNode(Node parent, double parent_cost, double[][] distances, int[] active_set, int index) {
        super(parent, parent_cost, distances, active_set, index);
    }

    public RoutingGamesNode createChild (double parent_cost, int[] active_set, int index){
        return new RoutingGamesNode(this, parent_cost, distances, active_set, index);
    }

    public double getLowerBound() {

        if(active_set.length == 1)
            return getParentCost() + distances[active_set[0]][0];

        if(active_set.length == 2) {
            double cost = distances[active_set[0]][0];
            if (active_set[0] == index)
                cost = distances[active_set[1]][0];
            return getParentCost() + distances[active_set[0]][active_set[1]] + cost;
        }

        double value = 0;

        for(int location : active_set) {
            double low = distances[location][0];
            for(int other: active_set) {
                if(other == location)
                    continue;
                double cost = distances[location][other];
                if(cost < low) {
                    low = cost;
                }
            }
            value += low;
        }

        return getParentCost() + value ;
    }

    public double getPathCost() {
        return distances[0][index] + getParentCost();
    }
}
