public class RideSharingNode extends Node{

    public RideSharingNode(Node parent, double parent_cost, double[][] distances, int[] active_set, int index) {
        super(parent, parent_cost, distances, active_set, index);
    }

    public RideSharingNode createChild (double parent_cost, int[] active_set, int index){
        return new RideSharingNode(this, parent_cost, distances, active_set, index);
    }

    public double getLowerBound() {
        double value = 0;
        if(active_set.length == 1)
            return getPathCost();
        if(active_set.length == 2)
            return getParentCost() + distances[active_set[0]][active_set[1]];
        double highest = 0;
        for(int location : active_set) {
            double low = Double.MAX_VALUE;
            for(int other: active_set) {
                if(other == location)
                    continue;
                double cost = distances[location][other];
                if(cost < low)
                    low = cost;
            }
            if(low>highest)
                highest=low;
            value += low;
        }
        return getParentCost() + value - highest;
    }

    public double getPathCost() {
        return getParentCost();
    }
}
