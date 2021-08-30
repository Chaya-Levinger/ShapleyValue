import java.util.HashSet;
import java.util.LinkedList;

public class RoutingGamesSolver extends Solver{

    public RoutingGamesSolver(double[][] Distances) {
        super(Distances);
        int[] active_set = new int[distances.length];
        for(int i = 0; i < active_set.length; i++)
            active_set[i] = i;
        root = new RoutingGamesNode(null, 0, distances, active_set, 0);
    }

    protected double GreedyInitialization(int startInd, int i, HashSet<Integer> location_set, LinkedList<Integer> greedy_best_path) {

        greedy_best_path.add(0, i);
        location_set.remove(i);

        if(location_set.isEmpty()){
            best_path = new int[greedy_best_path.size()];
            for(int j=0; j<greedy_best_path.size(); j++)
                best_path[greedy_best_path.size()-j-1] = greedy_best_path.get(j);
            return distances[startInd][i];
        }
        double lowest = Double.MAX_VALUE;
        int closest = 0;
        for(int location : location_set) {
            double cost = distances[i][location];
            if(cost < lowest) {
                lowest = cost;
                closest = location;
            }
        }
        return lowest + GreedyInitialization(startInd, closest, location_set, greedy_best_path);
    }

    protected RoutingGamesSolver newSolver(double[][] Distances){
        return new RoutingGamesSolver(Distances);
    }

}
