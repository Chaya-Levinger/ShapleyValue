import org.jgrapht.Graph;
import org.jgrapht.alg.tour.ChristofidesThreeHalvesApproxMetricTSP;
import org.jgrapht.graph.*;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;

public abstract  class Solver {
    double[][] distances;
    double best_cost;
    int[] best_path;
    Node root;

    public Solver(double[][] Distances) {
        distances = new double[Distances.length][Distances[0].length];
        for(int i = 0; i < distances.length; i++) {
            for(int j = 0; j < distances[0].length; j++)
                distances[i][j] = Distances[i][j];
        }
    }

    public double getCost() {
        return best_cost;
    }

    public int[] calculate() {
        HashSet<Integer> location_set = new HashSet<Integer>(distances.length);
        for(int i = 0; i < distances.length; i++)
            location_set.add(i);

        best_cost = GreedyInitialization(0,0, location_set, new LinkedList<>());

        traverse(root);

        return best_path;
    }

    public Integer[] ChristofidesGreedyCalculate() {
        Graph<Integer, DefaultWeightedEdge> graph = new DefaultUndirectedWeightedGraph<>(DefaultWeightedEdge.class);
        for(int i=0; i<distances.length; i++)
            graph.addVertex(i);
        for(int i=0; i<distances.length; i++)
            for(int j=i+1; j<distances.length; j++){
                graph.addEdge(i,j);
                graph.setEdgeWeight(i,j,distances[i][j]);
            }
        ChristofidesThreeHalvesApproxMetricTSP<Integer, DefaultWeightedEdge> Christo = new ChristofidesThreeHalvesApproxMetricTSP<>();
        List<Integer> lis = Christo.getTour(graph).getVertexList();
        lis.remove(lis.size()-1);
        Integer[] arr = new Integer[lis.size()];
        lis.toArray(arr);
        return arr;
    }

    public int[] NNGreedyCalculate(int startInd) {
        HashSet<Integer> location_set = new HashSet<Integer>(distances.length);
        for(int i = 0; i < distances.length; i++)
            location_set.add(i);

        best_cost = GreedyInitialization(startInd, startInd, location_set, new LinkedList<>());
        return best_path;
    }

    private void traverse(Node parent) {
        Node[] children = parent.generateChildren();

        for(Node child : children) {
            if(child.isTerminal()) {
                double cost = child.getPathCost();
                if(cost < best_cost) {
                    best_cost = cost;
                    best_path = child.getPath();
                }
            }
            else if(child.getLowerBound() <= best_cost) {
                traverse(child);
            }
        }
    }

    protected abstract double GreedyInitialization(int startInd, int i, HashSet<Integer> location_set, LinkedList<Integer> greedy_best_path);

    protected abstract Solver newSolver(double[][] Distances);
}