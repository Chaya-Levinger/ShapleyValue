public abstract class Node {
    public Node parent;
    protected double parent_cost;

    protected double[][] distances;
    protected int[] active_set;

    public int index;

    public Node(Node parent, double parent_cost, double[][] distances, int[] active_set, int index) {
        this.parent = parent;
        this.parent_cost = parent_cost;
        this.distances = distances;
        this.active_set = active_set;
        this.index = index;
    }

    public Node[] generateChildren() {
        Node[] children = new Node[active_set.length - 1];

        int[] new_set = new int[active_set.length - 1];
        int i = 0;
        for(int location : active_set) {
            if(location == index)
                continue;

            new_set[i] = location;
            i++;
        }

        for(int j = 0; j < children.length; j++)
            children[j] = this.createChild(distances[index][new_set[j]], new_set, new_set[j]);

        return children;
    }

    public boolean isTerminal() {
        return active_set.length == 1;
    }

    public int[] getPath() {
        int depth = distances.length - active_set.length + 1;
        int[] path = new int[depth];
        getPathIndex(path, depth - 1);
        return path;
    }

    public void getPathIndex(int[] path, int i) {
        path[i] = index;
        if (parent != null)
            parent.getPathIndex(path, i - 1);
    }

    public double getParentCost() {
        if (parent == null)
            return 0;

        return parent_cost + parent.getParentCost();
    }

    public abstract Node createChild (double parent_cost, int[] active_set, int index);

    public abstract double getLowerBound();

    public abstract double getPathCost();

}