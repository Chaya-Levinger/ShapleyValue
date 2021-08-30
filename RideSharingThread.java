import java.util.Arrays;

public class RideSharingThread extends myThread {

    RideSharingThread(String caption, boolean toulouseGraph, int n, int OriginIndex, int ind) {
        super(caption, toulouseGraph, n, OriginIndex, ind);
        solv = new RideSharingSolver(Distances);
    }

    protected double costOfByOrder(int[] ind) {
        double cost=0;
        for(int i=1; i<ind.length; i++)
            cost += Distances[ind[i-1]][ind[i]];
        return cost;
    }

    protected double[] SHAPO(){
        double[] shapo = new double[num];
        for(int i=1; i<=num; i++){
            shapo[i-1] += 1.0/i*Distances[0][i];
            for(int q=i+1; q<=num; q++) {
                shapo[i - 1] -= 1.0 / (q * (q - 1)) * Distances[0][q];
                shapo[i-1] += 1.0/((q-i)*(q-i+1))*Distances[i][q];
            }
            for(int p=1; p<i; p++) {
                shapo[i - 1] += 1.0 / ((i - p) * (i - p + 1)) * Distances[p][i];
                for (int q = i + 1; q <= num; q++)
                    shapo[i - 1] -= 2.0 / ((q - p) * (q - p + 1) * (q - p - 1)) * Distances[p][q];
            }
        }
        return shapo;
    }

    protected double[] Shortcut(double totalCost){
        double[] p = new double[num];
        double sum = 0;
        for(int i=1; i<num; i++)
            p[i-1] = Distances[i-1][i] + Distances[i][i+1] - Distances[i-1][i+1];
        p[num-1] = Distances[num-1][num];
        sum = Arrays.stream(p).sum();
        for(int i=0; i<num; i++)
            p[i] = (p[i]/sum)*totalCost;
        return p;
    }

}