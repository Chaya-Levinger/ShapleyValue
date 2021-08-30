import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.*;
import java.util.stream.*;

public abstract class myThread extends Thread {

    protected static int iterations = 100;
    protected int num;
    protected double[][] Distances;
    protected int originIndex;
    protected float sumReRoutTime = (float)0.;
    protected float sumTSPTime = (float)0.;
    protected boolean ToulouseGraph;
    protected Solver solv;
    public int index;
    public double totalCost;

    myThread(String caption, boolean toulouseGraph, int n, int OriginIndex, int ind){
        super(caption);
        num = n;
        this.originIndex = OriginIndex;
        Distances = new double[num+1][num+1];
        ToulouseGraph = toulouseGraph;
        index = ind;
    }

    public void run() {
        if(ToulouseGraph)
            System.out.println("start\tToulouse\t" + num + " passengrs\t" + originIndex);
        else System.out.println("start\tNewYork\t\t" + num + " passengrs\t" + originIndex);
        double[] shapley = new double[num];
        double[] SHAPO = new double[num];
        double[] ReRout = new double[num];
        double[] Short = new double[num];
        double[] depot = new double[num];
        double[] Appro1 = new double[num];
        double[] Appro2 = new double[num];

        double[] shapley_SHAPO = new double[5];// {sumPercent, sumABS, sumPow, sumMaxABS, RMSE}
        double[] shapley_Short = new double[5];// {sumPercent, sumABS, sumPow, sumMaxABS, RMSE}
        double[] shapley_depot = new double[5];// {sumPercent, sumABS, sumPow, sumMaxABS, RMSE}
        double[] shapley_ReRute = new double[5];// {sumPercent, sumABS, sumPow, sumMaxABS, RMSE}
        double[] shapley_Appro1 = new double[5];// {sumPercent, sumABS, sumPow, sumMaxABS, RMSE}
        double[] shapley_Appro2 = new double[5];// {sumPercent, sumABS, sumPow, sumMaxABS, RMSE}

        int[] inputIndexes;
        int[] initOrder;

        double sumShapley = 0;
        float sumShapleyTime = (float)0.;
        float sumShortTime = (float)0.;
        float sumDepotTime = (float)0.;
        float sumSHAPOTime = (float)0.;
        float sumAppro1Time = (float)0.;
        float sumAppro2Time = (float)0.;
        long StartTime, EndTime;

        for(int i=0; i<iterations; i++){
            if(ToulouseGraph)
                System.out.println("start\tToulouse\t" + num + " passengrs\t" + originIndex + "_" + index + "\t iteration " + i);
            else System.out.println("start\tNewYork\t\t" + num + " passengrs\t" + originIndex + "_" + index + "\t iteration " + i);
            inputIndexes = randomInput();
            initOrder = setInitOrder();
            StartTime = System.currentTimeMillis();
            shapley = Shapley(ReRout);
            EndTime = System.currentTimeMillis();
            sumShapleyTime += (float)(EndTime-StartTime)/1000;
            StartTime = System.currentTimeMillis();
            SHAPO = SHAPO();
            EndTime = System.currentTimeMillis();
            sumSHAPOTime += (float)(EndTime-StartTime)/1000;
            StartTime = System.currentTimeMillis();
            Short = Shortcut(totalCost);
            EndTime = System.currentTimeMillis();
            sumShortTime += (float)(EndTime-StartTime)/1000;
            StartTime = System.currentTimeMillis();
            depot = depot(totalCost);
            EndTime = System.currentTimeMillis();
            sumDepotTime += (float)(EndTime-StartTime)/1000;
            StartTime = System.currentTimeMillis();
            ReRout(ReRout, totalCost);
            EndTime = System.currentTimeMillis();
            sumReRoutTime += (float)(EndTime-StartTime)/1000;
            StartTime = System.currentTimeMillis();
            Appro1 = Appro1(totalCost);
            EndTime = System.currentTimeMillis();
            sumAppro1Time += (float)(EndTime-StartTime)/1000;
            StartTime = System.currentTimeMillis();
            Appro2 = Appro2(totalCost);
            EndTime = System.currentTimeMillis();
            sumAppro2Time += (float)(EndTime-StartTime)/1000;

            outputValues(i+1, inputIndexes, initOrder, shapley, SHAPO, Short, depot, ReRout, Appro1, Appro2);

            double[] forIterationBySHAPO = new double[2];// {MaxABS, RMSE}
            double[] forIterationShort = new double[2];// {MaxABS, RMSE}
            double[] forIterationDepot = new double[2];// {MaxABS, RMSE}
            double[] forIterationReRout = new double[2];// {MaxABS, RMSE}
            double[] forIterationAppro1 = new double[2];// {MaxABS, RMSE}
            double[] forIterationAppro2 = new double[2];// {MaxABS, RMSE}

            for(int j=0; j<num; j++)
            {
                sumShapley += shapley[j];
                calc(shapley_SHAPO, shapley[j], SHAPO[j]);
                calc(shapley_Short, shapley[j], Short[j]);
                calc(shapley_depot, shapley[j], depot[j]);
                calc(shapley_ReRute, shapley[j], ReRout[j]);
                calc(shapley_Appro1, shapley[j], Appro1[j]);
                calc(shapley_Appro2, shapley[j], Appro2[j]);

                calcForIteration(forIterationBySHAPO, shapley[j], SHAPO[j]);
                calcForIteration(forIterationShort, shapley[j], Short[j]);
                calcForIteration(forIterationDepot, shapley[j], depot[j]);
                calcForIteration(forIterationReRout, shapley[j], ReRout[j]);
                calcForIteration(forIterationAppro1, shapley[j], Appro1[j]);
                calcForIteration(forIterationAppro2, shapley[j], Appro2[j]);
            }

            shapley_SHAPO[3] += forIterationBySHAPO[0];
            shapley_Short[3] += forIterationShort[0];
            shapley_depot[3] += forIterationDepot[0];
            shapley_ReRute[3] += forIterationReRout[0];
            shapley_Appro1[3] += forIterationAppro1[0];
            shapley_Appro2[3] += forIterationAppro2[0];

            shapley_SHAPO[4] += Math.sqrt(forIterationBySHAPO[1]/num);
            shapley_Short[4] += Math.sqrt(forIterationShort[1]/num);
            shapley_depot[4] += Math.sqrt(forIterationDepot[1]/num);
            shapley_ReRute[4] += Math.sqrt(forIterationReRout[1]/num);
            shapley_Appro1[4] += Math.sqrt(forIterationAppro1[1]/num);
            shapley_Appro2[4] += Math.sqrt(forIterationAppro2[1]/num);

            output(i+1, sumShapley, sumShapleyTime, shapley_SHAPO, shapley_Short, shapley_depot, shapley_ReRute, shapley_Appro1, shapley_Appro2);
            outputTime(i+1, sumShapleyTime, sumShortTime, sumDepotTime,sumSHAPOTime, sumAppro1Time, sumAppro2Time);
        }
    }

    protected int[] randomInput(){
        double[][] cityDistances;
        if(ToulouseGraph)
            cityDistances = ShapleyValue.ToulouseDistances;
        else cityDistances = ShapleyValue.NewYorkDistances;
        int[] arr = new int [num+1];
        arr[0] = this.originIndex;
        Random rand = new Random();
        for(int i=1; i<num+1; i++){
            do {
                arr[i] = rand.nextInt(cityDistances.length);
            }while (arr[i]== this.originIndex);
        }
        for(int i=0; i<num+1; i++)
            for(int j=0; j<num+1; j++)
                Distances[i][j] = cityDistances[arr[i]][arr[j]];
        return arr;
    }

    protected int[] setInitOrder(){
        long StartTime, EndTime;
        StartTime = System.currentTimeMillis();
        int[] opt2 = Opt2();
        double opt2Cost = costOfByOrder(opt2);
        EndTime = System.currentTimeMillis();
        totalCost = opt2Cost;
        reorderMatrix(opt2);
        return opt2;
    }

    protected void reorderMatrix(int [] initOrder){
        double newMatrix[][] = new double[num+1][num+1];
        for(int i=0; i<num+1; i++)
            for(int j=0; j<num+1; j++)
                newMatrix[i][j] = Distances[initOrder[i]][initOrder[j]];
        for(int i=0; i<num+1; i++)
            for(int j=0; j<num+1; j++)
                Distances[i][j] = newMatrix[i][j];
    }

    protected int[] Opt2(){
        int opt2[] = new int[Distances.length];
        for(int i=0; i<opt2.length; i++)
            opt2[i]=i;
        int count =0;
        boolean flag = true;
        while(flag && count < Math.pow(num,2.)){
            flag = false;
            count++;
            outer:
            for(int i=0; i<opt2.length; i++){
                for(int j=i+2; j<opt2.length; j++){
                    int temp[] = new int[Distances.length];
                    for(int k=0; k<temp.length; k++)
                        if(k>i && k<=j)
                            temp[k] = opt2[j-(k-i-1)];
                        else temp[k] = opt2[k];
                    if(costOfByOrder(temp)<costOfByOrder(opt2)){
                        for (int k = 0; k < temp.length; k++)
                            opt2[k] = temp[k];
                        flag=true;
                        break outer;
                    }
                }
            }
        }
        while(opt2[0]!=0){
            int temp = opt2[0];
            for(int i=0; i<opt2.length-1; i++)
                opt2[i] = opt2[i+1];
            opt2[opt2.length-1] = temp;
        }
        return opt2;
    }

    protected int[] bestOrder (List<Integer> ind) {
        double[][] dist = new double[ind.size()][ind.size()];
        for (int i=0; i<ind.size(); i++){
            for(int j=0; j<ind.size(); j++)
                dist[i][j] = Distances[ind.get(i)][ind.get(j)];
        }
        Solver s = solv.newSolver(dist);
        s.calculate();
        return Arrays.stream(s.best_path).map(i->ind.get(i)).toArray();
    }

    protected double costOf(List<Integer> ind){
        int[] best = bestOrder(ind);
        return costOfByOrder(best);
    }

    protected double[] Shapley(double[] ReRout) {
        double[] p = new double[num];
        List<Integer> ind = new LinkedList();
        for (int i = 1; i <= num; i++)
            ind.add(i);
        List<List<Integer>> AllCoalitions = new LinkedList(allCoalitions(ind));
        List<Double> AllCosts = new LinkedList();
        for (int i = 0; i < AllCoalitions.size(); i++) {
            long StartTime, EndTime;
            List<Integer> temp = (AllCoalitions.get(i));
            StartTime = System.currentTimeMillis();
            temp.add(0,0);
            double cost = costOf(temp);
            EndTime = System.currentTimeMillis();
            if(temp.size()==num+1)
                cost=totalCost;
            AllCosts.add(cost);
            temp.remove(0);
            if(temp.size()==num-1){
                for(int pass=1; pass<=num; pass++)
                    if(!temp.contains(pass))
                        ReRout[pass-1] = cost;
                sumReRoutTime += (float)(EndTime-StartTime)/1000;
            }
        }
        double[] costsBySizes = new double[num];
        int[] countBySizes = new int[num];
        for (int i=0; i<num; i++){
            Arrays.fill(costsBySizes, 0);
            Arrays.fill(countBySizes, 0);
            for(int j=0; j<AllCoalitions.size(); j++){
                List<Integer> temp = new LinkedList (AllCoalitions.get(j));
                if(temp.contains(i+1)){
                    costsBySizes[temp.size()-1]+=AllCosts.get(j);
                    countBySizes[temp.size()-1]++;
                    if (temp.size() != 1) {
                        temp.remove(temp.indexOf(i+1));
                        costsBySizes[temp.size()]-= AllCosts.get(AllCoalitions.indexOf(temp));
                    }
                }
            }
            for(int t=0; t<costsBySizes.length; t++)
                p[i]+= (costsBySizes[t] / countBySizes[t]) / p.length;
        }
        return p;
    }

    protected List<List<Integer>> allCoalitions (List<Integer> ind){
        List<List<Integer>> AllCoalitions = new LinkedList();
        if(ind.size() == 1){
            AllCoalitions.add(new LinkedList<>(ind));
            return AllCoalitions;
        }
        AllCoalitions.add(ind);
        List<List<Integer>> MiniAllCoalitions = new LinkedList();
        for (int i=0; i<ind.size(); i++){
            List<Integer> temp = new LinkedList<>(ind);
            temp.remove(temp.get(i));
            MiniAllCoalitions = allCoalitions(temp);
            for (int j=0; j<MiniAllCoalitions.size(); j++) {
                if (!AllCoalitions.contains(MiniAllCoalitions.get(j)))
                    AllCoalitions.add(MiniAllCoalitions.get(j));
            }
        }
        return AllCoalitions;
    }

    protected abstract double costOfByOrder(int[] ind);
    protected abstract double[] SHAPO(); //our approximation.
    protected abstract double[] Shortcut(double totalCost);

    protected double[] depot(double totalCost){
        double[] p = new double[num];
        double sumPrivate = 0;
        for(int i=1; i<=num; i++)
            sumPrivate += Distances[0][i];
        for(int i=1; i<=num; i++)
            p[i-1] = (Distances[0][i]/sumPrivate)* totalCost;
        return p;
    }

    protected void ReRout(double[] reRout, double totalCost){
        double sum = Arrays.stream(reRout).sum();
        for(int i=0; i<num; i++)
            reRout[i] = (reRout[i]/sum)* totalCost;
    }

    protected double[] Appro1(double totalCost){
        double[] p = new double[num];
        Double[] SC = new Double[num-1];

        for(int i=1; i<=num; i++) {
            int ind = 0;
            for(int j=1; j<=num; j++){
                if (j==i) continue;
                SC[ind] = s(j,i);
                ind++;
            }
            p[i-1] = s(i,i) - SC(SC);
        }

        double sum = Arrays.stream(p).sum();
        for(int i=0; i<num; i++)
            p[i] = (p[i]/sum)* totalCost;
        return p;
    }

    private double s (int i, int j){
        return Distances[0][i] + Distances[0][j] - Distances[i][j];
    }

    private double SC(Double[] sc){
        double res = 0;
        Arrays.sort(sc, Collections.reverseOrder());
        for(int i=1; i<=sc.length; i++)
            res += sc[i-1]/(i*(i+1));
        return res;
    }

    private double PartSum(Double[] v){
        int m = v.length;
        Arrays.sort(v);
        double s=0;
        for(int i=1; i<m; i++)
            s += (v[m-i]-v[m-i-1])/i;
        s += v[0]/m;
        return s;
    }

    protected double[] Appro2(double totalCost){
        double[] CfVals = new double[num];
        int n = Distances[0].length;
        Double[][] DS = new Double[n-1][n-1];
        for(int i=0; i<n-1; i++)
            for(int j=0; j<n-1; j++)
                DS[i][j] = s(i+1,j+1);
        for (int i=0; i<n-1; i++){
            Double[] w = Arrays.copyOf(DS[i], DS[i].length);;
            double vali = PartSum(w);
            double [][] DM1 = new double[n-2][n-2];
            int ctneg = 0;
            for(int j=0; j<n-1; j++){
                if (j==i) continue;
                for(int k=j+1; k<n-1; k++) {
                    if (k == i) continue;
                    double d1 = Double.min(DS[i][j], DS[i][k]);
                    double d1M = Double.max(DS[i][j], DS[i][k]);
                    double d2 = DS[j][k];

                    int j1 = (j<i) ? j : j-1;
                    int k1 = (k<i) ? k : k-1;

                    DM1[k1][j1] = DM1[j1][k1] = d1;

                    if (d2 < d1){
                        ctneg ++;
                        DM1[k1][j1] = DM1[j1][k1] = d2;
                    }
                }
            }
            double[] col = new double[n-2];
            int ind = 0;
            for(int l=0; l<n-1; l++) {
                if (l == i) continue;
                col[ind] = DS[l][i];
                ind++;
            }

            int[] S = IntStream.range(0, col.length)
                    .boxed().sorted((a, b) -> col[a] == col[b] ? 0 : col[a]<col[b]? -1 :1)
                    .mapToInt(ele -> ele).toArray();
            Arrays.sort(col);

            double[][] copyDM1 = new double[DM1.length][DM1[0].length];
            for(int ind1 = 0; ind1<S.length; ind1++)
                for(int ind2 = 0; ind2<S.length; ind2++)
                    copyDM1[ind1][ind2] = DM1[S[ind1]][S[ind2]];
            double corr = ShValTFO2(col, copyDM1);
            CfVals[i] = DS[i][i] - corr;
        }

        double sum = Arrays.stream(CfVals).sum();
        for(int i=0; i<num; i++)
            CfVals[i] = (CfVals[i]/sum)* totalCost;
        return CfVals;
    }

    protected double ShValTFO2(double[] v,double[][]tab) {
        int n = v.length;
        double[] str = new double[n];
        double[][] tabcounts = new double[n][n];
        double eps = 0.000001;
        double s = 0.5 * DoubleStream.of(v).sum();
        for (int i = 0; i < n - 1; i++){
            for (int j = i + 1; j < n; j++) {
                double val = tab[i][j];
                int k=0;
                for(int l=0; l<n; l++){
                    if(l==i || l==j)
                        continue;
                    if ((tab[i][l] > val+eps || (Math.abs(tab[i][l]-val)<eps && l>j)) && (tab[l][j] > val+eps || (Math.abs(tab[l][j]-val)<eps && l>i)))
                        k++;
                }
                tabcounts[i][j] = tabcounts[j][i] = k;
                str[k+1]++;
            }
        }
        int n2 = n*(n-1)/2;
        Double [][] new_v = new Double[n2][2];
        int ct = 0;
        for(int i=0; i<n-1; i++)
            for(int j=i+1; j<n; j++){
                new_v[ct][0] = tabcounts[i][j];
                new_v[ct][1] = tab[i][j];
                ct++;
            }
        List<List<Double>> vs = Stream.of(new_v).map(a-> Arrays.asList(a)).collect(Collectors.toList());
        Collections.sort(vs,new Comparator<List<Double>>(){
            @Override
            public int compare(final List<Double> lhs,List<Double> rhs) {
                if(-lhs.get(0) > -rhs.get(0)) return 1;
                if(-lhs.get(0) < -rhs.get(0)) return -1;
                if(lhs.get(1) > rhs.get(1)) return 1;
                if(lhs.get(1) < rhs.get(1)) return -1;
                return 0;
            }
        });
        int ck = n-2;
        int cct = 1;
        for(int i=0; i<n2; i++){
            double val = vs.get(i).get(1);
            s -= 2*val/((ck+3)*(ck+2)*(ck+1));
            cct--;
            if (cct == 0) {
                ck--;
                cct = n-1-ck;
            }
        }
        return s;
    }

    protected void calc(double[] cal, double shapley, double n)
    {
        if (shapley == 0)
            System.out.println("problem");
        else
        {
            cal[0] += Math.abs(n - shapley) / shapley;
            cal[1] += Math.abs(n - shapley);
            cal[2] += Math.pow(n - shapley, 2.);
        }
    }

    protected void calcForIteration(double[] cal, double shapley, double n){
        if(Math.abs(n-shapley)>cal[0])
            cal[0] = Math.abs(n-shapley);
        cal[1] += Math.pow(n-shapley,2.);
    }

    protected void output(int i, double sumShapley, float sumShapleyTime, double[] shapley_SHAPO, double[] shapley_Short, double[] shapley_depot, double[] shapley_ReRout, double[] shapley_Appro1, double[] shapley_Appro2){
        try {
            String fileName;
            if(ToulouseGraph)
                fileName = "Toulouse_";
            else fileName = "NewYork_";
            fileName += this.originIndex + "_1_" + index + "_" + Integer.toString(num) + "Passengers.txt";
            File file = new File(fileName);
            PrintWriter writer = new PrintWriter(file, "UTF-8");
            int count = num * i;
            writer.println("iterations - " + i);
            writer.println("count - " + count);
            writer.println("sumShapley - " + sumShapley);
            writer.println("sumShapleyTime - " + sumShapleyTime);

            writer.println();
            writer.println("SHAPO - ");
            writer.println();
            writer.println("SumPrecent - " + shapley_SHAPO[0]);
            writer.println("SumABS - " + shapley_SHAPO[1]);
            writer.println("SumPow - " + shapley_SHAPO[2]);
            writer.println("SumMaxABS - " + shapley_SHAPO[3]);
            writer.println("RMSE - " + shapley_SHAPO[4]);

            writer.println();
            writer.println("Shortcut - ");
            writer.println();
            writer.println("SumPrecent - " + shapley_Short[0]);
            writer.println("SumABS - " + shapley_Short[1]);
            writer.println("SumPow - " + shapley_Short[2]);
            writer.println("SumMaxABS - " + shapley_Short[3]);
            writer.println("RMSE - " + shapley_Short[4]);

            writer.println();
            writer.println("depot - ");
            writer.println();
            writer.println("SumPrecent - " + shapley_depot[0]);
            writer.println("SumABS - " + shapley_depot[1]);
            writer.println("SumPow - " + shapley_depot[2]);
            writer.println("SumMaxABS - " + shapley_depot[3]);
            writer.println("RMSE - " + shapley_depot[4]);

            writer.println();
            writer.println("ReRout - ");
            writer.println();
            writer.println("SumPrecent - " + shapley_ReRout[0]);
            writer.println("SumABS - " + shapley_ReRout[1]);
            writer.println("SumPow - " + shapley_ReRout[2]);
            writer.println("SumMaxABS - " + shapley_ReRout[3]);
            writer.println("RMSE - " + shapley_ReRout[4]);

            writer.println();
            writer.println("Appro1 - ");
            writer.println();
            writer.println("SumPrecent - " + shapley_Appro1[0]);
            writer.println("SumABS - " + shapley_Appro1[1]);
            writer.println("SumPow - " + shapley_Appro1[2]);
            writer.println("SumMaxABS - " + shapley_Appro1[3]);
            writer.println("RMSE - " + shapley_Appro1[4]);

            writer.println();
            writer.println("Appro2 - ");
            writer.println();
            writer.println("SumPrecent - " + shapley_Appro2[0]);
            writer.println("SumABS - " + shapley_Appro2[1]);
            writer.println("SumPow - " + shapley_Appro2[2]);
            writer.println("SumMaxABS - " + shapley_Appro2[3]);
            writer.println("RMSE - " + shapley_Appro2[4]);

            writer.close();
        }
        catch (IOException e) {}
    }

    protected void outputTime(int i, float sumShapleyTime, float sumShortTime, float sumDepotTime, float sumSHAPOTime, float sumAppro1Time, float sumAppro2Time){
        try {
            String fileName;
            if(ToulouseGraph)
                fileName = "Toulouse_";
            else fileName = "NewYork_";
            fileName += this.originIndex + "_1_" + index + "_Time_" + Integer.toString(num) + "Passengers.txt";
            File file = new File(fileName);
            PrintWriter writer = new PrintWriter(file, "UTF-8");
            int count = num * i;
            writer.println("count - " + count);
            writer.println("SumShapleyTime - " + sumShapleyTime);
            writer.println("SumShortTime - " + sumShortTime);
            writer.println("SumDepotTime - " + sumDepotTime);
            writer.println("SumSHAPOTime - " + sumSHAPOTime);
            writer.println("SumReRoutTime - " + sumReRoutTime);
            writer.println("SumAppro1Time - " + sumAppro1Time);
            writer.println("SumAppro2Time - " + sumAppro2Time);
            writer.println("SumTSPTime - " + sumTSPTime);
        }
        catch (IOException e) {}
    }

    protected void outputValues(int i, int[] inputIndex, int[] initOrder, double[] Shapley, double[] SHAPO, double[] Short, double[] Depot, double[] ReRout, double[] Appro1, double[] Appro2) {
        try {
            String fileName = "Values_";
            if(ToulouseGraph)
                fileName += "Toulouse_";
            else fileName += "NewYork_";
            fileName += this.originIndex + "_" + index + "_" + Integer.toString(num) + "Passengers.txt";
            FileWriter file = new FileWriter(fileName, true);
            PrintWriter writer = new PrintWriter(file);
            writer.println("Location Indexes - " + Arrays.toString(inputIndex));
            writer.println("Initial Order - " + Arrays.toString(initOrder));
            writer.println("Shapley Values - " + Arrays.toString(Shapley));
            writer.println("SHAPO Values - " + Arrays.toString(SHAPO));
            writer.println("ShortCut Values - " + Arrays.toString(Short));
            writer.println("Depot Values - " + Arrays.toString(Depot));
            writer.println("Rerout Values - " + Arrays.toString(ReRout));
            writer.println("Appro1 Values - " + Arrays.toString(Appro1));
            writer.println("Appro2 Values - " + Arrays.toString(Appro2));
            writer.println();
            writer.println("-----------------------------------------------------------");
            writer.println();
            writer.close();
        } catch (IOException e) {}
    }
}