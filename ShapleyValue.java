import java.io.*;
import java.util.Arrays;

public class ShapleyValue {

    static double[][] ToulouseDistances;
    static double[][] NewYorkDistances;
    static int[] ToulouseOrigins;
    static int[] newYorkOrigins;
    public static int min = 3;
    public static int max = 14;
    private static boolean RoutingGames = true;

    public static void main(String [] args) {

        System.out.println("Reading Toulouse matrix...");
        ToulouseMatrix();
        System.out.println("Reading NewYork matrix...");
        NewYorkMatrix();
        origins();

        System.out.println("starting test Shapley...");
        int numOfOrigins = 11;
        int threads_foreach_origin = 1;
        int len = 2*numOfOrigins*threads_foreach_origin; //origins*graphs*threads_foreach_origin
        myThread[] Arr = new myThread[len];
        for (int n = min; n <= max; n++) {
            boolean Toulouse = true;
            int ind = 0;
            while (ind < len) {
                for (int j = 0; j < ToulouseOrigins.length; j++) {
                    for (int i = 0; i < threads_foreach_origin; i++) {
                        if (Toulouse) {
                            if (RoutingGames)
                                Arr[ind] = new RoutingGamesThread("ToulouseThread _" + ToulouseOrigins[j] + "_" + Integer.toString(n), true, n, ToulouseOrigins[j], i);
                            else
                                Arr[ind] = new RideSharingThread("ToulouseThread _" + ToulouseOrigins[j] + "_" + Integer.toString(n), true, n, ToulouseOrigins[j], i);
                        } else {
                            if (RoutingGames)
                                Arr[ind] = new RoutingGamesThread("NewYorkThread_" + newYorkOrigins[j] + "_" + Integer.toString(n), false, n, newYorkOrigins[j], i);
                            else
                                Arr[ind] = new RideSharingThread("NewYorkThread_" + newYorkOrigins[j] + "_" + Integer.toString(n), false, n, newYorkOrigins[j], i);
                        }
                        Arr[ind].start();
                        ind++;
                    }
                }
                Toulouse = false;
            }
            boolean flag;
            do{
                flag = true;
                for(myThread t : Arr)
                    if(t.isAlive())
                        flag = false;
            }
            while(!flag);
        }
    }

    private static void ToulouseMatrix(){
        int ToulouseSize = 40000;
        ToulouseDistances = new double[ToulouseSize][ToulouseSize];
        try {
            String fileName = new String("/home/azariaa/DNN/Chaya/ToulouseFloydWarshall_"+Integer.toString(ToulouseSize)+".txt");//should be "newToulouseFloyedWarshal
            BufferedReader reader = new BufferedReader(new FileReader(fileName));
            for(int i=0;i<ToulouseSize; i++) {
                String line = reader.readLine();
                String[] SL = line.split(", ");
                for(int j=0; j<ToulouseSize; j++)
                    ToulouseDistances[i][j] = Double.parseDouble(SL[j]);
            }
        }
        catch (IOException e) {}
    }
    private static void NewYorkMatrix(){
        int NewYorkSize = 8001;
        NewYorkDistances = new double[NewYorkSize][NewYorkSize];
        try {
            String fileName = new String("/home/azariaa/DNN/Chaya/NewYorkFloydWarshall.txt");
            BufferedReader reader = new BufferedReader(new FileReader(fileName));
            for(int i=0;i<NewYorkSize; i++) {
                String line = reader.readLine();
                String[] SL = line.split(", ");
                for(int j=0; j<NewYorkSize; j++)
                    NewYorkDistances[i][j] = Double.parseDouble(SL[j]);
            }
        }
        catch (IOException e) {}
    }
    private static void origins(){
        int ToulouseAirportIndex = 13994;
        int ToulouseRandomOrigin1 = 774;
        int ToulouseRandomOrigin2 = 39930;
        int ToulouseRandomOrigin3 = 19768;
        int ToulouseRandomOrigin4 = 5109;
        int ToulouseRandomOrigin5 = 4982;
        int ToulouseRandomOrigin6 = 20449;
        int ToulouseRandomOrigin7 = 19248;
        int ToulouseRandomOrigin8 = 16420;
        int ToulouseRandomOrigin9 = 877;
        int ToulouseRandomOrigin10 = 16249;


        int NewYorkTrainStation = 2878;
        int NewYorkRandomOrigin1 = 3096;
        int NewYorkRandomOrigin2 = 1050;
        int NewYorkRandomOrigin3 = 2932;
        int NewYorkRandomOrigin4 = 1151;
        int NewYorkRandomOrigin5 = 5883;
        int NewYorkRandomOrigin6 = 1816;
        int NewYorkRandomOrigin7 = 7636;
        int NewYorkRandomOrigin8 = 675;
        int NewYorkRandomOrigin9 = 3671;
        int NewYorkRandomOrigin10 = 1494;

        ToulouseOrigins = new int[]{ToulouseAirportIndex,
                ToulouseRandomOrigin1,
                ToulouseRandomOrigin2,
                ToulouseRandomOrigin3,
                ToulouseRandomOrigin4,
                ToulouseRandomOrigin5,
                ToulouseRandomOrigin6,
                ToulouseRandomOrigin7,
                ToulouseRandomOrigin8,
                ToulouseRandomOrigin9,
                ToulouseRandomOrigin10};

        newYorkOrigins = new int[]{NewYorkTrainStation,
                NewYorkRandomOrigin1,
                NewYorkRandomOrigin2,
                NewYorkRandomOrigin3,
                NewYorkRandomOrigin4,
                NewYorkRandomOrigin5,
                NewYorkRandomOrigin6,
                NewYorkRandomOrigin7,
                NewYorkRandomOrigin8,
                NewYorkRandomOrigin9,
                NewYorkRandomOrigin10};
    }
}