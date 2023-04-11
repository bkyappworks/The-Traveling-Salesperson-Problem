import java.io.File;
import java.io.FileNotFoundException;
import java.util.LinkedList;
import java.util.Scanner;

public class TSPoptimal {
    // build graph -> find mst with prim -> solve tsp
    public static Vertex[] crimeArray;
    public static double[][] distanceArray;
    public static void buildGraph(String crimeDataLocation, String start, String end) {
        StringBuilder countSize = new StringBuilder();
        // In Part 1, we will be using the (X, Y) pairs to compute the distance between each crime.
        File file = new File(crimeDataLocation);
        Scanner scanner = null;
        int index = 0;
        try {
            scanner = new Scanner(file);
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }
        scanner.nextLine();
        while (scanner.hasNext()) {
            String line = scanner.nextLine();
            // filter by date
            String date = line.split(",")[5];
            int month = Integer.parseInt(date.split("/")[0]);
            int day = Integer.parseInt(date.split("/")[1]);
            int year = Integer.parseInt(date.split("/")[2]);
            if (
                    year >= Integer.parseInt(start.split("/")[2]) &&
                            month >= Integer.parseInt(start.split("/")[0]) &&
                            day >= Integer.parseInt(start.split("/")[1]) &&
                            year <= Integer.parseInt(end.split("/")[2]) &&
                            month <= Integer.parseInt(end.split("/")[0]) &&
                            day <= Integer.parseInt(end.split("/")[1])

            ) {
                System.out.println(line);
                countSize.append(line).append("!");
            }
        }
        // The graph will be implemented as a one dimensional array of crime records
        String[] data = countSize.toString().split("!");
        crimeArray = new Vertex[data.length];
        for (int i = 0; i < data.length; i++) {
            double x = Double.parseDouble(data[i].split(",")[0]);
            double y = Double.parseDouble(data[i].split(",")[1]);
            Vertex v = new Vertex(index,x,y);
            crimeArray[i] = v;
            index++;
        }

        // along with a two-dimensional array of doubles – holding the computed distance between each pair
        distanceArray = new double[data.length][data.length];
        for (int i = 0; i < crimeArray.length; i++) {
            for (int j = 0; j < crimeArray.length; j++) {
                distanceArray[i][j] = calculateDistance(i,j);
            }
        }
    }
    private static double calculateDistance(int i, int j) {
        Vertex p1 = crimeArray[i];
        Vertex p2 = crimeArray[j];
        double distance = Math.pow((Math.pow(p1.getX()-p2.getX(),2) + Math.pow(p1.getY()-p2.getY(),2)),0.5);
        return distance;
    }
    private static double calculateDistance2(Vertex i, Vertex j) {
        Vertex p1 = i;
        Vertex p2 = j;
        double distance = Math.pow((Math.pow(p1.getX()-p2.getX(),2) + Math.pow(p1.getY()-p2.getY(),2)),0.5);
        return distance;
    }
    public static LinkedList<Vertex>[] primMST(Vertex root, Vertex[] crimeArray, double[][] distanceArray) {
        Vertex[] parentArray = new Vertex[crimeArray.length];
        Vertex[] queueArray = new Vertex[crimeArray.length];
        // populate keyArray
        double[] keyArray = new double[crimeArray.length];
        // fill with infinity
        double inf = Double.POSITIVE_INFINITY;
        for (int i = 0; i < keyArray.length; i++) {
            keyArray[i] = inf;
        }
        // populate queueArray
        for (int i = 0; i < crimeArray.length; i++) {
            addMinHeap(crimeArray[i],queueArray);
        }
        fixMinHeap(queueArray,keyArray);
        keyArray[root.getIndex()] = 0;
        parentArray[root.getIndex()] = null;
        while (queueArray.length != 0) {
            Vertex u = queueArray[0]; // u = ExtractMin(Q)
            queueArray = popMinHeap(queueArray);
            fixMinHeap(queueArray,keyArray);
            int uIndex = u.getIndex();
            StringBuilder adjacencySB = new StringBuilder();
            for (int i = 0; i < distanceArray[uIndex].length; i++) {
                if (distanceArray[uIndex][i] != 0) {
                    adjacencySB.append(i).append("!");
                }
            }
            String [] adjacencyIndex = adjacencySB.toString().split("!");
            Vertex[] adjacency = new Vertex[adjacencyIndex.length];
            for (int i = 0; i < adjacencyIndex.length; i++) {
                adjacency[i] = crimeArray[Integer.parseInt(adjacencyIndex[i])];
            }

            for (Vertex v: adjacency) { // foreach v in Adj[u]
                double distanceUV = calculateDistance2(u,v);
                if (contains(v,queueArray) && distanceUV < keyArray[v.getIndex()]) { // if v in Q && w(u,v) < key[v]
                    // then pi[v] = u
                    parentArray[v.getIndex()] = u;
                    // key[v] = w(u,v)
                    keyArray[v.getIndex()] = distanceUV;
                }
            }
        }
        LinkedList[] mst = new LinkedList[parentArray.length + 1];
        for (int i = 0; i < mst.length; i++) {
            mst[i] = new LinkedList<>();
        }
        for (int i = 0; i < parentArray.length; i++) {
            if (parentArray[i] == null) {
                // save root at the last index of mst array
                mst[mst.length-1].add(crimeArray[i]);
            } else {
                mst[parentArray[i].getIndex()].add(crimeArray[i]);
            }

        }
        return mst;
    }
    public static LinkedList<Vertex> ApproxTSPTour(Vertex root, Vertex[] crimeArray, double[][] distanceArray) {
        // 1. Select a vertex r e V[G] to be a root vertex
        // 2. Compute a minimum spanning tree T for G from root r using MST-Prim(G,c,r)
        LinkedList<Vertex>[] mst = primMST(root,crimeArray, distanceArray);
        // 3. Let L be the list of vertices visited in a preorder tree walk of T
        LinkedList<Vertex> L = new LinkedList<>();
        L.add(root);
        for (int i = 0; i < mst.length; i++) {
            preOrderWalk(root,mst,mst[i],L); // mst[i]:LL
        }
        // 4. Return the Hamiltonian cycle H that visits the vertices in the order L
        LinkedList<Vertex> H = L;
        for (int i = 0; i < H.size(); i++) {
            System.out.print(H.get(i).getIndex()+" ");
        }
        return H;
    }
    public static double cycleLength(Vertex root) {
        LinkedList<Vertex>[] mst = primMST(root,crimeArray, distanceArray);
        LinkedList<Vertex> L = new LinkedList<>();
        L.add(root);
        for (int i = 0; i < mst.length; i++) {
            preOrderWalk(root,mst,mst[i],L); // mst[i]:LL
        }
        LinkedList<Vertex> H = L;
        double length = 0;
        for (int i = 0; i < H.size() - 1; i++) {
            length += distanceArray[H.get(i).getIndex()][H.get(i+1).getIndex()];
        }
        double convertTomiles = 0.00018939;
        return length*convertTomiles;
    }
    public static double cycleLength2(LinkedList<Double> cycle) {
        cycle.set(cycle.size() - 1,cycle.get(0));
        double length = 0;
        for (int i = 0; i < cycle.size() - 1; i++) {
            double idx = cycle.get(i);
            int id = (int)idx;
            double next = cycle.get(i+1);
            int nextid = (int)next;
            length += distanceArray[id][nextid];
        }
        double convertTomiles = 0.00018939;
        return length*convertTomiles;
    }
    private static void preOrderWalk(Vertex root, LinkedList[] mst, LinkedList<Vertex> ll,LinkedList<Vertex> L) {
        if (ll.size() == 0 || ll == null) {
            return;
        }
        // walk through each LL
        for (Vertex v: ll) {
            if (!contains2(v,L) || v.equals(root)) {
                L.add(v);
            }
            preOrderWalk(root,mst,mst[v.getIndex()],L);
        }
    }
    private static boolean contains(Vertex v, Vertex[] queueArray) {
        for (Vertex vertex: queueArray) {
            if (vertex.equals(v)) {
                return true;
            }
        }
        return false;
    }
    private static boolean contains2(Vertex v, LinkedList<Vertex> ll) {
        for (int i = 0; i < ll.size(); i++) {
            if (ll.get(i).equals(v)) {
                return true;
            }
        }
        return false;
    }
    private static void fixMinHeap(Vertex[] queueArray, double[] keyArray) {
//        System.out.println("queueArray.length: "+queueArray.length);
        if (queueArray.length <= 2) {
            return;
        }
        // identify the last parent node -> compare with left & right child -> swap if needed
        int lastParentNodeIdx = (queueArray.length - 2)/2;
        for (int parent = lastParentNodeIdx; parent >= 0; parent--) {
            int leftChildIdx = 2 * parent + 1; // left child is stored at index 2k + 1
            int rightChildIdx = 2 * parent + 2;// right child at index 2k + 2
            double leftChild = keyArray[queueArray[leftChildIdx].getIndex()];
            if (rightChildIdx == queueArray.length) {
                // swap with left child
                Vertex tmp = queueArray[parent];
                queueArray[parent] = queueArray[leftChildIdx];
                queueArray[leftChildIdx] = tmp;
            } else {
                double rightChild = keyArray[queueArray[rightChildIdx].getIndex()];
                if (rightChild > leftChild) {
                    // swap with left child
                    Vertex tmp = queueArray[parent];
                    queueArray[parent] = queueArray[leftChildIdx];
                    queueArray[leftChildIdx] = tmp;
                }
                if (leftChild > rightChild) {
                    // swap with left child
                    Vertex tmp = queueArray[parent];
                    queueArray[parent] = queueArray[rightChildIdx];
                    queueArray[rightChildIdx] = tmp;
                }
            }
        }
    }
    private static void addMinHeap(Vertex vertex,Vertex[] queueArray) {
        int index = 0;
        for (int i = 0; i < queueArray.length; i++) {
            if (queueArray[i] == null) {
                index = i;
                break;
            }
        }
        queueArray[index] = vertex;
    }
    private static Vertex[] popMinHeap(Vertex[] queueArray) {
        Vertex[] newArray = new Vertex[queueArray.length - 1];
        for (int i = 1; i < queueArray.length; i++) {
            newArray[i-1] = queueArray[i];
        }
        return newArray;
    }
    // https://stackoverflow.com/questions/2920315/permutation-of-array
    public static LinkedList<Double> permute(Vertex root) {
        LinkedList<Double> winner = new LinkedList<>();
        //populate winner with a size of crimeArray.length + 1, last idx saves cycle length
        for (int i = 0; i <= crimeArray.length; i++) {
            winner.add(Double.POSITIVE_INFINITY); // size: 5
        }
        permute(crimeArray,0, root, winner);
        return winner;
    }
    public static void permute(Vertex[] arr, int k, Vertex root, LinkedList<Double> winner) {
        LinkedList<Double> cycle = new LinkedList<>();
        for (int i = k; i < arr.length; i++){
            swap(arr, i, k);
            permute(arr, k+1, root, winner);
            swap(arr, k, i);
        }
        if (k == arr.length -1 && arr[0].getIndex() == root.getIndex()){
            // populate cycle
            for (int i = 0; i <= arr.length; i++) {
                cycle.add(Double.NEGATIVE_INFINITY);
            }
            // update cycle idx 0 - arr.length - 1
            for (int i = 0; i < arr.length; i++) {
                cycle.set(i,Double.parseDouble(String.valueOf(arr[i].getIndex())));
            }
            // update cycle last idx
            cycle.set(cycle.size() - 1, cycleLength2(cycle));
            // update winner when there's a smaller cycle length
            if (cycle.get(cycle.size() - 1) < winner.get(winner.size() - 1)) {
                for (int i = 0; i < cycle.size(); i++) {
                    winner.set(i,cycle.get(i));
                }
            }
        }
    }
    private static void swap(Vertex[] arr, int i, int k) {
        Vertex tmp = arr[i];
        arr[i] = arr[k];
        arr[k] = tmp;
    }
    public static void main(String[] args) {
        Scanner scanner = new Scanner(System.in);
        // prompt the user for two dates
        System.out.println("Enter Start Date: ");
        String start = scanner.nextLine();
        System.out.println("Enter End Date: ");
        String end = scanner.nextLine();
        System.out.println();
        System.out.println("Crime records between "+start+" and "+end);
        buildGraph("CrimeLatLonXY1990.csv",start,end);
        System.out.println();
        System.out.println("Hamiltonian Cycle (not necessarily optimum):");
        Vertex root = crimeArray[0];
        ApproxTSPTour(root,crimeArray,distanceArray);
        System.out.println();
        System.out.println("Length Of cycle: "+cycleLength(root)+" miles");
        System.out.println();
        System.out.println();
        System.out.println("Looking at every permutation to find the optimal solution");
        System.out.println("The best permutation");
        LinkedList<Double> ans = permute(root);
        for (int i = 0; i < ans.size() - 1; i++) {
            double d = ans.get(i);
            System.out.print((int)d +" ");
        }
        double r = ans.get(0);
        System.out.print((int)r);
        System.out.println();
        System.out.println();
        System.out.println("Optimal Cycle length = "+ans.get(ans.size() - 1)+" miles");
    }
    static class Vertex {
        int index;
        double X;
        double Y;
        Vertex(int index, double x, double  y) {
            this.index = index;
            this.X = x;
            this.Y = y;
        }

        public int getIndex() {
            return index;
        }

        public double getX() {
            return X;
        }

        public double getY() {
            return Y;
        }

//        @Override
//        public String toString() {
//            return "Vertex{" +
//                    "index=" + index +
//                    ", X=" + X +
//                    ", Y=" + Y +
//                    '}';
//        }
@Override
public String toString() {
    return "Vertex{" +
            "index=" + index +
            '}';
}

        public boolean equals(Object x) {
            if (!(x instanceof Vertex))
                return false;
            Vertex key = ((Vertex)x);

            return this.getIndex() == key.getIndex()
                    && this.getX() == key.getX()
                    && this.getY() == key.getY();
        }
    }
}
