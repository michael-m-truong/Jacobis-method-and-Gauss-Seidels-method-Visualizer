import java.io.File;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Scanner;

public class Driver {

    private static int numOfEquations;
    private static double[][] testMatrix1 = {
        {2,-1,0,1},
        {-1,3,-1,8},
        {0,-1,2,-5}
    }; 
    private static double[][] testMatrix2 = {
        {5,-1,2,12},
        {3,8,-2,-25},
        {1,1,4,6}
    }; 
    private static double[][] userAugmentedMatrix;
    private static double userError; // default
    private static double[] startingSolution; 

    public static void main(String[] args) {
        Scanner input = new Scanner(System.in);
        System.out.println("Welcome to Jacobi's Method and Gauss-Seidel method visualizer!");
        System.out.println("1. Select from file");
        System.out.println("2. Manually enter matrix");
        System.out.println("3. Use test case matrix");
        System.out.print("Enter number choice: ");
        int numberChoice = Integer.parseInt(input.nextLine());

        if (numberChoice == 1) {                         // user input file
            List<double[]> listA = new ArrayList<>(); 
            System.out.print("Enter text file name with .txt : ");
            String userInputFile = input.nextLine();
            try {
                File myObj = new File(userInputFile);
                Scanner myReader = new Scanner(myObj);
                
                while (myReader.hasNextLine()) {
                    String data = myReader.nextLine();
                    if (data.replaceAll("\\s", "") == "") break;
                    String arr[] = data.split(" ");
                    arr = Arrays.copyOfRange(arr, 0, arr.length);
                    double [] arrA = Arrays.stream(arr).mapToDouble(Double::valueOf).toArray();
                    listA.add(arrA);
                }
                myReader.close();
            } catch (FileNotFoundException e) {
                System.out.println("An error occurred.");
            }
            userAugmentedMatrix = listA.toArray(new double[0][]);
            numOfEquations = userAugmentedMatrix.length;
            startingSolution = new double[numOfEquations];
            System.out.print("Enter desired error: ");
            userError = input.nextDouble();
            for (int i = 0; i < numOfEquations; i++) {
                System.out.print("Enter starting solution for x" + (i+1) + ": ");
                double startingVal = input.nextDouble();
                startingSolution[i] = startingVal;
            }
        }
        if (numberChoice == 2) {                                     // user input
            System.out.print("Enter number of equations: ");
            numOfEquations = Integer.parseInt(input.nextLine());
            userAugmentedMatrix = new double[numOfEquations][numOfEquations+1];
            startingSolution = new double[numOfEquations];
            for (int i = 0; i < numOfEquations; i++) {
                for (int j = 0; j < numOfEquations; j++) {
                    System.out.print("Enter value for a" + (i+1) + (j+1) + ": ");
                    int val = Integer.parseInt(input.nextLine());
                    userAugmentedMatrix[i][j] = val;
                }
                System.out.print("Enter value for b" + (i+1) + ": ");
                int val = Integer.parseInt(input.nextLine());
                userAugmentedMatrix[i][numOfEquations] = val;
            }
            System.out.print("Enter desired error: ");
            userError = input.nextDouble();
            for (int i = 0; i < numOfEquations; i++) {
                System.out.print("Enter starting solution for x" + (i+1) + ": ");
                double startingVal = input.nextDouble();
                startingSolution[i] = startingVal;
            }
        }
        if (numberChoice == 3) {                // test case
            userAugmentedMatrix = testMatrix2;
            numOfEquations = userAugmentedMatrix.length;
            startingSolution = new double[numOfEquations];
            System.out.print("Enter desired error: ");
            userError = input.nextDouble();
            for (int i = 0; i < numOfEquations; i++) {
                System.out.print("Enter starting solution for x" + (i+1) + ": ");
                double startingVal = input.nextDouble();
                startingSolution[i] = startingVal;
            }
        }
        System.out.println("");
        System.out.println("Starting matrix:");
        for (int i = 0; i < numOfEquations; i++) {
            System.out.println(Arrays.toString(userAugmentedMatrix[i]));
        }
        
        input.close();
        System.out.println();
        System.out.println("Jacobi's Method:");
        Jacobi(testMatrix2, Arrays.copyOf(startingSolution, startingSolution.length), userError);
        System.out.println("---------------------------------------------------------------");
        System.out.println("Gauss-Seidel's Method:");
        GaussSeidel(testMatrix2, Arrays.copyOf(startingSolution, startingSolution.length), userError);
    }


    public static double[] Jacobi(double[][] A, double[] x, double err) {
        int maxIteration = 50;  
        double prevL2Norm = 0;
		int n = x.length;
		double xTemp[] = new double[n];
        System.out.println(x[0]);
		for(int i=0; i<n; i++) 
			xTemp[i] = 0;
		
		for(int k=1; k<=maxIteration; k++)
		{
			double total = 0;
            //calc l2 norm of previous iteration
            for (int i = 0; i < A.length; i++) {  
                total += Math.pow(x[i], 2);
            }
            prevL2Norm = Math.sqrt(total);
			for(int i=0; i<n; i++)
			{
				double sum = 0.0;

				for(int j=0; j<n; j++)
				{
					if(i!=j)
						sum = sum + A[i][j]*x[j];
				}

				double temp = (A[i][n] - sum)/A[i][i];
                

				xTemp[i] = temp;
			}
			

			for(int i=0; i<n; i++)
			{
				x[i] = xTemp[i];
			}

            double sum = 0;
			System.out.print("\nIteration:"+k+": [");
			for(int i=0; i<n; i++)
			{
                sum += Math.pow(x[i], 2);
				System.out.print(String.format("%.9f", x[i])+" ");
                if (i == n-1) {
                    System.out.println("]");
                }
			}
            double L2Norm = Math.sqrt(sum);
            System.out.println("L2Norm: " + L2Norm);
            double absoluteError = Math.abs(L2Norm - prevL2Norm);
            System.out.println("Absolute Error: " + absoluteError);

			if(absoluteError < err) 
			{
				System.out.println("\nConverged");
				return x;
			}
		}
		System.out.println("Did not converge in 50 iterations. Relative error was not less than user defined error");
        return x;
    }

    public static double[] GaussSeidel(double[][] A, double[] x, double err) {
        int maxIteration = 50;  
        double prevL2Norm = 0;
		int n = x.length;
		for(int k=1; k<=maxIteration; k++)
		{
            double total = 0;
            //calc l2 norm of previous iteration
            for (int i = 0; i < A.length; i++) {  
                total += Math.pow(x[i], 2);
            }
            prevL2Norm = Math.sqrt(total);
			for(int i=0; i<n; i++)
			{
				double sum = 0.0;
				
				for(int j=0; j<n; j++)
				{
					if(i!=j)
						sum = sum + A[i][j]*x[j];
				}
				
				double temp = (A[i][n] - sum)/A[i][i];
								
				x[i] = temp;
			}

            double sum = 0;
			System.out.print("\nIteration"+k+": [");
			for(int i=0; i<n ;i++)
			{
                sum += Math.pow(x[i], 2);
				System.out.print(String.format("%.9f", x[i])+" ");
                if (i == n-1) {
                    System.out.println("]");
                }
			}
			double L2Norm = Math.sqrt(sum);
            System.out.println("L2Norm: " + L2Norm);
            double absoluteError = Math.abs(L2Norm - prevL2Norm);
            System.out.println("Absolute Error: " + absoluteError);

			if(absoluteError < err) 
			{
				System.out.println("\nConverged");
				return x;
			}
		}
		System.out.println("Did not converge in 50 iterations. Relative error was not less than user defined error");
        return x;
    }
}

