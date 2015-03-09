/**
 * CS350 Term Project
 * @authors Parker Berger, Todd Brochu, Seyed Nima Sajadpour
 * Portland State University
 * CS350, Winter 2015
 * 
 * This program implements two versions of Convex Hull:
 * QuickHull and brute force.  It collects performance
 * data on both.
 *
 */

import java.lang.management.*;
import java.lang.Math;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashSet;
import java.util.Set;
import java.util.TreeSet;
import java.io.*;
import java.text.SimpleDateFormat;
import java.util.Calendar;
import java.util.Scanner;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.Files;
import java.nio.charset.Charset;
import java.nio.charset.StandardCharsets;

public class Main {

	public static void main(String[] args) {
		int i = 0;
		int errOpt = 0;
		boolean quit = false; 
		int menuSelect = 0;
		final int numberOfRuns = 25; 
		final int base = 2;
		int exponent = 4;		//increment this up to 10 to collect performance data
		int arraySize;
		ArrayList<CartesianPoint> points = new ArrayList<CartesianPoint>();
		Set<CartesianPoint> bruteForcePoints = new HashSet<CartesianPoint>();
		Set<CartesianPoint> quickHullPoints = new HashSet<CartesianPoint>();
		Set<CartesianPoint> correctTestPoints = new HashSet<CartesianPoint>();		
		Scanner in = new Scanner(System.in);
	
		do {
			System.out.println("\n\n*********** Convex Hull ***********");
			System.out.println("Do you want to load random points or use points from a file?");
			System.out.println("	0. Quit");
			System.out.println("	1. Random Points");
			System.out.println("	2. Load Points From a File");	
			menuSelect = in.nextInt();

			switch (menuSelect) 
			{
			case 1:
				System.out.print("OK, I will load some random points for you...");
				arraySize = (int) java.lang.Math.pow(base,exponent); 
				seedArray(points, arraySize);		// Load random points into the ArrayList
				break;
			case 2:
				System.out.println("OK, Let me load up your file for you...");
				System.out.println("Your file must be in the format: <x> <y> (one entry per line)");
				System.out.print("    What is the name of your file?: ");				
				String nameOfFile = in.next();
				System.out.println();
				errOpt = seedArray(points, nameOfFile);
				if(errOpt == 1)
					System.out.println("Your file is formatted wrong!");
				else if(errOpt == 2)
					System.out.println("Your file name is wrong! (did you add \".txt\"?)");
				else
					break;
			case 0:
				quit = true;
				System.out.println("Adios!");
				System.exit(1);
				break;
			default:
				System.out.println("I Don't See That Option???");
				System.exit(1);
			}

			System.out.println("\n\nExecuting Brute Force Convex Hull " + numberOfRuns + " times.");
			//begin measuring the execution time
			long startUserTimeNano = System.nanoTime();
			for (i = 0; i < numberOfRuns; i++) 
			{
				//run brute force version 25 times
				bruteForcePoints = BruteForceConvexHull(points);
			}
			// stop the time and calculate difference
			long bruteForceTime = System.nanoTime() - startUserTimeNano;
		
			System.out.println("Executing QuickHull " + numberOfRuns + " times.");
			//begin measuring the execution time
			startUserTimeNano = System.nanoTime();
			for (i = 0; i < numberOfRuns; i++) 
			{
				//run quickhull algorithm 25 times
				quickHullPoints = QuickHull(points);
			}
			// stop the time and calculate difference
			long quickHullTime = System.nanoTime() - startUserTimeNano;
			
			System.out.println("Executing Control Algorithm " + numberOfRuns + " times.");
			//begin measuring the execution time
			startUserTimeNano = System.nanoTime();
			for (i = 0; i < numberOfRuns; i++) 
			{
				//run known correct algorithm 25 times
				correctTestPoints = testAlgorithm(points);
			}
			// stop the time and calculate difference
			long knownTestTime = System.nanoTime() - startUserTimeNano;
			
			System.out.println("Testing complete.  Writing to file.");
			writeFile(points, bruteForcePoints, quickHullPoints, correctTestPoints, bruteForceTime, quickHullTime, knownTestTime);

			testSuit(bruteForcePoints, quickHullPoints, correctTestPoints, bruteForceTime, quickHullTime, knownTestTime);
			clearArrays(points, bruteForcePoints, quickHullPoints, correctTestPoints);
			System.out.println("Done!");
		} while (!quit);

	}

	
	/*********************************************BRUTE FORCE*********************************/
	
	private static Set<CartesianPoint> BruteForceConvexHull(ArrayList<CartesianPoint> points) 
	{
		Set<CartesianPoint> hull = new TreeSet<CartesianPoint>();
		Set<Integer> sides = new HashSet<Integer>();
		int a, b, c;
		int position; //the result of ax+by for each evaluated point
		int j = 1;
		
		for (int i = 0; i < points.size()-1; ) 
		{
			//calculate a,b and c to perform geometry
			while (j < points.size()) {
				a = points.get(j).y - points.get(i).y;
				b = points.get(i).x - points.get(j).x;
				c = points.get(i).x * points.get(j).y - points.get(j).x * points.get(i).y;
				
				//evaluate relative position of all points not related to the current pair
				for (CartesianPoint cp : points) 
				{ 
					if (points.indexOf(cp) != i && points.indexOf(cp) != j) 
					{
						position = a*cp.x + b*cp.y;
						
						if (position < c) sides.add(-1); //if the evaluated point is less than c, then mark with -1
						else if (position == c) sides.add(0); //if the evaluated point is equal to c, then mark with 0
						else sides.add(1);  //if the evaluated point is greater than c, then mark with 1
					}		
				}
				
				//if all evaluated points are on the same side of the line thru pair (i, j),
				//then add the pair to hull
				if (sides.size() == 1 && !sides.contains(0)) 
				{ 
					hull.add(points.get(i));
					hull.add(points.get(j));
				}
				
				//clear the set for use on the next round of evaluations
				sides.clear();
				
				j++;
			}
			i++;
			j = i+1;			
		}
		return hull;
	}

	/*********************************************QUICK HULL*********************************/
	
	private static Set<CartesianPoint> QuickHull(ArrayList<CartesianPoint> points) 
		{
		Set<CartesianPoint> hull = new HashSet<CartesianPoint>();
		
		ArrayList<CartesianPoint> upperHull = new ArrayList<CartesianPoint>(); // Upper hull
		ArrayList<CartesianPoint> lowerHull = new ArrayList<CartesianPoint>(); // Lower hull
		
		int a, b, c;
		int position; //the result of ax+by for each evaluated point
		
		int minVal = points.get(0).x;	// Left most point
		int maxVal = points.get(0).x;	// Right mot point
		int indexMin = 0;		// Index of minimum point
		int indexMax = 0; 	// Index of maximum point
		
		// Find left and right most points
		for(int i = 0; i < points.size(); ++i) 
		{
			if(minVal > points.get(i).x)
			{
				minVal = points.get(i).x;
				indexMin = i;
			} else if(maxVal < points.get(i).x)
			{
				maxVal = points.get(i).x;
				indexMax = i;
			}
		}
		
		// Add left and right most points to the convex hull
		hull.add(points.get(indexMin));
		hull.add(points.get(indexMax));
		
		// Segment the points on the upper and lower hull
		a = points.get(indexMax).y - points.get(indexMin).y;
		b = points.get(indexMin).x - points.get(indexMax).x;
		c = points.get(indexMin).x * points.get(indexMax).y - points.get(indexMax).x * points.get(indexMin).y;
		
		// Find and load the upper and lower hull arrays
		for(int i = 0; i < points.size(); ++i) 
		{
			position = a*points.get(i).x + b*points.get(i).y;
			
			// If point is on the top side of the line add to the upper hull
			if (position < c) 
				upperHull.add(points.get(i));
			// If point is on the bottom side of the line add to the lower hull
			else if (position > c) 
				lowerHull.add(points.get(i));
			// Else do nothing because the point is on the line on not on the convex hull
		}

		FindHull(hull, upperHull, points.get(indexMin), points.get(indexMax));
		FindHull(hull, lowerHull, points.get(indexMax), points.get(indexMin));
		
		return hull;
	}
	
	private static void FindHull(Set<CartesianPoint> hull, ArrayList<CartesianPoint> points, CartesianPoint p1, CartesianPoint p2) 
	{
		ArrayList<CartesianPoint> S1 = new ArrayList<CartesianPoint>(); // Upper hull
		ArrayList<CartesianPoint> S2 = new ArrayList<CartesianPoint>(); // Lower hull
		
		int a1, b1, c1;
		int a2, b2, c2;
		int position1;
		int position2;
		
		int maxDistance = 0;	// Maximum distance from the line
		int indexMax = 0;			// Index of point with maximum distance from the line
		
		//If points is empty return
		if(points.isEmpty())
			return;
    if(points.size() == 1)
    {
		  hull.add(points.get(0));
		  return;
    }

		a1 = p2.y - p1.y;
		b1 = p1.x - p2.x;
		
		// Find farthest point, C, from p1p2
		for(int i = 0; i < points.size(); ++i) 
		{
			position1 = a1*points.get(i).x + b1*points.get(i).y;
			
			if(Math.abs(position1) > maxDistance)
			{
				maxDistance = Math.abs(position1);
				indexMax = i;
			}	
		}
		
		// Add point at indexMax to convexHull
		hull.add(points.get(indexMax));
		
		// Find S1 and S2. Exclude S2.
		CartesianPoint C = points.get(indexMax);
		
		// Convex Hull = to the left
		a1 = C.y - p1.y;
		b1 = p1.x - C.x;
		c1 = p1.x * C.y - C.x * p1.y;
		
		// Convex Hull = to the right
		a2 = p2.y - C.y;
		b2 = C.x - p2.x;
		c2 = C.x * p2.y - p2.x * C.y;
		
		// Find and add points that are outside of the p1, p2, C boundery
		for(int i = 0; i < points.size(); ++i) 
		{
			position1 = a1*points.get(i).x + b1*points.get(i).y;
			position2 = a2*points.get(i).x + b2*points.get(i).y;
			
			if(position1 < c1)
				S1.add(points.get(i));
			else if(position2 < c2)
				S2.add(points.get(i));
		}
		
		FindHull(hull, S1, p1, C);
		FindHull(hull, S2, C, p2);
	}
	
	/*******************************KNOWN CORRECT TEST ALGORITHM*******************************/

	 private static Set<CartesianPoint> testAlgorithm(ArrayList<CartesianPoint> points)
    {
				ArrayList<CartesianPoint> convexHull = new ArrayList<CartesianPoint>();
	
        int minPoint = -1, maxPoint = -1;
        int minX = Integer.MAX_VALUE;
        int maxX = Integer.MIN_VALUE;

        for (int i = 0; i < points.size(); i++)
        {
            if (points.get(i).x < minX)
            {
                minX = points.get(i).x;
                minPoint = i;
            }
            if (points.get(i).x > maxX)
            {
                maxX = points.get(i).x;
                maxPoint = i;
            }
        }

        CartesianPoint A = points.get(minPoint);
        CartesianPoint B = points.get(maxPoint);

        convexHull.add(A);
        convexHull.add(B);

        ArrayList<CartesianPoint> leftSet = new ArrayList<CartesianPoint>();
        ArrayList<CartesianPoint> rightSet = new ArrayList<CartesianPoint>();

        for (int i = 0; i < points.size(); i++)
        {
            CartesianPoint p = points.get(i);
            if (pointLocation(A, B, p) == -1)
                leftSet.add(p);
            else if (pointLocation(A, B, p) == 1)
                rightSet.add(p);
        }

        hullSet(A, B, rightSet, convexHull);
        hullSet(B, A, leftSet, convexHull);

				Set<CartesianPoint> setConvexHull = new HashSet<CartesianPoint>(convexHull);

				return setConvexHull;
    }

		private static int distance(CartesianPoint A, CartesianPoint B, CartesianPoint C)
    {
        int ABx = B.x - A.x;
        int ABy = B.y - A.y;
        int num = ABx * (A.y - C.y) - ABy * (A.x - C.x);

        if (num < 0)
            num = -num;
        return num;
    }

    private static void hullSet(CartesianPoint A, CartesianPoint B, ArrayList<CartesianPoint> set,
            ArrayList<CartesianPoint> hull)
    {

        int insertPosition = hull.indexOf(B);

        if (set.size() == 0)
            return;
        if (set.size() == 1)
        {
            CartesianPoint p = set.get(0);
            hull.add(insertPosition, p);
            return;
        }

        int dist = Integer.MIN_VALUE;
        int furthestPoint = -1;

        for (int i = 0; i < set.size(); i++)
        {
            CartesianPoint p = set.get(i);
            int distance = distance(A, B, p);
            if (distance > dist)
            {
                dist = distance;
                furthestPoint = i;
            }
        }

        CartesianPoint P = set.get(furthestPoint);
        hull.add(insertPosition, P);

        ArrayList<CartesianPoint> leftSetAP = new ArrayList<CartesianPoint>();

        for (int i = 0; i < set.size(); i++)
        {
            CartesianPoint M = set.get(i);
            if (pointLocation(A, P, M) == 1)
            {
                leftSetAP.add(M);
            }
        }

        ArrayList<CartesianPoint> leftSetPB = new ArrayList<CartesianPoint>();

        for (int i = 0; i < set.size(); i++)
        {
            CartesianPoint M = set.get(i);
            if (pointLocation(P, B, M) == 1)
            {
                leftSetPB.add(M);
            }
        }

        hullSet(A, P, leftSetAP, hull);
        hullSet(P, B, leftSetPB, hull);
    }

    private static int pointLocation(CartesianPoint A, CartesianPoint B, CartesianPoint P)
    {
        int cp1 = (B.x - A.x) * (P.y - A.y) - (B.y - A.y) * (P.x - A.x);

        if (cp1 > 0)
            return 1;
        else if (cp1 == 0)
            return 0;
        else
            return -1;
    }
    
    /*********************************************UTILITY FUNCTIONS*********************************/
	
	// Loads the array passed in as an argument with Cartesian Points of random numbers
	private static void seedArray(ArrayList<CartesianPoint> points, int size) 
	{
		for (int i = 0; i < size; i++ ) 
		{
			points.add(new CartesianPoint((int)(Math.random()*(200 + 1)) - 100, (int)(Math.random()*(200 + 1)) - 100 ));
		}
	}

	// Loads the array passed in as an argument with Cartesian Points from a file
	private static int seedArray(ArrayList<CartesianPoint> points, String nameOfFile)
	{
		Path path = Paths.get(nameOfFile);

		try (BufferedReader reader = Files.newBufferedReader(path, StandardCharsets.UTF_8))
		{
			String line = null;
			while((line = reader.readLine()) != null)
			{
				String[] split = line.split("\\s+");

				if(split.length > 2)	// File is formatted wrong!
					return 1;

				int x = Integer.parseInt(split[0]);
				int y = Integer.parseInt(split[1]);

				points.add(new CartesianPoint(x, y));
			}
		}catch (Exception e) {
			return 2;	// File name can not be found
		}
		
		return 3; // No error
	}

	// Prints the results to a log .txt file
	private static void writeFile(ArrayList<CartesianPoint> points, Set<CartesianPoint> bruteForcePoints, Set<CartesianPoint> quickHullPoints, Set<CartesianPoint> testAlgorithmPoints, long BruteForceTime, long QuickHullTime, long testAlgorithmTime) {

		BufferedWriter writer = null;		
		int arraySize = points.size();

		try {
			String timeStamp = new SimpleDateFormat("MM-dd-yyyy_HH.mm.ss").format(Calendar.getInstance().getTime());
			File logFile = new File("logs/" + timeStamp + ".txt");
			
			writer = new BufferedWriter(new FileWriter(logFile));

			writer.write("Cartesian Points: \n");
			for (int i = 0; i < arraySize; i++ ) {
				writer.write("(" + points.get(i).x + ", " + points.get(i).y + ")\n");
			}

			writer.write("\nToatal number of nodes: " + arraySize);
			writer.write("\nBruteForce:    Average execution time: " + BruteForceTime);
			writer.write("\nQuickHull:     Average execution time: " + QuickHullTime);
			writer.write("\nTestAlgorithm: Average execution time: " + testAlgorithmTime);

			writer.write ("\n\nBruteForce: Convex Hull Points: \n");
			for (CartesianPoint p : bruteForcePoints) 
			{
				writer.write("(" + p.x + ", " + p.y + ")\n");
			}

			writer.write ("\n\nQuickHull: Convex Hull Points: \n");
			for (CartesianPoint s : quickHullPoints) 
			{
				writer.write("(" + s.x + ", " + s.y + ")\n");
			}

			writer.write ("\n\nKnown Correct Test Algorithm: Convex Hull Points: \n");
			for (CartesianPoint t : testAlgorithmPoints) 
			{
				writer.write("(" + t.x + ", " + t.y + ")\n");
			}

 		}	catch (Exception e) {
			e.printStackTrace();
		} finally {
			try {
				writer.close();
			} catch (Exception e) {

			}
		}

	}

private static void testSuit(Set<CartesianPoint> bruteForcePoints, Set<CartesianPoint> quickHullPoints, 
								 Set<CartesianPoint> correctTestPoints, long BruteForceTime, 
								 long QuickHullTime, long testAlgorithmTime) 
	{
		
		ArrayList<CartesianPoint> brutelist = new ArrayList<CartesianPoint>();
		ArrayList<CartesianPoint> quicklist = new ArrayList<CartesianPoint>();
		ArrayList<CartesianPoint> testlist = new ArrayList<CartesianPoint>();

		// Sort the BruteForce
		for(CartesianPoint b : bruteForcePoints)
			brutelist.add(new CartesianPoint( b.x, b.y ));
		brutelist = bubbleSort(brutelist, brutelist.size()-1);

		// Sort the quickHullPoints
		for(CartesianPoint q : quickHullPoints)
			quicklist.add(new CartesianPoint( q.x, q.y ));
		quicklist = bubbleSort(quicklist, quicklist.size()-1);
				
		// Sort the testAlgorithm (known correct algorithm)
		for(CartesianPoint t : correctTestPoints)
			testlist.add(new CartesianPoint( t.x, t.y ));
		testlist = bubbleSort(testlist, testlist.size()-1);

		for(int i=0; i < testlist.size(); i++) {
			// Compare the contents of bruteForcePoints w/ quickHullPoints
			if ((testlist.get(i).x == brutelist.get(i).x) && (testlist.get(i).y == brutelist.get(i).y) ) {
				if((testlist.get(i).x == quicklist.get(i).x)&& (testlist.get(i).y == quicklist.get(i).y) ) {
					System.out.println(" Perfect!!!! ");
					System.out.println("(" + testlist.get(i).x + ", " + testlist.get(i).y + ")");
				}else{
					//Incorect quickhull
					System.out.println(" Incorrect QuickHull");
					System.out.println("(" + testlist.get(i).x + ", " + testlist.get(i).y + ")");
				}
			}else if((testlist.get(i).x == quicklist.get(i).x) && (testlist.get(i).y == quicklist.get(i).y)){
				//incorrect bruteforce but correct quickhull
				System.out.println(" Incorrect Bruteforce");
				System.out.println("(" + testlist.get(i).x + ", " + testlist.get(i).y + ")");
			}else{
				System.out.println(" Both implemented algorithm have the wrong result!! Sorry!!");
				System.out.println("(" + testlist.get(i).x + ", " + testlist.get(i).y + ")");
			}
		}
		System.out.println("\n\n Number of bruteforce: " + brutelist.size());
		System.out.println(" Number of QuickHull: " + quicklist.size());
		System.out.println(" Number of test: " + testlist.size());
		System.out.println(" BruteForce Time is: " + BruteForceTime);
		System.out.println(" QuickHull Time is: " + QuickHullTime);
		System.out.println(" Test Time is: " + testAlgorithmTime);

		if( QuickHullTime < BruteForceTime )
			System.out.println(" QuickHull is faster than BruteForce by: " + (BruteForceTime-QuickHullTime) );
		else 
			System.out.println(" QuickHull is faster than BruteForce by: " + (QuickHullTime-BruteForceTime) );
		return;
	}
	// Clears all the data in the arrays passed in as arguments
	private static void clearArrays(ArrayList<CartesianPoint> points, Set<CartesianPoint> bruteForcePoints, Set<CartesianPoint> quickHullPoints, Set<CartesianPoint> correctTestPoints) 
	{
		points.clear();
		bruteForcePoints.clear();
		quickHullPoints.clear();
		correctTestPoints.clear();
	}


	private static ArrayList<CartesianPoint> bubbleSort(ArrayList<CartesianPoint> list, int length) {

		for (int i = 0; i < length; i++) {
			for (int j = 0; j < length - i; j++)
	      	{
	      		if (list.get(j).x > list.get(j+1).x)
	         	{
		         		int tempx = list.get(j).x;
		              	int tempy = list.get(j).y;
		              	list.get(j).x = list.get(j+1).x;
		              	list.get(j).y = list.get(j+1).y;
		              	list.get(j+1).x = tempx;
		              	list.get(j+1).y = tempy;
	        	}else if(list.get(j).x == list.get(j+1).x) {
	        		if(list.get(j).y > list.get(j+1).y){
		         		int tempx = list.get(j).x;
		              	int tempy = list.get(j).y;
		              	list.get(j).x = list.get(j+1).x;
		              	list.get(j).y = list.get(j+1).y;
		              	list.get(j+1).x = tempx;
		              	list.get(j+1).y = tempy;
	              	}
	        	} 
	      	} 
	    }
	    return list;
	}


}
