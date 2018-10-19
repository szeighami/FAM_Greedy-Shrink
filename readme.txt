Readme (ARR)
=========================
This package contains the source code for the Greedy-Shrink algorithm. 
We have used "Boost library" (www.boost.org) in our program for random sampling of utility functions from a uniform distribution. Thus, the user should download this library and place it under the folder "ARR" (with the folder name "boost"). If the the library is in any other directory, replace the directory in -I option below to the desired directory. (i.e. if the boost file is in /home/myFiles/boost, then your compilation command should include -I /home/myFiles)

Usage Step
==========
a. Compilation
	g++ -o run -I . arr.cpp qsort.cpp
b. Execution
	./run PointsFile k n d N
Where PointsFile contains the all the points in the dataset, k is the size of the solution returned, n is the size of the dataset, d is the dimensionality of the dataset and N is the sample size. The algorithm assumes a uniform distribution of linear utility funcitons and samples N utility functions from that distribution.
In the output file, you can see the results of the algorithm. In the output points section, you can find selection set returned by the algorithm. The first integer shows the zero-based number of the point in the input file and the other numbers are the attributes of the point in each dimension. "Average number of points avoided" shows the number of points avoided (for which the average regret ratio was not calculated) in each iteration averaged over the total number of iterations using the optimization techniques and improvements described in section 5.4 of the paper. "Average number of users avoided" shows the number of utility functions avoided (when calculating the avergae regret ratio) in each iteration averaged over the total number of iterations using the optimization techniques.

===========
About dataset:
The dataset should contain n points in d dimensions. Each point must be written in a separate line with its dimensions being separated by tabs.

===========
Example (The hotel example from the paper):
The dataset contains 4 points in 3 dimensions. Use
	./run ExamplePoints.txt 2 4 3 10000

	./run ExamplePoints.txt 2 4 3 4 ExampleUtilityFunctions.txt

to select 2 points from the 4 points in the 3-dimensional dataset using a sample size of 10000 utility functions. The program outputs the results to the file result.txt as well. You will get an output similar to this:
