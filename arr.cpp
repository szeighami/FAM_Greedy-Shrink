#include<iostream>
#include<string.h>
#include<fstream>
#include<algorithm>
#include<math.h>
#include<vector>
#include <stdlib.h>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/random/random_device.hpp>
#include <lp/data_utility.h>
#include <lp/evaluate.h>
#include <time.h>
#include <stdio.h>

using namespace std;

#define veryLongInt long long int
int memoryUsage = 0;

#ifndef WIN32
#include <sys/resource.h>
#include <sys/times.h>
#endif


#ifndef WIN32
void calculateExecutionTime(struct rusage *myTimeStart, struct rusage *myTimeEnd, float *userTime, float *sysTime)
{
	(*userTime) =
		((float) (myTimeEnd->ru_utime.tv_sec  - myTimeStart->ru_utime.tv_sec)) +
		((float) (myTimeEnd->ru_utime.tv_usec - myTimeStart->ru_utime.tv_usec)) * 1e-6;
	(*sysTime) =
		((float) (myTimeEnd->ru_stime.tv_sec  - myTimeStart->ru_stime.tv_sec)) +
		((float) (myTimeEnd->ru_stime.tv_usec - myTimeStart->ru_stime.tv_usec)) * 1e-6;
	
}
#endif


struct Node{
	int index;
	double arr;
	Node* next;
};

veryLongInt getSamplingSize(int dim, int kValue, float epsilonValue, float deltaValue)
{
	
	float sOne = dim - 1;
	
	float returnValuefloat = 3 * log( 1.0 / deltaValue ) / (epsilonValue*epsilonValue);
	
	veryLongInt returnValue = floor(returnValuefloat);
	
	return returnValue;
}

float** utilFuncArray;
float** raw;
float* satInD;
float* satInS;
int bestPoint;
float maxUtility;
int userID;


int main(int argc, char *argv[])
{ 
	for (int runCounter = 0; runCounter < 1; runCounter++)
	{
		#ifndef WIN32
			float  userTime, sysTime;
			struct rusage myTime_start, myTime_end;
			float  userTime2, sysTime2;
			struct rusage myTime_start2, myTime_end2;
			// start time 
			getrusage(RUSAGE_SELF,&myTime_start);
		#endif

        boost::mt19937 rng(time(0));
		boost::uniform_real<float> u(0, 1);
		boost::variate_generator<boost::mt19937&, boost::uniform_real<float> > gen(rng, u);

		
        //initializations
		FILE * fp;
        char resFile[70];
        sprintf(resFile, "result.txt");
		fp = fopen (resFile, "a");
		
		int dim, size, kInput, sampleSizeInput;
		bool inputSample = false;
		string filename = argv[1];
		kInput = atoi(argv[2]);
		size = atoi(argv[3]);
		dim = atoi(argv[4]);
		sampleSizeInput = atoi(argv[5]);
		
		ifstream inFile;
		inFile.open(filename.c_str());
		ifstream inFileUtil;
		if (argc == 7)
		{
			inputSample = true;
		    string samplefilename = argv[6];
			inFileUtil.open(samplefilename .c_str());
		}

		
		raw = new float*[size];
		memoryUsage += sizeof(float*) * size;
		bool* isResult = new bool[size];
		memoryUsage += sizeof(bool*) * size;
		Node** bestPointsUsers = new Node*[size];
		memoryUsage += sizeof(Node*) * size;

        //Reading the points
		for(int i = 0; i < size; i++)
		{
			raw[i] = new float[dim];
			memoryUsage += sizeof(float) * dim;
			isResult[i] = false;
			for(int j = 0; j < dim; j++)
			{
				inFile >> raw[i][j];
			}
			bestPointsUsers[i] = NULL;
		}
		inFile.close();
		cout << "Size of the dataset: " << size << endl;
		cout << "dimensionality of dataset: " << dim << endl;


		if(kInput >= size)
		{
			char inputChar;
			cout<< "The value of k is equal or greater than the data size, continue? (y/n)"<<endl;
			cout << "Memory usage:" << memoryUsage << endl;
			if(inputChar == 'y')
			{
				cout << "Exit!" << endl;
				return 100;
			}
			kInput = size;
		}

		veryLongInt samplingSize = 0;
		if (sampleSizeInput == 0)
			samplingSize = getSamplingSize(dim, kInput, 0.0001 /*epilson*/, 0.1 /*delta*/);
		else
			samplingSize  = sampleSizeInput;
		
		cout << "Value of k: " << kInput <<endl;
		cout << "Sampling Size: " << samplingSize << endl;
		cout << "-----------------" << endl;
		
		fprintf(fp, "k = %d, Sampling Size = %d, input = %s--------------------------------------------------------------------\n", kInput, sampleSizeInput, filename.c_str());
		utilFuncArray = new float*[samplingSize];
		memoryUsage += sizeof(float*) * samplingSize;


		satInD = new float[samplingSize];
		memoryUsage += sizeof(float) * samplingSize;
		satInS = new float[samplingSize];
		memoryUsage += sizeof(float) * samplingSize;
		Node* headS = NULL;
		int solutionSize = 0;
		int prevBestPoint = -1;
		bestPoint = -1;
		maxUtility = -1;
        if (inputSample)
            cout << "reading sample" << endl;
        //reading the utility functions and finding best points
		for (int i = 0; i < samplingSize; i++)
		{
			utilFuncArray[i] = new float[dim];
			memoryUsage += sizeof(float) * dim;
			for(int d = 0; d < dim; d++)
			{
				if (inputSample)
					inFileUtil >> utilFuncArray[i][d];
				else
					utilFuncArray[i][d] = gen();
			}
			
			maxUtility = -1;
            for (int j = 0; j < size; j++)
            {
                float utility = 0;
                for(int d = 0; d < dim; d++)
                {
                    utility += utilFuncArray[i][d] * raw[j][d];
                }
                if (utility > maxUtility)
                {
                    maxUtility = utility;
                    bestPoint = j;
                    satInD[i] = utility;
                    satInS[i] = utility;
                }
            }

			if (!isResult[bestPoint])
			{
				isResult[bestPoint] = true;
				solutionSize++;
				if (headS == NULL)
				{
					headS = new Node;
					memoryUsage += sizeof(Node);
					headS->index = bestPoint;
					headS->next = NULL;
				}
				else
				{
					Node* temp = new Node;
					memoryUsage += sizeof(Node);
					temp->index = bestPoint;
					temp->next = headS;
					headS = temp;
				}
			}
			if (bestPointsUsers[bestPoint] == NULL)
			{
				bestPointsUsers[bestPoint] = new Node;
				memoryUsage += sizeof(Node);
				bestPointsUsers[bestPoint]->index = i;
				bestPointsUsers[bestPoint]->next = NULL;
			}
			else
			{
				Node* temp = new Node;
				memoryUsage += sizeof(Node);
				temp->index = i;
				temp->next = bestPointsUsers[bestPoint];
				bestPointsUsers[bestPoint] = temp;
			}
		}

	#ifndef WIN32
		// end time 
		getrusage(RUSAGE_SELF,&myTime_end);
		// output execution time 
		calculateExecutionTime(&myTime_start, &myTime_end, &userTime, &sysTime);

		fprintf(fp, "Preprocessing time : \n");
		fprintf(fp, "User time : %f seconds\n", userTime);
		fprintf(fp, "System time : %f seconds\n\n", sysTime);
		// start time 
		getrusage(RUSAGE_SELF,&myTime_start2);
	#endif

		if (solutionSize <= kInput)
		{
			cout << "Fewer than k points for all the uses!" << endl;
			fprintf(fp, "Fewer than k points for all the uses!\n");
			fprintf(fp, "memoryUsage (bytes): %d\n", memoryUsage);
			cout << "Average Regret Ratio: " <<  0 << endl;
			fprintf(fp, "Average Regret Ratio: 0");
			fprintf(fp, "\n\n");
			fclose(fp);
            continue;
		}
        //Running the main algorithm
		cout << "Finding the solution" << endl;
		float minArr = 2;
		float sumRegretRatio = 0;
		double lowerBound = -1;
		int iterationsAvoided = 0;
		int usersAvoided = 0;
		int usersIterations = 0;
		for (int i = solutionSize; i > kInput; i--)
		{
			minArr = 2;
			int worstPoint = -1;
			Node* worstPointPointer = NULL;
			int iterationCount = 0;
			Node* worstPointPrev = NULL;
			Node* prev = NULL;
			for (Node* iter = headS; iter != NULL; prev = iter, iter = iter->next)
			{
				iterationCount++;
				if (iter->arr < 0)
				{	
					worstPoint = iter->index;
					worstPointPointer = iter;
					worstPointPrev = prev;
					break;	
				}
				if (lowerBound != -1 && iter->arr > lowerBound)
					break;

				int n = iter->index;
				int usersIterated = 0;
				float tempSumRegretRatio = sumRegretRatio;
				
				for (Node* usersIterator = bestPointsUsers[n]; usersIterator != NULL; usersIterator = usersIterator->next)
				{
					usersIterated++;
					int N = usersIterator->index;
					tempSumRegretRatio -= (satInD[N] - satInS[N]) / satInD[N];
		
					float maxUtility = -1;
					for (Node* iter2 = headS; iter2 != NULL; iter2 = iter2->next)
					{
						int j = iter2->index;
						if (j == n)
							continue;

						float utility = 0;
						for(int d = 0; d < dim; d++)
						{
							utility += utilFuncArray[N][d] * raw[j][d];
						}
						if (utility > maxUtility)
						{
							maxUtility = utility;
						}
					}
					tempSumRegretRatio += (satInD[N] - maxUtility) / satInD[N];
				}
				usersIterations++;
				usersAvoided += samplingSize - usersIterated;
				float arr = tempSumRegretRatio / samplingSize;
				iter->arr = arr;

				if (arr < minArr)
				{
					minArr = arr;
					worstPoint = n;
					worstPointPointer = iter;
					worstPointPrev = prev;
				}
			}
			iterationsAvoided += i - iterationCount;
		
			if (worstPointPrev !=NULL)
			{
				worstPointPrev->next = worstPointPointer->next;
			}
			else
			{
				headS = worstPointPointer->next;
			}
			delete worstPointPointer;

			for (Node* usersIterator = bestPointsUsers[worstPoint]; usersIterator != NULL;)
			{
				int N = usersIterator->index;
				sumRegretRatio -= (satInD[N] - satInS[N]) / satInD[N];
				int bestPoint = -1;
				float maxUtility = -1;			
				for (Node* iter2 = headS; iter2 != NULL; iter2 = iter2->next)
				{
					int j = iter2->index;
					float utility = 0;
					for(int d = 0; d < dim; d++)
					{
						utility += utilFuncArray[N][d] * raw[j][d];
					}
					if (utility > maxUtility)
					{
						maxUtility = utility;
						bestPoint = j;
						satInS[N] = utility;
					}
				}
				if(bestPointsUsers[bestPoint] == NULL)
				{
					bestPointsUsers[bestPoint] = new Node;
					bestPointsUsers[bestPoint]->index = N;
					bestPointsUsers[bestPoint]->next = NULL;
				}
				else
				{
					Node* temp = new Node;
					temp->index = N;
					temp->next = bestPointsUsers[bestPoint];
					bestPointsUsers[bestPoint] = temp;
				}
				sumRegretRatio += (satInD[N] - satInS[N]) / satInD[N];
				Node* temp = usersIterator;
				usersIterator = usersIterator->next;
				delete temp;
			}

			bestPointsUsers[worstPoint] = NULL;
			for (Node* iter1 = headS; iter1 != NULL; iter1 = iter1->next)
			{
				for (Node* iter2 = iter1; iter2 != NULL; iter2 = iter2->next)
				{
					if (iter1->arr > iter2->arr)
					{
						double temp = iter1->arr;
						iter1->arr = iter2->arr;
						iter2->arr = temp;

						int temp2 = iter1->index;
						iter1->index = iter2->index;
						iter2->index = temp2;
					}
				}
			}

            int counter = 0;
			Node* iter = headS;
			Node* prev1 = NULL;
            for (; iter != NULL && counter <= 1; prev1 = iter, iter = iter->next, counter++);
			if (iter != NULL)
			{
				int n = iter->index;
				float tempSumRegretRatio = sumRegretRatio;
				for (Node* usersIterator = bestPointsUsers[n]; usersIterator != NULL; usersIterator = usersIterator->next)
				{
					int N = usersIterator->index;
					tempSumRegretRatio -= (satInD[N] - satInS[N]) / satInD[N];
		
					float maxUtility = -1;
					for (Node* iter2 = headS; iter2 != NULL; iter2 = iter2->next)
					{
						int j = iter2->index;
						if (j == n)
							continue;

						float utility = 0;
						for(int d = 0; d < dim; d++)
						{
							utility += utilFuncArray[N][d] * raw[j][d];
						}
						if (utility > maxUtility)
						{
							maxUtility = utility;
						}
					}
					tempSumRegretRatio += (satInD[N] - maxUtility) / satInD[N];
				}
				float arr = tempSumRegretRatio / samplingSize;
				lowerBound = arr;
				iter->arr = arr;
			
				if (prev1 != NULL)
					prev1->next = iter->next;
				
				iter->next = headS;
				headS = iter;
				
				for (; iter->next!= NULL && iter->next->arr < iter->arr; iter = iter->next)
				{	
					Node* iter2 = iter->next;
					double temp = iter->arr;
					iter->arr = iter2->arr;
					iter2->arr = temp;

					int temp2 = iter->index;
					iter->index = iter2->index;
					iter2->index = temp2;
					
				}
			}
		}

	#ifndef WIN32
		// end time 
		getrusage(RUSAGE_SELF,&myTime_end2);
		// output execution time 
		calculateExecutionTime(&myTime_start2, &myTime_end2, &userTime2, &sysTime2);

		fprintf(fp, "Query time : \n");
		fprintf(fp, "User time : %f seconds\n", userTime2);
		fprintf(fp, "System time : %f seconds\n\n", sysTime2);
	#endif
		cout << "memoryUsage (bytes):" << memoryUsage << endl;
		fprintf(fp, "memoryUsage (bytes): %d \n\n", memoryUsage);
		cout << "Average Regret Ratio: " <<  minArr << endl;
        
        //calculating more statistics
        for (int i = 0; i < samplingSize; i++)
        {
            for(int d = 0; d < dim; d++)
            {
                if (inputSample)
					inFileUtil >> utilFuncArray[i][d];
				else
                    utilFuncArray[i][d] = gen();
            }
            float maxUtility  = 0;
            for (int j = 0; j < size; j++)
            {
                float utility = 0;
                for(int d = 0; d < dim; d++)
                {
                    utility += utilFuncArray[i][d] * raw[j][d];
                }
                if (utility > maxUtility)
                {
                    maxUtility = utility;
                    satInD[i] = utility;
                }
            }
        }

        float sumVariance = 0;
        float sampleMRR = 0;
        vector<float> user_rr;
        for (int user = 0; user < samplingSize; user++)
        {
            float maxUtil = 0;
            for (Node* iter = headS; iter != NULL; iter = iter->next)
            {
                int i = iter->index;
                float utility = 0;
                for(int d = 0; d < dim; d++)
                    utility += utilFuncArray[user][d] * raw[i][d];
                if (utility > maxUtil)
                    maxUtil = utility;
            }
            satInS[user] = maxUtil;
            float rr = (satInD[user]-maxUtil)/satInD[user];
            user_rr.push_back(rr);
            if (rr > sampleMRR)
                sampleMRR = rr;
            float diff = (rr - minArr);
            sumVariance += diff*diff;
        }
        float variance = sumVariance/samplingSize;
        float sd = sqrt(variance);
        int notWithinSD = 0;

        for (int user = 0; user < samplingSize; user++)
        {
            float rr = (satInD[user]-satInS[user])/satInD[user];
            if (rr > minArr + sd || rr < minArr - sd)
                notWithinSD++;
        }

        std::sort(user_rr.begin(), user_rr.end());
		float avgIterAvoided = iterationsAvoided * 1.0 / (solutionSize - kInput);
		float avgUsersAvoided = usersAvoided * 1.0 / (usersIterations);
		fprintf(fp, "Average Regret Ratio: %f\n", minArr);
		fprintf(fp, "Sample sd: %f\n",sd);
		fprintf(fp, "No within SD: %d\n", notWithinSD);
		fprintf(fp, "70\%: %f\n", user_rr[(int)(samplingSize*0.7)]);
		fprintf(fp, "80\%: %f\n", user_rr[(int)(samplingSize*0.8)]);
		fprintf(fp, "90\%: %f\n", user_rr[(int)(samplingSize*0.9)]);
		fprintf(fp, "95\%: %f\n", user_rr[(int)(samplingSize*0.95)]);
		fprintf(fp, "99\%: %f\n", user_rr[(int)(samplingSize*0.99)]);
		fprintf(fp, "Sample mrr: %f\n", sampleMRR);
		fprintf(fp, "Average number of points avoided: %f\n", avgIterAvoided);
		fprintf(fp, "Average number of users avoided: %f\n", avgUsersAvoided);
		fprintf(fp, "\n\n");
		fclose(fp);
	}
}
