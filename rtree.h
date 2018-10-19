#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <math.h>
#ifndef WIN32
#include <sys/resource.h>
#include <sys/times.h>
#include <unistd.h>
#endif
#include "qsort.h"

#define DEBUG 1


#define CONFIG_FILE	"rtree.config"
#define SAVE_RTREE_FILE "save_rtree_file"

#define FALSE    	0
#define TRUE     	1

#define RANGE_SEARCH 0
#define kNN_SEARCH 1
#define CHOICE kNN_SEARCH

#define ASC_NUM  	48
#define NO_ID	 	-1
#define FOUND		1
#define NOT_FOUND 	0

#define ROOT  0   
#define LEAF  1
#define NODE  2

//#define INFINITY  FLT_MAX
#define INFINITY  1E+37 //FLT_MAX
#define UNDEFINED -3  // for id of entries in PR

#define R_FLOAT
//#define R_TYPE int
#define R_TYPE float

/* Global variable ******************
m: min. number entries of each node;
M: max. number entries of each node;
dim: dimension of the rtree data.
*************************************/

typedef struct rtree_info_s
{
	int m, M, dim, reinsert_p, no_histogram;
	int extra_level;
} rtree_info;

typedef struct node { 
	R_TYPE *a;
	R_TYPE *b;
	int id;
	int attribute;
	int vacancy;
	struct node *parent;
	struct node **ptr; 
}   node_type;
	
typedef struct NN {	
		double dist;
		int oid;
		struct node *pointer; 
		int level;
		struct NN *next; 
} NN_type;
		
typedef struct BranchArray {	
	double min;
	node_type *node;
} ABL;


typedef struct config { 
	int dim;
	int m;
	int M;
	int reinsert_p;
	int no_histogram;
	//char nodefile[FILENAME_MAX];
	//char rootfile[FILENAME_MAX];
	char queryfile[FILENAME_MAX];
	char positionfile[FILENAME_MAX]; 
	char save_tree_file[FILENAME_MAX]; 
}   config_type;



struct setNode_s
{
	int noOfNode;
	node_type **elt;
	int *level;
};

typedef struct setNode_s setNode;


void overflow(node_type *over_node, int over_level, int old_level, node_type
			  *extra_node, node_type *root, rtree_info *aInfo);

void initialize(config_type *config, char *configFile, rtree_info *aInfo);
void tree_node_allocate(node_type **node, rtree_info *aInfo);
void NN_freeChain(NN_type *aNN);

void k_NN_search(node_type *root, R_TYPE *query, int k, NN_type **returnResult, rtree_info *aInfo);
