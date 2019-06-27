/* Engineering Computation Assignment 2
   Kelvin Liao			16/05/2018			*/
		
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define Y_AXIS_LIMIT 60
#define X_AXIS_LIMIT 70
#define MAX_ELS	52

#define NUM_OF_HEADER_LINES 1
#define NUM_OF_DATA_PER_LINE 5
#define MM_TO_LITERS 1e-6 
#define MULTIPLES_OF_TEN 10
#define MINIMUM_STRESS_TO_DIE 1

/* define it is -1 as the values from 0 to the maximum number of trees are 
	possible whereas -1 number of trees dead is impossible */
#define NO_TREES_DEAD -1			

#define STAGE_3 3
#define STAGE_4 4

typedef struct {
	char	label;
	double	xloc;
	double	yloc;
	int		liters;
	double	rootrad;
}	tree_t;

/* Function declarations organized according to stages */
/**************************** Stage 1 *******************************/
int		mygetchar(void);
void	scan_file(tree_t tree_data[], int *num_data_lines);
double	sum(tree_t tree_data[], int num_data_lines);
void	print_stage_1(tree_t tree_data[], int num_data_lines);
/********************************************************************/



/**************************** Stage 2 *******************************/
void 	print_stage_2(tree_t tree_data[], int num_data_lines);
void	trees_in_conflict(tree_t tree_data[], int num_data_lines, 
						  int current_tree); 
double	calc_distance(double x1, double x2, double y1, double y2);
/********************************************************************/



/**************************** Stage 3 *******************************/
void	print_graph_with_axes(tree_t tree_data[], int num_data_lines, 
							  int stage);
void 	print_x_axis(int stage);
void	loop_x_axis(tree_t tree_data[], int num_data_lines, int ctr_pt_y);
int 	tree_found_in_point(tree_t tree_data[], double ctr_pt_y, 
							double ctr_pt_x, int current_tree);
void	print_tree(tree_t tree_data[], int num_data_lines, 
				   double ctr_pt_x, double ctr_pt_y, int found_tree_index);
void	print_y_axis(int y, int stage);
void	print_graph(tree_t tree_data[], int num_data_lines, 
					double ctr_pt_x, double ctr_pt_y, int found_tree_index, 
					int found_tree);
int		find_tree_on_pt(tree_t tree_data[], int num_data_lines,
				double ctr_pt_y, double ctr_pt_x, int *found_tree_index);
int		compare_tree_distances(tree_t tree_data[], int num_data_lines, 
							   double ctr_pt_x, double ctr_pt_y, 
							   double tree_on_pt, int found_tree_index);
void	x_axis_tags(int stage);
void	y_axis_scale(int stage, int y);
/********************************************************************/



/**************************** Stage 4 *******************************/
void	read_annual_rainfall(int argc, char* argv[], double *annual_rainfall);
double	calc_water_for_survival(tree_t tree_data[], int current_tree, 
								int space_count[]);
double	stress_factor(tree_t tree_data[], int current_tree, 
					  int space_count[], double annual_rainfall);
void	compute_stress_factor_for_all_trees(tree_t tree_data[], 
											int remaining_trees, 
											int space_count[], 
											double annual_rainfall, 
											double stress_factors[]);
void	check_dying_trees(tree_t tree_data[], int num_data_lines, 
						  int space_count[], double stress_factors[], 
						  double annual_rainfall, int *remaining_trees);
void	kill_trees(tree_t tree_data[], int tree_to_die, int *remaining_trees, 
		   		   double highest_stress_factor);
void	initialize_array(int A[], int n);
void	count_spaces_occupied_by_trees(tree_t tree_data[], int num_data_lines,
							   		   double ctr_pt_y, double ctr_pt_x, 
							   		   int space_count[]);
void	tree_occupied_space(int space_count[], tree_t tree_data[], 
							int num_data_lines);
void	print_stage_4(tree_t tree_data[], int num_data_lines, 
					  int space_count[], double stress_factors[], 
					  double annual_rainfall,  int stage);
/********************************************************************/


/* Main Function */
int 
main (int argc, char* argv[]) {
	
	/* declarations */
	tree_t	tree_data[MAX_ELS];
	int 	num_data_lines = 0;
	
	/* scan the data into the struct named tree_data */
	scan_file(tree_data, &num_data_lines);
	
	/* print stage 1 */
	print_stage_1(tree_data, num_data_lines);
	
	/* print stage 2 */
	print_stage_2(tree_data, num_data_lines);
	
	/* print stage 3 */
	print_graph_with_axes(tree_data, num_data_lines, STAGE_3);
	
	double annual_rainfall = 0.0;
	/* if there is an input from the command line, read it */
	read_annual_rainfall(argc, argv, &annual_rainfall);
	
	/* space_count stores the number of spaces occupied by each tree */
	int space_count[MAX_ELS] = {0};
	
	/* stress_factors stores the values of the stress factor for each tree */
	double stress_factors[MAX_ELS] = {0};
	
	/* prints stage 4 */
	print_stage_4(tree_data, num_data_lines, space_count, 
				  stress_factors, annual_rainfall, STAGE_4);
	
	/* Finish! */
	return 0;	
}



/********************************************************************/
/* Code by Professor Alistair Moffat to accommodate for both 
   windows and Linux computers */
int	
mygetchar(void) {
	int c;
	while ((c=getchar())=='\r') {
	}
	return c;
}
/********************************************************************/



/********************************************************************/
/* Code written by Kelvin Liao
	which is found from assignment 1, Engineering Computation */

/* scans the file into 5 arrays */
void
scan_file(tree_t tree_data[], int *num_data_lines) {
	/* scan 5 values from 1 line */
	int nvals_read = NUM_OF_DATA_PER_LINE;

	int heading_lines = 0;
	int c;
	
	/* 	Look for the new line character and after one new line character,
		the rest of the file will be the desired data */
	while ((c = mygetchar()) != EOF) {
		if (c == '\n') {
			heading_lines ++;	
		}
		
		if (heading_lines == NUM_OF_HEADER_LINES) {
			break;
		}
	}

	while (nvals_read != EOF) {
		nvals_read = scanf("%c%lf%lf%d%lf", 
			&tree_data[*num_data_lines].label,
			&tree_data[*num_data_lines].xloc,
			&tree_data[*num_data_lines].yloc,
			&tree_data[*num_data_lines].liters,
			&tree_data[*num_data_lines].rootrad);
		
		/* count the number of rows */
		if (nvals_read == NUM_OF_DATA_PER_LINE) {
			(*num_data_lines)++;
		}
	}
} 
/********************************************************************/
		


/********************************************************************/
/* counts the total liters of water required for a year */
double 
sum(tree_t tree_data[], int num_data_lines) {
	int i;
	double sum = 0;
	
	for (i=0; i < num_data_lines; i++) {
		sum = sum + tree_data[i].liters;
	}
	return sum*(MM_TO_LITERS);
}
/********************************************************************/



/********************************************************************/
/* print stage 1 */
void
print_stage_1(tree_t tree_data[], int num_data_lines) {
	printf("\n");
	
	printf("S1: total data lines   = %5d trees\n", num_data_lines);
	
	printf("S1: total water needed = %4.3f megaliters per year\n", 
		   sum(tree_data, num_data_lines));
	
	printf("\n");
}
/********************************************************************/



/********************************************************************/
/* calculates the distance between two points */
double
calc_distance(double x1, double x2, double y1, double y2) {
	return sqrt((x1 - x2)*(x1 - x2) + (y1 - y2)*(y1 - y2));
}
/********************************************************************/



/********************************************************************/
/* checks the trees that are in conflict with the current tree */
void
trees_in_conflict(tree_t tree_data[], int num_data_lines, int current_tree) {
	int next;
	
	/* loops through all the trees*/
	for (next=0; next < num_data_lines; next++) {
		/* if we are comparing the same tree, skip it */
		if (next == current_tree) {
			next++;	
		}
		
		/* compares the distances between two trees and sees if they are 
			sharing any spots */
		if (calc_distance(tree_data[current_tree].xloc, tree_data[next].xloc, 
						  tree_data[current_tree].yloc, tree_data[next].yloc) 
			<=
			(tree_data[current_tree].rootrad + tree_data[next].rootrad)) {
			printf(" ");
			printf("%c", tree_data[next].label);
			}
	}
}
/********************************************************************/



/********************************************************************/
/* prints stage 2 */
void
print_stage_2(tree_t tree_data[], int num_data_lines) {
	int i;

	/* loops through all the trees */
	for (i=0; i < num_data_lines; i++) {
		printf("S2: tree %c is in conflict with", tree_data[i].label);
		
		/* prints the trees that are in conflict with the current tree i */
		trees_in_conflict(tree_data, num_data_lines, i);
		printf("\n");
	}
	printf("\n");
}
/********************************************************************/



/********************************************************************/
/* prints the y axis according to whether it is in stage 3 or stage 4 */	
void
print_y_axis(int y, int stage) {
	
	if (y == Y_AXIS_LIMIT) {
		y_axis_scale(stage, y);
	}
	else if ((y%MULTIPLES_OF_TEN == 0) || (y==0)) {
		y_axis_scale(stage, y);
	} 
	else {
		if (stage == STAGE_3) {
			printf("S3:    |");	
		} else {
			printf("S4:    |");	
		}
	}
}
/********************************************************************/



/********************************************************************/
/* prints the y axis scales as stage 3 or stage 4 */
void
y_axis_scale(int stage, int y) {
	if (stage == STAGE_3) {
		printf("S3: %2d +", y);	
	} else {
		printf("S4: %2d +", y);	
	}
}
/********************************************************************/



/********************************************************************/
/* prints the x axis according to whether it is stage 3 or stage 4 */
void 
print_x_axis(int stage) {
	int x;
	int next_point;
	x_axis_tags(stage);
	
	/* prints the "+----------+" sequence */
	for (x=0; x <= X_AXIS_LIMIT; x++) {
		if ( (x==0) || (x%MULTIPLES_OF_TEN == 0)) {
			printf("+");	
		} else {
			printf("-");	
		}
	}
	
	printf("\n");
	x_axis_tags(stage);
	
	/* prints the x label scales i.e. "0    10    20" and so on */
	for (x=0; x <= X_AXIS_LIMIT; x++) {
		next_point = x + 1;
		
		if (x==0) {
			printf("%-d", x);	
		}
		/* if the next point is a multiple of 10, print 10 to align the number
		   with the previous +------+ sequence */
		else if (next_point%MULTIPLES_OF_TEN == 0) {
			printf("%-d", next_point);	
		} 
		/* if the next point is not a multiple of 10, print spaces */
		else if (x%MULTIPLES_OF_TEN != 0) {
			printf(" ");
		}
	}
	printf("\n");
}
/********************************************************************/



/********************************************************************/
/* prints the x axis tags as stage 3 or stage 4 */
void
x_axis_tags(int stage) {
	if (stage == STAGE_3) {
		printf("S3:     ");
	} else {
		printf("S4:     ");
	}
}
/********************************************************************/



/********************************************************************/
/* check if any tree can reach the current point of the graph */	
int 
tree_found_in_point(tree_t tree_data[], double ctr_pt_y, 
					double ctr_pt_x, int current_tree) {

	/* calculate the distance between the center of the tree and the point */
	double distance = calc_distance(ctr_pt_x,tree_data[current_tree].xloc,
									ctr_pt_y,tree_data[current_tree].yloc);
	
	/* if the distance is smaller than the tree's radius, then it is in the 
		catchment zone */
	if (distance <= tree_data[current_tree].rootrad) {
		return 1;	
	}
	return 0;
}
/********************************************************************/



/********************************************************************/
/* If a tree is found, it returns 1 and the index of the tree via pointer
	If a tree is not found, it returns 0 */
int
find_tree_on_pt(tree_t tree_data[], int num_data_lines,
				double ctr_pt_y, double ctr_pt_x, int *found_tree_index) {
	int found_tree = 0;
	int i;	
		
	/* loop_trees */
	for (i=0; i < num_data_lines; i++) {
		found_tree = tree_found_in_point(tree_data, ctr_pt_y, 
											ctr_pt_x, i);
		
		/* if a tree is found, get the index of the tree */
		if (found_tree == 1) {
			*found_tree_index = i;
			break;
		}
	}
	
	return found_tree;
}
/********************************************************************/



/********************************************************************/
/* Prints a tree at the point when called */
void
print_tree(tree_t tree_data[], int num_data_lines, 
			double ctr_pt_x, double ctr_pt_y, int found_tree_index) {
	
	/* intialize distance */
	double tree_on_pt = calc_distance(ctr_pt_x,
										  tree_data[found_tree_index].xloc,
										  ctr_pt_y,
										  tree_data[found_tree_index].yloc);
	
	/* Find all the trees that exist at that point and the 
	   lowest distance wins the spot */
	found_tree_index = compare_tree_distances(tree_data, num_data_lines,
							ctr_pt_x, ctr_pt_y, 
							tree_on_pt, found_tree_index);
	
	/* prints the tree label */
	printf("%c", tree_data[found_tree_index].label);
}
/********************************************************************/



/********************************************************************/
/* compare the distances between all trees and the point to determine 
   which tree gets the spot if there is an overlap */
int
compare_tree_distances(tree_t tree_data[], int num_data_lines,
						double ctr_pt_x, 
						   double ctr_pt_y, double tree_on_pt, 
						   int found_tree_index) { 
	double other_trees;
	int i;
	
	/* calculates the distance between the tree and the point */
	tree_on_pt = calc_distance(ctr_pt_x, tree_data[found_tree_index].xloc,
							   ctr_pt_y, tree_data[found_tree_index].yloc);
	
	/* loop trees */
	for (i=0; i < num_data_lines; i++) {
		
		/* get the distance from other trees to this point */
		other_trees = calc_distance(ctr_pt_x,tree_data[i].xloc,
											  ctr_pt_y,tree_data[i].yloc);
		
		/* if a tree cannot reach the point, skip it*/
		if (other_trees > tree_data[i].rootrad) {
			continue;
		}
		
		/* if a tree can reach the point AND it is closer to the point*/
		if (other_trees < tree_on_pt) {
			/* change the index to change which tree gets printed */
			found_tree_index = i;
			tree_on_pt = other_trees;
		}
	}	
	return found_tree_index;				   
}
/********************************************************************/



/********************************************************************/
/* prints the whole graph with the trees and spaces */
void
print_graph(tree_t tree_data[], int num_data_lines, 
			double ctr_pt_x, double ctr_pt_y, int found_tree_index, 
			int found_tree) {

	/* there is a tree! */
	if (found_tree) {
		print_tree(tree_data, num_data_lines, ctr_pt_x, ctr_pt_y,
					found_tree_index); 
	}
	
	/* print a blank if there is no tree */
	else {
		printf(" ");	
	}
}
/********************************************************************/



/********************************************************************/
/* loops across the x axis and prints the trees or spaces */
void
loop_x_axis(tree_t tree_data[], int num_data_lines, int ctr_pt_y) {
	
	int found_tree_index = 0;
	int found_tree, x;
	double ctr_pt_x; /* center point x axis value */
	
	/* loop x axis */
	for (x=0; x < X_AXIS_LIMIT; x++) {
		found_tree = 0;
		/* add 0.5 as we want the center of the point */
		ctr_pt_x = x + 0.5;
		
		found_tree = find_tree_on_pt(tree_data, num_data_lines, ctr_pt_y, 
									 ctr_pt_x, &found_tree_index);
		
		print_graph(tree_data, num_data_lines, ctr_pt_x, ctr_pt_y,
					found_tree_index, found_tree);
	}
}
/********************************************************************/




/********************************************************************/
void
print_graph_with_axes(tree_t tree_data[], int num_data_lines, int stage) {
	int y;
	double ctr_pt_y; /* center point y axis value */
	
	/* loop y axis, move down by -2 as one space on the terminal 
		corresponds to 2 meters*/
	for (y= Y_AXIS_LIMIT; y >= 0; y=y-2) {
		print_y_axis(y, stage);
		/* add 1 as we want the center of the point */
		ctr_pt_y = y + 1;
		
		loop_x_axis(tree_data, num_data_lines, ctr_pt_y);
		
		/* next line on y axis*/
		printf("\n");
	}
	print_x_axis(stage);
}
/********************************************************************/



/***************************************************************************/
/* prints the new graph for stage 4 with the 
   annual rainfall taken into account 			*/
void
print_stage_4(tree_t tree_data[], int num_data_lines, int space_count[],
			  double stress_factors[], double annual_rainfall,  int stage) {

	int remaining_trees = num_data_lines;
	
	check_dying_trees(tree_data, num_data_lines, space_count, 
					  stress_factors, annual_rainfall, &remaining_trees);
	
	printf("\n");
	print_graph_with_axes(tree_data, remaining_trees, STAGE_4);
	printf("\n");
}
/***************************************************************************/



/***************************************************************************/
/* reads the annual rainfall from the commandline if it is provided */
void
read_annual_rainfall(int argc, char* argv[], double *annual_rainfall) {
	/* if the annual rainfall is provided in the commandline, read it */
	if (argc > 1) {
		*annual_rainfall = atof(argv[1]);	
	} 
	
	/* print the annual rainfall value that was read */
	printf("\n");
	printf("S4: rainfall amount = %4.1f\n", *annual_rainfall);
}
/***************************************************************************/



/***************************************************************************/
/* calculates the amount of water in mm needed by a tree to survive */
double
calc_water_for_survival(tree_t tree_data[], int current_tree, 
						int space_count[]) {

	double liters_per_year = tree_data[current_tree].liters;
	double occupied_space_by_tree = space_count[current_tree];
	
	/* multiplied by 2 as each cell space on the graph corresponds to 
		2 square meters */
	occupied_space_by_tree = occupied_space_by_tree*2; 
	
	return liters_per_year/occupied_space_by_tree;	
}
/***************************************************************************/




/***************************************************************************/
/* computes the stress factor for the given tree */
double
stress_factor(tree_t tree_data[], int current_tree, 
			  int space_count[], double annual_rainfall) {
	
	double required_for_survival_rainfall = 
			calc_water_for_survival(tree_data,current_tree, space_count);
	
	return (required_for_survival_rainfall/annual_rainfall);
}
/***************************************************************************/



/***************************************************************************/
/* computes the stress factors of each tree */
void
compute_stress_factor_for_all_trees(tree_t tree_data[], int remaining_trees, 
									int space_count[], double annual_rainfall, 
									double stress_factors[]) {
	int i;
	for (i=0; i < remaining_trees; i++) {
		tree_occupied_space(space_count, tree_data, remaining_trees);
		
		stress_factors[i] = stress_factor(tree_data, i, 
										  space_count, annual_rainfall);
	}
}
/***************************************************************************/



/***************************************************************************/
/* removes the tree from tree_data effectively "killing it" */
void
kill_trees(tree_t tree_data[], int tree_to_die, int *remaining_trees, 
		   double highest_stress_factor) {
	
	int n;
	/* if there exists a tree that is going to die, kill it */
	if (tree_to_die != NO_TREES_DEAD) {
		printf("S4: tree %c has stress factor %3.2f and dies next\n",
				tree_data[tree_to_die].label, highest_stress_factor);
		
		/* no need to iterate all as 1 will be taken out anyway*/
		/* delete the tree that will die */
		for (n = tree_to_die; n < *remaining_trees-1; n++) {
			tree_data[n].label = tree_data[n+1].label;
			tree_data[n].xloc = tree_data[n+1].xloc;
			tree_data[n].yloc = tree_data[n+1].yloc;
			tree_data[n].liters = tree_data[n+1].liters;
			tree_data[n].rootrad = tree_data[n+1].rootrad;
		}
		/* decrease as one tree has been removed */
		(*remaining_trees)--;
	}
}
/***************************************************************************/



/***************************************************************************/
/* looks for the trees that would die depending on the annual rainfall */
void
check_dying_trees(tree_t tree_data[], int num_data_lines, int space_count[],
				  double stress_factors[], double annual_rainfall, 
				  int *remaining_trees) {
		
	int m,k, tree_to_die;
	double highest_stress_factor;
	
	/* iterate every tree as one may die allowing others to live */
	for (m=0; m < num_data_lines; m++) {
		tree_to_die = -1;
		
		/* get all the stress factors of each tree */
		compute_stress_factor_for_all_trees(tree_data, *remaining_trees, 
											space_count, annual_rainfall, 
											stress_factors);
		
		/* initialize variable */
		highest_stress_factor = MINIMUM_STRESS_TO_DIE;
		
		/* look for the tree with the highest stress factor */
		for (k=0; k < *remaining_trees; k++) {
			/* if the stress factor of this tree is greater than 1 AND 
			   the stress factor is greater than the other trees */
			if ( (stress_factors[k] > MINIMUM_STRESS_TO_DIE) && 
				 (stress_factors[k] > highest_stress_factor)) {
			
				highest_stress_factor = stress_factors[k];
				tree_to_die = k;
			}
		}
		/* kill the tree with the highest stress factor */
		kill_trees(tree_data, tree_to_die, remaining_trees, 
				   highest_stress_factor);
	}
}
/***************************************************************************/



/***************************************************************************/
/* initializes all values in an array to 0 */
void
initialize_array(int A[], int n) {
	int i;
	/* initialize each time the function is called */
	for (i=0; i < n; i++) {
		A[i] = 0;
	}
}
/***************************************************************************/



/***************************************************************************/
/* determines which tree is at the point specified and record the space
	occupied by that tree in the array "space_count" */
void
count_spaces_occupied_by_trees(tree_t tree_data[], int num_data_lines,
							   double ctr_pt_y, double ctr_pt_x, 
							   int space_count[]) {

	int found_tree_index = 0;
	int found_tree;
	double tree_on_pt;
	
	/* found_tree_index corresponds to the first tree found on the point */
	found_tree = find_tree_on_pt(tree_data, num_data_lines, ctr_pt_y, 
								 ctr_pt_x, &found_tree_index);
	
	/* if you found a tree, the found_tree_index would be that tree */
	if (found_tree == 1) {
				
		/* first tree on that point*/
		tree_on_pt = calc_distance(ctr_pt_x,
								  tree_data[found_tree_index].xloc,
								  ctr_pt_y,
								  tree_data[found_tree_index].yloc);

		/*Compare with other trees and the winning tree gets the spot 
		  If only 1 tree on that spot, found_tree_index does not change */
		found_tree_index = compare_tree_distances(tree_data, 
										num_data_lines,
										ctr_pt_x, ctr_pt_y, 
										tree_on_pt, found_tree_index);
		
		/* counts the spaces occupied for the tree on that spot */
		(space_count[found_tree_index]) ++;
	}
}
/***************************************************************************/



/***************************************************************************/
/* goes through the entire graph and counts the space occupied by each tree */
void
tree_occupied_space(int space_count[], tree_t tree_data[], 
					int num_data_lines) {

	int x,y;
	double ctr_pt_y, ctr_pt_x; 
	
	initialize_array(space_count, MAX_ELS);
	
	/* loop y axis */
	for (y= Y_AXIS_LIMIT; y >= 0; y=y-2) {
		ctr_pt_y = y + 1;
		
		/* loop x axis */
		for (x=0; x < X_AXIS_LIMIT; x++) {
			ctr_pt_x = x + 0.5;
			
			/* determines which tree occupies the space and 
				records the space occupied by that tree		*/
			count_spaces_occupied_by_trees(tree_data, num_data_lines,
										   ctr_pt_y, ctr_pt_x, space_count);
		}
	}	
}
/***************************************************************************/



/* ======================= Programming is fun! ========================*/