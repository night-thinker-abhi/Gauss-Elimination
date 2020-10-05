#include <sys/time.h>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include "math.h"


#define SUCCESS 0
#define ERROR -1

#define EPSILON 0.000001

void read_matrix_size_from_file(char *filename, int *rows, int *columns);
double ** read_user_matrix_from_file(char *filename, int rows, int columns, int rank, int nprocs);
double ** allocate_matrix(int, int);
void free_matrix(double **matrix, int rows);
void divide_by_max(double **, int, int);
void input_clicking_probabilities(double **, int, int, double *);
void write_clicking_probabilities_to_file(double *cp, int rows);
void print_best_acceptance_threshold(double *, int);
void print_matrix(double **, int, int);

int main (int argc, char** argv)
{
    if (argc != 2)
    {
        printf("please provide a user matrix!\n");
        return ERROR;
    }

    
    int rank, nprocs;
    MPI_Init( &argc, &argv);
    MPI_Comm_rank( MPI_COMM_WORLD, &rank);
    MPI_Comm_size( MPI_COMM_WORLD, &nprocs);
    clock_t begin1, end1, begin2, end2;

    
    int rows, columns;
    double **A;
    double *cp;
    if (rank == 0)
    {
        printf("malloc cp process %d\n", rank); 
        read_matrix_size_from_file(argv[1], &rows, &columns);
        cp = malloc(columns * sizeof(double));
    }
   

    MPI_Bcast(&rows, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&columns, 1, MPI_INT, 0, MPI_COMM_WORLD);


   
    //int rpp = rows / rank;
    //int rrpp = rows % rank;
    printf("WE reached here process %d\n", rank);

    int i,j,k;
    int *map = malloc(sizeof(int) * rows);

    int stop = 0;
    while(1)
    {
        for(i=0;i<nprocs; i++)
        {
            if (rank == i)
            {
                printf("Reading for process %d\n", rank);
                A = read_user_matrix_from_file(argv[1], rows, columns, rank, nprocs);
                stop = 1;
            }   
        }
        if (stop)
        {
            break;
        }
    }




    MPI_Barrier(MPI_COMM_WORLD);

    for (i = 0; i < rows; i++)
    {
        map[i] = i % nprocs;
    }
    printf("Process %d is here\n", rank);

    for(k=0; k<rows; k++)
    {
        //printf("Broadcast recieved/sent in process %d from %d\n", rank, map[k]);
        MPI_Bcast (&A[k][k],rows-k,MPI_DOUBLE,map[k],MPI_COMM_WORLD);
        printf("process %d, broadcast b[%d] as %G\n", rank,k,A[k][columns-1] );
        MPI_Bcast (&A[k][columns-1],1,MPI_DOUBLE,map[k],MPI_COMM_WORLD);
        double pivot;
        for(i= k+1; i<rows; i++) 
    {
      if(map[i] == rank)
      {
          pivot = A[i][k]/A[k][k];
      }
    }               
    for(i= k+1; i<rows; i++) 
    {       
      if(map[i] == rank)
      {
        for(j=0;j<rows;j++)
        {
            A[i][j]=A[i][j]-( pivot * A[k][j] );
        }
        A[i][columns-1]= A[i][columns-1]-( pivot * A[k][columns-1] );
        printf("b[%d] is %G\n", i, A[i][columns-1]);
      }
    }
    }

    MPI_Barrier(MPI_COMM_WORLD);
    printf("Process %d finished\n", rank);

    if (rank == 0)
    {
            
        int row,row2;
        for (row = rows-1; row >= 0; row--) {
            A[row][columns-1] = A[row][columns-1] / A[row][row];
            printf("divide %G by %G\n",A[row][columns-1] , A[row][row]);
            A[row][row] = 1;
            for (row2 = row-1; row2 >= 0; row2--) {
                A[row2][columns-1] += A[row2][row]*A[row][columns-1];
                A[row2][row] = 0;
            }
        }
    }

    if (rank == 0)
    {
        print_matrix(A, rows, columns);
        divide_by_max(A, rows, columns);
        
        input_clicking_probabilities(A, rows, columns, cp);
        print_best_acceptance_threshold(cp, rows);
        write_clicking_probabilities_to_file(cp, rows);
        printf("Free cp process %d\n", rank);
        free(cp);
        printf("Freed cp\n");
    }

  
    free_matrix(A, rows);
    MPI_Finalize();
    return SUCCESS;

}


void read_matrix_size_from_file(char *filename, int *rows, int *columns)
{
    FILE *file;
    file = fopen(filename, "r");

  
    *rows = 1;
    *columns = 1;
    char c;
    int columns_known = 0;
    while(!feof(file)) {
        c = fgetc(file);
        if (c == ' ') {
            if (!columns_known) (*columns)++;
        } 

        if (c == '\n') {
            (*rows)++;
            columns_known = 1;
            continue;
        }
    }

    printf("There are %d rows and %d columns\n", *rows, *columns);

    fclose(file);
}

double ** read_user_matrix_from_file(char *filename, int rows, int columns, int rank, int nprocs)
{
    FILE *file;
    file = fopen(filename, "r");

   
    //rewind(file);
    printf("Rank is %d\n", rank);
    printf("Nprocess is %d\n", nprocs);
    printf("Rows is %d and columns is %d\n", rows, columns);
    double **matrix = allocate_matrix(rows, columns);
    int i,j;
    for (i = 0; i < rows; i++)
    {
        printf("process %d, %d mod %d = %d\n", world_rank, i, nprocs, i% nprocs);
        if (world_rank == i % nprocs)
        {
            for (j = 0; j < columns; j++)
            {
                fscanf(file, "%lf", &matrix[i][j]);
            }
        }else {
            for ( j= 0; j < columns; j++)
            {
                fscanf(file, "%lf", &buf[j]);
                matrix[i][j] = 0;
            }
        }

    }

    return matrix;
}

double ** allocate_matrix(int rows, int columns)
{
    double ** matrix = (double **) malloc(sizeof(double *) * rows);
    int i;
    for (i = 0; i < rows; i++)
    {
        matrix[i] = (double *) malloc(sizeof(double) * columns);
    }

    return matrix;
}

void free_matrix(double **matrix, int rows)
{
    int i;
    for (i = 0; i < rows; i++)
    {
        free(matrix[i]);
    }
    free(matrix);
}

void input_clicking_probabilities(double **matrix, int rows, int columns, double *cp) {
    int row;
    for (row = 0; row < rows; row++) {
        cp[row] = matrix[row][columns-1];
    }
}

void write_clicking_probabilities_to_file(double *cp, int rows)
{
   
    FILE *output_file;
    int row;
    output_file = fopen("clicking_probabilities.txt","w");
    for (row = 0; row < rows; row++) {
        fprintf(output_file, "%lf\n", cp[row]);
    }

    fclose(output_file);
}

void print_matrix(double **matrix, int rows, int columns)
{
    FILE *output_file;
    int row, column;
    output_file = fopen("row_reduced_matrix.txt","w");
    for (row = 0; row < rows; row++) {
        for (column = 0; column < columns; column++) {
            fprintf(output_file,"%lf ",matrix[row][column]);
        }
        fprintf(output_file,"\n");
    }   
    fclose(output_file);
}

void print_best_acceptance_threshold(double *cp, int rows) {

}

void divide_by_max(double **matrix, int rows, int columns) {
    double max = 0; 
    int row, column;

  
    for (row = 0; row < rows; row++) {      
        if (max < fabs(matrix[row][columns-1])) max = fabs(matrix[row][columns-1]);
    }

    for (row = 0; row < rows; row++) {
       
        if (equals(max,0)) {
            matrix[row][columns-1] = 0;
        } else {
            matrix[row][columns-1] = fabs (matrix[row][columns-1]) / max;
        }
    }
}

int equals(double a, double b) {
    if (fabs(a-b) < EPSILON) return 1;
    else return 0;
}