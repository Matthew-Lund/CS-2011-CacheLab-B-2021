/* 
 * trans.c - Matrix transpose B = A^T
 *Team: mtlund
 Member(s): Matthew Lund
 Email: mtlund@wpi.edu
 * Each transpose function must have a prototype of the form:
 * void trans(int M, int N, int A[N][M], int B[M][N]);
 *
 * A transpose function is evaluated by counting the number of misses
 * on a 1KB direct mapped cache with a block size of 32 bytes.
 */ 
#include <stdio.h>
#include "cachelab.h"

int is_transpose(int M, int N, int A[N][M], int B[M][N]);

/* 
 * transpose_submit - This is the solution transpose function that you
 *     will be graded on for Part B of the assignment. Do not change
 *     the description string "Transpose submission", as the driver
 *     searches for that string to identify the transpose function to
 *     be graded. 
 */
char transpose_submit_desc[] = "Transpose submission";
void transpose_submit(int M, int N, int A[N][M], int B[M][N])
{
	//Variables Used
	int chunkSize; //var. for size of each chunk/block; Used in every iteration
	int chunkRow , chunkCol;//iterating over chunks, in outer loops
	int r , c; //iterating through each chunk, used primarily in inner loops
	int tmp = 0, d = 0; //d = diagonal and tmp for temporary values
	int v0, v1, v2, v3, v4;//Variables used in the N == 64 case
	
	if(N == 32)
	{
	//Using chunkSize = 8 in this case. Only N == 32 in condtion because
	// matrix transpose can only happen for any a*b and c*a where a needs
	// to be the same
	//Blocking is utilized here; 4 levels of loops used. 2 outer loops for across chunks and 2 inner loops for in each chunk
		chunkSize = 8;
		for(chunkCol = 0; chunkCol < N; chunkCol += 8)
		{
			for(chunkRow = 0; chunkRow < M; chunkRow += 8)
			{
				for(r = chunkRow; r < chunkRow + 8; r++)
				{
					for(c = chunkCol; c < chunkCol + 8; c++)
					{
						//if rows and columns not equal
						if (r != c)
						{
							B[c][r] = A[r][c];
						}
						//if they are equal
						else
						{
						
							tmp = A[r][c];
							d = r;
						}
					}
					//Don't move elements on diags. because transposing a square matrix
					if(chunkRow == chunkCol)
					{
						B[d][d] = tmp;
					}
				}
			}
		}
	}
	
	
	else if (N == 64)
	{
	//Using chunkSize = 4; 2 levels of loops used
	//Assigning elements in each row individually; Reduces misses
		chunkSize = 4;
		for(r = 0; r < N; r+= chunkSize)
		{
			for(c = 0; c < M; c += chunkSize)
			{
				//Elements in A[r][], A[r+2][] and A[r+1][] assigned to vars for use in loop
				//This is because we can only modify matrix B and not matrix A.
				
				v0 = A[r][c];
				v1 = A[r+1][c];
				v2 = A[r+2][c];
				v3 = A[r+2][c+1];
				v4 = A[r+2][c+2];
				//B[c+3][] assigned
				B[c+3][r] = A[r][c+3];
				B[c+3][r+1] = A[r+1][c+3];
				B[c+3][r+2] = A[r+2][c+3];
				//B[c+2][] assigned
				B[c+2][r] = A[r][c+2];
				B[c+2][r+1] = A[r+1][c+2];
				B[c+2][r+2] = v4;
				v4 = A[r+1][c+1];
				//B[c+1][] assigned
				B[c+1][r] = A[r][c+1];
				B[c+1][r+1] = v4;
				B[c+1][r+2] = v3;
				//B[c][] assigned
				B[c][r] = v0;
				B[c][r+1] = v1;
				B[c][r+2] = v2;
				//Row A[r+3][] assigned to left-most elements in B
				B[c][r+3] = A[r+3][c];
				B[c+1][r+3] = A[r+3][c+1];
				B[c+2][r+3] = A[r+3][c+2];
				v0 = A[r+3][c+3];
				//B[c+3][] assigned
				B[c+3][r+3] = v0;
			}
		}
	}
	
	
	else
	{
	//Case for random matrix size; Use chunkSize = 16
	//2 levels of loops to iterate over chunks in column major iteration and 2 levels to go through chunks.
		chunkSize = 16;
		for(chunkCol = 0; chunkCol < M; chunkCol += chunkSize)
		{
			for(chunkRow = 0; chunkRow < N; chunkRow += chunkSize)
			{
			//B/c sizes can be odd, not all chunks are square. Certain case: if (chunkRow + 16 > N) is invalid access
			//Also do regular check for i < N and j < M
				for(r = chunkRow; (r < N) && (r < chunkRow + chunkSize); r++)
				{
					for(c = chunkCol; (c < M) && (c < chunkCol + chunkSize); c++)
					{
						//If rows != columns
						if(r != c)
						{
							B[c][r] = A[r][c];
						}
						//If they are equal
						else
						{
							tmp = A[r][c];
							d = r;
						}	
					}
					//Row and Col. # are same in chunks , assign diag. elements
					if(chunkRow == chunkCol)
					{
						B[d][d] = tmp;					
					}
				}
			}
		}
	}
}

/* 
 * You can define additional transpose functions below. We've defined
 * a simple one below to help you get started. 
 */ 

/* 
 * trans - A simple baseline transpose function, not optimized for the cache.
 */
char trans_desc[] = "Simple row-wise scan transpose";
void trans(int M, int N, int A[N][M], int B[M][N])
{
    int i, j, tmp;

    for (i = 0; i < N; i++) {
        for (j = 0; j < M; j++) {
            tmp = A[i][j];
            B[j][i] = tmp;
        }
    }    

}

/*
 * registerFunctions - This function registers your transpose
 *     functions with the driver.  At runtime, the driver will
 *     evaluate each of the registered functions and summarize their
 *     performance. This is a handy way to experiment with different
 *     transpose strategies.
 */
void registerFunctions()
{
    /* Register your solution function */
    registerTransFunction(transpose_submit, transpose_submit_desc); 

    /* Register any additional transpose functions */
    registerTransFunction(trans, trans_desc); 

}

/* 
 * is_transpose - This helper function checks if B is the transpose of
 *     A. You can check the correctness of your transpose by calling
 *     it before returning from the transpose function.
 */
int is_transpose(int M, int N, int A[N][M], int B[M][N])
{
    int i, j;

    for (i = 0; i < N; i++) {
        for (j = 0; j < M; ++j) {
            if (A[i][j] != B[j][i]) {
                return 0;
            }
        }
    }
    return 1;
}

