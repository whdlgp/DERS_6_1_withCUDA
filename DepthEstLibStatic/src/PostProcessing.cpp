#include "PostProcessing.h"

#include <stdio.h>
#include <string.h>

#define mmax(x,y) (x>y ? (x) : (y))
#define mmin(x,y) (x<y ? (x) : (y))

#define DILATION				0
#define EROSION					1
#define CLOSING					2
#define OPENING					3
#define OPENING_RECON			4  
#define CLOSING_RECON			5  
#define CLOSING_OPENING_RECON	6
#define OPENING_CLOSING_RECON	7

#define ITER	200

void mirroring(unsigned char **in, int width, int height, int h, int w)
{
    int i, j, k;
	
	if (h == 0);
	else
	{
		for (i=0,k=2*h-1; i<h; i++,k--)
			for (j=w; j<w+width; j++)
			{
				in[i][j] = in[k][j];
			}
		for (i=height+h,k=height+h-1; i<height+2*h; i++,k--)
			for (j=w; j<w+width; j++)
			{
				in[i][j] = in[k][j];
			}
		for (j=0,k=2*w-1; j<w; j++,k--)
			for (i=0; i<2*h+height; i++)
			{
				in[i][j] = in[i][k];
			}
		for (j=width+w,k=width+w-1; j<width+2*w; j++,k--)
			for (i=0; i<2*h+height;i++)
			{
				in[i][j] = in[i][k];
			}
    }
}


void entry_of_morphology(unsigned char** input, int flag, int sesize, int width, int height)
{  
	int i,j;
	int offset;  
	unsigned char** inY;
	unsigned char** in;
	unsigned char** out;
	unsigned char** ref;

	offset = (sesize-1)/2;

	inY = (unsigned char **)calloc(height, sizeof(unsigned char *));
	for (i=0; i<height; i++) inY[i] = (unsigned char *)calloc(width, sizeof(unsigned char));

	for (i=0; i<height; i++)
		for (j=0; j<width; j++)
		{
			inY[i][j] = input[i][j];
		}

	in  = (unsigned char **)calloc(height+2*offset, sizeof(unsigned char *));
	ref = (unsigned char **)calloc(height+2*offset, sizeof(unsigned char *));
	out = (unsigned char **)calloc(height+2*offset, sizeof(unsigned char *));

	for (i=0; i<height+2*offset; i++)
	{
		in[i]  = (unsigned char *)calloc(width+2*offset, sizeof(unsigned char));
		ref[i] = (unsigned char *)calloc(width+2*offset, sizeof(unsigned char));
		out[i] = (unsigned char *)calloc(width+2*offset, sizeof(unsigned char));
	}

	for (i=offset; i<height+offset; i++)
		for (j=offset; j<width+offset; j++)
		{
			in[i][j] = inY[i-offset][j-offset];
		}

	mirroring(in, width, height, offset, offset);

	for (i=0; i<height+2*offset; i++)
		for (j=0; j<width+2*offset; j++)
		{
			ref[i][j] = in[i][j];
			out[i][j] = 0;
		}

	if (flag == CLOSING_OPENING_RECON)			// closingopening by reconstruction
	{  
		closing_recon(in, out, ref, offset, width, height);
		for (i=0; i<height+2*offset; i++)
			for (j=0; j<width+2*offset; j++)
			{
				ref[i][j] = in[i][j];
			}

		opening_recon(in, out, ref, offset, width, height);
		for (i=offset; i<height+offset; i++)
			for (j=offset; j<width+offset; j++)
			{
				inY[i-offset][j-offset]=in[i][j];
			}
	}

	else if (flag == EROSION)					// erosion
	{  
		erosion(in, out, offset, offset, width, height);
		for (i=offset; i<height+offset; i++)
			for (j=offset; j<width+offset; j++)
			{
				inY[i-offset][j-offset] = out[i][j];
			}
	}

	else if (flag == DILATION)					// dilation
	{  
		dilation(in, out, offset, offset, width, height);
		for (i=offset; i<height+offset; i++)
			for (j=offset; j<width+offset; j++)
			{
				inY[i-offset][j-offset] = out[i][j];
			}
	}

	else if (flag == OPENING_CLOSING_RECON)		// opening_closing by reconstruction
	{  
		opening_recon(in, out, ref, offset, width, height);
		for (i=0; i<height+2*offset; i++)
			for (j=0; j<width+2*offset; j++)
			{
				ref[i][j] = in[i][j];
			}
			
		closing_recon(in, out, ref, offset, width, height);
		for (i=offset; i<height+offset; i++)
			for (j=offset; j<width+offset; j++)
			{
				inY[i-offset][j-offset] = in[i][j];
			}
	}

	else if (flag == CLOSING_RECON)				// closing by reconstruction
	{  
		closing_recon(in, out, ref, offset, width, height);
		for (i=offset; i<height+offset; i++)
			for (j=offset; j<width+offset; j++)
			{
				inY[i-offset][j-offset] = in[i][j];
			}
	}

	else if (flag == OPENING_RECON)				// opening by reconstruction
	{  
		opening_recon(in, out, ref, offset, width, height);
		for (i=offset; i<height+offset; i++)
			for (j=offset; j<width+offset; j++)
			{
				inY[i-offset][j-offset] = in[i][j];
			}
	}

	else if (flag == OPENING)					// opening
	{  
		opening(in, out, offset, offset, width, height);
		for (i=offset; i<height+offset; i++)
			for (j=offset; j<width+offset; j++)
			{
				inY[i-offset][j-offset] = out[i][j];
			}
	}

	else if (flag == CLOSING)					// closing
	{  
		closing(in, out, offset, offset, width, height);
		for (i=offset; i<height+offset; i++)
			for (j=offset; j<width+offset; j++)
			{
				inY[i-offset][j-offset] = out[i][j];
			}
	}

	for (i=0; i<height; i++)
		for (j=0; j<width; j++)
		{
			input[i][j] = inY[i][j];
		}

	for (i=0; i<height; i++) free(inY[i]);
	for (i=0; i<height+2*offset; i++)
	{
		free(in[i]);	
		free(ref[i]);	
		free(out[i]);
	}
	free(inY);	
	free(in);	
	free(ref);	
	free(out);
}


void closing_recon(unsigned char** in, unsigned char** out, unsigned char** ref, int se, int width, int height)
{
	int offset;
	int i, j, k;
	
	offset = se;
	dilation(in, out, se, offset, width, height);

	for (i=0; i<height+2*offset; i++)
		memcpy(in[i], out[i], width+2*offset);

	for (k=0; k<ITER; k++)
	{
		erosion(in, out, 1, offset, width, height); 
		for (i=0; i<height+2*offset; i++)
			for (j=0; j<width+2*offset; j++)
			{
				in[i][j] = mmax(ref[i][j], out[i][j]);
			}
	}
}

void opening_recon(unsigned char** in, unsigned char** out, unsigned char** ref, int se, int width, int height)
{
	int offset;
	int i, j, k;
	
	offset = se;
	erosion(in, out, se, offset, width, height);

	for (i=0; i<height+2*offset; i++)
		for (j=0; j<width+2*offset; j++)
		{
			in[i][j] = out[i][j];
		}

	for (k=0; k<se*2; k++)
	{
		dilation(in, out, 1, offset, width, height); 
		for (i=0; i<height+2*offset; i++)
			for (j=0; j<width+2*offset; j++)
			{
				in[i][j] = mmin(ref[i][j], out[i][j]);
			}
	}
}

void closing(unsigned char** in, unsigned char** out, int se, int offset, int width, int height)
{
	int i, j;

	dilation(in, out, se, offset, width, height);
	for (i=0; i<height+2*se; i++)
		for (j=0; j<width+2*se; j++)
		{
			in[i][j] = out[i][j];
		}
	erosion(in, out, se, offset, width, height);
}

void opening(unsigned char** in, unsigned char** out, int se, int offset, int width, int height)
{
	int i, j;

	erosion(in, out, se, offset, width, height);
	for (i=0; i<height+2*se; i++)
		for (j=0; j<width+2*se; j++)
		{
			in[i][j] = out[i][j];
		}
	dilation(in, out, se, offset, width, height);
}

void dilation(unsigned char** in, unsigned char** out, int se, int offset, int width, int height)
{
	int i, j, x, y, max;

	for (i=offset; i<height+offset; i++)
		for (j=offset; j<width+offset; j++)
		{
			max = 0;
			for (x=i-se; x<=i+se; x++)
				for(y=j-se; y<=j+se; y++)
				{
					if (in[x][y] > max) max = in[x][y];
				}
			out[i][j] = max;
		}
}

void erosion(unsigned char** in, unsigned char** out, int se, int offset, int width, int height)
{
	int i, j, x, y, min;

	for (i=offset; i<height+offset; i++)
		for (j=offset; j<width+offset; j++)
		{
			min = 10000;
			for (x=i-se; x<=i+se; x++)
				for(y=j-se; y<=j+se; y++)
				{
					if (in[x][y] < min) min = in[x][y];
				}
			out[i][j] = min;
		}
}


void FullSearch(unsigned char **Prev, unsigned char **Curr, int *motion_x, int *motion_y, int height, int width, int blocksize)
{
	int count = 0;
	int i, j, m, n, y, x;
	int diff, sum, min;
	int SearchRange = blocksize * 2;

	for (j=0; j<height; j+=blocksize)
		for (i=0; i<width; i+=blocksize)
		{
			min = 10000;
			sum = 0;
			for (y=0; y<blocksize; y++)
				for (x=0; x<blocksize; x++)
				{
					diff = abs(Curr[j+y][i+x] - Prev[j+y][i+x]);
					sum += diff;
				}
			if(sum > blocksize*blocksize) {
				for (n=-SearchRange; n<SearchRange; n++)
					for (m=-SearchRange; m<SearchRange; m++)
					{
						if ( (i+m>=0) && (i+m<width-blocksize) && (j+n>=0) && (j+n<height-blocksize) )
	
	
	
	
	
						{
							sum = 0;
							for (y=0; y<blocksize; y++)
								for (x=0; x<blocksize; x++)
	
	
	
	
	
	
	
								{
									diff = abs(Curr[j+y+n][i+x+m] - Prev[j+y][i+x]);
									sum += diff;
	
								}
							if (sum < min)
							{
								min = sum;
								motion_x[count] = m;
								motion_y[count] = n;
							}
						}
					}
			} else {
				motion_x[count] = 0;
                motion_y[count] = 0;
			}
			count++;
		}
}


void MakeBlockMap(int **tempP, int **tempB, int **PretempB, int new_height, int new_width, int precision, int BSize)
{
	int i, j, x, y;

	// ------------- Make MovingObjects Block -------------
	for (j=0; j<new_height; j+=BSize) 
		for (i=0; i<new_width; i+=BSize*precision) 
		{
			int count = 0;
			for (y=j; y<j+BSize; y++) 
				for (x=i; x<i+BSize*precision; x++) 
				{
					if (tempP[y][x] == 255) count++;
				}
			
			for (y=j; y<j+BSize; y++) 
				for (x=i; x<i+BSize*precision; x++) 
				{
					if (PretempB[j][i] == 255) 
					{
						if (count > (int)(0.25*(double)BSize*(double)precision)) tempB[y][x] = 255;
						else tempB[y][x] = 0;
					} 
					else 
					{
						if (count > (int)(0.5 *(double)BSize*(double)precision)) tempB[y][x] = 255;
						else tempB[y][x] = 0;
					}
				}
		}

	// ------------- Hole Filling of BlockMap --------------
	for (j=0; j<new_height; j+=BSize) 
		for (i=0; i<new_width; i+=BSize*precision) 
		{
			int count = 0;
			if ((i - BSize*precision >= 0) && (i + BSize*precision < new_width))
			{
				if (tempB[j][i - BSize*precision] == 255) count++;
				if (tempB[j][i + BSize*precision] == 255) count++;
			}
			else count++;

			if ((j - BSize >= 0) && (j + BSize < new_height))
			{
				if (tempB[j - BSize][i] == 255) count++;
				if (tempB[j + BSize][i] == 255) count++;
			}
			else count++;

			if (count >= 3)
			{
				for (y=j; y<j+BSize; y++)
					for (x=i; x<i+BSize*precision; x++)
					{
						tempB[y][x] = 255;
					}
			}
		}

	for (j=0; j<new_height; j++)
		for (i=0; i<new_width; i++)
		{
			PretempB[j][i] = tempB[j][i];
		}

	// ------------- Expansion of BlockMap --------------
	for (j=0; j<new_height; j+=BSize)
		for (i=0; i<new_width; i+=BSize*precision)
		{
			if ((i - BSize*precision >= 0) && (i + BSize*precision < new_width))
			{
				if (PretempB[j][i - BSize*precision] == 0 && PretempB[j][i] == 255) 
				{
					for (y=j; y<j+BSize; y++)
					{
						for (x=i-BSize*precision; x<i; x++)
						{
							tempB[y][x] = 255;
						}
					}
				}
				if (PretempB[j][i + BSize*precision] == 0 && PretempB[j][i] == 255) 
				{
					for (y=j; y<j+BSize; y++)
					{
						for (x=i+BSize*precision; x<i+2*BSize*precision-1; x++)
						{
							tempB[y][x] = 255;
						}
					}
				}
			}

			if ((j - BSize >= 0) && (j + BSize < new_height))
			{
				if (PretempB[j - BSize][i] == 0 && PretempB[j][i] == 255) 
				{
					for (y=j-BSize; y<j; y++)
					{
						for (x=i; x<i+BSize*precision; x++)
						{
							tempB[y][x] = 255;
						}
					}
				}
				if (PretempB[j+BSize][i] == 0 && PretempB[j][i] == 255) 
				{
					for (y=j+BSize; y<j+2*BSize-1; y++)
					{
						for (x=i; x<i+BSize*precision; x++)
						{
							tempB[y][x] = 255;
						}
					}
				}
			}
		}

}