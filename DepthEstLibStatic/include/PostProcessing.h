#include "stdlib.h"
#include "math.h"
#include "version.h"

template <class T>
void quick_sort(T *data, int left, int right )
{
	int i,j;
	T x,y;

	i = left;
	j = right;

	x = data[(left+right)/2];
	
	do {
		while(data[i]<x && i<right) i++;
		while(x<data[j] && j>left) j--;

		if(i<=j) {
			y = data[i];
			data[i] = data[j];
			data[j] = y;
			i++;
			j--;
		}
	}while(i<=j);

	if(left<j) quick_sort(data,left,j);
	if(i<right) quick_sort(data,i,right);
}

void KSWAP(unsigned char *data, int x, int y);

template <class T>
void median(int height, int width, int masksize, T **src_image)
{
	int i, j, m, n, cnt; 
	
  T var;
  T *dst_image;
	dst_image = (T *)calloc(height*width, sizeof(T));
	memset(dst_image, 255, height*width*sizeof(T));
	T *mask = (T *)calloc(masksize*masksize, sizeof(T));

	for(i=0; i<height-masksize+((masksize-1)/2); i++)
	{
		for(j=0; j<width-masksize+((masksize-1)/2); j++)
		{
			cnt = -1;
	
			if( (i+masksize > height) || (j+masksize > width) ) 
				continue;
			for(m=0; m<masksize; m++)
			{
				for(n=0; n<masksize; n++)
				{
					mask[++cnt] = /*(unsigned char)*/src_image[m+i][n+j];
				}
			}

			quick_sort( mask, 0, cnt );
			var = mask[cnt/2];
			
			dst_image[(i+((masksize-1)/2))*width+j+((masksize-1)/2)] = var;
		}
	}

	for(i=0; i<height; i++)
	{
		for(j=0; j<width; j++)
		{
                        if(dst_image[i*width+j] != ((1<<(sizeof(T)*8))-1)) 
				src_image[i][j] = dst_image[i*width+j];
		}
	}

	free(mask);
	free(dst_image);
}



void mirroring(unsigned char **in, int width, int height, int h, int w);
void entry_of_morphology(unsigned char** input, int flag, int sesize, int width, int height);
void closing_recon(unsigned char** in, unsigned char** out, unsigned char** ref, int se, int width, int height);
void opening_recon(unsigned char** in, unsigned char** out, unsigned char** ref, int se, int width, int height);
void closing(unsigned char** in, unsigned char** out, int se, int offset, int width, int height);
void opening(unsigned char** in, unsigned char** out, int se, int offset, int width, int height);
void dilation(unsigned char** in, unsigned char** out, int se, int offset, int width, int height);
void erosion(unsigned char** in, unsigned char** out, int se, int offset, int width, int height);
void FullSearch(unsigned char **Prev, unsigned char **Curr, int *motion_x, int *motion_y, int height, int width, int blocksize);
void MakeBlockMap(int **tempP, int **tempB, int **PretempB, int new_height, int new_width, int precision, int BSize);