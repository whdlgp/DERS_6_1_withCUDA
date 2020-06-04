/********************************************************************************************
* Implementing Graph Cuts on CUDA using algorithm given in CVGPU '08                       **
* paper "CUDA Cuts: Fast Graph Cuts on GPUs"                                               **
*                                                                                          **
* Copyright (c) 2008 International Institute of Information Technology.                    **
* All rights reserved.                                                                     **
*                                                                                          **
* Permission to use, copy, modify and distribute this software and its documentation for   **
* educational purpose is hereby granted without fee, provided that the above copyright     **
* notice and this permission notice appear in all copies of this software and that you do  **
* not sell the software.                                                                   **
*                                                                                          **
* THE SOFTWARE IS PROVIDED "AS IS" AND WITHOUT WARRANTY OF ANY KIND,EXPRESS, IMPLIED OR    **
* OTHERWISE.                                                                               **
*                                                                                          **
* Created By Vibhav Vineet.                                                                **
********************************************************************************************/

#include "CudaCuts.h"

#include <algorithm>

using namespace std;

typedef struct inputfile_t
{
	char* filename;
	int width;
	int height;
	int numOfLabels;
	int* dataterm_error;
	int* smoothness_table;
	int* hcue;
	int* vcue;
}inputfile_t;

void loadFile(inputfile_t& input_file);
void writePGM(const char *filename, CudaCuts& cuts);

int main(int argc, char* argv[])
{
	if(argc != 2)
	{
		printf("usage: <datafilename.txt>");
	}
	char* dataFile = argv[1];

	// Read input file
	inputfile_t input_file;
	input_file.filename = dataFile;
	loadFile(input_file);

	// Set-up graph with input
	CudaCuts cuts(input_file.width, input_file.height, input_file.numOfLabels
				, input_file.dataterm_error, input_file.smoothness_table
				, input_file.hcue, input_file.vcue);

	// Generate input labels, 0 to (numOfLabels)
	std::vector<int> labels(input_file.numOfLabels);
	int n=0;
	generate(labels.begin(), labels.end(), [&n] { return n++;});
	cuts.run(labels);
	
	//save output file
	const char* output_name = "result_sponge/flower_cuda_test.pgm";
	writePGM(output_name, cuts);

	return 0;
}

void writePGM(const char* filename, CudaCuts& cuts)
{
	int** out_pixel_values = (int**)malloc(sizeof(int*)*cuts.height);

	for (int i = 0; i < cuts.height; i++)
	{
		out_pixel_values[i] = (int*)malloc(sizeof(int)* cuts.width);
		for (int j = 0; j < cuts.width; j++) {
			out_pixel_values[i][j] = 0;
		}
	}
	for (int i = 0; i < cuts.graph_size1; i++)
	{

		int row = i / cuts.width1, col = i % cuts.width1;

		if (row >= 0 && col >= 0 && row <= cuts.height - 1 && col <= cuts.width - 1)
			out_pixel_values[row][col] = (int)(float(cuts.pixelLabel[row*cuts.width + col])/(cuts.num_Labels-1) * 255);
	}
	FILE* fp = fopen(filename, "w");

	fprintf(fp, "%c", 'P');
	fprintf(fp, "%c", '2');
	fprintf(fp, "%c", '\n');
	fprintf(fp, "%d %c %d %c ", cuts.width, ' ', cuts.height, '\n');
	fprintf(fp, "%d %c", 255, '\n');

	for (int i = 0; i<cuts.height; i++)
	{
		for (int j = 0; j<cuts.width; j++)
		{
			fprintf(fp, "%d\n", out_pixel_values[i][j]);
		}
	}
	fclose(fp);
	for (int i = 0; i < cuts.height; i++)
		free(out_pixel_values[i]);
	free(out_pixel_values);
}

void loadFile(inputfile_t& input_file)
{
	printf("enterd\n");
	int &width = input_file.width;
	int &height = input_file.height;
	int &nLabels = input_file.numOfLabels;
	
	int *&dataCostArray = input_file.dataterm_error;
	int *&smoothCostArray = input_file.smoothness_table;
	int *&hCue = input_file.hcue;
	int *&vCue = input_file.vcue;

	char* filename = input_file.filename;
	FILE *fp = fopen(filename, "r");

	fscanf(fp, "%d %d %d", &width, &height, &nLabels);

	int i, n, x, y;
	int gt;
	for (i = 0; i < width * height; i++)
		fscanf(fp, "%d", &gt);

	dataCostArray = (int*)malloc(sizeof(int)* width * height * nLabels);

	for (int c = 0; c < nLabels; c++) {
		n = c;
		for (i = 0; i < width * height; i++) {
			fscanf(fp, "%d", &dataCostArray[n]);
			n += nLabels;
		}
	}

	hCue = (int*)malloc(sizeof(int)* width * height);
	vCue = (int*)malloc(sizeof(int)* width * height);

	n = 0;
	for (y = 0; y < height; y++) {
		for (x = 0; x < width - 1; x++) {
			fscanf(fp, "%d", &hCue[n++]);
		}
		hCue[n++] = 0;
	}

	n = 0;
	for (y = 0; y < height - 1; y++) {
		for (x = 0; x < width; x++) {
			fscanf(fp, "%d", &vCue[n++]);
		}
	}
	for (x = 0; x < width; x++) {
		vCue[n++] = 0;
	}

	fclose(fp);
	smoothCostArray = (int*)malloc(sizeof(int)*nLabels * nLabels);
	for(int i = 0; i < nLabels; i++)
	{
		for(int j = 0; j < nLabels; j++)
		{
			smoothCostArray[i*nLabels + j] = abs(i-j);
			//smoothCostArray[i*nLabels + j] = (i == j) ? 0 : 255;
		}
	}
}