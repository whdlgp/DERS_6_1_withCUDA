
#ifndef CLIP3
#define CLIP3(x,min,max) ( (x)<(min)?(min):((x)>(max)?(max):(x)) )
#endif

#ifndef TURN3
#define TURN3(x,min,max) ( (x)<(min)?(-x): ((x)>(max)?(2*(max)-(x)):(x)) )
#endif

#ifndef BYTE
#define BYTE unsigned char
#endif

void HorizontalLinearFilter_1D_half(BYTE *in, BYTE *out, int width, int padding_size)
{
	int i, left, right;
	int max_width = width + padding_size;
	int width_minus1 = width-1;

	for(i=-padding_size; i<max_width; i++)
	{
		left = CLIP3(i, 0, width_minus1);
		right = CLIP3(i+1, 0, width_minus1);
		out[ (i+padding_size)<<1   ] =  in[left];
		out[((i+padding_size)<<1)+1] = (in[left]+in[right]+1)>>1;
	}
}

void HorizontalLinearFilter_2D_half(BYTE **in, BYTE **out, int width, int height, int padding_size)
{
	int i, left, right;
	int max_width = width + padding_size;
	int width_minus1 = width-1;
	int h, int_pel, half_pel;

	for(i=-padding_size; i<max_width; i++)
	{
		left = CLIP3(i, 0, width_minus1);
		right = CLIP3(i+1, 0, width_minus1);
		int_pel = (i+padding_size)<<1;
		half_pel = int_pel+1;
		for(h=0; h<height; h++)
		{
			out[h][int_pel] = in[h][left];
			out[h][half_pel] = (in[h][left]+in[h][right]+1)>>1;
		}
	}
}

void HorizontalLinearFilter_1D_qpel(BYTE *in, BYTE *out, int width, int padding_size)
{
	int i, left, right;
	int max_width = width + padding_size;
	int width_minus1 = width-1;

	for(i=-padding_size; i<max_width; i++)
	{
		left = CLIP3(i, 0, width_minus1);
		right = CLIP3(i+1, 0, width_minus1);
		out[ (i+padding_size)<<2   ] =  in[left];
		out[((i+padding_size)<<2)+1] = (in[left]*3 + in[right]   + 2)>>2;
		out[((i+padding_size)<<2)+2] = (in[left]   + in[right]   + 1)>>1;
		out[((i+padding_size)<<2)+3] = (in[left]   + in[right]*3 + 2)>>2;
	}
}

void HorizontalLinearFilter_2D_qpel(BYTE **in, BYTE **out, int width, int height, int padding_size)
{
	int i, left, right;
	int max_width = width + padding_size;
	int width_minus1 = width-1;
	int h, pel[4];

	for(i=-padding_size; i<max_width; i++)
	{
		left = CLIP3(i, 0, width_minus1);
		right = CLIP3(i+1, 0, width_minus1);
		pel[0] = (i+padding_size)<<2;
		pel[1] = pel[0]+1;
		pel[2] = pel[1]+1;
		pel[3] = pel[2]+1;
		for(h=0; h<height; h++)
		{
			out[h][pel[0]] =  in[h][left];
			out[h][pel[1]] = (in[h][left]*3 + in[h][right]   + 2)>>2;
			out[h][pel[2]] = (in[h][left]   + in[h][right]   + 1)>>1;
			out[h][pel[3]] = (in[h][left]   + in[h][right]*3 + 2)>>2;
		}
	}
}

void HorizontalAVCFilter_1D_half(BYTE *in, BYTE *out, int width, int padding_size)
{
	int i, pel[6];
	int max_width = width + padding_size;
	int width_minus1 = width-1;

	for(i=-padding_size; i<max_width; i++)
	{
		pel[0] = CLIP3(i-2, 0, width_minus1);
		pel[1] = CLIP3(i-1, 0, width_minus1);
		pel[2] = CLIP3(i  , 0, width_minus1);
		pel[3] = CLIP3(i+1, 0, width_minus1);
		pel[4] = CLIP3(i+2, 0, width_minus1);
		pel[5] = CLIP3(i+3, 0, width_minus1);

		out[ (i+padding_size)<<1   ] =  in[pel[2]];
		out[((i+padding_size)<<1)+1] = CLIP3( (20*(in[pel[2]]+in[pel[3]]) - 5*(in[pel[1]]+in[pel[4]]) + (in[pel[0]]+in[pel[5]]) +16)>>5, 0, 255 );
	}
}

void HorizontalAVCFilter_2D_half(BYTE **in, BYTE **out, int width, int height, int padding_size)
{
	int i, pel[6];
	int max_width = width + padding_size;
	int width_minus1 = width-1;
	int h, int_pel, half_pel;

	for(i=-padding_size; i<max_width; i++)
	{
		pel[0] = CLIP3(i-2, 0, width_minus1);
		pel[1] = CLIP3(i-1, 0, width_minus1);
		pel[2] = CLIP3(i  , 0, width_minus1);
		pel[3] = CLIP3(i+1, 0, width_minus1);
		pel[4] = CLIP3(i+2, 0, width_minus1);
		pel[5] = CLIP3(i+3, 0, width_minus1);

		int_pel = (i+padding_size)<<1;
		half_pel = int_pel+1;
		for(h=0; h<height; h++)
		{
			out[h][int_pel] = in[h][pel[2]];
			out[h][half_pel] = CLIP3( (20*(in[h][pel[2]]+in[h][pel[3]]) - 5*(in[h][pel[1]]+in[h][pel[4]]) + (in[h][pel[0]]+in[h][pel[5]]) +16)>>5, 0, 255 );
		}
	}
}

void HorizontalAVCFilter_1D_qpel(BYTE *in, BYTE *out, int width, int padding_size)
{
	int i, pel[6];
	int max_width = width + padding_size;
	int width_minus1 = width-1;

	for(i=-padding_size; i<max_width; i++)
	{
		pel[0] = CLIP3(i-2, 0, width_minus1);
		pel[1] = CLIP3(i-1, 0, width_minus1);
		pel[2] = CLIP3(i  , 0, width_minus1);
		pel[3] = CLIP3(i+1, 0, width_minus1);
		pel[4] = CLIP3(i+2, 0, width_minus1);
		pel[5] = CLIP3(i+3, 0, width_minus1);

		out[ (i+padding_size)<<2   ] = in[pel[2]];
		out[((i+padding_size)<<2)+2] = CLIP3( (20*(in[pel[2]]+in[pel[3]]) - 5*(in[pel[1]]+in[pel[4]]) + (in[pel[0]]+in[pel[5]]) +16)>>5, 0, 255 );
		out[((i+padding_size)<<2)+1] = (in[pel[2]] + out[((i+padding_size)<<2)+2] + 1)>>1;
		out[((i+padding_size)<<2)+3] = (in[pel[3]] + out[((i+padding_size)<<2)+2] + 1)>>1;
	}
}

void HorizontalAVCFilter_2D_qpel(BYTE **in, BYTE **out, int width, int height, int padding_size)
{
	int i;
	int max_width = width + padding_size;
	int width_minus1 = width-1;
	int h, ipel[6], opel[4];

	for(i=-padding_size; i<max_width; i++)
	{
		ipel[0] = CLIP3(i-2, 0, width_minus1);
		ipel[1] = CLIP3(i-1, 0, width_minus1);
		ipel[2] = CLIP3(i  , 0, width_minus1);
		ipel[3] = CLIP3(i+1, 0, width_minus1);
		ipel[4] = CLIP3(i+2, 0, width_minus1);
		ipel[5] = CLIP3(i+3, 0, width_minus1);

		opel[0] = (i+padding_size)<<2;
		opel[1] = opel[0]+1;
		opel[2] = opel[1]+1;
		opel[3] = opel[2]+1;
		for(h=0; h<height; h++)
		{
			out[h][opel[0]] =  in[h][ipel[2]];
			out[h][opel[2]] = CLIP3( (20*(in[h][ipel[2]]+in[h][ipel[3]]) - 5*(in[h][ipel[1]]+in[h][ipel[4]]) + (in[h][ipel[0]]+in[h][ipel[5]]) +16)>>5, 0, 255 );
			out[h][opel[1]] = (in[h][ipel[2]] + out[h][opel[2]] + 1)>>1;
			out[h][opel[3]] = (in[h][ipel[3]] + out[h][opel[2]] + 1)>>1;
		}
	}
}

void Horizontal6tapFilter_1D_half(BYTE *in, BYTE *out, int width, int padding_size, int coeff[6])
{
	int i, ipel[6], opel[2];
	int max_width = width + padding_size;
	int width_minus1 = width-1;
	int sum_coeff, offset;

	for(i=sum_coeff=0; i<6; i++)
		sum_coeff += coeff[i];
	offset = sum_coeff>>1;

	for(i=-padding_size; i<max_width; i++)
	{
		ipel[0] = CLIP3(i-2, 0, width_minus1);
		ipel[1] = CLIP3(i-1, 0, width_minus1);
		ipel[2] = CLIP3(i  , 0, width_minus1);
		ipel[3] = CLIP3(i+1, 0, width_minus1);
		ipel[4] = CLIP3(i+2, 0, width_minus1);
		ipel[5] = CLIP3(i+3, 0, width_minus1);
		opel[0] = (i+padding_size)<<1;
		opel[1] = opel[0]+1;

		out[opel[0]] = in[ipel[2]];
		out[opel[1]] = CLIP3( (in[ipel[0]]*coeff[0] + in[ipel[1]]*coeff[1] + in[ipel[2]]*coeff[2] + 
													 in[ipel[3]]*coeff[3] + in[ipel[4]]*coeff[4] + in[ipel[5]]*coeff[5] + offset)/sum_coeff, 0, 255);
	}
}

void Horizontal6tapFilter_2D_half(BYTE **in, BYTE **out, int width, int height, int padding_size, int coeff[6])
{
	int i, ipel[6], opel[2];
	int max_width = width + padding_size;
	int width_minus1 = width-1;
	int h;
	int sum_coeff, offset;

	for(i=sum_coeff=0; i<6; i++)
		sum_coeff += coeff[i];
	offset = sum_coeff>>1;

	for(i=-padding_size; i<max_width; i++)
	{
		ipel[0] = CLIP3(i-2, 0, width_minus1);
		ipel[1] = CLIP3(i-1, 0, width_minus1);
		ipel[2] = CLIP3(i  , 0, width_minus1);
		ipel[3] = CLIP3(i+1, 0, width_minus1);
		ipel[4] = CLIP3(i+2, 0, width_minus1);
		ipel[5] = CLIP3(i+3, 0, width_minus1);

		opel[0] = (i+padding_size)<<1;
		opel[1] = opel[0]+1;

		for(h=0; h<height; h++)
		{
			out[h][opel[0]] = in[h][ipel[2]];
			out[h][opel[1]] = CLIP3( (in[h][ipel[0]]*coeff[0] + in[h][ipel[1]]*coeff[1] + in[h][ipel[2]]*coeff[2] + 
															 in[h][ipel[3]]*coeff[3] + in[h][ipel[4]]*coeff[4] + in[h][ipel[5]]*coeff[5] + offset)/sum_coeff, 0, 255);
		}
	}
}

void HorizontalCubicFilter_1D_half(BYTE *in, BYTE *out, int width, int padding_size)
{
	int i, ipel[4], opel[2];
	int max_width = width + padding_size;
	int width_minus1 = width-1;

	for(i=-padding_size; i<max_width; i++)
	{
//		ipel[0] = CLIP3(i-1, 0, width_minus1);
//		ipel[1] = CLIP3(i  , 0, width_minus1);
//		ipel[2] = CLIP3(i+1, 0, width_minus1);
//		ipel[3] = CLIP3(i+2, 0, width_minus1);
		ipel[0] = TURN3(i-1, 0, width_minus1);
		ipel[1] = TURN3(i  , 0, width_minus1);
		ipel[2] = TURN3(i+1, 0, width_minus1);
		ipel[3] = TURN3(i+2, 0, width_minus1);

		opel[0] = (i+padding_size)<<1;
		opel[1] = opel[0]+1;

		out[opel[0]] = in[ipel[1]];
		out[opel[1]] = CLIP3( (5*(in[ipel[1]] +    in[ipel[2]])-   in[ipel[0]] -   in[ipel[3]] +  4)>>3, 0, 255 );
	}
}

void HorizontalCubicFilter_2D_half(BYTE **in, BYTE **out, int width, int height, int padding_size)
{
	int i, ipel[4], opel[2];
	int max_width = width + padding_size;
	int width_minus1 = width-1;
	int h;

	for(i=-padding_size; i<max_width; i++)
	{
//		ipel[0] = CLIP3(i-1, 0, width_minus1);
//		ipel[1] = CLIP3(i  , 0, width_minus1);
//		ipel[2] = CLIP3(i+1, 0, width_minus1);
//		ipel[3] = CLIP3(i+2, 0, width_minus1);
		ipel[0] = TURN3(i-1, 0, width_minus1);
		ipel[1] = TURN3(i  , 0, width_minus1);
		ipel[2] = TURN3(i+1, 0, width_minus1);
		ipel[3] = TURN3(i+2, 0, width_minus1);

		opel[0] = (i+padding_size)<<1;
		opel[1] = opel[0]+1;

		for(h=0; h<height; h++)
		{
			out[h][opel[0]] = in[h][ipel[1]];
			out[h][opel[1]] = CLIP3( (5*(in[h][ipel[1]] +    in[h][ipel[2]])-   in[h][ipel[0]] -   in[h][ipel[3]] +  4)>>3, 0, 255 );
		}
	}
}

void HorizontalCubicFilter_1D_qpel(BYTE *in, BYTE *out, int width, int padding_size)
{
	int i, ipel[4], opel[4];
	int max_width = width + padding_size;
	int width_minus1 = width-1;

	for(i=-padding_size; i<max_width; i++)
	{
//		ipel[0] = CLIP3(i-1, 0, width_minus1);
//		ipel[1] = CLIP3(i  , 0, width_minus1);
//		ipel[2] = CLIP3(i+1, 0, width_minus1);
//		ipel[3] = CLIP3(i+2, 0, width_minus1);
		ipel[0] = TURN3(i-1, 0, width_minus1);
		ipel[1] = TURN3(i  , 0, width_minus1);
		ipel[2] = TURN3(i+1, 0, width_minus1);
		ipel[3] = TURN3(i+2, 0, width_minus1);

		opel[0] = (i+padding_size)<<2;
		opel[1] = opel[0]+1;
		opel[2] = opel[0]+2;
		opel[3] = opel[0]+3;

		out[opel[0]] = in[ipel[1]];
		out[opel[1]] = CLIP3( (57*in[ipel[1]] + 19*in[ipel[2]] - 9*in[ipel[0]] - 3*in[ipel[3]] + 32)>>6, 0, 255 );
		out[opel[2]] = CLIP3( (5*(in[ipel[1]] +    in[ipel[2]])-   in[ipel[0]] -   in[ipel[3]] +  4)>>3, 0, 255 );
		out[opel[3]] = CLIP3( (19*in[ipel[1]] + 57*in[ipel[2]] - 3*in[ipel[0]] - 9*in[ipel[3]] + 32)>>6, 0, 255 );
	}
}

void HorizontalCubicFilter_2D_qpel(BYTE **in, BYTE **out, int width, int height, int padding_size)
{
	int i, ipel[4], opel[4];
	int max_width = width + padding_size;
	int width_minus1 = width-1;
	int h;

	for(i=-padding_size; i<max_width; i++)
	{
//		ipel[0] = CLIP3(i-1, 0, width_minus1);
//		ipel[1] = CLIP3(i  , 0, width_minus1);
//		ipel[2] = CLIP3(i+1, 0, width_minus1);
//		ipel[3] = CLIP3(i+2, 0, width_minus1);
		ipel[0] = TURN3(i-1, 0, width_minus1);
		ipel[1] = TURN3(i  , 0, width_minus1);
		ipel[2] = TURN3(i+1, 0, width_minus1);
		ipel[3] = TURN3(i+2, 0, width_minus1);

		opel[0] = (i+padding_size)<<2;
		opel[1] = opel[0]+1;
		opel[2] = opel[0]+2;
		opel[3] = opel[0]+3;

		for(h=0; h<height; h++)
		{
			out[h][opel[0]] = in[h][ipel[1]];
			out[h][opel[1]] = CLIP3( (57*in[h][ipel[1]] + 19*in[h][ipel[2]] - 9*in[h][ipel[0]] - 3*in[h][ipel[3]] + 32)>>6, 0, 255 );
			out[h][opel[2]] = CLIP3( (5*(in[h][ipel[1]] +    in[h][ipel[2]])-   in[h][ipel[0]] -   in[h][ipel[3]] +  4)>>3, 0, 255 );
			out[h][opel[3]] = CLIP3( (19*in[h][ipel[1]] + 57*in[h][ipel[2]] - 3*in[h][ipel[0]] - 9*in[h][ipel[3]] + 32)>>6, 0, 255 );
		}
	}
}

void DummyFilter_1D(BYTE *in, BYTE *out, int width, int padding_size)
{
	return;
}

void DummyFilter_2D(BYTE **in, BYTE **out, int width, int height, int padding_size)
{
  int h,i;
  int max_width = width + padding_size;
  for(h=0;h<height;h++)
  {
    for(i=-padding_size; i<max_width; i++)
	  {	 
      out[h][i] = in[h][i];
    }
  }
	return;
}

void VerticalLinearFilter_2D_half(BYTE **in, BYTE **out, int width, int height, int padding_size)
{
	int i, top, bottom;
	int max_width = width + padding_size;
  int height_minus1 = height-1;
	int h, int_pel, half_pel;

  for(h=0; h<height; h++)
  {
    top = CLIP3(h,0,height_minus1);
    bottom = CLIP3(h+1,0,height_minus1);
    int_pel = (h)<<1;
	  half_pel = int_pel+1;
	
	  for(i=-padding_size; i<max_width; i++)
	  {	 
      out[int_pel][i] = in[top][i];
			out[half_pel][i] = (in[top][i]+in[bottom][i]+1)>>1;
	  }
  }
}

void VerticalLinearFilter_2D_qpel(BYTE **in, BYTE **out, int width, int height, int padding_size)
{
	int i, top, bottom;
	int max_width = width + padding_size;
	int height_minus1 = height-1;
	int h, pel[4];

  for(h=0; h<height; h++)
	{
		top = CLIP3(h, 0, height_minus1);
		bottom = CLIP3(h+1, 0, height_minus1);
		pel[0] = (h)<<2;
		pel[1] = pel[0]+1;
		pel[2] = pel[1]+1;
		pel[3] = pel[2]+1;
		for(i=-padding_size; i<max_width; i++)
		{
			out[pel[0]][i] =  in[top][i];
			out[pel[1]][i] = (in[top][i]*3 + in[bottom][i]   + 2)>>2;
			out[pel[2]][i] = (in[top][i]   + in[bottom][i]   + 1)>>1;
			out[pel[3]][i] = (in[top][i]   + in[bottom][i]*3 + 2)>>2;
		}
	}
}

void VerticalAVCFilter_2D_half(BYTE **in, BYTE **out, int width, int height, int padding_size)
{
	int i, pel[6];
	int max_width = width + padding_size;
	int height_minus1 = height-1;
	int h, int_pel, half_pel;

  for(h=0; h<height; h++)
	{
    pel[0] = CLIP3(h-2, 0, height_minus1);
    pel[1] = CLIP3(h-1, 0, height_minus1);
		pel[2] = CLIP3(h  , 0, height_minus1);
		pel[3] = CLIP3(h+1, 0, height_minus1);
		pel[4] = CLIP3(h+2, 0, height_minus1);
		pel[5] = CLIP3(h+3, 0, height_minus1);

		int_pel = (h)<<1;
		half_pel = int_pel+1;
    for(i=-padding_size; i<max_width; i++)
		{
			out[int_pel][i] = in[pel[2]][i];
			out[half_pel][i] = CLIP3( (20*(in[pel[2]][i]+in[pel[3]][i]) - 5*(in[pel[1]][i]+in[pel[4]][i]) + (in[pel[0]][i]+in[pel[5]][i]) +16)>>5, 0, 255 );
		}
	}
}

void VerticalAVCFilter_2D_qpel(BYTE **in, BYTE **out, int width, int height, int padding_size)
{
	int i;
	int max_width = width + padding_size;
	int height_minus1 = height-1;
	int h, ipel[6], opel[4];

	for(h=0; h<height; h++)
	{
		ipel[0] = CLIP3(h-2, 0, height_minus1);
		ipel[1] = CLIP3(h-1, 0, height_minus1);
		ipel[2] = CLIP3(h  , 0, height_minus1);
		ipel[3] = CLIP3(h+1, 0, height_minus1);
		ipel[4] = CLIP3(h+2, 0, height_minus1);
		ipel[5] = CLIP3(h+3, 0, height_minus1);

		opel[0] = (h)<<2;
		opel[1] = opel[0]+1;
		opel[2] = opel[1]+1;
		opel[3] = opel[2]+1;
		
    for(i=-padding_size; i<max_width; i++)
		{
			out[opel[0]][i] =  in[ipel[2]][i];
			out[opel[2]][i] = CLIP3( (20*(in[ipel[2]][i]+in[ipel[3]][i]) - 5*(in[ipel[1]][i]+in[ipel[4]][i]) + (in[ipel[0]][i]+in[ipel[5]][i]) +16)>>5, 0, 255 );
			out[opel[1]][i] = (in[ipel[2]][i] + out[opel[2]][i] + 1)>>1;
			out[opel[3]][i] = (in[ipel[3]][i] + out[opel[2]][i] + 1)>>1;
		}
	}
}


void Vertical6tapFilter_2D_half(BYTE **in, BYTE **out, int width, int height, int padding_size, int coeff[6])
{
	int i, ipel[6], opel[2];
	int max_width = width + padding_size;
	int height_minus1 = height-1;
	int h;
	int sum_coeff, offset;

	for(i=sum_coeff=0; i<6; i++)
		sum_coeff += coeff[i];
	offset = sum_coeff>>1;

	for(h=0; h<height; h++)
	{
		ipel[0] = CLIP3(h-2, 0, height_minus1);
		ipel[1] = CLIP3(h-1, 0, height_minus1);
		ipel[2] = CLIP3(h  , 0, height_minus1);
		ipel[3] = CLIP3(h+1, 0, height_minus1);
		ipel[4] = CLIP3(h+2, 0, height_minus1);
		ipel[5] = CLIP3(h+3, 0, height_minus1);

		opel[0] = (h)<<1;
		opel[1] = opel[0]+1;

    for(i=-padding_size; i<max_width; i++)
		{
			out[opel[0]][i] = in[ipel[2]][i];
			out[opel[1]][i] = CLIP3( (in[ipel[0]][i]*coeff[0] + in[ipel[1]][i]*coeff[1] + in[ipel[2]][i]*coeff[2] + 
															 in[ipel[3]][i]*coeff[3] + in[ipel[4]][i]*coeff[4] + in[ipel[5]][i]*coeff[5] + offset)/sum_coeff, 0, 255);
		}
	}
}

void VerticalCubicFilter_2D_half(BYTE **in, BYTE **out, int width, int height, int padding_size)
{
	int i, ipel[4], opel[2];
	int max_width = width + padding_size;
	int height_minus1 = height-1;
	int h;

	for(h=0; h<height; h++)
	{
//		ipel[0] = CLIP3(i-1, 0, width_minus1);
//		ipel[1] = CLIP3(i  , 0, width_minus1);
//		ipel[2] = CLIP3(i+1, 0, width_minus1);
//		ipel[3] = CLIP3(i+2, 0, width_minus1);
    ipel[0] = TURN3(h-1, 0, height_minus1);
    ipel[1] = TURN3(h  , 0, height_minus1);
    ipel[2] = TURN3(h+1, 0, height_minus1);
    ipel[3] = TURN3(h+2, 0, height_minus1);

    opel[0] = (h)<<1;
		opel[1] = opel[0]+1;

    for(i=-padding_size; i<max_width; i++)
		{
			out[opel[0]][i] = in[ipel[1]][i];
			out[opel[1]][i] = CLIP3( (5*(in[ipel[1]][i] +    in[ipel[2]][i])-   in[ipel[0]][i] -   in[ipel[3]][i] +  4)>>3, 0, 255 );
		}
	}
}

void VerticalCubicFilter_2D_qpel(BYTE **in, BYTE **out, int width, int height, int padding_size)
{
	int i, ipel[4], opel[4];
	int max_width = width + padding_size;
	int height_minus1 = height-1;
	int h;

	for(h=0; h<height; h++)
	{
//		ipel[0] = CLIP3(i-1, 0, width_minus1);
//		ipel[1] = CLIP3(i  , 0, width_minus1);
//		ipel[2] = CLIP3(i+1, 0, width_minus1);
//		ipel[3] = CLIP3(i+2, 0, width_minus1);
		ipel[0] = TURN3(h-1, 0, height_minus1);
		ipel[1] = TURN3(h  , 0, height_minus1);
		ipel[2] = TURN3(h+1, 0, height_minus1);
		ipel[3] = TURN3(h+2, 0, height_minus1);

		opel[0] = (h)<<2;
		opel[1] = opel[0]+1;
		opel[2] = opel[0]+2;
		opel[3] = opel[0]+3;

    for(i=-padding_size; i<max_width; i++)
		{
			out[opel[0]][i] = in[ipel[1]][i];
			out[opel[1]][i] = CLIP3( (57*in[ipel[1]][i] + 19*in[ipel[2]][i] - 9*in[ipel[0]][i] - 3*in[ipel[3]][i] + 32)>>6, 0, 255 );
			out[opel[2]][i] = CLIP3( (5*(in[ipel[1]][i] +    in[ipel[2]][i])-   in[ipel[0]][i] -   in[ipel[3]][i] +  4)>>3, 0, 255 );
			out[opel[3]][i] = CLIP3( (19*in[ipel[1]][i] + 57*in[ipel[2]][i] - 3*in[ipel[0]][i] - 9*in[ipel[3]][i] + 32)>>6, 0, 255 );
		}
	}
}
