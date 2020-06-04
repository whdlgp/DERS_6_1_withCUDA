#ifndef __INCLUDE_UPSAMPLE_H__
#define __INCLUDE_UPSAMPLE_H__


#ifndef BYTE
#define BYTE unsigned char
#endif

void HorizontalLinearFilter_1D_half(BYTE *in, BYTE *out, int width, int padding_size);
void HorizontalLinearFilter_2D_half(BYTE **in, BYTE **out, int width, int height, int padding_size);
void HorizontalLinearFilter_1D_qpel(BYTE *in, BYTE *out, int width, int padding_size);
void HorizontalLinearFilter_2D_qpel(BYTE **in, BYTE **out, int width, int height, int padding_size);
void HorizontalAVCFilter_1D_half(BYTE *in, BYTE *out, int width, int padding_size);
void HorizontalAVCFilter_2D_half(BYTE **in, BYTE **out, int width, int height, int padding_size);
void HorizontalAVCFilter_1D_qpel(BYTE *in, BYTE *out, int width, int padding_size);
void HorizontalAVCFilter_2D_qpel(BYTE **in, BYTE **out, int width, int height, int padding_size);
void Horizontal6tapFilter_1D_half(BYTE *in, BYTE *out, int width, int padding_size, int coeff[6]);
void Horizontal6tapFilter_2D_half(BYTE **in, BYTE **out, int width, int height, int padding_size, int coeff[6]);
void HorizontalCubicFilter_1D_half(BYTE *in, BYTE *out, int width, int padding_size);
void HorizontalCubicFilter_2D_half(BYTE **in, BYTE **out, int width, int height, int padding_size);
void HorizontalCubicFilter_1D_qpel(BYTE *in, BYTE *out, int width, int padding_size);
void HorizontalCubicFilter_2D_qpel(BYTE **in, BYTE **out, int width, int height, int padding_size);

//void VerticalLinearFilter_1D_half(BYTE *in, BYTE *out, int width, int padding_size);
void VerticalLinearFilter_2D_half(BYTE **in, BYTE **out, int width, int height, int padding_size);
//void VerticalLinearFilter_1D_qpel(BYTE *in, BYTE *out, int width, int padding_size);
void VerticalLinearFilter_2D_qpel(BYTE **in, BYTE **out, int width, int height, int padding_size);
//void VerticalAVCFilter_1D_half(BYTE *in, BYTE *out, int width, int padding_size);
void VerticalAVCFilter_2D_half(BYTE **in, BYTE **out, int width, int height, int padding_size);
//void VerticalAVCFilter_1D_qpel(BYTE *in, BYTE *out, int width, int padding_size);
void VerticalAVCFilter_2D_qpel(BYTE **in, BYTE **out, int width, int height, int padding_size);
//void Vertical6tapFilter_1D_half(BYTE *in, BYTE *out, int width, int padding_size, int coeff[6]);
void Vertical6tapFilter_2D_half(BYTE **in, BYTE **out, int width, int height, int padding_size, int coeff[6]);
//void VerticalCubicFilter_1D_half(BYTE *in, BYTE *out, int width, int padding_size);
void VerticalCubicFilter_2D_half(BYTE **in, BYTE **out, int width, int height, int padding_size);
//void VerticalCubicFilter_1D_qpel(BYTE *in, BYTE *out, int width, int padding_size);
void VerticalCubicFilter_2D_qpel(BYTE **in, BYTE **out, int width, int height, int padding_size);

void DummyFilter_1D(BYTE *in, BYTE *out, int width, int padding_size);
void DummyFilter_2D(BYTE **in, BYTE **out, int width, int height, int padding_size);

#endif