/* Aramis' Bitmap Lib - Copyright (c) 2017 Aramis Hornung Moraes */

#define AHMBMP_VERSION "0.9.6"

#include <stdio.h>
#include "ahmbmp.h"
#include "ahmycbcr.h"
#include "ahmlux.h"
#include "ahmskin.h"
#include "ahmfilters.h"
#include "string.h"
#include "math.h"

/* debug stuff for MSVC */
#ifdef WIN32
#ifdef _MSC_VER
#define _CRTDBG_MAP_ALLOC
#include <stdlib.h>
#include <crtdbg.h>
#endif
#endif
#ifdef WIN32
#include "windows.h"
#endif



// return the mask's relative position considering as reference the origin/target pixel and its size in layers
int local_sp(int global_p, int taget_p, int k)
{
	int start_edge = (taget_p - k<0 ? 0 : taget_p - k);
	return ((global_p - start_edge) >(k * 2) ? (k * 2) : (global_p - start_edge));
}


/* int x = pixel x position
   int y = pixel y position
   int K = amoeba maximumlayere size
   int L = amoeba maximum variance threshold*/
void amoeba_process_pixel(ahm_bitmap *img, ahm_bitmap *trgt_img, int x, int y, int K, int L)
{
	const unsigned int maxWidth = img->Width;
	const unsigned int maxHeight = img->Height;
	ahm_bitmap *amoeba_mask = create_ahmBitmap((K*2)+1, (K*2) + 1);
	int i, j; i = j = 0;

	// media adjacente ao pixel alvo (M1)
	/*
	int M1 = 0;
	
	for(j = (y-1 < 0? 0 : y-1); (j <= y+1 && j <= maxHeight); ++j)
	{
		for(i = (x-1 < 0? 0 : x-1); (i <= x+1 && i <= maxWidth); ++i) // iterate the row from left to right
		{
			if(!(i == x && j == y)) // sum everything except the very target pixel
			{
				M1 += get_pixel(img, i, j).r_;
				set_pixel(amoeba_mask, i, j, 255, 255, 255); // the core is always marked
			}
		}
	}
	M1 = M1/8;
	*/

	set_pixel(amoeba_mask, K, K, 255, 255, 255); // the core is always marked

	// draw the amoeba mask
	for(int k = 1; k < K; ++k) // iterate amoeba layers
	{
		// for each line in the window
		for(j = (y-k<0?0:y-k); j <= (y+k>maxHeight?maxHeight:y+k); ++j)
		{

			// if its the window cap or botton, process the entire row, 
			if(j == 0 || j == y-k || j == y+k || j == maxHeight)
			{
				for(i = (x-k<0?0:x-k); i <= (x+k>maxWidth?maxWidth:x+k); ++i) // iterate the row from left to right
				{

					if(j == 0 || j == y-k) // if its cap, analyze adjacent bottom pixels
					{
						int x1;
						x1 = i;
						if(i == (x-k)) // far left in top, read inner corn of the structure
							x1 = i+1;
						else if (i == (x+k)) // far right in top, read inner corn of the structure
							x1 = i-1;
						

						//ABS(Pd[n]-Pd[n-1])<=ABS(M1-L)
						int Pdn = get_pixel(img, i, j).r_;
						//int Pdn1 = get_pixel(img, x1, j+1).r_;
						int Pdn1 = get_pixel(img, x, y).r_;
						int r_val =  get_pixel(amoeba_mask, local_sp(x1,x,K), local_sp(j + 1, y, K)).r_;

						if(r_val == 255 &&
						abs(Pdn-Pdn1)<=abs(L)) {
							set_pixel(amoeba_mask, local_sp(i, x, K), local_sp(j, y, K), 255, 255, 255); // mark it
						}
					}
					if(j == maxHeight|| j == y+k) // if its bottom, analyze adjacent upper pixels
					{
						int x1;
						x1 = i;
						if(i == (x-k)) // far left in bottom, read inner corn of the structure
							x1 = i+1;
						else if (i == (x+k)) // far right in bottom, read inner corn of the structure
							x1 = i-1;
						

						int Pdn = get_pixel(img, i, j).r_;
						//int Pdn1 = get_pixel(img, x1, j-1).r_;
						int Pdn1 = get_pixel(img, x, y).r_;
						int r_val =  get_pixel(amoeba_mask, local_sp(x1, x, K), local_sp(j - 1, y, K)).r_;

						if(r_val == 255 &&
						abs(Pdn-Pdn1)<=abs(L)) {
							set_pixel(amoeba_mask, local_sp(i, x, K), local_sp(j, y, K), 255, 255, 255); // mark it
						}
					}
				}
			}
			else // otherwise, only its far left and far right side pixels
			{
				// left
				{
					i = (x-k<0?0:x-k); // far left pixel in the layer-window

					//ABS(Pd[n]-Pd[n-1])<=ABS(M1-L)
					int Pdn = get_pixel(img, i, j).r_;
					//int Pdn1 = get_pixel(img, i+1, j).r_;
					int Pdn1 = get_pixel(img, x, y).r_;
					int r_val =  get_pixel(amoeba_mask, local_sp(i+1, x, K), local_sp(j, y, K)).r_;

					if(r_val == 255 &&
					abs(Pdn-Pdn1)<=abs(L)) {
						set_pixel(amoeba_mask, local_sp(i, x, K), local_sp(j, y, K), 255, 255, 255); // mark it
					}

				}
				//right
				{
					i = (x+k>maxWidth?maxWidth:x+k); // far left pixel in the layer-window

					//ABS(Pd[n]-Pd[n-1])<=ABS(M1-L)
					int Pdn = get_pixel(img, i, j).r_;
					//int Pdn1 = get_pixel(img, i-1, j).r_;
					int Pdn1 = get_pixel(img, x, y).r_;
					int r_val =  get_pixel(amoeba_mask, local_sp(i-1, x, K), local_sp(j, y, K)).r_;

					if(r_val == 255 &&
					abs(Pdn-Pdn1)<=abs(L)) {
						set_pixel(amoeba_mask, local_sp(i, x, K), local_sp(j, y, K), 255, 255, 255); // mark it
					}
				}
			} // top/bottom right/left layer sides
		} // window lines
	} // k layers

	
	// by now, we have our amoeba mask, a bitmap wehre with the exact same
	// size of the original image where the white pixels represent the amoeba
	// for the pixel P<x, y>,
	// we can either calculate mean or median and apply the filter only
	// considering those white pixels (amoeba)
	unsigned int pixel_count = 0;
	unsigned int mean_color = 0;
	for(i = 0; i < maxWidth; ++i)
	{
		for(j = 0; j < maxHeight; ++j)
		{
			if(get_pixel(amoeba_mask, local_sp(i, x, K), local_sp(j, y, K)).r_ == 255)
			{
				++pixel_count;
				mean_color += (unsigned int)get_pixel(img, i, j).r_;
			}
		}
	}
	if(pixel_count)
		mean_color = mean_color/pixel_count;


	set_pixel(trgt_img, x, y, mean_color, mean_color, mean_color);

	destroy_ahmBitmap(amoeba_mask);
}

void main_process_image(char *file_path)
{
	ahm_bitmap *myBmp;
	ahm_bitmap *ycbcr;
	ahm_bitmap *ycbcr_y;
	ahm_bitmap *ycbcr_cb;
	ahm_bitmap *ycbcr_cr;
	ahm_bitmap *trgt_img;
	char file_output_name[4096];

	/* iterators */
	unsigned int i = 0;
	unsigned int j = 0;

    myBmp = create_bmp_from_file(file_path);
    if(myBmp == 0)
    {
        destroy_ahmBitmap(myBmp);
		printf("Error: File is not supported. Only 24bit bitmaps can be processed.\n");
		return;
    }
	trgt_img=create_ahmBitmap(myBmp->Width,myBmp->Height);
	/* allocate image structures */
	ycbcr=create_ahmBitmap(myBmp->Width,myBmp->Height);
	ycbcr_y = create_ahmBitmap(myBmp->Width, myBmp->Height);
	ycbcr_cb = create_ahmBitmap(myBmp->Width, myBmp->Height);
	ycbcr_cr = create_ahmBitmap(myBmp->Width, myBmp->Height);

	/* start processing */
	create_ycbcr(myBmp, ycbcr);
	create_ycbcr_channels(myBmp, ycbcr_y, ycbcr_cb, ycbcr_cr);
	//ycbcr_to_rgb(ycbcr,myBmp); // testing

	for(i = 0; i < ycbcr_y->Width; ++i)
	{
		for(j = 0; j < ycbcr_y->Height; ++j)
		{
			printf("%f%%\n", (float)((float)i/(float)ycbcr->Width)*100.0f);
			printf("%f%%", (float)((float)j/(float)ycbcr->Height)*100.0f);
			amoeba_process_pixel(ycbcr_y, trgt_img, i, j, 4, 20);
		}
	}

	save_bmp(trgt_img, "amoeba_test.bmp"); // after all, restore the processed image to colors


	destroy_ahmBitmap(ycbcr);
	destroy_ahmBitmap(ycbcr_y);
	destroy_ahmBitmap(ycbcr_cb);
	destroy_ahmBitmap(ycbcr_cr);
	destroy_ahmBitmap(myBmp);
}

int main(int argc, char *argv[])
{
#ifdef _MSC_VER
#ifdef WIN32
	_CrtSetDbgFlag(_CRTDBG_ALLOC_MEM_DF | _CRTDBG_LEAK_CHECK_DF);
#endif

#endif
 

	if(argc > 1)
        main_process_image(argv[1]);
	else
	{
		printf("Usage: parameter 1 = 24bits Bitmap file\n");
		fflush(stdin);
		getchar();
	}

#ifdef _MSC_VER
#ifdef WIN32
	_CrtSetReportMode(_CRT_ERROR, _CRTDBG_MODE_DEBUG);
#endif
#endif
	return 0;
}
/*
#if defined(_WIN32)
int APIENTRY _WinMain(HINSTANCE hInstance, HINSTANCE hPevInstance, LPSTR lpCmdLine, int nCmdShow)
{
	return main( __argc, __argv);
}
#endif*/
