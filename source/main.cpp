/*	Aramis' Bitmap Lib
copyright (c) 2017 by Aramis Hornung Moraes
This file is part of the ahmBitmap project
read ahmbmp.h for copying conditions.
*/


#define AHMBMP_VERSION "0.9.8"

#include <stdio.h>
#include "ahmbmp.h"
#include "ahmycbcr.h"
#include "ahmlux.h"
#include "ahmskin.h"
#include "ahmfilters.h"
#include "string.h"
#include "math.h"
#include "time.h"

/* debug stuff for MSVC */
#ifdef _WIN32_debug
#ifdef _MSC_VER
#define _CRTDBG_MAP_ALLOC
#include <stdlib.h>
#include <crtdbg.h>
#endif
#endif
#ifdef _WIN32
#include "windows.h"
#endif

// return the mask's relative position considering as reference the origin/target pixel and its size in layers
int local_sp(int global_p, int target_p, int k)
{
	//int start_edge = (target_p - k<0 ? 0 : target_p - k);
	//return ((global_p - start_edge) >((k * 2)) ? ((k * 2)) : (global_p - start_edge));

	int start_edge = global_p - k;
	return target_p - start_edge;
}


/* int x = pixel x position
   int y = pixel y position
   int K = amoeba maximumlayere size
   int L = amoeba maximum variance threshold*/
void amoeba_process_pixel(ahm_bitmap *img, ahm_bitmap *trgt_img, int x, int y, int K, int L)
{
	const unsigned int maxWidth = img->Width;
	const unsigned int maxHeight = img->Height;
	ahm_bitmap *amoeba_mask = create_ahmBitmap((K*2) + 1, (K*2) + 1);
	int i, j; i = j = 0;

	set_pixel(amoeba_mask, K, K, 255, 255, 255); // the core is always marked

	// draw the amoeba mask
	for(int k = 1; k < K; ++k) // iterate amoeba layers
	{
		// for each line in the window
		for(j = (y-k<0?0:y-k); j <= (y+k>maxHeight-1?maxHeight-1:y+k); ++j)
		{
			// if its the window cap or botton, process the entire row, 
			if(j == 0 || j == y-k || j == y+k || j == maxHeight-1)
			{
				for(i = (x-k<0?0:x-k); i <= (x+k>maxWidth-1?maxWidth-1:x+k); ++i) // iterate the row from left to right
				{
					if(j == 0 || j == y-k) // if its cap, analyze adjacent bottom pixels
					{
						int x1;
						x1 = i;
						if(i == (x-k) || i == 0) // far left in top, read inner corn of the structure
							x1 = i+1;
						else if (i == (x+k) || i == maxWidth - 1) // far right in top, read inner corn of the structure
							x1 = i-1;

						int Pdn = get_pixel(img, i, j).r_;
						int Pdn1 = get_pixel(img, x, y).r_;
						int r_val =  get_pixel(amoeba_mask, local_sp(x1,x,K), local_sp(((j + 1)>maxHeight-1? maxHeight - 1:j+1), y, K)).r_;

						if(r_val == 255 &&
						abs(Pdn-Pdn1)<=abs(L)) {
							set_pixel(amoeba_mask, local_sp(i, x, K), local_sp(j, y, K), 255, 255, 255); // mark it
						}
					}
					if(j == maxHeight-1|| j == y+k) // if its bottom, analyze adjacent upper pixels
					{
						int x1;
						x1 = i;
						if(i == (x-k) || i == 0) // far left in bottom, read inner corn of the structure
							x1 = i+1;
						else if (i == (x+k) || i == maxWidth - 1) // far right in bottom, read inner corn of the structure
							x1 = i-1;
						

						int Pdn = get_pixel(img, i, j).r_;
						//int Pdn1 = get_pixel(img, x1, j-1).r_;
						int Pdn1 = get_pixel(img, x, y).r_;
						int r_val =  get_pixel(amoeba_mask, local_sp(x1, x, K), local_sp((j - 1<0?0:j-1), y, K)).r_;

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

					int Pdn = get_pixel(img, i, j).r_;
					int Pdn1 = get_pixel(img, x, y).r_;
					int r_val =  get_pixel(amoeba_mask, local_sp(((i+1)<0?0:i+1), x, K), local_sp(j, y, K)).r_;

					if(r_val == 255 &&
					abs(Pdn-Pdn1)<=abs(L)) {
						set_pixel(amoeba_mask, local_sp(i, x, K), local_sp(j, y, K), 255, 255, 255); // mark it
					}

				}
				//right
				{
					i = (x+k>maxWidth-1?maxWidth-1:x+k); // far left pixel in the layer-window

					int Pdn = get_pixel(img, i, j).r_;
					int Pdn1 = get_pixel(img, x, y).r_;
					int r_val =  get_pixel(amoeba_mask, local_sp((i-1>maxWidth -1? maxWidth - 1 :i-1), x, K), local_sp(j, y, K)).r_;

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

void main_process_image(char *file_path, int k, int l)
{
	ahm_bitmap *myBmp;
	ahm_bitmap *ycbcr;
	ahm_bitmap *ycbcr_y;
	ahm_bitmap *ycbcr_cb;
	ahm_bitmap *ycbcr_cr;
	ahm_bitmap *trgt_img;
	char *input_file_name;
	char file_output_name[4096];

	

	/* iterators */
	unsigned int i = 0;
	unsigned int j = 0;

    myBmp = create_bmp_from_file(file_path);
	input_file_name = strtok(file_path, ".");

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
	
	clock_t start_t, end_t, total_t;
	int iaux = 10;
	start_t = clock();
	printf("Estimating time...");
	int seq = 0;
	for(i = 0; i < ycbcr_y->Width; ++i)
	{
		for(j = 0; j < ycbcr_y->Height; ++j)
		{
			amoeba_process_pixel(ycbcr_y, trgt_img, i, j, k, l);
			//sprintf(file_output_name, "%s_prog_%d.bmp\0", input_file_name, seq);
			//++seq;
			//save_bmp(trgt_img, file_output_name); // after all, restore the processed image to colors
		}
		if (iaux == 0)
		{
			end_t = clock();
			total_t = (float)((((end_t - start_t)/10.0f) / CLOCKS_PER_SEC)*(ycbcr_y->Width - i));
#ifdef _WIN32
			system("cls");
#endif
			printf("Processing image K = %d L = %d> %d%% completed\nTime remaining: %d minut%s and %d second%s\n",
				k, l,
													(int)(((float)i/(float)ycbcr_y->Width)*100.0f),
													total_t/60,
													((total_t / 60)>1?"es":"e"),
													total_t%60,
													((total_t % 60)>1 ? "s" : ""));
			start_t = clock();
			iaux = 10;
		}
		else --iaux;
		
	}
	//sprintf(file_output_name, "%s_Filtered_Y_Channel.bmp\0", input_file_name);
	//save_bmp(trgt_img, file_output_name); // after all, restore the processed image to colors

	printf("\nWWConverting image to RGB...\n");

	for (i = 0; i < ycbcr->Width; ++i)
	{
		for (j = 0; j < ycbcr->Height; ++j)
		{
			set_pixel(ycbcr, i, j, get_pixel(trgt_img, i, j).r_,
				get_pixel(ycbcr, i, j).g_,
				get_pixel(ycbcr, i, j).b_);
		}

	}

	ycbcr_to_rgb(ycbcr, myBmp);
	sprintf(file_output_name, "%s_mean_amoeba_filtered_K%d_L%d.bmp\0", input_file_name, k, l);
	save_bmp(myBmp, file_output_name); // after all, restore the processed image to colors

	printf("\nDone!");

	destroy_ahmBitmap(trgt_img);
	destroy_ahmBitmap(ycbcr_cr);
	destroy_ahmBitmap(ycbcr_cb);
	destroy_ahmBitmap(ycbcr_y);
	destroy_ahmBitmap(ycbcr);
	destroy_ahmBitmap(myBmp);
}

int main(int argc, char *argv[])
{
#ifdef _WIN32_debug
#ifdef _MSC_VER
#ifdef WIN32
	_CrtSetDbgFlag(_CRTDBG_ALLOC_MEM_DF | _CRTDBG_LEAK_CHECK_DF);
#endif
#endif
#endif
	int k = (argc > 2? atoi(argv[2]): 4);
	int l = (argc > 3 ? atoi(argv[3]) : 24);
	if(argc > 1)
        main_process_image(argv[1], k, l);
	else
	{
		printf("Usage:\nparameter 1 = 24bits Bitmap file\nparameter 2 = Amoeba maximum possible size (radius) in pixels (default %d)\nparameter 3 = variance threshold (default %d)\n", k, l);
		fflush(stdin);
		getchar();
	}
#ifdef _WIN32_debug
#ifdef _MSC_VER
#ifdef WIN32
	_CrtSetReportMode(_CRT_ERROR, _CRTDBG_MODE_DEBUG);
#endif
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
