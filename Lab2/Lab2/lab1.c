#define _CRT_SECURE_NO_WARNINGS
#include <malloc.h>
#include <stdio.h> 
#include <stdint.h> 
#include <conio.h>
#include <math.h>
#pragma pack(push, 2) 
struct 
{ 
  uint16_t signature;
  uint32_t fileSize;
  uint16_t foo1;
  uint16_t foo2;
  uint32_t fileHeadSize;
  uint32_t headerSize; 
  uint32_t width;
  uint32_t height; 
  uint16_t colorDimensions;
  uint16_t density;
  uint32_t zipType;
  uint32_t massiveLength;
  uint32_t horResolution;
  uint32_t verResolution;
  uint32_t colorsQuantity;
  uint32_t shadesQuantity;
} fileStruck, fileStruck2; 
#pragma pack(pop) 
#define SQRT_2PI 2.506628274631

int winsize(double);
double Gauss(double, double);
double *GKernel_1D(int, double);

int main(int argc, char** argv) 
{ 
	char *InputName, *OutputName,*src, *dest;
		FILE *InputFile, *OutputFile;
		int offset;

	InputName = "in.bmp"; OutputName = "out.bmp";
	if (0 != fopen_s(&InputFile, InputName, "rb"))
	{
		printf("File not found\n");
		_getch();
		return 1;
	}
	fread_s(&fileStruck, sizeof(fileStruck), sizeof(fileStruck), 1, InputFile);
	fileStruck2 = fileStruck;
	if (fileStruck.signature != 19778 || fileStruck.density != 24)
    { 
		printf("Wrong file type\n"); 
		_getch(); 
		return 1; 
    } 
	src = (char*)malloc(fileStruck.massiveLength);
	dest = (char*)malloc(fileStruck.massiveLength);
    fread_s(src, fileStruck.massiveLength, fileStruck.massiveLength, 1, InputFile); 
    fclose(InputFile);
	
	if ((fileStruck.width * 3) % 4 == 0) { offset = 0; }
		else { offset = 4 - (fileStruck.width * 3) % 4; }
	
	Gblur(src, dest, fileStruck.width*3 + offset, fileStruck.height, 4.0); 
	
	fopen_s(&OutputFile, OutputName, "wb");
	fwrite(&fileStruck, sizeof(fileStruck), 1, OutputFile); 
    fwrite(dest, fileStruck.massiveLength, 1, OutputFile); 
    fclose(OutputFile);
    free(src); 

    return 0; 
} 

int Gblur(char* src, char* dest, int width, int height, double sigma) // sigma - радиус размытия
{
	if (!src || !dest) return -1;

	int	row, col, col_r, col_g, col_b, winsize, halfsize, k, offset, count = 0, rows, count1, count2, count3;
	double  row_g, row_b, row_r, col_all;
	unsigned char  r_r, r_b, r_g, c_all;
	unsigned char *tmp = NULL;
	double *kernel; 

		winsize = win_size(sigma); 
		kernel = GKernel_1D(winsize, sigma); //  вычислил окно 1 раз для данного радиуса размытия
		winsize *= 3; // размер окна в байтах
		halfsize = winsize / 2;

	if ((tmp = (char*)calloc(width * height, sizeof(char))) == NULL) // промежуточное изображение, размытое по строкам
		return -1;
	// вычисления по строкам
	for (row = 0; row < height; row++)
	{ 
		col_r = 0;
		col_g = 1;
		col_b = 2;
		for (rows = 0; rows < width; rows+=3)
		{
			row_r = row_g = row_b = 0.0;
			count1 = count2 = count3 = 0;

			for (k = 1; k < winsize; k += 3)
			{
				if ((k + col_r - halfsize >= 0) && (k + col_r - halfsize < width))
				{
					r_r = *(src + row * width + col_r + k - halfsize);
					row_r += (int)(r_r)* kernel[count1];
					count1++;
				}

				if ((k + col_g - halfsize >= 0) && (k + col_g - halfsize < width))
				{
					r_g = *(src + row * width + col_g + k - halfsize);
					row_g += (int)(r_g)* kernel[count2];
					count2++;
				}

				if ((k + col_b - halfsize >= 0) && (k + col_b - halfsize < width))
				{
					r_b = *(src + row * width + col_b + k - halfsize);
					row_b += (int)(r_b)* kernel[count3];
					count3++;
				}
			}

			*(tmp + row * width + col_r) = (unsigned char)(ceil(row_r));
			*(tmp + row * width + col_g) = (unsigned char)(ceil(row_g));
			*(tmp + row * width + col_b) = (unsigned char)(ceil(row_b));
			col_r += 3;
			col_g += 3;
			col_b += 3; 
		}
		}

	// вычисления по столбцам
	winsize /= 3;
	halfsize = winsize / 2;
	for (col = 0; col < width; col++)
		for (row = 0; row < height; row++)
		{
		col_all = 0.0;
		for (k = 0; k < winsize; k++) 
			if ((k + row - halfsize >= 0) && (k + row - halfsize < height))
			{
			c_all = *(tmp + (row + k - halfsize) * width + col);
			col_all += ((int)c_all) * kernel[k];
			}
		*(dest + row * width + col) = (unsigned char)(ceil(col_all));
		}
		
	free(tmp);
	free(kernel);

	return 0;
}

double* GKernel_1D(int win_size, double sigma)
{
	int wincenter, x;
	double *kern, sum = 0.0;
	wincenter = win_size / 2;
	kern = (double*)calloc(win_size, sizeof(double));

	for (x = 0; x < wincenter + 1; x++)
	{
		kern[wincenter - x] = kern[wincenter + x] = Gauss(sigma, x);
		sum += kern[wincenter - x] + ((x != 0) ? kern[wincenter + x] : 0.0);
	}
	for (x = 0; x < win_size; x++)
		kern[x] /= sum;

	return kern;
}

double Gauss(double sigma, double x)
{
	return exp(-(x * x) / (2.0 * sigma * sigma)) / (sigma * SQRT_2PI);
}

int win_size(double sigma)
{
	double SIGMA_FACTOR = 3; // пр-ло 3-х сигм
	return (1 + (((int)ceil(SIGMA_FACTOR * sigma)) * 2)); // + 1 чтоб был центр ядра
}

