#include <stdio.h>                          /* Sobel.c */
#include <math.h>

/*

   ANSWER: 
      High Threshold: 110
      Low Threshold: 40
      CL: 
         gcc -o sobel sobel.c -lm
         ./sobel garb34.pgm output_magnitude.pgm output_high.pgm output_low.pgm 110 40

*/


int pic[256][256];
int outpicx[256][256];
int outpicy[256][256];
int maskx[3][3] = {{-1,0,1},{-2,0,2},{-1,0,1}};
int masky[3][3] = {{1,2,1},{0,0,0},{-1,-2,-1}};
double ival[256][256],maxival;

int main(argc,argv) int argc; char **argv; {
      int i,j,p,q,mr,sum1,sum2;
      double hi_threshold, lo_threshold;
      FILE *fo1, *fo2, *fo3, *fp1, *fopen();
      char *foobar;


      // reading input file
      argc--; argv++;
      foobar = *argv;
      fp1=fopen(foobar,"rb");

      // scanning in and discarding pgm header
      char header[3];
      int width, height, maxval;
      fscanf(fp1, "%s %d %d %d%*c", header, &width, &height, &maxval);



      // magnitude output file
      argc--; argv++;
      foobar = *argv;
      fo1=fopen(foobar,"wb");

      // threshold high output file
      argc--; argv++;
      foobar = *argv;
      fo2=fopen(foobar,"wb");

      // threshold low output file
      argc--; argv++;
      foobar = *argv;
      fo3=fopen(foobar,"wb");

      
      
      // reading in threshold 1
      argc--; argv++;
      foobar = *argv;
      hi_threshold = atof(foobar);

      // reading in threshold 2
      argc--; argv++;
      foobar = *argv;
      lo_threshold = atof(foobar);



      // Write PGM header to output file
      fprintf(fo1, "P5\n");
      fprintf(fo1, "%d %d\n", 256, 256);
      fprintf(fo1, "255\n");

      fprintf(fo2, "P5\n");
      fprintf(fo2, "%d %d\n", 256, 256);
      fprintf(fo2, "255\n");

      fprintf(fo3, "P5\n");
      fprintf(fo3, "%d %d\n", 256, 256);
      fprintf(fo3, "255\n");


      for (i=0;i<256;i++)
      { for (j=0;j<256;j++)
         {
            pic[i][j]  =  getc (fp1);
            pic[i][j]  &= 0377;
         }
      }
        
      mr = 1;
      for (i=mr;i<256-mr;i++)
      { for (j=mr;j<256-mr;j++)
         {
            sum1 = 0;
            sum2 = 0;
            for (p=-mr;p<=mr;p++)
            {
               for (q=-mr;q<=mr;q++)
               {
                  sum1 += pic[i+p][j+q] * maskx[p+mr][q+mr];
                  sum2 += pic[i+p][j+q] * masky[p+mr][q+mr];
               }
            }
            outpicx[i][j] = sum1;
            outpicy[i][j] = sum2;
         }
      }

      maxival = 0;
      for (i=mr;i<256-mr;i++)
      { for (j=mr;j<256-mr;j++)
         {
            ival[i][j]=sqrt((double)((outpicx[i][j]*outpicx[i][j]) + (outpicy[i][j]*outpicy[i][j])));
            if (ival[i][j] > maxival)
               maxival = ival[i][j];

         }
      }



      for (i=0;i<256;i++)
         { for (j=0;j<256;j++)
            {
               ival[i][j] = (ival[i][j] / maxival) * 255;            
               fprintf(fo1,"%c",(char)((int)(ival[i][j])));
               
               // High threshold output
               if (ival[i][j] > hi_threshold) {
                  fprintf(fo2, "%c", (char)255);
               } else {
                  fprintf(fo2, "%c", (char)0);
               }

               // Low threshold output
               if (ival[i][j] > lo_threshold) {
                  fprintf(fo3, "%c", (char)255);
               } else {
                  fprintf(fo3, "%c", (char)0);
               }

            }
         }

}
