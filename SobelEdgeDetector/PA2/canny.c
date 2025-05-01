#include <stdio.h>                  /*  Marr-Hildreth.c  (or marrh.c) */
#include <math.h>
#define  PICSIZE 256
#define  MAXMASK 100
#define HISTSIZE 256


/*

   ANSWER: 
      High Threshold: 110
      Low Threshold: 40
      CL: 
         gcc -o canny canny.c -lm
         ./canny garb34.pgm magnitude.pgm peaks.pgm final.pgm 1.0 100 35
         ./canny garb34.pgm magnitude.pgm peaks.pgm final.pgm 1.0 0.023

*/


int    pic[PICSIZE][PICSIZE];
double outpic1[PICSIZE][PICSIZE];
double outpic2[PICSIZE][PICSIZE];
int    edgeflag[PICSIZE][PICSIZE];
double maskx[MAXMASK][MAXMASK];
double masky[MAXMASK][MAXMASK];
double convx[PICSIZE][PICSIZE];
double convy[PICSIZE][PICSIZE];
double ival[PICSIZE][PICSIZE]; // our magnitude array
int cand[PICSIZE][PICSIZE];  
int peaks[PICSIZE][PICSIZE];  // New array for storing peaks
int final[PICSIZE][PICSIZE];  
int histogram[HISTSIZE];


int main(argc,argv) int argc; char **argv;{
      int     i,j,p,q,s,t,mr,centx,centy;
      double  maskval,sum1,sum2,sig,maxival,minival,maxval,hi_threshold, lo_threshold, percent;
      FILE    *fo1,*fo2,*fo3,*fp1, *fopen();
      char    *foobar;

      argc--; argv++;
      foobar = *argv;
      fp1=fopen(foobar,"rb");

      // scanning in and discarding pgm header
      char header[3];
      int width, height, temp_maxval;
      fscanf(fp1, "%s %d %d %d%*c", header, &width, &height, &temp_maxval);


      argc--; argv++;
      foobar = *argv;
      fo1=fopen(foobar,"wb");

      argc--; argv++;
      foobar = *argv;
      fo2=fopen(foobar,"wb");

      argc--; argv++;
      foobar = *argv;
      fo3=fopen(foobar,"wb");

      argc--; argv++;
      foobar = *argv;
      sig = atof(foobar);
      
      /*
      // reading in threshold 1
      argc--; argv++;
      foobar = *argv;
      hi_threshold = atof(foobar);

      // reading in threshold 2
      argc--; argv++;
      foobar = *argv;
      lo_threshold = atof(foobar);
      */
      
      
      // reading in percent
      argc--; argv++;
      foobar = *argv;
      percent = atof(foobar);
      


      // Write PGM header to output file
      fprintf(fo1, "P5\n");
      fprintf(fo1, "%d %d\n", PICSIZE, PICSIZE);
      fprintf(fo1, "255\n");

      // Write PGM header to output file
      fprintf(fo2, "P5\n");
      fprintf(fo2, "%d %d\n", PICSIZE, PICSIZE);
      fprintf(fo2, "255\n");

      // Write PGM header to output file
      fprintf(fo3, "P5\n");
      fprintf(fo3, "%d %d\n", PICSIZE, PICSIZE);
      fprintf(fo3, "255\n");
      


      mr = (int)(sig * 3);
      centx = (MAXMASK / 2);
      centy = (MAXMASK / 2);

      for (i=0;i<256;i++)
      { for (j=0;j<256;j++)
         {
         pic[i][j]  =  getc (fp1);
         }
      }

      // gradient 
      for (p=-mr;p<=mr;p++)
      {  for (q=-mr;q<=mr;q++)
         {
            maskval = (q* (exp(-1*(((p*p)+(q*q))/(2*(sig*sig))))));
            (maskx[p+centy][q+centx]) = maskval;

            maskval = (p* (exp(-1*(((p*p)+(q*q))/(2*(sig*sig))))));
            (masky[p+centy][q+centx]) = maskval;
         }
      }

      // convolution
      for (i=mr;i<=255-mr;i++)
      { for (j=mr;j<=255-mr;j++)
         {
            sum1 = 0.0;
            sum2 = 0.0;
            for (p=-mr;p<=mr;p++)
            {
               for (q=-mr;q<=mr;q++)
               {
                  sum1 += pic[i+p][j+q] * maskx[p+centy][q+centx];
                  sum2 += pic[i+p][j+q] * masky[p+centy][q+centx];
               }
            }
            outpic1[i][j] = sum1;
            convx[i][j] = sum1;
            outpic2[i][j] = sum2;
            convy[i][j] = sum2;
         }
      }

      // magnitude
      maxival = 0;
      for (i=mr;i<256-mr;i++)
      { for (j=mr;j<256-mr;j++)
         {
            ival[i][j]=sqrt((double)((outpic1[i][j]*outpic1[i][j]) + (outpic2[i][j]*outpic2[i][j])));
            if (ival[i][j] > maxival)
               maxival = ival[i][j];
         }
      }
      
      // peaks
      for (i=mr;i<256-mr;i++) {
         for (j=mr;j<256-mr;j++) {

            if ((convx[i][j])==0.0){
               convx[i][j]=0.00001;
            }
            double slope = convy[i][j]/convx[i][j];
            // horizontal 
            if((slope <= 0.4142) && (slope > -0.4142)) {
               if((ival[i][j] > ival[i][j-1]) && (ival[i][j] > ival[i][j+1])) {
                  cand[i][j] = 255;
               }
            }
            // Diagonal edges (positive slope)
            else if((slope <= 2.4142) && (slope > 0.4142)) {
               if((ival[i][j] > ival[i-1][j-1]) && (ival[i][j] > ival[i+1][j+1])) {
                  cand[i][j] = 255;
               }
            }
            // Diagonal (negative slope)
            else if((slope <= -0.4142) && (slope > -2.4142)) {
               if((ival[i][j] > ival[i+1][j-1]) && (ival[i][j] > ival[i-1][j+1])) {
                  cand[i][j] = 255;
               }
            }
            // vertical
            else {
               if((ival[i][j] > ival[i-1][j]) && (ival[i][j] > ival[i+1][j])) {
                  cand[i][j] = 255;
               }
            }
         }
      }

      
      // copying over peaks array since real peaks array changed during thresholding
      for (i=0; i<PICSIZE; i++) {
         for (j=0; j<PICSIZE; j++) {
             peaks[i][j] = cand[i][j];
         }
      }
          
      
      // code to automatically get high and low thresholds     
      // Compute histogram of scaled magnitudes
      for(i=mr; i<256-mr; i++) {
            for(j=mr; j<256-mr; j++) {
               // Note: ival has already been scaled to 0-255 range
               histogram[(int)ival[i][j]]++;
            }
      }
      // Calculate cutoff based on image size
      int cutoff = (int)(percent * (PICSIZE-2*mr) * (PICSIZE-2*mr));
      int areaOfTops = 0;
      int hi;

      // Find threshold by accumulating histogram from top
      for(i=HISTSIZE-1; i>=0; i--) {
         areaOfTops += histogram[i];
         if(areaOfTops > cutoff) {
            break;
         }
      }
      hi_threshold = i;  // This is our high threshold
      lo_threshold = hi_threshold * 0.35;  // Low threshold is 35% of high threshold
      
     

      
      // 1st pass, all above high
      for (i=mr;i<256-mr;i++) {
         for (j=mr;j<256-mr;j++) {
            if (cand[i][j]==255) {
               if (ival[i][j]>hi_threshold){
                  cand[i][j]=0;
                  final[i][j]=255;
               }
               else if (ival[i][j]<lo_threshold){
                  cand[i][j]=0;
                  final[i][j]=0;
               }
            }
         }
      }
      // second pass, between low and high but connected to high
      int moretodo = 1;
      while (moretodo==1) {
         moretodo = 0;
         for (i=mr;i<256-mr;i++) {
            for (j=mr;j<256-mr;j++) {
               if (cand[i][j]==255 && final[i][j]!=255){
                  for (int p=-1; p<=1; p++) {
                     for (int q=-1; q<=1; q++) {
                        if (final[i+p][j+q]==255){
                           cand[i][j]=0;
                           final[i][j]=255;
                           moretodo = 1;
                        }
                     }
                  }
               }

            }
         }
      }  
      

      for (i=0;i<256;i++) 
         { for (j=0;j<256;j++)
            {  
               // magnitude output file
               ival[i][j] = (ival[i][j] / maxival) * 255;            
               fprintf(fo1,"%c",(char)((int)(ival[i][j])));

               // peaks output file
               fprintf(fo2,"%c",(char)((int)(peaks[i][j])));
               
               // final output file
               fprintf(fo3,"%c",(char)((int)(final[i][j])));

            }
         }
  }

