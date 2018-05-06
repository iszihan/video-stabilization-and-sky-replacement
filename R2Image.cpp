//Source file for image class



// Include files 

#include "R2/R2.h"
#include "R2Pixel.h"
#include "R2Image.h"
#include "svd.h"
#include "iostream"
#include "algorithm"
#include "cmath"
#include "vector"

using namespace std;


////////////////////////////////////////////////////////////////////////
// Constructors/Destructors
////////////////////////////////////////////////////////////////////////


R2Image::
R2Image(void)
  : pixels(NULL),
    npixels(0),
    width(0), 
    height(0)
{
}



R2Image::
R2Image(const char *filename)
  : pixels(NULL),
    npixels(0),
    width(0), 
    height(0)
{
  // Read image
  Read(filename);
}



R2Image::
R2Image(int width, int height)
  : pixels(NULL),
    npixels(width * height),
    width(width), 
    height(height)
{
  // Allocate pixels
  pixels = new R2Pixel [ npixels ];
  assert(pixels);
}



R2Image::
R2Image(int width, int height, const R2Pixel *p)
  : pixels(NULL),
    npixels(width * height),
    width(width), 
    height(height)
{
  // Allocate pixels
  pixels = new R2Pixel [ npixels ];
  assert(pixels);

  // Copy pixels 
  for (int i = 0; i < npixels; i++) 
    pixels[i] = p[i];
}



R2Image::
R2Image(const R2Image& image)
  : pixels(NULL),
    npixels(image.npixels),
    width(image.width), 
    height(image.height)
    
{
  // Allocate pixels
  pixels = new R2Pixel [ npixels ];
  assert(pixels);

  // Copy pixels 
  for (int i = 0; i < npixels; i++) 
    pixels[i] = image.pixels[i];
}



R2Image::
~R2Image(void)
{
  // Free image pixels
  if (pixels) delete [] pixels;
}



R2Image& R2Image::
operator=(const R2Image& image)
{
  // Delete previous pixels
  if (pixels) { delete [] pixels; pixels = NULL; }

  // Reset width and height
  npixels = image.npixels;
  width = image.width;
  height = image.height;

  // Allocate new pixels
  pixels = new R2Pixel [ npixels ];
  assert(pixels);

  // Copy pixels 
  for (int i = 0; i < npixels; i++) 
    pixels[i] = image.pixels[i];

  // Return image
  return *this;
}


void R2Image::
svdTest()
{
  int n=4;
  int n2=2*n;
  int coord [6][n];
  //pass input for the coordinates of point correspondences
  for(int i=1;i<=n;i++){
    cout<<"Please enter the x coordinate of your "<<i<<"th point in image 1.";
    cin >> coord[1][i];
    cout<<"Please enter the y coordinate of your "<<i<<"th point in image 1.";
    cin >> coord[2][i];
    cout<<"Please enter the x coordinate of your "<<i<<"th point in image 2.";
    cin >> coord[4][i];
    cout<<"Please enter the y coordinate of your "<<i<<"th point in image 2.";
    cin >> coord[5][i];
    coord[3][i]=1;
    coord[6][i]=1;
  }

  // build the 2nx9 matrix of equations
  double** linEquations = dmatrix(1,n2,1,9);
  for(int j=1;j<=n2;j+=2){//calculate the two equations for each point correspondence
    int x=(j+1)/2;
    linEquations[j][1] = 0;
    linEquations[j][2] = 0;
    linEquations[j][3] = 0;
    linEquations[j][4] = -coord[1][x];
    linEquations[j][5] = -coord[2][x];
    linEquations[j][6] = -coord[3][x];
    linEquations[j][7] = coord[1][x]*coord[5][x];
    linEquations[j][8] = coord[2][x]*coord[5][x];
    linEquations[j][9] = coord[3][x]*coord[5][x];

    linEquations[j+1][1] = coord[1][x];
    linEquations[j+1][2] = coord[2][x];
    linEquations[j+1][3] = coord[3][x];
    linEquations[j+1][4] = 0;
    linEquations[j+1][5] = 0;
    linEquations[j+1][6] = 0;
    linEquations[j+1][7] = -coord[1][x]*coord[4][x];
    linEquations[j+1][8] = -coord[2][x]*coord[4][x];
    linEquations[j+1][9] = -coord[3][x]*coord[4][x];
  }


    // compute the SVD
    double singularValues[10]; // 1..9
    double** nullspaceMatrix = dmatrix(1,9,1,9);
    svdcmp(linEquations, n2, 9, singularValues, nullspaceMatrix);

    // get the result
    printf("\n Singular values: %f, %f, %f, %f, %f, %f, %f, %f, %f\n",singularValues[1],singularValues[2],singularValues[3],singularValues[4],singularValues[5],singularValues[6],singularValues[7],singularValues[8],singularValues[9]);

    // find the smallest singular value:
    int smallestIndex = 1;
    for(int i=2;i<10;i++) if(singularValues[i]<singularValues[smallestIndex]) smallestIndex=i;

    // solution is the nullspace of the matrix, which is the column in V corresponding to the smallest singular value (which should be 0)
    printf("Singular vector that gives H: %f, %f, %f, %f, %f, %f, %f, %f, %f\n",nullspaceMatrix[1][smallestIndex],nullspaceMatrix[2][smallestIndex],nullspaceMatrix[3][smallestIndex],nullspaceMatrix[4][smallestIndex],nullspaceMatrix[5][smallestIndex],nullspaceMatrix[6][smallestIndex],nullspaceMatrix[7][smallestIndex],nullspaceMatrix[8][smallestIndex],nullspaceMatrix[9][smallestIndex]);
  
    return;	
}


double** dlt(double** N)
{
 
  int n = 4;
  int n2 = 2*n;

  // build the 2nx9 matrix of equations
  double** linEquations = dmatrix(1,n2,1,9);
  for(int j=1;j<=n2;j+=2){//calculate the two equations for each point correspondence
    int x=(j+1)/2;
    linEquations[j][1] = 0;
    linEquations[j][2] = 0;
    linEquations[j][3] = 0;
    linEquations[j][4] = -N[1][x];
    linEquations[j][5] = -N[2][x];
    linEquations[j][6] = -N[3][x];
    linEquations[j][7] = N[1][x]*N[5][x];
    linEquations[j][8] = N[2][x]*N[5][x];
    linEquations[j][9] = N[3][x]*N[5][x];

    linEquations[j+1][1] = N[1][x];
    linEquations[j+1][2] = N[2][x];
    linEquations[j+1][3] = N[3][x];
    linEquations[j+1][4] = 0;
    linEquations[j+1][5] = 0;
    linEquations[j+1][6] = 0;
    linEquations[j+1][7] = -N[1][x]*N[4][x];
    linEquations[j+1][8] = -N[2][x]*N[4][x];
    linEquations[j+1][9] = -N[3][x]*N[4][x];
  }


    // compute the SVD
    double singularValues[10]; // 1..9
    double** nullspaceMatrix = dmatrix(1,9,1,9);
    svdcmp(linEquations, n2, 9, singularValues, nullspaceMatrix);

    // get the result
    //printf("\n Singular values: %f, %f, %f, %f, %f, %f, %f, %f, %f\n",singularValues[1],singularValues[2],singularValues[3],singularValues[4],singularValues[5],singularValues[6],singularValues[7],singularValues[8],singularValues[9]);

    // find the smallest singular value:
    int smallestIndex = 1;
    for(int i=2;i<10;i++) if(singularValues[i]<singularValues[smallestIndex]) smallestIndex=i;

    // solution is the nullspace of the matrix, which is the column in V corresponding to the smallest singular value (which should be 0)
    //printf("Singular vector that gives H: %f, %f, %f, %f, %f, %f, %f, %f, %f\n",nullspaceMatrix[1][smallestIndex],nullspaceMatrix[2][smallestIndex],nullspaceMatrix[3][smallestIndex],nullspaceMatrix[4][smallestIndex],nullspaceMatrix[5][smallestIndex],nullspaceMatrix[6][smallestIndex],nullspaceMatrix[7][smallestIndex],nullspaceMatrix[8][smallestIndex],nullspaceMatrix[9][smallestIndex]);
    
    double** H = dmatrix(1,3,1,3);
    H[1][1]=nullspaceMatrix[1][smallestIndex];
    H[1][2]=nullspaceMatrix[2][smallestIndex];
    H[1][3]=nullspaceMatrix[3][smallestIndex];
    H[2][1]=nullspaceMatrix[4][smallestIndex];
    H[2][2]=nullspaceMatrix[5][smallestIndex];
    H[2][3]=nullspaceMatrix[6][smallestIndex];
    H[3][1]=nullspaceMatrix[7][smallestIndex];
    H[3][2]=nullspaceMatrix[8][smallestIndex];
    H[3][3]=nullspaceMatrix[9][smallestIndex];
    return H;	
}



////////////////////////////////////////////////////////////////////////
// Image processing functions
// YOU IMPLEMENT THE FUNCTIONS IN THIS SECTION
////////////////////////////////////////////////////////////////////////

// Per-pixel Operations ////////////////////////////////////////////////

void R2Image::
Translate(double dx, double dy){

  R2Image Temp(*this);
  
  for( int i=0; i<width;i++){
    for(int j=0; j<height;j++){
      if(i-dx<width && i-dx>0 && j-dy<height && j-dy>0)
	{
	  Pixel(i,j)=Temp.Pixel(i-dx,j-dy);
	}
      else{
	Pixel(i,j)=R2black_pixel;
      }		
    }
  }
}


void R2Image::
Brighten(double factor)
{
  // Brighten the image by multiplying each pixel component by the factor.
  // This is implemented for you as an example of how to access and set pixels
  for (int i = 0; i < width; i++) {
    for (int j = 0;  j < height; j++) {
      Pixel(i,j) *= factor;
      Pixel(i,j).Clamp();
    }
  }
}

void R2Image::
SobelX(void)
{

  R2Image Temp(*this); 
  // R2Pixel n(0.5,0.5,0.5,1.0);//normalizer pixel
  //sobelX kernel
  double sobel_x [3][3];
  sobel_x[0][0] = 1; sobel_x[0][1] = 0; sobel_x[0][2] = -1; 
  sobel_x[1][0] = 2; sobel_x[1][1] = 0; sobel_x[1][2] = -2;  
  sobel_x[2][0] = 1; sobel_x[2][1] = 0; sobel_x[2][2] = -1;

  //copy the original image to the temp image
  for (int i = 0; i < width; i++) {
    for (int j = 0;  j < height; j++) {
      double r= Pixel(i,j).Red();
      double g= Pixel(i,j).Green();
      double b= Pixel(i,j).Blue();

      Temp.Pixel(i,j).SetRed(r);
      Temp.Pixel(i,j).SetGreen(g);
      Temp.Pixel(i,j).SetBlue(b);
    }
  }
 

  for (int i = 1; i < width-1; i++) {
    for (int j = 1;  j < height-1; j++) {
      Pixel(i,j).SetRed(0);
      Pixel(i,j).SetGreen(0);
      Pixel(i,j).SetBlue(0);
      for (int lx=-1; lx <= 1; lx++){
	for (int ly= -1; ly <= 1; ly++){
	  Pixel(i,j) += Temp.Pixel(i+lx,j+ly) * sobel_x[lx+1][ly+1];
	}
      }
      //Pixel(i,j) += n; 
      //Pixel(i,j).Clamp();	  
    }
  }
 
}

void R2Image::
SobelY(void)
{
  
  R2Image Temp(*this); 
  //R2Pixel n(0.5,0.5,0.5,1.0);//normalizer
  
  //copy the original image to the temporary image
  for (int i = 0; i < width; i++) {
    for (int j = 0;  j < height; j++) {
      double r= Pixel(i,j).Red();
      double g= Pixel(i,j).Green();
      double b= Pixel(i,j).Blue();

      Temp.Pixel(i,j).SetRed(r);
      Temp.Pixel(i,j).SetGreen(g);
      Temp.Pixel(i,j).SetBlue(b);
    }
  }
  //sobelY kernel
  double sobel_y [3][3];
  sobel_y[0][0] = 1; sobel_y[0][1] = 2; sobel_y[0][2] = 1; 
  sobel_y[1][0] = 0; sobel_y[1][1] = 0; sobel_y[1][2] = 0;  
  sobel_y[2][0] = -1; sobel_y[2][1] = -2; sobel_y[2][2] = -1;

  for (int i = 1; i < width-1; i++) {
    for (int j = 1;  j < height-1; j++) {
      Pixel(i,j).SetRed(0);
      Pixel(i,j).SetGreen(0);
      Pixel(i,j).SetBlue(0);
      for (int lx= -1; lx <= 1; lx++){
        for (int ly= -1; ly <= 1; ly++){
          Pixel(i,j) += Temp.Pixel(i+lx,j+ly) * sobel_y[lx+1][ly+1];
	}
      }
      //Pixel(i,j) += n;
      // Pixel(i,j).Clamp();
    }
  }
}

//median helper funciton
bool pCompare(R2Pixel n, R2Pixel m)
{
  return (n.Red()+n.Green()+n.Blue())<(m.Red()+m.Green()+m.Blue());
}

void R2Image::
MedianFilter(void)
{
  R2Image Temp(*this);
  R2Pixel k[25];
  
  for(int i=0;i<width;i++){
    for(int j=0;j<height;j++){
      for(int lx=-2;lx<=2;lx++){
	for(int ly=-2;ly<=2;ly++){
	  int index = 5*(lx+2)+(ly+2);
	  k[index]=Pixel(i+lx,j+ly);
	}
      }
      sort(k,k+25,pCompare);
      Temp.Pixel(i,j)=k[12];
    }
  }

  for(int i=0;i<width;i++){
    for(int j=0;j<height;j++){
      Pixel(i,j)=Temp.Pixel(i,j);
    }
  }
}

double pDiff(R2Pixel n, R2Pixel m){
  R2Pixel dif = n-m;
  double vd = sqrt(dif.Red()*dif.Red()+dif.Blue()*dif.Blue()+dif.Green()*dif.Green());
  return vd;
}

double weight(double n){
  int sigma = 2;
  double m = 2*sigma*sigma;
  double denom = sqrt(m*M_PI);
  double nom = exp((-n*n)/m);
  return nom/denom;
}

void R2Image::
BilateralFilter(void)
{
 
  R2Image Temp(*this);
  int r = 6;
  float dt; //distance between pixels
  float cd; //color difference between pixels
  float cw; //color weight
  float w; //weight assigned to that pixel
  double sum=0.0; //store sum of weights for normalization

  for (int i = r; i < width-r; i++) {
    for (int j = r;  j < height-r; j++) {
      Temp.Pixel(i,j).SetRed(0);
      Temp.Pixel(i,j).SetGreen(0);
      Temp.Pixel(i,j).SetBlue(0);
      sum = 0.0;
      for(int ly=-r;ly<=r;ly++){
	for (int lx=-r;lx<=r;lx++){
	  
	  dt = sqrt(lx*lx+ly*ly);
	  cd = pDiff(Pixel(i,j),Pixel(i+lx,j+ly));
	  float cd_n = 1-cd;
	  if (cd_n<0.7){
	    cw = 0;
	  }else if(cd_n>0.7){
	    cw = 1;
	  }else{
	    cw = cd_n;
	  }
	  
	  w = weight(dt)*cw;
	  Temp.Pixel(i,j) += Pixel(i+lx,j+ly)*w;
	  sum += w;
	    
	}
      }
      Temp.Pixel(i,j) /= sum;
      Temp.Pixel(i,j).Clamp();
     
    }
  }
  

  for(int i=0; i<width; i++){
    for(int j =0; j<height;j++){
      Pixel(i,j)=Temp.Pixel(i,j);

    }
  }


}

 
void R2Image::
LoG(void)
{
  // Apply the LoG operator to the image
  
  // FILL IN IMPLEMENTATION HERE (REMOVE PRINT STATEMENT WHEN DONE)
  fprintf(stderr, "LoG() not implemented\n");
}



void R2Image::
ChangeSaturation(double factor)
{

  for (int i = 0; i < width; i++) {
    for (int j = 0;  j < height; j++) {
      double r = Pixel(i,j).Red();
      double b = Pixel(i,j).Blue();
      double g = Pixel(i,j).Green();
      double l = 0.3*r + 0.6*g + 0.1*b; //grayscale
      Pixel(i,j).SetRed(l+factor*(r-l)); 
      Pixel(i,j).SetBlue(l+factor*(b-l));
      Pixel(i,j).SetGreen(l+factor*(g-l));
      Pixel(i,j).Clamp();
    }
  }

}


// Linear filtering ////////////////////////////////////////////////
void R2Image::
Blur(double sigma)
{
  R2Image Temp(*this);
  int kernel_w = 6*sigma+1;
  float kernel [kernel_w][1];
  
  double sum=0.0;//for later precomputation of kernel values
  int r = 3*sigma;
  float m = 2*sigma*sigma;
  float denom = sqrt(m*M_PI);
  
  for (int i=0; i<kernel_w; i++){
    float n = i-r;
    float nom = exp((-n*n)/m);
    kernel[i][0] = nom/denom;
    sum += kernel[i][0];
  }
  
  for (int i=0; i<kernel_w; i++){
    kernel[i][0]/=sum; //normalize the kernel
  }

  for (int i = 0; i < width; i++) {
    for (int j = 0;  j < height; j++) {
      Temp.Pixel(i,j).SetRed(0);
      Temp.Pixel(i,j).SetGreen(0);
      Temp.Pixel(i,j).SetBlue(0);
      
      for(int ly=-r;ly<=r;ly++){
	if(j+ly<0){ //edge processing
	  Temp.Pixel(i,j) += Pixel(i,-j-ly)*kernel[ly+r][0];
	}else{
	  if(j+ly>=height){ //edge processing
	    Temp.Pixel(i,j) += Pixel(i,2*height-1-j-ly)*kernel[ly+r][0];
	  }else{
	    Temp.Pixel(i,j) += Pixel(i,j+ly)*kernel[ly+r][0];
	  }
	}
      }
    }
  }

  for (int i = 0; i < width; i++) {
    for (int j = 0;  j < height; j++) {
      Pixel(i,j).SetRed(0);
      Pixel(i,j).SetGreen(0);
      Pixel(i,j).SetBlue(0);
      
      for(int lx=-r;lx<=r;lx++){
	if(i+lx<0){ //edge
	  Pixel(i,j) += Temp.Pixel(-i-lx,j)*kernel[lx+r][0];
	}else{
	  if(i+lx>=width){//edge
	    Pixel(i,j) += Temp.Pixel(2*width-1-i-lx,j)*kernel[lx+r][0];
	  }
	  else{
	    Pixel(i,j) += Temp.Pixel(i+lx,j)*kernel[lx+r][0];
	  }
	}
      }
      
      //Pixel(i,j).Clamp();
    }
  }
  
 }


void R2Image::
Harris(double sigma)
{
  R2Image Temp1(*this);
  R2Image Temp2(*this);
  R2Image Temp3(*this);
  R2Pixel n(0.5,0.5,0.5,1);
  Temp1.SobelX();
  Temp2.SobelY();

  for(int i =0;i<width;i++){
    for(int j=0;j<height;j++){
      Temp3.Pixel(i,j)=Temp1.Pixel(i,j)*Temp2.Pixel(i,j);
      Temp1.Pixel(i,j)=Temp1.Pixel(i,j)*Temp1.Pixel(i,j);
      Temp2.Pixel(i,j)=Temp2.Pixel(i,j)*Temp2.Pixel(i,j);
    }
  }

  Temp1.Blur(sigma);
  Temp2.Blur(sigma);
  Temp3.Blur(sigma);

  for(int i =0;i<width;i++){
    for(int j=0;j<height;j++){
      Pixel(i,j)=Temp1.Pixel(i,j)*Temp2.Pixel(i,j)-Temp3.Pixel(i,j)*Temp3.Pixel(i,j)-0.04*((Temp1.Pixel(i,j)+Temp2.Pixel(i,j))*(Temp1.Pixel(i,j)+Temp2.Pixel(i,j)));
      Pixel(i,j)+=n;
      //Pixel(i,j).Clamp();
      
    }
  }
}

//helper functions
double pValue(R2Pixel n){

  return n.Red()+n.Blue()+n.Green();
  
}

void R2Image::
featureDetect(){
  
  R2Image Temp(*this);
  R2Pixel red(255,0,0,1);
  R2Pixel black(0,0,0,1);
  
  Temp.Harris(2);
  double max;
  int max_x;
  int max_y;
  
  for(int n=1;n<=150;n++){
    max = pValue(black);
    
    for(int i=10;i<width-10;i++){
      for(int j=10;j<height-10;j++){
	double v = pValue(Temp.Pixel(i,j));
	if (v > max){
	  max = v;
	  max_x = i;
	  max_y = j;
	}
      }
    }
    
    for(int dx = -20;dx<=20;dx++){
      for(int dy=-20;dy<=20;dy++){
	Temp.Pixel(max_x+dx,max_y+dy)=black;
      }
    }
    
    for(int k=-5;k<=5;k++){
      for (int l=-5; l<=5; l++){
	Pixel(max_x+k,max_y+l)=red;
     
      }
    }
  }

}


void R2Image::
Sharpen()
{

R2Image Temp(*this); 

for (int i = 0; i < width; i++) {
    for (int j = 0;  j < height; j++) {
	double r= Pixel(i,j).Red();
        double g= Pixel(i,j).Green();
        double b= Pixel(i,j).Blue();

	Temp.Pixel(i,j).SetRed(r);
        Temp.Pixel(i,j).SetGreen(g);
        Temp.Pixel(i,j).SetBlue(b);
	}
}

double sharp[3][3];
sharp[0][0] = 0; sharp[0][1] = -1; sharp[0][2] = 0; 
sharp[1][0] = -1; sharp[1][1] = 4; sharp[1][2] = -1;  
sharp[2][0] = 0; sharp[2][1] = -1; sharp[2][2] = 0;

for (int i=1; i<width-1; i++){
   for (int j=1; j<height-1; j++){
	for (int lx=-1; lx<=1; lx++){
           for (int ly=-1; ly<=1; ly++){
	       Pixel(i,j) += Temp.Pixel(i+lx,j+ly) * sharp[lx+1][-ly+1];
           }
        }
    Pixel(i,j).Clamp();
   }
}
 
}

void R2Image::
line(int x0, int x1, int y0, int y1, float r, float g, float b)
{
  if(x0>x1)
    {
      int x=y1;
      y1=y0;
      y0=x;

      x=x1;
      x1=x0;
      x0=x;
      }
  int deltax = x1 - x0;
  int deltay = y1 - y0;
  float error = 0;
  float deltaerr = 0.0;
  if(deltax!=0) deltaerr =fabs(float(float(deltay) / deltax));    // Assume deltax != 0 (line is not vertical),
  // note that this division needs to be done in a way that preserves the fractional part
  int y = y0;
  for(int x=x0;x<=x1;x++)
    {
      Pixel(x,y).Reset(r,g,b,1.0);
      error = error + deltaerr;
      if(error>=0.5)
	{
	  if(deltay>0) y = y + 1;
	  else y = y - 1;

	  error = error - 1.0;
	}
    }
  /* if(x0>3 && x0<width-3 && y0>3 && y0<height-3)
    {
      for(int x=x0-3;x<=x0+3;x++)
	{
	  for(int y=y0-3;y<=y0+3;y++)
	    {
	      Pixel(x,y).Reset(r,g,b,1.0);
	    }
	}
	}*/
}

//calculate the SSD between two pixels
float SSD(R2Pixel n, R2Pixel m){
  R2Pixel r=n-m;
  float result = r.Red()*r.Red()+r.Green()*r.Green()+r.Blue()*r.Blue();
  return result;

}

void R2Image::
blendOtherImageTranslated(R2Image * otherImage)
{
 
  R2Image Temp1(*this);//for feature search
  R2Image Temp2(*this);//keep original image
  R2Pixel black(0,0,0,1);
  R2Pixel red(255,0,0,1);
  
  Temp1.Harris(1.5);//without the clamp
  
  int feature[4][150];//store the coordinates of the features on image1 and the matching features in the other image
  memset(feature,0,sizeof(int)*4*150);//set all array values to zero
  float max;
  int max_x,max_y;
  
  int w = width*0.1;//search window half width
  int h = height*0.1;//search window half height
  int bx,by;
  float ssd,mssd;
  int mx,my;
  int wstart, wend, hstart, hend;

  //feature search on image 1 and store the coordinate in the feature array
  for(int n=0;n<150;n++){
    max = pValue(black);
    
    for(int i=10;i<width-10;i++){
      for(int j=10;j<height-10;j++){
	 
	if (pValue(Temp1.Pixel(i,j)) >= max){
	  max = pValue(Temp1.Pixel(i,j));
	  max_x = i;
	  max_y = j;
	}
      }
    }
    
    feature[0][n]=max_x;
    feature[1][n]=max_y;
    
    for(int dx=-15; dx<=15;dx++){
      for(int dy=-15; dy<=15;dy++){
	Temp1.Pixel(feature[0][n]+dx,feature[1][n]+dy)=black;//black out the surrounding area
	
      }
    }
    /*
    for(int k1=-5;k1<=5;k1++){
      for(int k2=-5;k2<=5;k2++){
	Pixel(feature[0][n]+k1,feature[1][n]+k2)=red;//draw filled boxes around the found features on image 1
      }
    }
    */
  }
  
  //do local search for each feature on the other image
  for(int l=0;l<150;l++){
  
    bx=feature[0][l];
    by=feature[1][l];
  
    wstart = bx-w;//x coordinate of the search window's left edge
    if(wstart<3){
      wstart=3;
    }
  
    wend = bx+w;//x coordinate of the search window's right edge
    if(wend>width-4){
      wend=width-4;
    }
  
    hstart =by-h;//y coordinate of the search window's bottom edge
    if(hstart<3){
      hstart=3;
    }
  
    hend = by+h;//y coordinate of the search window's top edge
    if(hend>height-4){
      hend=height-4;
    }
    
    mssd=10000000000000; //initiate the minimum ssd to a large value
  
    for(int sw=wstart;sw<wend;sw++){
      for(int sh=hstart;sh<hend;sh++){
	ssd=0;//reset the ssd value for each pair of pixels to be compared
      
	for(int m=-3;m<=3;m++){//for each pair of pixels, compare the ssd of a 7*7 window
	  for(int n=-3;n<=3;n++){
	    float dif = SSD(Temp2.Pixel(bx+m,by+n),otherImage->Pixel(sw+m,sh+n));
	    ssd += dif;
	  }
	}
      
	if(ssd<mssd)
	  {//set the minimum ssd to the current smallest ssd and store the x,y coordinates with mx, my
	    mssd=ssd;
	    mx=sw;
	    my=sh;
	  }
      
      }
    }
  
    feature[2][l]=mx;//store the best location within the search window
    feature[3][l]=my;
  
  }
 
  for(int p=0;p<150;p++){
    
   for(int k1=-5;k1<=5;k1++){
      for(int k2=-5;k2<=5;k2++){
	Pixel(feature[0][p]+k1,feature[1][p]+k2)=red;//draw filled boxes around the found features on image 1
      }
   }
    
    this->line(feature[0][p],feature[2][p],feature[1][p],feature[3][p],255,0,0);//draw the matching lines
   
  }
 
}


void R2Image::
Ransac(R2Image * otherImage)
{
 
  R2Image Temp1(*this);//for feature search
  R2Image Temp2(*this);//keep original image
  R2Pixel black(0,0,0,1);
  R2Pixel red(255,0,0,1);
  R2Pixel green(0,255,0,1);
  
  Temp1.Harris(2);
  
  int feature[5][150];//store the coordinates of the features on image1 and the matching features in the other image
  memset(feature,0,sizeof(int)*5*150);//set all array values to zero
  float max;
  int max_x,max_y;
  
  int w = width*0.1;//search window half width
  int h = height*0.1;//search window half height
  int bx,by;
  float ssd,mssd;
  int mx,my;
  int wstart, wend, hstart, hend;

  //feature search on image 1 and store the coordinate in the feature array
  for(int n=0;n<150;n++){
    max = pValue(Temp1.Pixel(10,10));
    
    for(int i=10;i<width-10;i++){
      for(int j=10;j<height-10;j++){
	 
	if (pValue(Temp1.Pixel(i,j)) >= max){
	  max = pValue(Temp1.Pixel(i,j));
	  max_x = i;
	  max_y = j;
	}
      }
    }
    
    for(int dx = -15;dx<=15;dx++){
      for(int dy=-15;dy<=15;dy++){
	Temp1.Pixel(max_x+dx,max_y+dy)=black;//black out the surrounding area
	
      }
    }
    feature[0][n]=max_x;
    feature[1][n]=max_y;
    
    for(int k1=-5;k1<=5;k1++){
      for(int k2=-5;k2<=5;k2++){
	Pixel(max_x+k1,max_y+k2)=green;//draw filled boxes around the found features on image 1
      }
    }
  }
  
  //do local search for each feature on the other image
  for(int l=0;l<150;l++){
  
    bx=feature[0][l];
    by=feature[1][l];
  
    wstart = bx-w;//x coordinate of the search window's left edge
    if(wstart<3){
      wstart=3;
    }
  
    wend = bx+w;//x coordinate of the search window's right edge
    if(wend>width-4){
      wend=width-4;
    }
  
    hstart =by-h;//y coordinate of the search window's bottom edge
    if(hstart<3){
      hstart=3;
    }
  
    hend = by+h;//y coordinate of the search window's top edge
    if(hend>height-4){
      hend=height-4;
    }
    
    mssd=10000000000000; //initiate the minimum ssd to a large value
  
    for(int sw=wstart;sw<wend;sw++){
      for(int sh=hstart;sh<hend;sh++){
	ssd=0;//reset the ssd value for each pair of pixels to be compared
      
	for(int m=-3;m<=3;m++){//for each pair of pixels, compare the ssd of a 7*7 window
	  for(int n=-3;n<=3;n++){
	    float dif = SSD(Temp2.Pixel(bx+m,by+n),otherImage->Pixel(sw+m,sh+n));
	    ssd += dif;
	  }
	}
      
	if(ssd<mssd)
	  {//set the minimum ssd to the current smallest ssd and store the x,y coordinates with mx, my
	    mssd=ssd;
	    mx=sw;
	    my=sh;
	  }
      
      }
    }
  
    feature[2][l]=mx;//store the best location within the search window
    feature[3][l]=my;
  
  }

  int r = 0;
  int x1,y1,x2,y2,target;
  int inCount;//cross product and count of inliers
  int thr=0;
  float sim,base;

 
  while(r<150){
    
    int v = rand() % 150; //randomly choose a motion vector
    x1 = feature[2][v]-feature[0][v];
    y1 = feature[3][v]-feature[1][v];
    inCount = 0;
    
    for(int z=0;z<150;z++){//loop through all the motion vectors to count the inliers
      x2 = feature[2][z]-feature[0][z];
      y2 = feature[3][z]-feature[1][z];
      sim =sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2));//calculate the distance between two motion vectors
    
      if(sim<=6){
	inCount++;
      }
    }

    if(inCount>thr){
      target = v;
      thr = inCount;
    }
    r++;
  }
  
  
  //result[0] = feature[2][target]-feature[0][target]; //store dx,dy to result;
  //result[1] = feature[3][target]-feature[1][target];

  
  x1 = feature[2][target]-feature[0][target]; //store dx,dy to result;
  y1 = feature[3][target]-feature[1][target];


  
  //cout<<"target:"<<target;
  
  for(int z=0;z<150;z++){//redraw the translation line
    x2 = feature[2][z]-feature[0][z];
    y2 = feature[3][z]-feature[1][z];
    sim =sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2));
    cout<<sim<<' ';
    if(sim<=6){
      this->line(feature[0][z],feature[2][z],feature[1][z],feature[3][z],0,255,0);//draw the matching lines
    }
    else{
      for(int h1=-5;h1<=5;h1++){
	for(int h2=-5;h2<=5;h2++){
	  Pixel(feature[0][z]+h1,feature[1][z]+h2)=red;//draw filled boxes around these outliers
	}
      }
      this->line(feature[0][z],feature[2][z],feature[1][z],feature[3][z],255,0,0);//draw the matching lines
    }
  }
  
}



void R2Image::
RansacT(R2Image * otherImage, double result[2])
{
 
  R2Image Temp1(*this);//for feature search
  R2Image Temp2(*this);//keep original image
  R2Pixel black(0,0,0,1);
  R2Pixel red(255,0,0,1);
  R2Pixel green(0,255,0,1);
  
  Temp1.Harris(2);
  
  int feature[5][150];//store the coordinates of the features on image1 and the matching features in the other image
  memset(feature,0,sizeof(int)*5*150);//set all array values to zero
  float max;
  int max_x,max_y;
  
  int w = width*0.1;//search window half width
  int h = height*0.1;//search window half height
  int bx,by;
  float ssd,mssd;
  int mx,my;
  int wstart, wend, hstart, hend;

  //feature search on image 1 and store the coordinate in the feature array
  for(int n=0;n<150;n++){
    max = pValue(Temp1.Pixel(10,10));
    
    for(int i=10;i<width-10;i++){
      for(int j=10;j<height-10;j++){
	 
	if (pValue(Temp1.Pixel(i,j)) >= max){
	  max = pValue(Temp1.Pixel(i,j));
	  max_x = i;
	  max_y = j;
	}
      }
    }
    
    for(int dx = -15;dx<=15;dx++){
      for(int dy=-15;dy<=15;dy++){
	Temp1.Pixel(max_x+dx,max_y+dy)=black;//black out the surrounding area
	
      }
    }
    feature[0][n]=max_x;
    feature[1][n]=max_y;
    
    for(int k1=-5;k1<=5;k1++){
      for(int k2=-5;k2<=5;k2++){
	Pixel(max_x+k1,max_y+k2)=green;//draw filled boxes around the found features on image 1
      }
    }
  }
  
  //do local search for each feature on the other image
  for(int l=0;l<150;l++){
  
    bx=feature[0][l];
    by=feature[1][l];
  
    wstart = bx-w;//x coordinate of the search window's left edge
    if(wstart<3){
      wstart=3;
    }
  
    wend = bx+w;//x coordinate of the search window's right edge
    if(wend>width-4){
      wend=width-4;
    }
  
    hstart =by-h;//y coordinate of the search window's bottom edge
    if(hstart<3){
      hstart=3;
    }
  
    hend = by+h;//y coordinate of the search window's top edge
    if(hend>height-4){
      hend=height-4;
    }
    
    mssd=10000000000000; //initiate the minimum ssd to a large value
  
    for(int sw=wstart;sw<wend;sw++){
      for(int sh=hstart;sh<hend;sh++){
	ssd=0;//reset the ssd value for each pair of pixels to be compared
      
	for(int m=-3;m<=3;m++){//for each pair of pixels, compare the ssd of a 7*7 window
	  for(int n=-3;n<=3;n++){
	    float dif = SSD(Temp2.Pixel(bx+m,by+n),otherImage->Pixel(sw+m,sh+n));
	    ssd += dif;
	  }
	}
      
	if(ssd<mssd)
	  {//set the minimum ssd to the current smallest ssd and store the x,y coordinates with mx, my
	    mssd=ssd;
	    mx=sw;
	    my=sh;
	  }
      
      }
    }
  
    feature[2][l]=mx;//store the best location within the search window
    feature[3][l]=my;
  
  }

  int r = 0;
  int x1,y1,x2,y2,target;
  int inCount;//cross product and count of inliers
  int thr=0;
  float sim,base;

 
  while(r<150){
    
    int v = rand() % 150; //randomly choose a motion vector
    x1 = feature[2][v]-feature[0][v];
    y1 = feature[3][v]-feature[1][v];
    inCount = 0;
    
    for(int z=0;z<150;z++){//loop through all the motion vectors to count the inliers
      x2 = feature[2][z]-feature[0][z];
      y2 = feature[3][z]-feature[1][z];
      sim =sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2));//calculate the distance between two motion vectors
    
      if(sim<=6){
	inCount++;
      }
    }

    if(inCount>thr){
      target = v;
      thr = inCount;
    }
    r++;
  }
  
  
  result[0] = feature[2][target]-feature[0][target]; //store dx,dy to result;
  result[1] = feature[3][target]-feature[1][target];

  /*
  x1 = feature[2][target]-feature[0][target]; //store dx,dy to result;
  y1 = feature[3][target]-feature[1][target];
  */
  
  //cout<<"target:"<<target;
  /*
  for(int z=0;z<150;z++){//redraw the translation line
    x2 = feature[2][z]-feature[0][z];
    y2 = feature[3][z]-feature[1][z];
    sim =sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2));
    cout<<sim<<' ';
    if(sim<=6){
      this->line(feature[0][z],feature[2][z],feature[1][z],feature[3][z],0,255,0);//draw the matching lines
    }
    else{
      for(int h1=-5;h1<=5;h1++){
	for(int h2=-5;h2<=5;h2++){
	  Pixel(feature[0][z]+h1,feature[1][z]+h2)=red;//draw filled boxes around these outliers
	}
      }
      this->line(feature[0][z],feature[2][z],feature[1][z],feature[3][z],255,0,0);//draw the matching lines
    }
  }
  */
  
}


void R2Image::
RansacP(R2Image * otherImage)
{
 
  R2Image Temp1(*this);//for feature search
  R2Image Temp2(*this);//keep original image
  R2Pixel black(0,0,0,1);
  R2Pixel red(255,0,0,1);
  R2Pixel green(0,255,0,1);
  
  Temp1.Harris(1.5);
  
  int feature[4][150];//store the coordinates of the features on image1 and the matching features in the other image
  memset(feature,0,sizeof(int)*4*150);//set all array values to zero
  float max;
  int max_x,max_y;
  
  int w = width*0.1;//search window half width
  int h = height*0.1;//search window half height
  int bx,by;
  float ssd,mssd;
  int mx,my;
  int wstart, wend, hstart, hend;

  //feature search on image 1 and store the coordinate in the feature array
  for(int n=0;n<150;n++){
    max = pValue(black);
    
    for(int i=10;i<width-10;i++){
      for(int j=10;j<height-10;j++){
	 
	if (pValue(Temp1.Pixel(i,j)) > max){
	  max = pValue(Temp1.Pixel(i,j));
	  max_x = i;
	  max_y = j;
	}
      }
    }
    
    for(int dx = -20;dx<=20;dx++){
      for(int dy=-20;dy<=20;dy++){
	Temp1.Pixel(max_x+dx,max_y+dy)=black;//black out the surrounding area
	
      }
    }
    feature[0][n]=max_x;
    feature[1][n]=max_y;
    
    for(int k1=-5;k1<=5;k1++){
      for(int k2=-5;k2<=5;k2++){
	Pixel(max_x+k1,max_y+k2)=green;//draw filled boxes around the found features on image 1
      }
    }
  }
  
  //do local search for each feature on the other image
  for(int l=0;l<150;l++){
  
    bx=feature[0][l];
    by=feature[1][l];
  
    wstart = bx-w;//x coordinate of the search window's left edge
    if(wstart<3){
      wstart=3;
    }
  
    wend = bx+w;//x coordinate of the search window's right edge
    if(wend>width-4){
      wend=width-4;
    }
  
    hstart =by-h;//y coordinate of the search window's bottom edge
    if(hstart<3){
      hstart=3;
    }
  
    hend = by+h;//y coordinate of the search window's top edge
    if(hend>height-4){
      hend=height-4;
    }
    
    mssd=10000000000000; //initiate the minimum ssd to a large value
  
    for(int sw=wstart;sw<wend;sw++){
      for(int sh=hstart;sh<hend;sh++){
	ssd=0;//reset the ssd value for each pair of pixels to be compared
      
	for(int m=-3;m<=3;m++){//for each pair of pixels, compare the ssd of a 7*7 window
	  for(int n=-3;n<=3;n++){
	    float dif = SSD(Temp2.Pixel(bx+m,by+n),otherImage->Pixel(sw+m,sh+n));
	    ssd += dif;
	  }
	}
      
	if(ssd<mssd)
	  {//set the minimum ssd to the current smallest ssd and store the x,y coordinates with mx, my
	    mssd=ssd;
	    mx=sw;
	    my=sh;
	  }
      
      }
    }
  
    feature[2][l]=mx;//store the best location within the search window
    feature[3][l]=my;
  
  }

  
  
  int r = 0;
  double x1,y1,w1,x2,y2,w2,x2_,y2_,w2_;
  int inCount;//count of inliers
  int thr=0;//threshold of inliers
  double sim;//similarity
  int v [1][4];//matrix to store the randomly selected four point correspondences
  double** coord=dmatrix(1,6,1,4);//coordinate matrix for the four point correspondences
  double** target=dmatrix(1,3,1,3);//to store the optimal matrix
 
  while(r<250){
    
    v[0][0] = rand() % 150; //randomly choose 4 distinct motion vectors
    v[0][1] = rand() % 150;
    if(v[0][1]==v[0][0]){
      while(v[0][1]==v[0][0]){//make sure v2 is not the same as v1
	v[0][1] = rand() % 150;
      }
    }
    v[0][2] = rand() % 150;
    if(v[0][2]==v[0][1] |v[0][2]==v[0][0]){
      while(v[0][2]==v[0][1]|v[0][2]==v[0][0]){//make sure v3 is different
	v[0][2] = rand()%150;
      }
    }
    v[0][3] = rand() % 150;
    if(v[0][3]==v[0][2]|v[0][3]==v[0][1]|v[0][3]==v[0][0]){//make sure v4 is different
      while(v[0][3]==v[0][2]|v[0][3]==v[0][1]|v[0][3]==v[0][0]){
	v[0][3] = rand()%150;
      }
      
    }
    
    for(int y=1;y<=4;y++){//copy the randomly selected point correspondences to the the coordinate matrix
      for(int w=1;w<=2;w++){
	int tv = v[0][y-1];
	coord[w][y]=feature[w-1][tv];
	coord[3][y]=1;
	coord[6][y]=1;
      }
    }

    for(int y=1;y<=4;y++){
      for(int w=4;w<=5;w++){
	int tv = v[0][y-1];
	coord[w][y]=feature[w-2][tv];
      }
    }
   

    //pass down the four point correspondences to DLT and solve for H
    double** H= dmatrix(1,3,1,3);
    H = dlt(coord);   

    inCount = 0;
    for(int z=0;z<150;z++){//loop through all the correspondences to count the inliers
      
      x1 = feature[0][z];//v1's coordinates 
      y1 = feature[1][z];
      w1 = 1;
      
      x2 = feature[2][z];//v2's coordinates
      y2 = feature[3][z];
      w2 = 1;

      x2_=x1*H[1][1]+y1*H[1][2]+H[1][3];//Hv1's coordinates
      y2_=x1*H[2][1]+y1*H[2][2]+H[2][3];
      w2_=x1*H[3][1]+y1*H[3][2]+H[3][3];
      //factor to same scale
      x2_ /= w2_;
      y2_ /= w2_;
      //cout<<"x2,x2_:"<<x2<<","<<x2_<<' ';
      
      sim =sqrt((x2_-x2)*(x2_-x2)+(y2_-y2)*(y2_-y2));//calculate the distance between Hv1 and v2;
      //cout<<sim<<"&";
      
      if(sim<=6){
	inCount++;
      }
    }
    
    if(inCount>thr){//if the counted inliers is greater than the current threshold, store the current matrix to the target matrix
      for(int a=1; a<=3; a++){
	for(int b=1; b<=3; b++){
	  target[a][b]=H[a][b];
	}
      }
      thr = inCount;//set the threshold to the new standard
    }
    r++;
  }
  
  
  for(int z=0;z<150;z++){//redraw the translation line
    x1 = feature[0][z];//v1's coordinates 
    y1 = feature[1][z];
    w1 = 1;
      
    x2 = feature[2][z];//v2's coordinates
    y2 = feature[3][z];
    w2 = 1;
    
    x2_=x1*target[1][1]+y1*target[1][2]+target[1][3];//Hv1's coordinates
    y2_=x1*target[2][1]+y1*target[2][2]+target[2][3];
    w2_=x1*target[3][1]+y1*target[3][2]+target[3][3];
    x2_ /= w2_;
    y2_ /= w2_;
      
    sim =sqrt((x2_-x2)*(x2_-x2)+(y2_-y2)*(y2_-y2));//calculate the distance between Hv1 and v2;
    //cout<<sim<<' ';
    if(sim<=6){
      this->line(feature[0][z],feature[2][z],feature[1][z],feature[3][z],0,255,0);//draw the matching lines
    }
    else{
      for(int h1=-5;h1<=5;h1++){
	for(int h2=-5;h2<=5;h2++){
	  Pixel(feature[0][z]+h1,feature[1][z]+h2)=red;//draw filled boxes around these outliers
	}
      }
      this->line(feature[0][z],feature[2][z],feature[1][z],feature[3][z],255,0,0);//draw the matching lines
    }
  }
  
}



void R2Image::
blendOtherImageHomography(R2Image * otherImage)
{
 
  R2Image Temp1(*this);//for feature search
  R2Image Temp2(*this);//keep original image
  R2Pixel black(0,0,0,1);
  R2Pixel red(255,0,0,1);
  R2Pixel green(0,255,0,1);
  
  Temp1.Harris(1.5);
  
  int feature[4][150];//store the coordinates of the features on image1 and the matching features in the other image
  memset(feature,0,sizeof(int)*4*150);//set all array values to zero
  float max;
  int max_x,max_y;
  
  int w = width*0.1;//search window half width
  int h = height*0.1;//search window half height
  int bx,by;
  float ssd,mssd;
  int mx,my;
  int wstart, wend, hstart, hend;

  //feature search on image 1 and store the coordinate in the feature array
  for(int n=0;n<150;n++){
    max = pValue(black);
    
    for(int i=10;i<width-10;i++){
      for(int j=10;j<height-10;j++){
	 
	if (pValue(Temp1.Pixel(i,j)) > max){
	  max = pValue(Temp1.Pixel(i,j));
	  max_x = i;
	  max_y = j;
	}
      }
    }
    
    for(int dx = -20;dx<=20;dx++){
      for(int dy=-20;dy<=20;dy++){
	Temp1.Pixel(max_x+dx,max_y+dy)=black;//black out the surrounding area
	
      }
    }
    feature[0][n]=max_x;
    feature[1][n]=max_y;
    
  }
  
  //do local search for each feature on the other image
  for(int l=0;l<150;l++){
  
    bx=feature[0][l];
    by=feature[1][l];
  
    wstart = bx-w;//x coordinate of the search window's left edge
    if(wstart<3){
      wstart=3;
    }
  
    wend = bx+w;//x coordinate of the search window's right edge
    if(wend>width-4){
      wend=width-4;
    }
  
    hstart =by-h;//y coordinate of the search window's bottom edge
    if(hstart<3){
      hstart=3;
    }
  
    hend = by+h;//y coordinate of the search window's top edge
    if(hend>height-4){
      hend=height-4;
    }
    
    mssd=10000000000000; //initiate the minimum ssd to a large value
  
    for(int sw=wstart;sw<wend;sw++){
      for(int sh=hstart;sh<hend;sh++){
	ssd=0;//reset the ssd value for each pair of pixels to be compared
      
	for(int m=-3;m<=3;m++){//for each pair of pixels, compare the ssd of a 7*7 window
	  for(int n=-3;n<=3;n++){
	    float dif = SSD(Temp2.Pixel(bx+m,by+n),otherImage->Pixel(sw+m,sh+n));
	    ssd += dif;
	  }
	}
      
	if(ssd<mssd)
	  {//set the minimum ssd to the current smallest ssd and store the x,y coordinates with mx, my
	    mssd=ssd;
	    mx=sw;
	    my=sh;
	  }
      
      }
    }
  
    feature[2][l]=mx;//store the best location in image 2 within the search window
    feature[3][l]=my;
  
  }
  
  int r = 0;//loop counter
  double x1,y1,w1,x2,y2,w2,x2_,y2_,w2_;
  int inCount;//count of inliers
  int thr=0;//threshold of inliers
  double sim;//similarity
  int v [1][4];//matrix to store the index of the randomly selected four point correspondences
  double** coord=dmatrix(1,6,1,4);//coordinate matrix for the four point correspondences
  double** target=dmatrix(1,3,1,3);//to store the optimal matrix
 
  while(r<250){
    
    v[0][0] = rand() % 150; //randomly choose 4 distinct motion vectors
    v[0][1] = rand() % 150;
    if(v[0][1]==v[0][0]){
      while(v[0][1]==v[0][0]){//make sure v2 is not the same as v1
	v[0][1] = rand() % 150;
      }
    }
    v[0][2] = rand() % 150;
    if(v[0][2]==v[0][1] |v[0][2]==v[0][0]){
      while(v[0][2]==v[0][1]|v[0][2]==v[0][0]){//make sure v3 is different
	v[0][2] = rand()%150;
      }
    }
    v[0][3] = rand() % 150;
    if(v[0][3]==v[0][2]|v[0][3]==v[0][1]|v[0][3]==v[0][0]){//make sure v4 is different
      while(v[0][3]==v[0][2]|v[0][3]==v[0][1]|v[0][3]==v[0][0]){
	v[0][3] = rand()%150;
      }
      
    }
    
    for(int y=1;y<=4;y++){//copy the randomly selected point correspondences to the the coordinate matrix
      for(int w=1;w<=2;w++){
	int tv = v[0][y-1];
	coord[w][y]=feature[w-1][tv];
	coord[3][y]=1;
	coord[6][y]=1;
      }
    }

    for(int y=1;y<=4;y++){
      for(int w=4;w<=5;w++){
	int tv = v[0][y-1];
	coord[w][y]=feature[w-2][tv];
      }
    }
   

    //pass down the four point correspondences to DLT and solve for H
    double** H= dmatrix(1,3,1,3);
    H = dlt(coord);   

    inCount = 0;
    for(int z=0;z<150;z++){//loop through all the correspondences to count the inliers
      
      x1 = feature[0][z];//v1's coordinates 
      y1 = feature[1][z];
      w1 = 1;
      
      x2 = feature[2][z];//v2's coordinates
      y2 = feature[3][z];
      w2 = 1;

      x2_=x1*H[1][1]+y1*H[1][2]+H[1][3];//Hv1's coordinates
      y2_=x1*H[2][1]+y1*H[2][2]+H[2][3];
      w2_=x1*H[3][1]+y1*H[3][2]+H[3][3];
      //factor to same scale
      x2_ /= w2_;
      y2_ /= w2_;
      
      
      sim =sqrt((x2_-x2)*(x2_-x2)+(y2_-y2)*(y2_-y2));//calculate the distance between Hv1 and v2;
      
      
      if(sim<=6){
	inCount++;
      }
    }
    
    
    if(inCount>thr){//if the counted inliers is greater than the current threshold, store the current matrix to the target matrix
      for(int a=1; a<=3; a++){
	for(int b=1; b<=3; b++){
	  target[a][b]=H[a][b];
	}
      }
      thr = inCount;//set the threshold to the new standard
    }
    r++;
  }
  
  //below we apply the same dlt algorithm to the robust features found through RANSAC above:
  //iteratively randomly selec 4 point correspondences, estimate the H matrix, count the inliers, choose the optimal matrix.
  
  int rfeature[4][thr];//array to store the robust feature locations found above
  int c=0;
  for(int z=0;z<150;z++){//copy the robust correspondence points to the rfeature array
    
    x1 = feature[0][z];//v1's coordinates 
    y1 = feature[1][z];
    w1 = 1;
      
    x2 = feature[2][z];//v2's coordinates
    y2 = feature[3][z];
    w2 = 1;
    
    x2_=x1*target[1][1]+y1*target[1][2]+target[1][3];//Hv1's coordinates
    y2_=x1*target[2][1]+y1*target[2][2]+target[2][3];
    w2_=x1*target[3][1]+y1*target[3][2]+target[3][3];
    x2_ /= w2_;
    y2_ /= w2_;
      
    sim =sqrt((x2_-x2)*(x2_-x2)+(y2_-y2)*(y2_-y2));//calculate the distance between Hv1 and v2;
   
    if(sim<=6 && c<thr){
      rfeature[0][c]=feature[0][z];
      rfeature[1][c]=feature[1][z];
      rfeature[2][c]=feature[2][z];
      rfeature[3][c]=feature[3][z];
      c++;
    }
  }
 

  r=0;//reinitiate
  int rthr=0;//the new threshold
  
  while(r<100){
    
    v[0][0] = rand() % thr; //randomly choose 4 distinct motion vectors
    v[0][1] = rand() % thr;
    if(v[0][1]==v[0][0]){
      while(v[0][1]==v[0][0]){//make sure v2 is not the same as v1
	v[0][1] = rand() % thr;
      }
    }
    v[0][2] = rand() % thr;
    if(v[0][2]==v[0][1] |v[0][2]==v[0][0]){
      while(v[0][2]==v[0][1]|v[0][2]==v[0][0]){//make sure v3 is different
	v[0][2] = rand()%thr;
      }
    }
    v[0][3] = rand() % thr;
    if(v[0][3]==v[0][2]|v[0][3]==v[0][1]|v[0][3]==v[0][0]){//make sure v4 is different
      while(v[0][3]==v[0][2]|v[0][3]==v[0][1]|v[0][3]==v[0][0]){
	v[0][3] = rand()% thr;
      }
    }
    
    for(int y=1;y<=4;y++){//copy the randomly selected point correspondences to the the coordinate matrix
      for(int w=1;w<=2;w++){
	int tv = v[0][y-1];
	coord[w][y]=rfeature[w-1][tv];
	coord[3][y]=1;
	coord[6][y]=1;
      }
    }

    for(int y=1;y<=4;y++){
      for(int w=4;w<=5;w++){
	int tv = v[0][y-1];
	coord[w][y]=rfeature[w-2][tv];
      }
    }
   
    //pass down the four point correspondences to DLT and solve for H
    double** H = dmatrix(1,3,1,3);
    H = dlt(coord);   

    inCount = 0;
    for(int z=0;z<thr;z++){//loop through all the robust correspondences to count the inliers
      
      x1 = rfeature[0][z];//v1's coordinates 
      y1 = rfeature[1][z];
      w1 = 1;
      
      x2 = rfeature[2][z];//v2's coordinates
      y2 = rfeature[3][z];
      w2 = 1;

      x2_=x1*H[1][1]+y1*H[1][2]+H[1][3];//Hv1's coordinates
      y2_=x1*H[2][1]+y1*H[2][2]+H[2][3];
      w2_=x1*H[3][1]+y1*H[3][2]+H[3][3];
      //factor to same scale
      x2_ /= w2_;
      y2_ /= w2_;
      sim =sqrt((x2_-x2)*(x2_-x2)+(y2_-y2)*(y2_-y2));//calculate the distance between Hv1 and v2;
     
      if(sim<=4){
	inCount++;
      }
    }
    
    if(inCount>rthr){//if the counted inliers is greater than the current threshold, store the current matrix to the target matrix
      for(int a=1; a<=3; a++){
	for(int b=1; b<=3; b++){
	  target[a][b]=H[a][b];
	}
      }
      rthr = inCount;//set the threshold to the new standard
    }
    //cout<<"robust threshold up to now: "<<rthr<<endl;
    //cout<<"r="<<r<<";"<<endl;
    
    if(r==99){
      break;
    }
    else if(r<99){
      r++;
    }
  }
  
  double u_,v_,w_; //warped coordinates
  int u1,v1;
  for(int u=0;u<width;u++){
    for(int v=0;v<height;v++){
    
      u_=u*target[1][1]+v*target[1][2]+target[1][3];//Hx1's coordinates
      v_=u*target[2][1]+v*target[2][2]+target[2][3];
      w_=u*target[3][1]+v*target[3][2]+target[3][3];
     
      u_/=w_;
      v_/=w_;
      if(u_>=0 && u_<width && v_>=0 && v_<height){
	u1 = int(u_+0.5);
	v1 = int(v_+0.5);
	R2Pixel newP = otherImage->Pixel(u1,v1); //blend 50-50
	newP *= 0.5;
	
	Pixel(u,v) *= 0.5;
	Pixel(u,v) += newP;
	//Pixel(u,v)=otherImage->Pixel(u1,v1); //just the warped image
      } 

    }
  }
}


////////////////////////////////////////////////////////////////////////
// Video-Processing functions 
////////////////////////////////////////////////////////////////////////




////////////////////////////////////////////////////////////////////////
// I/O Functions
////////////////////////////////////////////////////////////////////////

int R2Image::
Read(const char *filename)
{
  // Initialize everything
  if (pixels) { delete [] pixels; pixels = NULL; }
  npixels = width = height = 0;

  // Parse input filename extension
  char *input_extension;
  if (!(input_extension = (char*)strrchr(filename, '.'))) {
    fprintf(stderr, "Input file has no extension (e.g., .jpg).\n");
    return 0;
  }
  
  // Read file of appropriate type
  if (!strncmp(input_extension, ".bmp", 4)) return ReadBMP(filename);
  else if (!strncmp(input_extension, ".ppm", 4)) return ReadPPM(filename);
  else if (!strncmp(input_extension, ".jpg", 4)) return ReadJPEG(filename);
  else if (!strncmp(input_extension, ".jpeg", 5)) return ReadJPEG(filename);
  
  // Should never get here
  fprintf(stderr, "Unrecognized image file extension");
  return 0;
}



int R2Image::
Write(const char *filename) const
{
  // Parse input filename extension
  char *input_extension;
  if (!(input_extension = (char*)strrchr(filename, '.'))) {
    fprintf(stderr, "Input file has no extension (e.g., .jpg).\n");
    return 0;
  }
  
  // Write file of appropriate type
  if (!strncmp(input_extension, ".bmp", 4)) return WriteBMP(filename);
  else if (!strncmp(input_extension, ".ppm", 4)) return WritePPM(filename, 1);
  else if (!strncmp(input_extension, ".jpg", 5)) return WriteJPEG(filename);
  else if (!strncmp(input_extension, ".jpeg", 5)) return WriteJPEG(filename);

  // Should never get here
  fprintf(stderr, "Unrecognized image file extension");
  return 0;
}



////////////////////////////////////////////////////////////////////////
// BMP I/O
////////////////////////////////////////////////////////////////////////

#if (RN_OS == RN_LINUX) && !WIN32

typedef struct tagBITMAPFILEHEADER {
  unsigned short int bfType;
  unsigned int bfSize;
  unsigned short int bfReserved1;
  unsigned short int bfReserved2;
  unsigned int bfOffBits;
} BITMAPFILEHEADER;

typedef struct tagBITMAPINFOHEADER {
  unsigned int biSize;
  int biWidth;
  int biHeight;
  unsigned short int biPlanes;
  unsigned short int biBitCount;
  unsigned int biCompression;
  unsigned int biSizeImage;
  int biXPelsPerMeter;
  int biYPelsPerMeter;
  unsigned int biClrUsed;
  unsigned int biClrImportant;
} BITMAPINFOHEADER;

typedef struct tagRGBTRIPLE {
  unsigned char rgbtBlue;
  unsigned char rgbtGreen;
  unsigned char rgbtRed;
} RGBTRIPLE;

typedef struct tagRGBQUAD {
  unsigned char rgbBlue;
  unsigned char rgbGreen;
  unsigned char rgbRed;
  unsigned char rgbReserved;
} RGBQUAD;

#endif

#define BI_RGB        0L
#define BI_RLE8       1L
#define BI_RLE4       2L
#define BI_BITFIELDS  3L

#define BMP_BF_TYPE 0x4D42 /* word BM */
#define BMP_BF_OFF_BITS 54 /* 14 for file header + 40 for info header (not sizeof(), but packed size) */
#define BMP_BI_SIZE 40 /* packed size of info header */


static unsigned short int WordReadLE(FILE *fp)
{
  // Read a unsigned short int from a file in little endian format 
  unsigned short int lsb, msb;
  lsb = getc(fp);
  msb = getc(fp);
  return (msb << 8) | lsb;
}



static void WordWriteLE(unsigned short int x, FILE *fp)
{
  // Write a unsigned short int to a file in little endian format
  unsigned char lsb = (unsigned char) (x & 0x00FF); putc(lsb, fp); 
  unsigned char msb = (unsigned char) (x >> 8); putc(msb, fp);
}



static unsigned int DWordReadLE(FILE *fp)
{
  // Read a unsigned int word from a file in little endian format 
  unsigned int b1 = getc(fp);
  unsigned int b2 = getc(fp);
  unsigned int b3 = getc(fp);
  unsigned int b4 = getc(fp);
  return (b4 << 24) | (b3 << 16) | (b2 << 8) | b1;
}



static void DWordWriteLE(unsigned int x, FILE *fp)
{
  // Write a unsigned int to a file in little endian format 
  unsigned char b1 = (x & 0x000000FF); putc(b1, fp);
  unsigned char b2 = ((x >> 8) & 0x000000FF); putc(b2, fp);
  unsigned char b3 = ((x >> 16) & 0x000000FF); putc(b3, fp);
  unsigned char b4 = ((x >> 24) & 0x000000FF); putc(b4, fp);
}



static int LongReadLE(FILE *fp)
{
  // Read a int word from a file in little endian format 
  int b1 = getc(fp);
  int b2 = getc(fp);
  int b3 = getc(fp);
  int b4 = getc(fp);
  return (b4 << 24) | (b3 << 16) | (b2 << 8) | b1;
}



static void LongWriteLE(int x, FILE *fp)
{
  // Write a int to a file in little endian format 
  char b1 = (x & 0x000000FF); putc(b1, fp);
  char b2 = ((x >> 8) & 0x000000FF); putc(b2, fp);
  char b3 = ((x >> 16) & 0x000000FF); putc(b3, fp);
  char b4 = ((x >> 24) & 0x000000FF); putc(b4, fp);
}



int R2Image::
ReadBMP(const char *filename)
{
  // Open file
  FILE *fp = fopen(filename, "rb");
  if (!fp) {
    fprintf(stderr, "Unable to open image file: %s\n", filename);
    return 0;
  }

  /* Read file header */
  BITMAPFILEHEADER bmfh;
  bmfh.bfType = WordReadLE(fp);
  bmfh.bfSize = DWordReadLE(fp);
  bmfh.bfReserved1 = WordReadLE(fp);
  bmfh.bfReserved2 = WordReadLE(fp);
  bmfh.bfOffBits = DWordReadLE(fp);
  
  /* Check file header */
  assert(bmfh.bfType == BMP_BF_TYPE);
  /* ignore bmfh.bfSize */
  /* ignore bmfh.bfReserved1 */
  /* ignore bmfh.bfReserved2 */
  assert(bmfh.bfOffBits == BMP_BF_OFF_BITS);
  
  /* Read info header */
  BITMAPINFOHEADER bmih;
  bmih.biSize = DWordReadLE(fp);
  bmih.biWidth = LongReadLE(fp);
  bmih.biHeight = LongReadLE(fp);
  bmih.biPlanes = WordReadLE(fp);
  bmih.biBitCount = WordReadLE(fp);
  bmih.biCompression = DWordReadLE(fp);
  bmih.biSizeImage = DWordReadLE(fp);
  bmih.biXPelsPerMeter = LongReadLE(fp);
  bmih.biYPelsPerMeter = LongReadLE(fp);
  bmih.biClrUsed = DWordReadLE(fp);
  bmih.biClrImportant = DWordReadLE(fp);
  
  // Check info header 
  assert(bmih.biSize == BMP_BI_SIZE);
  assert(bmih.biWidth > 0);
  assert(bmih.biHeight > 0);
  assert(bmih.biPlanes == 1);
  assert(bmih.biBitCount == 24);  /* RGB */
  assert(bmih.biCompression == BI_RGB);   /* RGB */
  int lineLength = bmih.biWidth * 3;  /* RGB */
  if ((lineLength % 4) != 0) lineLength = (lineLength / 4 + 1) * 4;
  assert(bmih.biSizeImage == (unsigned int) lineLength * (unsigned int) bmih.biHeight);

  // Assign width, height, and number of pixels
  width = bmih.biWidth;
  height = bmih.biHeight;
  npixels = width * height;

  // Allocate unsigned char buffer for reading pixels
  int rowsize = 3 * width;
  if ((rowsize % 4) != 0) rowsize = (rowsize / 4 + 1) * 4;
  int nbytes = bmih.biSizeImage;
  unsigned char *buffer = new unsigned char [nbytes];
  if (!buffer) {
    fprintf(stderr, "Unable to allocate temporary memory for BMP file");
    fclose(fp);
    return 0;
  }

  // Read buffer 
  fseek(fp, (long) bmfh.bfOffBits, SEEK_SET);
  if (fread(buffer, 1, bmih.biSizeImage, fp) != bmih.biSizeImage) {
    fprintf(stderr, "Error while reading BMP file %s", filename);
    return 0;
  }

  // Close file
  fclose(fp);

  // Allocate pixels for image
  pixels = new R2Pixel [ width * height ];
  if (!pixels) {
    fprintf(stderr, "Unable to allocate memory for BMP file");
    fclose(fp);
    return 0;
  }

  // Assign pixels
  for (int j = 0; j < height; j++) {
    unsigned char *p = &buffer[j * rowsize];
    for (int i = 0; i < width; i++) {
      double b = (double) *(p++) / 255;
      double g = (double) *(p++) / 255;
      double r = (double) *(p++) / 255;
      R2Pixel pixel(r, g, b, 1);
      SetPixel(i, j, pixel);
    }
  }

  // Free unsigned char buffer for reading pixels
  delete [] buffer;

  // Return success
  return 1;
}



int R2Image::
WriteBMP(const char *filename) const
{
  // Open file
  FILE *fp = fopen(filename, "wb");
  if (!fp) {
    fprintf(stderr, "Unable to open image file: %s\n", filename);
    return 0;
  }

  // Compute number of bytes in row
  int rowsize = 3 * width;
  if ((rowsize % 4) != 0) rowsize = (rowsize / 4 + 1) * 4;

  // Write file header 
  BITMAPFILEHEADER bmfh;
  bmfh.bfType = BMP_BF_TYPE;
  bmfh.bfSize = BMP_BF_OFF_BITS + rowsize * height;
  bmfh.bfReserved1 = 0;
  bmfh.bfReserved2 = 0;
  bmfh.bfOffBits = BMP_BF_OFF_BITS;
  WordWriteLE(bmfh.bfType, fp);
  DWordWriteLE(bmfh.bfSize, fp);
  WordWriteLE(bmfh.bfReserved1, fp);
  WordWriteLE(bmfh.bfReserved2, fp);
  DWordWriteLE(bmfh.bfOffBits, fp);

  // Write info header 
  BITMAPINFOHEADER bmih;
  bmih.biSize = BMP_BI_SIZE;
  bmih.biWidth = width;
  bmih.biHeight = height;
  bmih.biPlanes = 1;
  bmih.biBitCount = 24;       /* RGB */
  bmih.biCompression = BI_RGB;    /* RGB */
  bmih.biSizeImage = rowsize * (unsigned int) bmih.biHeight;  /* RGB */
  bmih.biXPelsPerMeter = 2925;
  bmih.biYPelsPerMeter = 2925;
  bmih.biClrUsed = 0;
  bmih.biClrImportant = 0;
  DWordWriteLE(bmih.biSize, fp);
  LongWriteLE(bmih.biWidth, fp);
  LongWriteLE(bmih.biHeight, fp);
  WordWriteLE(bmih.biPlanes, fp);
  WordWriteLE(bmih.biBitCount, fp);
  DWordWriteLE(bmih.biCompression, fp);
  DWordWriteLE(bmih.biSizeImage, fp);
  LongWriteLE(bmih.biXPelsPerMeter, fp);
  LongWriteLE(bmih.biYPelsPerMeter, fp);
  DWordWriteLE(bmih.biClrUsed, fp);
  DWordWriteLE(bmih.biClrImportant, fp);

  // Write image, swapping blue and red in each pixel
  int pad = rowsize - width * 3;
  for (int j = 0; j < height; j++) {
    for (int i = 0; i < width; i++) {
      const R2Pixel& pixel = (*this)[i][j];
      double r = 255.0 * pixel.Red();
      double g = 255.0 * pixel.Green();
      double b = 255.0 * pixel.Blue();
      if (r >= 255) r = 255;
      if (g >= 255) g = 255;
      if (b >= 255) b = 255;
      fputc((unsigned char) b, fp);
      fputc((unsigned char) g, fp);
      fputc((unsigned char) r, fp);
    }

    // Pad row
    for (int i = 0; i < pad; i++) fputc(0, fp);
  }
  
  // Close file
  fclose(fp);

  // Return success
  return 1;  
}



////////////////////////////////////////////////////////////////////////
// PPM I/O
////////////////////////////////////////////////////////////////////////

int R2Image::
ReadPPM(const char *filename)
{
  // Open file
  FILE *fp = fopen(filename, "rb");
  if (!fp) {
    fprintf(stderr, "Unable to open image file: %s\n", filename);
    return 0;
  }

  // Read PPM file magic identifier
  char buffer[128];
  if (!fgets(buffer, 128, fp)) {
    fprintf(stderr, "Unable to read magic id in PPM file");
    fclose(fp);
    return 0;
  }

  // skip comments
  int c = getc(fp);
  while (c == '#') {
    while (c != '\n') c = getc(fp);
    c = getc(fp);
  }
  ungetc(c, fp);

  // Read width and height
  if (fscanf(fp, "%d%d", &width, &height) != 2) {
    fprintf(stderr, "Unable to read width and height in PPM file");
    fclose(fp);
    return 0;
  }
	
  // Read max value
  double max_value;
  if (fscanf(fp, "%lf", &max_value) != 1) {
    fprintf(stderr, "Unable to read max_value in PPM file");
    fclose(fp);
    return 0;
  }
	
  // Allocate image pixels
  pixels = new R2Pixel [ width * height ];
  if (!pixels) {
    fprintf(stderr, "Unable to allocate memory for PPM file");
    fclose(fp);
    return 0;
  }

  // Check if raw or ascii file
  if (!strcmp(buffer, "P6\n")) {
    // Read up to one character of whitespace (\n) after max_value
    int c = getc(fp);
    if (!isspace(c)) putc(c, fp);

    // Read raw image data 
    // First ppm pixel is top-left, so read in opposite scan-line order
    for (int j = height-1; j >= 0; j--) {
      for (int i = 0; i < width; i++) {
        double r = (double) getc(fp) / max_value;
        double g = (double) getc(fp) / max_value;
        double b = (double) getc(fp) / max_value;
        R2Pixel pixel(r, g, b, 1);
        SetPixel(i, j, pixel);
      }
    }
  }
  else {
    // Read asci image data 
    // First ppm pixel is top-left, so read in opposite scan-line order
    for (int j = height-1; j >= 0; j--) {
      for (int i = 0; i < width; i++) {
	// Read pixel values
	int red, green, blue;
	if (fscanf(fp, "%d%d%d", &red, &green, &blue) != 3) {
	  fprintf(stderr, "Unable to read data at (%d,%d) in PPM file", i, j);
	  fclose(fp);
	  return 0;
	}

	// Assign pixel values
	double r = (double) red / max_value;
	double g = (double) green / max_value;
	double b = (double) blue / max_value;
        R2Pixel pixel(r, g, b, 1);
        SetPixel(i, j, pixel);
      }
    }
  }

  // Close file
  fclose(fp);

  // Return success
  return 1;
}



int R2Image::
WritePPM(const char *filename, int ascii) const
{
  // Check type
  if (ascii) {
    // Open file
    FILE *fp = fopen(filename, "w");
    if (!fp) {
      fprintf(stderr, "Unable to open image file: %s\n", filename);
      return 0;
    }

    // Print PPM image file 
    // First ppm pixel is top-left, so write in opposite scan-line order
    fprintf(fp, "P3\n");
    fprintf(fp, "%d %d\n", width, height);
    fprintf(fp, "255\n");
    for (int j = height-1; j >= 0 ; j--) {
      for (int i = 0; i < width; i++) {
        const R2Pixel& p = (*this)[i][j];
        int r = (int) (255 * p.Red());
        int g = (int) (255 * p.Green());
        int b = (int) (255 * p.Blue());
        fprintf(fp, "%-3d %-3d %-3d  ", r, g, b);
        if (((i+1) % 4) == 0) fprintf(fp, "\n");
      }
      if ((width % 4) != 0) fprintf(fp, "\n");
    }
    fprintf(fp, "\n");

    // Close file
    fclose(fp);
  }
  else {
    // Open file
    FILE *fp = fopen(filename, "wb");
    if (!fp) {
      fprintf(stderr, "Unable to open image file: %s\n", filename);
      return 0;
    }
    
    // Print PPM image file 
    // First ppm pixel is top-left, so write in opposite scan-line order
    fprintf(fp, "P6\n");
    fprintf(fp, "%d %d\n", width, height);
    fprintf(fp, "255\n");
    for (int j = height-1; j >= 0 ; j--) {
      for (int i = 0; i < width; i++) {
        const R2Pixel& p = (*this)[i][j];
        int r = (int) (255 * p.Red());
        int g = (int) (255 * p.Green());
        int b = (int) (255 * p.Blue());
        fprintf(fp, "%c%c%c", r, g, b);
      }
    }
    
    // Close file
    fclose(fp);
  }

  // Return success
  return 1;  
}



////////////////////////////////////////////////////////////////////////
// JPEG I/O
////////////////////////////////////////////////////////////////////////


// #define USE_JPEG
#ifdef USE_JPEG
  extern "C" { 
#   define XMD_H // Otherwise, a conflict with INT32
#   undef FAR // Otherwise, a conflict with windows.h
#   include "jpeg/jpeglib.h"
  };
#endif



int R2Image::
ReadJPEG(const char *filename)
{
#ifdef USE_JPEG
  // Open file
  FILE *fp = fopen(filename, "rb");
  if (!fp) {
    fprintf(stderr, "Unable to open image file: %s\n", filename);
    return 0;
  }

  // Initialize decompression info
  struct jpeg_decompress_struct cinfo;
  struct jpeg_error_mgr jerr;
  cinfo.err = jpeg_std_error(&jerr);
  jpeg_create_decompress(&cinfo);
  jpeg_stdio_src(&cinfo, fp);
  jpeg_read_header(&cinfo, TRUE);
  jpeg_start_decompress(&cinfo);

  // Remember image attributes
  width = cinfo.output_width;
  height = cinfo.output_height;
  npixels = width * height;
  int ncomponents = cinfo.output_components;

  // Allocate pixels for image
  pixels = new R2Pixel [ npixels ];
  if (!pixels) {
    fprintf(stderr, "Unable to allocate memory for BMP file");
    fclose(fp);
    return 0;
  }

  // Allocate unsigned char buffer for reading image
  int rowsize = ncomponents * width;
  if ((rowsize % 4) != 0) rowsize = (rowsize / 4 + 1) * 4;
  int nbytes = rowsize * height;
  unsigned char *buffer = new unsigned char [nbytes];
  if (!buffer) {
    fprintf(stderr, "Unable to allocate temporary memory for JPEG file");
    fclose(fp);
    return 0;
  }

  // Read scan lines 
  // First jpeg pixel is top-left, so read pixels in opposite scan-line order
  while (cinfo.output_scanline < cinfo.output_height) {
    int scanline = cinfo.output_height - cinfo.output_scanline - 1;
    unsigned char *row_pointer = &buffer[scanline * rowsize];
    jpeg_read_scanlines(&cinfo, &row_pointer, 1);
  }

  // Free everything
  jpeg_finish_decompress(&cinfo);
  jpeg_destroy_decompress(&cinfo);

  // Close file
  fclose(fp);

  // Assign pixels
  for (int j = 0; j < height; j++) {
    unsigned char *p = &buffer[j * rowsize];
    for (int i = 0; i < width; i++) {
      double r, g, b, a;
      if (ncomponents == 1) {
        r = g = b = (double) *(p++) / 255;
        a = 1;
      }
      else if (ncomponents == 1) {
        r = g = b = (double) *(p++) / 255;
        a = 1;
        p++;
      }
      else if (ncomponents == 3) {
        r = (double) *(p++) / 255;
        g = (double) *(p++) / 255;
        b = (double) *(p++) / 255;
        a = 1;
      }
      else if (ncomponents == 4) {
        r = (double) *(p++) / 255;
        g = (double) *(p++) / 255;
        b = (double) *(p++) / 255;
        a = (double) *(p++) / 255;
      }
      else {
        fprintf(stderr, "Unrecognized number of components in jpeg image: %d\n", ncomponents);
        return 0;
      }
      R2Pixel pixel(r, g, b, a);
      SetPixel(i, j, pixel);
    }
  }

  // Free unsigned char buffer for reading pixels
  delete [] buffer;

  // Return success
  return 1;
#else
  fprintf(stderr, "JPEG not supported");
  return 0;
#endif
}


	

int R2Image::
WriteJPEG(const char *filename) const
{
#ifdef USE_JPEG
  // Open file
  FILE *fp = fopen(filename, "wb");
  if (!fp) {
    fprintf(stderr, "Unable to open image file: %s\n", filename);
    return 0;
  }

  // Initialize compression info
  struct jpeg_compress_struct cinfo;
  struct jpeg_error_mgr jerr;
  cinfo.err = jpeg_std_error(&jerr);
  jpeg_create_compress(&cinfo);
  jpeg_stdio_dest(&cinfo, fp);
  cinfo.image_width = width; 	/* image width and height, in pixels */
  cinfo.image_height = height;
  cinfo.input_components = 3;		/* # of color components per pixel */
  cinfo.in_color_space = JCS_RGB; 	/* colorspace of input image */
  cinfo.dct_method = JDCT_ISLOW;
  jpeg_set_defaults(&cinfo);
  cinfo.optimize_coding = TRUE;
  jpeg_set_quality(&cinfo, 100, TRUE);
  jpeg_start_compress(&cinfo, TRUE);
	
  // Allocate unsigned char buffer for reading image
  int rowsize = 3 * width;
  if ((rowsize % 4) != 0) rowsize = (rowsize / 4 + 1) * 4;
  int nbytes = rowsize * height;
  unsigned char *buffer = new unsigned char [nbytes];
  if (!buffer) {
    fprintf(stderr, "Unable to allocate temporary memory for JPEG file");
    fclose(fp);
    return 0;
  }

  // Fill buffer with pixels
  for (int j = 0; j < height; j++) {
    unsigned char *p = &buffer[j * rowsize];
    for (int i = 0; i < width; i++) {
      const R2Pixel& pixel = (*this)[i][j];
      int r = (int) (255 * pixel.Red());
      int g = (int) (255 * pixel.Green());
      int b = (int) (255 * pixel.Blue());
      if (r > 255) r = 255;
      if (g > 255) g = 255;
      if (b > 255) b = 255;
      *(p++) = r;
      *(p++) = g;
      *(p++) = b;
    }
  }



  // Output scan lines
  // First jpeg pixel is top-left, so write in opposite scan-line order
  while (cinfo.next_scanline < cinfo.image_height) {
    int scanline = cinfo.image_height - cinfo.next_scanline - 1;
    unsigned char *row_pointer = &buffer[scanline * rowsize];
    jpeg_write_scanlines(&cinfo, &row_pointer, 1);
  }

  // Free everything
  jpeg_finish_compress(&cinfo);
  jpeg_destroy_compress(&cinfo);

  // Close file
  fclose(fp);

  // Free unsigned char buffer for reading pixels
  delete [] buffer;

  // Return number of bytes written
  return 1;
#else
  fprintf(stderr, "JPEG not supported");
  return 0;
#endif
}






