#include <math.h>
#include <ctype.h>
#include <stdarg.h>
#include <stdint.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <string.h>
#include <time.h>
#include <vector>
#include <limits>
#include<algorithm>


#define SQR(A) ((A)*(A))
//#include <algorithm>
//#include "peakfinder8_my.h"
//#include "peakfinder9.cpp"
//#include "peakfinder8.cpp"
const int MASK_GOOD=1;
const int MASK_BAD=0;
const float MinVal = 1e-10;
const int loopsNum = 5; // number of loops to make
const int NumberOfRandomS = 100;

extern "C"{
float medianCutoff(float arr[], size_t n, float cutoff)  // cutoff from 0 to 1
{ if (n<1) return 0;
  const float maxval = 1e30;
  int upperlim = rint(cutoff*(n+0.01));
//  if (upperlim<1) upperlim=1;
  float _minim=arr[0];
  for (int jj=0; jj<upperlim; jj++)
  { _minim = maxval;
    int cmin = 0;
    for (int j=0; j<n; j++)
      if (arr[j]<_minim)
      { _minim = arr[j];
        cmin = j;
      }
    arr[cmin] = maxval;
  }
  return _minim;
}

    int SubLocalBG(float* inpAr, int stx, int enx, int sty, int eny, int fNumX, int radX, int radY, float badVal, float* smoothedAr)
{
  //float* inpAr = (float*)inpAr;
  float* arr1 = new float[(2*radX+1)*(2*radY+1)];
  int numX = enx-stx;
  int numY = eny-sty;
  size_t numComp = (size_t)numX*(size_t)numY;
  float* arrOld = new float[numComp];
  for (int xi=0; xi<numX; xi++)
    for (int yi=0; yi<numY; yi++)
      arrOld[xi+numX*yi] = inpAr[stx+xi+fNumX*(sty+yi)];

  if (smoothedAr != NULL)
    for (size_t i=0; i<numComp; i++)
      smoothedAr[i] = badVal;
//  for (size_t i=0; i<numComp; i++)
//    arrOld[i] = inpAr[i];

  for (int yi=0; yi<numY; yi++)
    for (int xi=0; xi<numX; xi++)
    { int nuel = 0;
      if (arrOld[xi+numX*yi]<badVal+1) continue;
      for (int bxi=-radX; bxi<=radX; bxi++)
        for (int byi=-radY; byi<=radY; byi++)
          if (xi+bxi>=0 && xi+bxi<numX && yi+byi>=0 && yi+byi<numY)
            if (arrOld[xi+bxi+numX*(yi+byi)]>badVal+1)
            { arr1[nuel] = arrOld[xi+bxi+numX*(yi+byi)];
              nuel++;
            }
      if (nuel>0)
      { //+float theval = quick_select(arr1, nuel);
        float theval = medianCutoff(arr1, nuel, 0.5);
        if (smoothedAr != NULL)
          smoothedAr[stx+xi+fNumX*(sty+yi)] = theval;
        else
          inpAr[stx+xi+fNumX*(sty+yi)] = arrOld[xi+numX*yi] - theval;
      }
    }

  delete arr1;
  delete arrOld;
  return 13;
}

void SubtractBgLoop(float* data, long pix_nn, int *pix_r, float hitfinderMinSNR, float ADCthresh, float bstpReg, float* radialcurve)
{
//! pix_r - array of radii of each pixel
//! hitfinderMinSNR - ?  it is =8 in cheetah.ini
//! ADCthresh -?  Some minimal threshold
  const float MinVal = 1e-10;
  const int loopsNum = 5; // number of loops to make

  // Apply mask (multiply data by 0 to ignore regions - this makes data below threshold for peak finding)
//mask  for(long i=0;i<pix_nn;i++) temp[i] *= mask[i];

  // Determine noise and offset as a funciton of radius
  float fminr, fmaxr;
  long lminr, lmaxr;
  fminr = 1e9;
  fmaxr = -1e9;

  // Figure out radius bounds
  for(long i=0;i<pix_nn;i++){
    if (pix_r[i] > fmaxr)
      fmaxr = pix_r[i];
    if (pix_r[i] < fminr)
      fminr = pix_r[i];
  }
  lmaxr = (int)ceil(fmaxr)+1;
  lminr = 0;
  //printf("%f, %f\n", fminr, fmaxr);
  // Allocate and zero arrays
  float* rsigma = (float*) calloc (lmaxr, sizeof(float));
  float* roffset = (float*) calloc (lmaxr, sizeof(float));
  long* rcount = (long*) calloc (lmaxr, sizeof(long));
  float* lowthreshold = (float*) calloc (lmaxr, sizeof(float));
  float* rthreshold = (float*) calloc (lmaxr, sizeof(float));

  for(long i=0; i<lmaxr; i++)
  { rthreshold[i] = 1e9;
    lowthreshold[i] = -1e9;
  }

  // Compute sigma and average of data values at each radius
  // From this, compute the ADC threshold to be applied at each radius
  // Iterate a few times to reduce the effect of positive outliers (ie: peaks)
  for(long counter=0; counter<loopsNum; counter++) {
    for(long i=0; i<lmaxr; i++) {
      roffset[i] = 0;
      rsigma[i] = 0;
      rcount[i] = 0;
    }
    for(long i=0;i<pix_nn;i++){
//mask      if(mask[i] != 0)
      if (data[i]>bstpReg+1)
      { long thisr = pix_r[i];
        if(data[i] < rthreshold[thisr] && data[i] > lowthreshold[thisr]) {
          roffset[thisr] += data[i];
          rsigma[thisr] += (data[i]*data[i]);
          rcount[thisr] += 1;
        }
      }
    }
    for(long i=0; i<lmaxr; i++) {
      if(rcount[i] == 0) {
        roffset[i] = 0;
        rsigma[i] = 0;
        rthreshold[i] = 1e9;
        lowthreshold[i] = -1e9;
      }
      else {
        float thisoffset = roffset[i]/(float)rcount[i];        //!!!!! this is average. I need array for each i at last loop
        roffset[i] = thisoffset;
        float someval = rsigma[i]/(float)rcount[i] - thisoffset*thisoffset;
        if (someval>MinVal) rsigma[i] = sqrt(someval);
        else rsigma[i] = 0;
//        rsigma[i] = sqrt(rsigma[i]/rcount[i] - thisoffset*thisoffset);;
        rthreshold[i] = roffset[i] + hitfinderMinSNR*rsigma[i];
        lowthreshold[i] = roffset[i] - hitfinderMinSNR*rsigma[i];
        if(rthreshold[i] < ADCthresh) rthreshold[i] = ADCthresh;
         //rthreshold[i] = ADCthresh;  // For testing
      }
    }
  }

//  FILE* afile = fopen("antaver.txt","wt");
//  for(long i=0; i<lmaxr; i++)
//    fprintf(afile,"%0.2f\n",roffset[i]);
//  fclose(afile);
  
  if (radialcurve != NULL){
    for(long i=0;i<lmaxr;i++)
      radialcurve[i] = roffset[i];}
  else{
    for(long i=0;i<pix_nn;i++)
      if (data[i]>bstpReg+1)
        data[i] -= roffset[pix_r[i]];
  }
  free (rsigma);
  free (roffset);
  free (rcount);
  free (rthreshold);
  free (lowthreshold);

  return;
}

int roundCpp(double xd)
{ //double _xd = xd + 103079215104.5;
  //return ((int*)&_xd)[0] >> 16;
  double t = ((xd) + 6755399441055744.0);
  return *((int *)(&t));
}

int roundMy(float x)
{ double x1 = x;
  return roundCpp(x1);
  //  return (int)round(x);
};

bool BuildRadialArray(size_t numEl, float *cox, float* coy, float istep,
                  int* pix_r, int* maxRad, float* pixelsR)
{
  // determine number of components in radial array

  *maxRad = 0;
  for (size_t i=0; i<numEl; i++)
  { float distR = sqrt(SQR(cox[i])+SQR(coy[i]))*istep;
    if (pixelsR != NULL) pixelsR[i] = distR;
    size_t _ind = roundMy(distR);
    if (_ind > *maxRad) *maxRad = _ind;
  }
  *maxRad += 1;
  
  // forming the array
  for (int i=0; i<numEl; i++)
  { int _ind = roundMy(sqrt(SQR(cox[i])+SQR(coy[i]))*istep);
  //should never happen
    if (_ind>=*maxRad)
      {printf("In radial averaging components are too far!\n"); continue;}
    pix_r[i] = _ind;
    //printf("%d\n", pix_r[i]);
  }
  return true;
}

bool MedianFilter1D(float* inpAr, int numX, int radX, float badVal)
{
  float* arr1 = new float[(2*radX+1)];
  float* arrOld = new float[numX];
  for (int i=0; i<numX; i++)
    arrOld[i] = inpAr[i];

  for (int xi=0; xi<numX; xi++)
    { int nuel = 0;
      for (int bxi=-radX; bxi<=radX; bxi++)
          if (xi+bxi>=0 && xi+bxi<numX)
            if (arrOld[xi+bxi]>badVal+1)
            { arr1[nuel] = arrOld[xi+bxi];
              nuel++;
            }
      if (nuel>0)
        inpAr[xi] = medianCutoff(arr1, nuel, 0.5);
      else inpAr[xi] = badVal;
    }
  delete[] arr1;
  delete[] arrOld;
  return true;
}

float PolarisationFactorDet(float detx, float dety, float detz, float degree)
{ if ((fabs(detx)<MinVal && fabs(dety)<MinVal)) return 1.;
//  float phi = atan2(dety,detx);
//  float teta = atan2(sqrt(detx*detx+dety*dety),detz);
//  float pol = degree*(1.-SQR(sin(teta)*cos(phi))) + (1.-degree)*(1.-SQR(sin(teta)*sin(phi)));

  float pdist2i = 1/(detx*detx + dety*dety + detz*detz);
//+  float pol = degree*(1.-SQR(detx)*pdist2i) + (1.-degree)*(1.-SQR(dety)*pdist2i);
  float pol = 1 - (SQR(dety)*(1.-degree) + SQR(detx)*degree)*pdist2i;

//  float pol = degree*(1.-SQR(sin(teta)*cos(phi))) ;//+ (1.-degree)*(1.-SQR(cos(teta)*sin(phi)));
//CF  pol = 1. - 2.*(1 - SQR(sin(phi)*sin(teta))) + 1 + SQR(cos(teta));
//old  float pol = 1.-SQR(sin(teta)*cos(phi));
  if (fabs(pol)<MinVal) pol=1;
//  else pol = 1/pol;
  return pol;
}

void CorrectPolarisation(float* im, int numcomp, float* cox, float* coy, float detdist, float badReg, float poldegree)
{ if (fabs(detdist)>MinVal)
    for (int i=0; i<numcomp; i++)
      if (im[i]>badReg+1) im[i] /= PolarisationFactorDet(cox[i],coy[i],detdist,poldegree);
}

bool MakePolarisationArray(float* pol, int numcomp, float* cox, float* coy, float detdist, float poldegree)
{ //if (fabs(detdist)<MinVal || pol==NULL) return false;
  for (int i=0; i<numcomp; i++)
    pol[i] += PolarisationFactorDet(cox[i],coy[i],detdist,poldegree);
  return true;
}



int MaskRingsSimple(float* data, char *mask, int *pix_r, int numPo, float badVal, float ringDiff, int smF)
//int MaskRingsSimple(float* data, char *mask, int numPoRad, int** radShells, int* radialNumComp, float badVal, float ringDiff, int smF)
{
  // find max rad
  int numPoRad = 0;
  for (int pix=0; pix<numPo; pix++)
    if (pix_r[pix] > numPoRad)
      numPoRad = pix_r[pix];
  numPoRad++;

  int smoothF = smF;
  if (smoothF<MinVal) smoothF = numPoRad/30; //???
  float difference = ringDiff;//0.2;

  // make radial average curve
  float* radaver = new float[numPoRad];
  int* irad = new int[numPoRad];
  for (int ring=0; ring<numPoRad; ring++)
  { radaver[ring] = 0.;
    irad[ring] = 0;
  }

  //calculate radial averaged curve
  for (int pix=0; pix<numPo; pix++)
  { int ring = pix_r[pix];
    if (mask != NULL)
    { if (mask[pix]!=MASK_GOOD) continue;
    } //else
    if (data[pix]<badVal+1) continue;
    radaver[ring] += data[pix];
    irad[ring]++;
  }
  for (int ring=0; ring<numPoRad; ring++)
    if (irad[ring]>0) radaver[ring] /= (float)irad[ring];
    else radaver[ring] = badVal-1;

  // smooth it
  float* smradaver = new float[numPoRad];
  for (int ring=0; ring<numPoRad; ring++)
    smradaver[ring] = radaver[ring];
  printf("%f\t%f\t\n",smradaver[100],smradaver[200]);
  MedianFilter1D(smradaver, numPoRad, smoothF, badVal);
  printf("%f\t%f\t\n",smradaver[100],smradaver[200]);

  // mark bad rings in irad
  for (int ring=0; ring<numPoRad; ring++)
    if (irad[ring]>0)
      if (fabs(smradaver[ring])>MinVal)
        if (fabs((smradaver[ring]-radaver[ring])/smradaver[ring])>difference)
          irad[ring] = -1;

  // mask all rings where difference between original and smoothed is too big
  for (int pix=0; pix<numPo; pix++)
    if (irad[pix_r[pix]]<0)
      if (mask == NULL) 
        data[pix] = badVal;
      else mask[pix] = MASK_BAD;

  delete[] irad;
  delete[] radaver;
  delete[] smradaver;
}

int maxmin_dbl(int *ad, float *max, float *min, long n)
{
    if (!ad || !max || !min) return 0;  /* validate parameters */

    *max = (float)LLONG_MIN;   /* initialize max/min to sufficiently */
    *min = (float)LLONG_MAX;   /* large negative/positive values.    */

    register long i;
    //for (i = 0; i < n; ++i) {
    //    if (ad[i] > *max) *max = ad[i];  /* test for new max */
    //    if (ad[i] < *min) *min = ad[i];  /* test for new min */
    //}

    for (i = n - 1; i>=0; --i)
    {
      if (ad[i] > *max) *max = ad[i];  /* test for new max */
      if (ad[i] < *min) *min = ad[i];  /* test for new min */
    }

    return 1;
}

void SubtractBgLoopNEW(float* data, long pix_nn, int *pix_r, float hitfinderMinSNR, float ADCthresh, float bstpReg, float* radialcurve)
{
//! pix_r - array of radii of each pixel
//! hitfinderMinSNR - ?  it is =8 in cheetah.ini
//! ADCthresh -?  Some minimal threshold
  //const float MinVal = 1e-10;
  //const int loopsNum = 5; // number of loops to make
  // Determine noise and offset as a funciton of radius
  float fminr, fmaxr;
  long lminr, lmaxr;
  fminr = 1e9;
  fmaxr = -1e9;
  
  // Figure out radius bounds
  if (not(maxmin_dbl(pix_r, &fmaxr, &fminr, pix_nn))){ /* get max/min from array */
      fprintf (stderr, "error: maxmin_dbl failed.\n");
  }
  //printf("%f, %f\n", fminr, fmaxr);
  lmaxr = (int)ceil(fmaxr)+1;
  lminr = 0;

  // Allocate and zero arrays
  //float* rsigma = (float*) calloc (lmaxr, sizeof(float));
  //float* roffset = (float*) calloc (lmaxr, sizeof(float));
  //long* rcount = (long*) calloc (lmaxr, sizeof(long));
  //float* lowthreshold = (float*) calloc (lmaxr, sizeof(float));
  //float* rthreshold = (float*) calloc (lmaxr, sizeof(float)); // calloc also set data to zero

  long lxmarFloatLength = lmaxr*sizeof(float);
  long lxmarLongLength = lmaxr*sizeof(long);

  float* rsigma = (float*) malloc (lxmarFloatLength);
  float* roffset = (float*) malloc (lxmarFloatLength);
  long* rcount = (long*) malloc (lxmarLongLength);
  float* lowthreshold = (float*) malloc (lxmarFloatLength);
  float* rthreshold = (float*) malloc (lxmarFloatLength); // calloc also set data to zero

  memset(rsigma, 0, lxmarFloatLength);
  memset(roffset, 0, lxmarFloatLength);
  memset(rthreshold, 0, lxmarFloatLength);
  memset(lowthreshold, 0, lxmarFloatLength);
  memset(rcount, 0, lxmarLongLength);



  //for(register long i=0; i<lmaxr; ++i)
  //{ 
  //  rthreshold[i] = 1e9;
  //  lowthreshold[i] = -1e9;
  //}

  for(register long i= lmaxr-1; i>=0; --i)
  { 
    rthreshold[i] = 1e9;
    lowthreshold[i] = -1e9;
  }


  // Compute sigma and average of data values at each radius
  // From this, compute the ADC threshold to be applied at each radius
  // Iterate a few times to reduce the effect of positive outliers (ie: peaks)
  
  //for(register long counter=0; counter<loopsNum; ++counter)
  for(register long counter=loopsNum-1; counter>=0; --counter) 
  {
    //for(register long i=0; i<pix_nn; ++i)
    for(register long i=pix_nn-1; i>=0; --i)
    {
      if (data[i]>bstpReg+1)
      { 
        long thisr = pix_r[i];
        if(data[i] < rthreshold[thisr] && data[i] > lowthreshold[thisr]) 
        {
          roffset[thisr] += data[i];
          rsigma[thisr] += (data[i]*data[i]);
          rcount[thisr] += 1;
        }
      }
    }

    //for(register long i=0; i<lmaxr; ++i) 
    for(register long i=lmaxr-1; i>=0; --i) 
    {
      if(rcount[i] == 0) 
      {
        roffset[i] = 0;
        rsigma[i] = 0;
        rthreshold[i] = 1e9;
        lowthreshold[i] = -1e9;
      }
      else  
      {
        float thisoffset = roffset[i]/(float)rcount[i];        //!!!!! this is average. I need array for each i at last loop
        roffset[i] = thisoffset;
        float someval = rsigma[i]/(float)rcount[i] - thisoffset*thisoffset;
        
        if (someval>MinVal) 
          {
            rsigma[i] = sqrt(someval);
          }
        else
          {
            rsigma[i] = 0;
          }
        float hitRsigma = hitfinderMinSNR*rsigma[i]; 
        rthreshold[i] = roffset[i] + hitRsigma;
        lowthreshold[i] = roffset[i] - hitRsigma;

        if(rthreshold[i] < ADCthresh) 
        {
          rthreshold[i] = ADCthresh;
        }
      }
    }
  }

  if (radialcurve != NULL){
    //for(register long i=0; i<lmaxr; ++i)
    for(register long i=lmaxr-1; i>=0; --i)
    {
      radialcurve[i] = roffset[i];
    }
  }
  else
  {
    //for(register long i=0;i<pix_nn;++i)
    for(register long i=pix_nn-1; i>=0; --i)
    {
      if (data[i]>bstpReg+1)
      {
        data[i] -= roffset[pix_r[i]];
      }
    }
  }
  free (rsigma);
  free (roffset);
  free (rcount);
  free (rthreshold);
  free (lowthreshold);

  return;
};

void Vector_SubtractBgLoopNEW(float* data, long pix_nn, int* pix_r, float hitfinderMinSNR,
                              float ADCthresh, float bstpReg, float* radialcurve)
{
  // Determine noise and offset as a funciton of radius
  const float MinVal = 1e-10;
  const int loopsNum = 5; // number of loops to make
  // Figure out radius bounds
  float fminr = 1e9;
  float fmaxr = 0;
  for (register long i = 0; i < pix_nn; ++i) {
    if (pix_r[i] < fminr) {
        fminr = pix_r[i];
    }
    if (pix_r[i] > fmaxr) {
        fmaxr = pix_r[i];
    }
  }


  //float fminr = *min_element(pix_r.begin(), pix_r.end());
  //float fmaxr  = *max_element(pix_r.begin(), pix_r.end());
 
  //printf("%f, %f\n",fminr, fmaxr);

  long lminr, lmaxr;

  lmaxr = (int)ceil(fmaxr)+1;
  lminr = 0;

  // Allocate and zero arrays

  long lxmarFloatLength = lmaxr*sizeof(float);
  long lxmarLongLength = lmaxr*sizeof(long);

  
  std::vector<float> rsigma(lxmarFloatLength, 0.);
  std::vector<float> roffset(lxmarFloatLength, 0.);
  std::vector<long> rcount(lxmarLongLength, 0.);
  std::vector<float> lowthreshold(lxmarFloatLength, -1e9);
  std::vector<float> rthreshold(lxmarFloatLength, 1e9);


  // Compute sigma and average of data values at each radius
  // From this, compute the ADC threshold to be applied at each radius
  // Iterate a few times to reduce the effect of positive outliers (ie: peaks)
  
  for(register long counter=0; counter<loopsNum; ++counter)
  {
    for(register long i=0; i<pix_nn; ++i)
    {
      if (data[i]>bstpReg+1)
      { 
        long thisr = pix_r[i];
        if(data[i] < rthreshold[thisr] && data[i] > lowthreshold[thisr]) 
        {
          roffset[thisr] += data[i];
          rsigma[thisr] += (data[i]*data[i]);
          rcount[thisr] += 1;
        }
      }
    }

    for(register long i=0; i<lmaxr; ++i) 
    {

      if(rcount[i] == 0) 
      {
        roffset[i] = 0;
        rsigma[i] = 0;
        rthreshold[i] = 1e9;
        lowthreshold[i] = -1e9;
      }
      else  
      {
        float thisoffset = roffset[i]/(float)rcount[i];        //!!!!! this is average. I need array for each i at last loop
        roffset[i] = thisoffset;
        float someval = rsigma[i]/(float)rcount[i] - thisoffset*thisoffset;

        rsigma[i] = someval > MinVal ? sqrt(someval) : 0;

        float hitRsigma = hitfinderMinSNR*rsigma[i]; 
        rthreshold[i] = roffset[i] + hitRsigma;
        lowthreshold[i] = roffset[i] - hitRsigma;

        if(rthreshold[i] < ADCthresh) 
        {
          rthreshold[i] = ADCthresh;
        }
      }
    }
  }

  if (radialcurve != NULL){
    for(register long i=0; i<lmaxr; ++i)
    {
      radialcurve[i] = roffset[i];
    }
  }
  else
  {
    for(register long i=0;i<pix_nn; ++i)
    {
      if (data[i]>bstpReg+1)
      {
        data[i] -= roffset[pix_r[i]];
      }
    }
  }
  //free (rsigma);
  //free (roffset);
  //free (rcount);
  //free (rthreshold);
  //free (lowthreshold);

  //rsigma.clear();
  //roffset.clear();
  //rcount.clear();
  //rthreshold.clear();
  //lowthreshold.clear();



  return;
};

bool MergeSeveralShells(size_t numEl, int maxShellsMerged, int minNumInShell, int* pix_r,
                        int* maxRad, int* numPerShell1, int*** shellElem1)
{
  if (maxShellsMerged <= 1) return false;
    int** shellElLoc = *shellElem1;

    int* numPerShell2 = new int[*maxRad];      // new one - now use this
    for (int i=0; i<*maxRad; i++)
      numPerShell2[i] = 0;
    int** shellElem2 = new int*[*maxRad];
    int numCombined = 1;
    int maxRadNew = 0;
    for (int i=0; i<*maxRad; i+=numCombined)
    { //numCombined = 1;  //?
      int curInShell = 0;//numPerShell1[i];
      int j;
      for (j=0; j<maxShellsMerged; j++)
      { if (i+j>=*maxRad)
          break;
        if (curInShell<minNumInShell)
          curInShell += numPerShell1[i+j];
        else
          break;
      }
      numCombined = j; // doesn't look right
      if (curInShell<1) continue;

      numPerShell2[maxRadNew] = curInShell;
      shellElem2[maxRadNew] = new int[numPerShell2[maxRadNew]];
      int kN=0;
      for (j=0; j<numCombined; j++)
        for (int k=0; k<numPerShell1[i+j]; k++)
        { shellElem2[maxRadNew][kN] = shellElLoc[i+j][k];
//        { shellElem2[maxRadNew][kN] = *shellElem1[i+j][k];
          kN++;
          if (kN>numPerShell2[maxRadNew])
          { printf("ups\n");
            break;
          }
        }
    // here kN must be = numPerShell[maxRadNew] = curInShell
      maxRadNew++;
    }
    for (int i=0; i<*maxRad; i++)
      if (shellElLoc[i] != NULL)
        delete shellElLoc[i];
    delete shellElLoc;

    //output
    for (int i=0; i<maxRadNew; i++) numPerShell1[i] = numPerShell2[i];
    delete numPerShell2;
    *maxRad = maxRadNew;
    *shellElem1 = shellElem2;

    // change pix_r to new bins
    if (pix_r != NULL)
    { for (int i=0; i<numEl; i++) pix_r[i] = -1;
      for (int i=0; i<*maxRad; i++)
        for (int k=0; k<numPerShell1[i]; k++)
          pix_r[shellElem2[i][k]] = i;
    }

  return true;
};



void RandomSelected(size_t numEl, int maxRad, int *pix_r){
  printf("-----Start RandomSelected-----\n");
  srand(time(NULL));
  int lenRandomIndex = NumberOfRandomS*sizeof(int);
  int* RandomIndex = (int*) malloc(lenRandomIndex);
  memset(RandomIndex, 0, lenRandomIndex);
  
  int lenNumPerShell1 = (maxRad) * sizeof(int);
  int* numPerShell1 =  (int*) malloc(lenNumPerShell1);
  memset(numPerShell1, 0, lenNumPerShell1);

  for (int i=0; i<numEl; i++){
    numPerShell1[pix_r[i]]++;
  }

  for(int r=0; r<maxRad; r++){
    int lenShellR = numPerShell1[r];
    if(lenShellR <= NumberOfRandomS){
      //printf("We skip %i because it has %i pixels in shell\n",r, lenShellR);
      continue;
    }

    int lenCurrentR = lenShellR*sizeof(int);
    int* currentR = (int *) malloc(lenCurrentR);
    memset(currentR, 0, lenCurrentR);
    
    int cInd = 0;
    for(size_t j=0; j<numEl; j++){
      if(pix_r[j] == r){
        currentR[cInd] = j;
        cInd++;
      }
    }

    //generate random numbers:
    int value[NumberOfRandomS];
    for (int i=0;i<NumberOfRandomS;i++)
    {
        int check; //variable to check or number is already used
        size_t pick_index; //variable to store the number in
        do
        {
        
        pick_index = rand() % lenShellR;
        //check or number is already used:
        check=1;
        for (int j=0;j<i;j++)
            if (pick_index == value[j]) //if number is already used
            {
                check=0; //set check to false
                break; //no need to check the other elements of value[]
            }
        } while (check == 0); //loop until new, unique number is found
        value[i]=pick_index; //store the generated number in the array
        RandomIndex[i] = currentR[pick_index];
    }
    for(size_t k=0; k < lenShellR; k++){
      int flag = 0;
      for (size_t q = 0; q < NumberOfRandomS; q++) { 
        if(RandomIndex[q] == currentR[k]){
          flag = 1;
        }
      }
      if(flag != 1){
        pix_r[currentR[k]] = -1; 
      }
    }

  }
  printf("-----End RandomSelected-----\n");
  return;
};


void RandomSelected2Arr(size_t numEl, int maxRad, int *pix_r, int **NonNegativeInd)
{
  srand(time(NULL));
  //int lenRandomIndex = NumberOfRandomS*sizeof(int);
  //int* RandomIndex = (int*) malloc(lenRandomIndex);
  //memset(RandomIndex, 0, lenRandomIndex);

  int RandomIndex [NumberOfRandomS];
  
  //int lenNumPerShell1 = (maxRad) * sizeof(int);
  //int* numPerShell1 =  (int*) malloc(lenNumPerShell1);
  //memset(numPerShell1, 0, lenNumPerShell1);
  
  int numPerShell1 [maxRad];
  for (int i=0; i<maxRad; ++i){
    numPerShell1[i]=0;
  }

  printf("I am in function\n");

  
  //Calculate the number of each pix_r per shell
  for (int i=0; i<numEl; ++i){
    numPerShell1[pix_r[i]]++;
  }
  
  int NlenNonNegativeIndexes = 0;  //number of nonnegative values in pix_r

  //Main part for random selection of NumberOfRandomS items 
  //for each pix_r
  for(int r=0; r<maxRad; ++r){
    int lenShellR = numPerShell1[r];

    //if number of items for this r is less than should be
    //selected, skip it. It means that we selected all items 
    //for this r
    if(lenShellR <= NumberOfRandomS){
      continue;
    }

    //int lenCurrentR = lenShellR*sizeof(int);
    //int* currentR = (int *) malloc(lenCurrentR); // array of indexes for this r
    //memset(currentR, 0, lenCurrentR);

    int currentR [lenShellR];
    
    //filling currentR array with all indexes for this r
    int cInd = 0;
    for(size_t j=0; j<numEl; ++j){
      if(pix_r[j] == r){
        currentR[cInd] = j;
        cInd++;
      }
    }

    //generate random indexes without repetiotion that should be selected from currentR
    //this indexes help us to save r value in these positions and others indexes for this r
    //set to -1
    int value[NumberOfRandomS];
    for (int i=0;i<NumberOfRandomS;++i)
    {
        int check; //variable to check or index is already used for this r
        size_t pick_index; //variable to store the random index in
        do
        {
        
        pick_index = rand() % lenShellR;
        //check or index is already used for this r:
        check=1;
        for (int j=0;j<i;++j)
            if (pick_index == value[j]) //if index is already used
            {
                check=0; //set check to false
                break; //no need to check the other elements of value[]
            }
        } while (check == 0); //loop until new, unique index is found
        value[i]=pick_index; //store the generated index in the array
        RandomIndex[i] = currentR[pick_index];
    }

    //set all positions for each r that are not on random selected to -1
    for(size_t k=0; k < lenShellR; ++k)
    {
      int flag = 0; // flag will be 1 if this index for this r in RandomIndex
      for (size_t q = 0; q < NumberOfRandomS; ++q) 
      { 
        if(RandomIndex[q] == currentR[k])
        {
          flag = 1; //this index is found
          break;
        }
      }
      if(flag != 1)
      {
        //index for this r not in RandomIndex, so set this index for this r to -1
        pix_r[currentR[k]] = -1;
      }
    }
  }

  
  for(size_t i =0; i < numEl; ++i)
  {
      NlenNonNegativeIndexes += (pix_r[i] > 0? 1:0);
  }

  printf("NlenNonNegativeIndexes %i\n", NlenNonNegativeIndexes);

  //int lenNonNegativeIndexes = (NlenNonNegativeIndexes) * sizeof(int);
  //int* NonNegativeIndexes =  (int*) malloc(lenNonNegativeIndexes);
  //memset(NonNegativeIndexes, 0, lenNonNegativeIndexes);

  int NonNegativeIndexes [NlenNonNegativeIndexes];

  for(size_t i=0;i<NlenNonNegativeIndexes;++i){
    NonNegativeIndexes[i] = 0;
  }
  
  size_t ind = 0;
  for(size_t i=0; i<numEl;++i){
    if(pix_r[i] > 0){
      NonNegativeIndexes[ind] = i;
      ind+=1;
    }
  }

  for(size_t i=0; i<NlenNonNegativeIndexes; i++){
    printf("NonNegativeIndexes[i] is %i\n", NonNegativeIndexes[i]);
  }
  
  *NonNegativeInd = NonNegativeIndexes;
  return;
};

//int ModRandomSelected2Arr(int numEl, int maxRad, int *pix_r, int **NonNegativeInd)
int ModRandomSelected2Arr(int numEl, int maxRad, int *pix_r)
{
  srand(time(NULL));

  //int RandomIndex [NumberOfRandomS];
  int lenRandomIndex = NumberOfRandomS*sizeof(int);
  int* RandomIndex = (int*) malloc(lenRandomIndex);
  memset(RandomIndex, 0, lenRandomIndex);
  
  
  //int numPerShell1 [maxRad];
  //for (int i=0; i<maxRad; ++i){
  //  numPerShell1[i]=0;
  //}

  int lenNumPerShell1 = (maxRad) * sizeof(int);
  int* numPerShell1 =  (int*) malloc(lenNumPerShell1);
  memset(numPerShell1, 0, lenNumPerShell1);

  printf("I am in function\n");

  
  //Calculate the number of each pix_r per shell
  for (int i=0; i<numEl; ++i){
    numPerShell1[pix_r[i]]++;
  }
  
  int NlenNonNegativeIndexes = 0;  //number of nonnegative values in pix_r

  //Main part for random selection of NumberOfRandomS items 
  //for each pix_r
  for(int r=0; r<maxRad; ++r){
    int lenShellR = numPerShell1[r];

    //if number of items for this r is less than should be
    //selected, skip it. It means that we selected all items 
    //for this r
    if(lenShellR <= NumberOfRandomS){
      continue;
    }

    //int currentR [lenShellR];
    int lenCurrentR = lenShellR*sizeof(int);
    int* currentR = (int *) valloc(lenCurrentR); // array of indexes for this r
    //memset(currentR, 0, lenCurrentR);
    
    //filling currentR array with all indexes for this r
    int cInd = 0;
    for(int j=0; j<numEl; ++j){
      if(pix_r[j] == r){
        currentR[cInd] = j;
        cInd++;
      }
    }

    //generate random indexes without repetiotion that should be selected from currentR
    //this indexes help us to save r value in these positions and others indexes for this r
    //set to -1
    int value[NumberOfRandomS];
    for (int i=0;i<NumberOfRandomS;++i)
    {
        int check; //variable to check or index is already used for this r
        int pick_index; //variable to store the random index in
        do
        {
        
        pick_index = rand() % lenShellR;
        //check or index is already used for this r:
        check=1;
        for (int j=0;j<i;++j)
            if (pick_index == value[j]) //if index is already used
            {
                check=0; //set check to false
                break; //no need to check the other elements of value[]
            }
        } while (check == 0); //loop until new, unique index is found
        value[i]=pick_index; //store the generated index in the array
        RandomIndex[i] = currentR[pick_index];
    }

    //set all positions for each r that are not on random selected to -1
    for(int k=0; k < lenShellR; ++k)
    {
      int flag = 0; // flag will be 1 if this index for this r in RandomIndex
      for (int q = 0; q < NumberOfRandomS; ++q) 
      { 
        if(RandomIndex[q] == currentR[k])
        {
          flag = 1; //this index is found
          break;
        }
      }
      if(flag != 1)
      {
        //index for this r not in RandomIndex, so set this index for this r to -1
        pix_r[currentR[k]] = -1;
      }
    }
  }

  
  for(int i =0; i < numEl; i++)
  {
      if(pix_r[i] > 0){
        NlenNonNegativeIndexes += 1;
      }
  }

  //int NonNegativeIndexes [NlenNonNegativeIndexes];

  printf("NlenNonNegativeIndexes in this function is  %i\n", NlenNonNegativeIndexes);
  //for(int i=0;i<NlenNonNegativeIndexes;++i){
  //  NonNegativeIndexes[i] = 0;
  //}

  int lenNlenNonNegativeIndexes = NlenNonNegativeIndexes*sizeof(int);
  int* NonNegativeIndexes = (int *) malloc(lenNlenNonNegativeIndexes);
  memset(NonNegativeIndexes, 0, lenNlenNonNegativeIndexes);

  //*NonNegativeInd = (int *) valloc(lenNlenNonNegativeIndexes);
  
  int ind = 0;
  for(int i=0; i<numEl;i++){
    if(pix_r[i] > 0){
      NonNegativeIndexes[ind] = i;
      //*NonNegativeInd[ind] = i;
      //printf("NonNegativeInd[%i] is  %i\n",ind, *NonNegativeInd[ind]);

      ind+=1;
    }
  }

  for(int i =0; i < NlenNonNegativeIndexes; i++)
  {
      printf("NonNegativeIndexes[%i] is  %i\n",i, NonNegativeIndexes[i]);
  }


  //*NonNegativeInd = NonNegativeIndexes;
  return NlenNonNegativeIndexes;
};

}