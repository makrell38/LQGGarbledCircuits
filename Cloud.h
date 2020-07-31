#ifndef CLOUD_H
#define CLOUD_H
#include "fixedPoint.h"

class Cloud{
        public:
        fixedPoint **gamma3;
        int sizegamma3[2];
        fixedPoint **gamma2;
        int sizegamma2[2];
        fixedPoint **gamma1;
        int sizegamma1[2];
        fixedPoint **ur;
        int sizeur[2];
        fixedPoint **xr;
        int sizexr[2];
        fixedPoint **L;
        int sizeL[2];
        fixedPoint **K;
        int sizeK[2];
        fixedPoint **uTilder;
        int sizeuTilder[2];
        fixedPoint **xgamma;
        int sizexgamma[2];
        fixedPoint **xHatk;
        int sizexHatk[2];
        fixedPoint **uk;
        int sizeuk[2];

        Cloud(){}

        void getInputs(fixedPoint **L, int *sizeL, fixedPoint **K, int * sizeK, fixedPoint **gamma3, int *sizegamma3, fixedPoint **gamma2, int *sizegamma2, fixedPoint **gamma1, int *sizegamma1, fixedPoint **xr, int *sizexr, fixedPoint **ur, int *sizeur, fixedPoint **xk, int *sizexk){
                this->sizeL[0] = sizeL[0];
                this->sizeL[1] = sizeL[1];
                this->L = new fixedPoint*[this->sizeL[0]];
                for(int i = 0; i < this->sizeL[0]; i++){
                        this->L[i] = new fixedPoint[this->sizeL[1]];
                }
                for(int i = 0; i < this->sizeL[0]; i++){
                        for(int j = 0; j < this->sizeL[1]; j++){
                                this->L[i][j] = L[i][j];
                        }
                }

                this->sizeK[0] = sizeK[0];
                this->sizeK[1] = sizeK[1];
                this->K = new fixedPoint*[this->sizeK[0]];
                for(int i = 0; i < this->sizeK[0]; i++){
                        this->K[i] = new fixedPoint[this->sizeK[1]];
                }
                for(int i = 0; i < this->sizeK[0]; i++){
                        for(int j = 0; j < this->sizeK[1]; j++){
                                this->K[i][j] = K[i][j];
                        }
                }

                this->sizegamma3[0] = sizegamma3[0];
                this->sizegamma3[1] = sizegamma3[1];
                this->gamma3 = new fixedPoint*[this->sizegamma3[0]];
                for(int i = 0; i < this->sizegamma3[0]; i++){
                        this->gamma3[i] = new fixedPoint[this->sizegamma3[1]];
                }
                for(int i = 0; i < this->sizegamma3[0]; i++){
                        for(int j = 0; j < this->sizegamma3[1]; j++){
                                this->gamma3[i][j] = gamma3[i][j];
                        }
                }

                this->sizegamma2[0] = sizegamma2[0];
                this->sizegamma2[1] = sizegamma2[1];
                this->gamma2 = new fixedPoint*[this->sizegamma2[0]];
                for(int i = 0; i < this->sizegamma2[0]; i++){
                        this->gamma2[i] = new fixedPoint[this->sizegamma2[1]];
                }
                for(int i = 0; i < this->sizegamma2[0]; i++){
                        for(int j = 0; j < this->sizegamma2[1]; j++){
                                this->gamma2[i][j] = gamma2[i][j];
                        }
                }

                this->sizegamma1[0] = sizegamma1[0];
                this->sizegamma1[1] = sizegamma1[1];
                this->gamma1 = new fixedPoint*[this->sizegamma1[0]];
                for(int i = 0; i < this->sizegamma1[0]; i++){
                        this->gamma1[i] = new fixedPoint[this->sizegamma1[1]];
                }
                for(int i = 0; i < this->sizegamma1[0]; i++){
                        for(int j = 0; j < this->sizegamma1[1]; j++){
                                this->gamma1[i][j] = gamma1[i][j];
                        }
                }

                this->sizexr[0] = sizexr[0];
                this->sizexr[1] = sizexr[1];
                this->xr = new fixedPoint*[this->sizexr[0]];
                for(int i = 0; i < this->sizexr[0]; i++){
                        this->xr[i] = new fixedPoint[this->sizexr[1]];
                }
                for(int i = 0; i < this->sizexr[0]; i++){
                        for(int j = 0; j < this->sizexr[1]; j++){
                                this->xr[i][j] = xr[i][j];
                        }
                }

                this->sizeur[0] = sizeur[0];
                this->sizeur[1] = sizeur[1];
                this->ur = new fixedPoint*[this->sizeur[0]];
                for(int i = 0; i < this->sizeur[0]; i++){
                        this->ur[i] = new fixedPoint[this->sizeur[1]];
                }
                for(int i = 0; i < this->sizeur[0]; i++){
                        for(int j = 0; j < this->sizeur[1]; j++){
                                this->ur[i][j] = ur[i][j];
                        }
                }

                this->sizexHatk[0] = sizexk[0];
                this->sizexHatk[1] = sizexk[1];
                this->xHatk = new fixedPoint*[this->sizexHatk[0]];
                for(int i = 0; i < this->sizexHatk[0]; i++){
                        this->xHatk[i] = new fixedPoint[this->sizexHatk[1]];
                }
                for(int i = 0; i < this->sizexHatk[0]; i++){
                        for(int j = 0; j < this->sizexHatk[1]; j++){
                                this->xHatk[i][j] = xk[i][j];
                        }
                }

                this->sizeuk[0] = sizeur[0];
                this->sizeuk[1] = sizeur[1];
                this->uk = new fixedPoint*[this->sizeuk[0]];
                for(int i = 0; i < this->sizeuk[0]; i++){
                        this->uk[i] = new fixedPoint[this->sizeuk[1]];
                }

        }

        void matrixMul(fixedPoint **A, fixedPoint **B, fixedPoint **ret, int *ASize, int *BSize){
                fixedPoint zero(0,24,24,PUBLIC);
                for(int i = 0; i < ASize[0]; i++){
                        for(int j = 0; j < BSize[1]; j++){
                                ret[i][j] = zero;
                                for(int k = 0; k < ASize[1]; k++){
                                        ret[i][j] = ret[i][j] + (A[i][k] * B[k][j]);
                                }
                        }
                }
        }


        void matrixVecMul(fixedPoint **A, fixedPoint **B, fixedPoint *ret, int *size){
                fixedPoint zero(0,24,24, PUBLIC);
                for(int i = 0; i < size[0]; i++){
                        ret[i] = zero;
                        for(int j = 0; j < size[1]; j++){
                                ret[i] = ret[i] + ((A[i][j]) * (B[j][0]));
                        }
                }

        }


        void computeConstants(){
                this->sizeuTilder[0] = sizeur[0];
                this->sizeuTilder[1] = sizeur[1];
                this->uTilder = new fixedPoint*[this->sizeuTilder[0]];
                for(int i = 0; i < this->sizeuTilder[0]; i++){
                        this->uTilder[i] = new fixedPoint[this->sizeuTilder[1]];
                }

                this->sizexgamma[0] = sizexr[0];
                this->sizexgamma[1] = sizexr[1];
                this->xgamma = new fixedPoint*[this->sizexgamma[0]];
                for(int i = 0; i < this->sizexgamma[0]; i++){
                        this->xgamma[i] = new fixedPoint[this->sizexgamma[1]];
                }

                fixedPoint **Kxr = new fixedPoint*[this->sizeK[0]];
                for(int i = 0; i < this->sizeK[0]; i++){
                        Kxr[i] = new fixedPoint[this->sizexr[1]];
                }
                matrixMul(this->K, this->xr, Kxr, this->sizeK, this->sizexr);
                for(int i = 0; i < this->sizeuTilder[0]; i++){
                        for(int j = 0; j < this->sizeuTilder[1]; j++){
                                this->uTilder[i][j] =  this->ur[i][j] + Kxr[i][j];
                        }
                }
                fixedPoint **gamma2xr = new fixedPoint*[this->sizegamma2[0]];
                for(int i = 0; i < this->sizegamma2[0]; i++){
                        gamma2xr[i] = new fixedPoint[this->sizexr[1]];
                }
                fixedPoint **gamma3ur = new fixedPoint*[this->sizegamma3[0]];
                for(int i = 0; i < this->sizegamma3[0]; i++){
                        gamma3ur[i] = new fixedPoint[this->sizeur[1]];
                }


                matrixMul(this->gamma2, this->xr, gamma2xr, this->sizegamma2, this->sizexr);
                matrixMul(this->gamma3, this->ur, gamma3ur, this->sizegamma3, this->sizeur);
                for(int i = 0; i < this->sizexgamma[0]; i++){
                        for(int j = 0; j < this->sizexgamma[1]; j++){
                                this->xgamma[i][j] = gamma3ur[i][j] + gamma2xr[i][j];
                        }
                }

        }

//Compute estimate xHat
        void computexHat(fixedPoint **zk){
                fixedPoint *xHatGamma = new fixedPoint[this->sizegamma1[0]];
                fixedPoint *Lzk = new fixedPoint[this->sizegamma1[0]];
                matrixVecMul(this->gamma1, this->xHatk, xHatGamma, this->sizegamma1);
                matrixVecMul(this->L, zk, Lzk, this->sizegamma1);
                for(int i = 0; i < this->sizegamma1[0]; i++){
                        xHatk[i][0] = xHatGamma[i] + Lzk[i] +xgamma[i][0];
                }
                delete[] xHatGamma;
                delete[] Lzk;
        }

//Computes control input uk
        void computeuk(){
                fixedPoint *kxHat = new fixedPoint[this->sizeK[0]];
                matrixVecMul(this->K, this->xHatk, kxHat, this->sizeK);
                for(int i = 0; i < this->sizeK[0]; i++){
                        this->uk[i][0] =  this->uTilder[i][0] - kxHat[i];
                }
                delete[] kxHat;
        }

};




#endif
