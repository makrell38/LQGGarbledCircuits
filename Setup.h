#ifndef SETUP_H
#define SETUP_H
#include "fixedPoint.h"

class Setup{
        public:
        fixedPoint **L;
        int sizeL[2];
        fixedPoint **K;
        int sizeK[2];
        fixedPoint **gamma3;
        int sizegamma3[2];
        fixedPoint **gamma2;
        int sizegamma2[2];
        fixedPoint **gamma1;
        int sizegamma1[2];
        fixedPoint **A;
        int sizeA[2];
        fixedPoint **B;
        int sizeB[2];
        fixedPoint **C;
        int sizeC[2];
        fixedPoint **LC;
        int sizeLC[2];
        fixedPoint **LCB;
        int sizeLCB[2];
        fixedPoint **LCA;
        int sizeLCA[2];

        Setup(){}

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


        void getFileSize(string input, int *size){
                int rowSize = 0, colSize = 0;
                fstream file;
                string inputLine, stringinput;
                file.open(input, ios::in);
                if (!file.is_open()){
                        cout<<"ERROR: file not opened"<<endl;
                        return;
                }
                while (getline(file, inputLine)){
                        rowSize++;
                        stringstream line(inputLine);
                        if(rowSize == 1){
                                while(getline(line, stringinput, ',')){
                                        colSize++;
                                }
                        }
                }
                size[0] = rowSize;
                size[1] = colSize;
                file.close();
        }

        void readFile(double **data, string input, int *size){
                int i = 0;
                int j = 0;
                fstream file;
                string inputLine, stringinput;
                file.open(input, ios::in);
                while(getline(file, inputLine)){
                        j = 0;
                        stringstream line(inputLine);
                        while(getline(line, stringinput, ',')){
                                data[i][j] = stod(stringinput);
                                j++;
                        }
                        i++;
                }
        }

        void inputData(){

                double **dataK;
                getFileSize("Data/K.txt", this->sizeK);
                dataK = new double*[this->sizeK[0]];
                for(int i = 0; i < this->sizeK[0]; i++){
                        dataK[i] = new double[this->sizeK[1]];
                }

                readFile(dataK, "Data/K.txt", this->sizeK);
                this->K = new fixedPoint*[this->sizeK[0]];
                for(int i = 0; i < this->sizeK[0]; i++){
                        this->K[i] = new fixedPoint[this->sizeK[1]];
                }

                double **dataL;
                getFileSize("Data/L.txt", this->sizeL);
                dataL = new double*[this->sizeL[0]];
                for(int i = 0; i < this->sizeL[0]; i++){
                        dataL[i] = new double[this->sizeL[1]];
                }

                readFile(dataL, "Data/L.txt", this->sizeL);
                this->L = new fixedPoint*[this->sizeL[0]];
                for(int i = 0; i < this->sizeL[0]; i++){
                        this->L[i] = new fixedPoint[this->sizeL[1]];
                }

                double **dataA;
                getFileSize("Data/A.txt", this->sizeA);
                dataA = new double*[this->sizeA[0]];
                for(int i = 0; i < this->sizeA[0]; i++){
                        dataA[i] = new double[this->sizeA[1]];
                }

                readFile(dataA, "Data/A.txt", this->sizeA);
                this->A = new fixedPoint*[this->sizeA[0]];
                for(int i = 0; i < this->sizeA[0]; i++){
                        this->A[i] = new fixedPoint[this->sizeA[1]];
                }


                double **dataB;
                getFileSize("Data/B.txt", this->sizeB);
                dataB = new double*[this->sizeB[0]];
                for(int i = 0; i < this->sizeB[0]; i++){
                        dataB[i] = new double[this->sizeB[1]];
                }

                readFile(dataB, "Data/B.txt", this->sizeB);
                this->B = new fixedPoint*[this->sizeB[0]];
                for(int i = 0; i < this->sizeB[0]; i++){
                        this->B[i] = new fixedPoint[this->sizeB[1]];
                }

                double **dataC;
                getFileSize("Data/C.txt", this->sizeC);
                dataC = new double*[this->sizeC[0]];
                for(int i = 0; i < this->sizeC[0]; i++){
                        dataC[i] = new double[this->sizeC[1]];
                }

                readFile(dataC, "Data/C.txt", this->sizeC);
                this->C = new fixedPoint*[this->sizeC[0]];
                for(int i = 0; i < this->sizeC[0]; i++){
                        this->C[i] = new fixedPoint[this->sizeC[1]];
                }

                //negative because K.txt holds the nagation of the control gain K
                for(int i = 0; i < this->sizeK[0]; i++){
                        for(int j = 0; j < this->sizeK[1]; j++){
                                this->K[i][j] = fixedPoint(-(dataK[i][j]), 24, 24, ALICE);
                        }
                }

                for(int i = 0; i < this->sizeL[0]; i++){
                        for(int j = 0; j < this->sizeL[1]; j++){
                                this->L[i][j] = fixedPoint(dataL[i][j], 24, 24, ALICE);
                        }
                }

                for(int i = 0; i < this->sizeA[0]; i++){
                        for(int j = 0; j < this->sizeA[1]; j++){
                                this->A[i][j] = fixedPoint(dataA[i][j], 24, 24, ALICE);
                        }
                }

                for(int i = 0; i < this->sizeB[0]; i++){
                        for(int j = 0; j < this->sizeB[1]; j++){
                                this->B[i][j] = fixedPoint(dataB[i][j], 24, 24, ALICE);
                        }
                }

                for(int i = 0; i < this->sizeC[0]; i++){
                        for(int j = 0; j < this->sizeC[1]; j++){
                                this->C[i][j] = fixedPoint(dataC[i][j], 24, 24, ALICE);
                        }
                }

                this->sizegamma3[0] = this->sizeB[0];
                this->sizegamma3[1] = this->sizeB[1];
                this->gamma3 = new fixedPoint*[this->sizegamma3[0]];
                for(int i = 0; i < this->sizegamma3[0]; i++){
                        this->gamma3[i] = new fixedPoint[this->sizegamma3[1]];
                }

                this->sizeLC[0] = this->sizeL[0];
                this->sizeLC[1] = this->sizeC[1];
                this->LC = new fixedPoint*[this->sizeLC[0]];
                for(int i = 0; i < this->sizeLC[1]; i++){
                        this->LC[i] = new fixedPoint[this->sizeLC[1]];
                }

                this->LCB = new fixedPoint*[this->sizeLC[0]];
                for(int i = 0; i < this->sizeLC[0]; i++){
                        this->LCB[i] = new fixedPoint[this->sizeB[1]];
                }

                this->sizegamma2[0] = this->sizegamma3[0];
                this->sizegamma2[1] = this->sizeK[1];
                this->gamma2 = new fixedPoint*[this->sizegamma2[0]];
                for(int i = 0; i < this->sizegamma2[0]; i++){
                        this->gamma2[i] = new fixedPoint[this->sizegamma2[1]];
                }

                this->sizegamma1[0] = this->sizeA[0];
                this->sizegamma1[1] = this->sizeA[1];
                this->gamma1 = new fixedPoint*[this->sizegamma1[0]];
                for(int i = 0; i < this->sizegamma1[0]; i++){
                        this->gamma1[i] = new fixedPoint[this->sizegamma1[1]];
                }

                this->LCA = new fixedPoint*[this->sizeA[0]];
                for(int i = 0; i < this->sizeA[0]; i++){
                        this->LCA[i] = new fixedPoint[this->sizeA[1]];
                }
        }

        void computeConstants(){
                //compute gamma3
                matrixMul(this->L,this->C,this->LC, this->sizeL, this->sizeC);
                matrixMul( this->LC, this->B, this->LCB, this->sizeLC, this->sizeB);
                for(int i = 0; i < this->sizegamma3[0]; i++){
                        for(int j = 0; j < this->sizegamma3[1]; j++){
                                this->gamma3[i][j] = this->B[i][j] - this->LCB[i][j];
                        }
                }

                //compute gamma2
                matrixMul(this->gamma3, this->K, this->gamma2, this->sizegamma3, this->sizeK);

                //compute gamma1
                matrixMul(this->LC, this->A, this->LCA, this->sizeLC, this->sizeA);
                for(int i = 0; i < this->sizegamma1[0]; i++){
                        for(int j = 0; j < this->sizegamma1[0]; j++){
                                this->gamma1[i][j] = (this->A[i][j] - this->LCA[i][j]) - this->gamma2[i][j];
                        }
                }

        }

};

#endif
