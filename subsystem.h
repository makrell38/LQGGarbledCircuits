#ifndef SUBSYSTEM_H
#define SUBSYSTEM_H
#include "fixedPoint.h"

class subSystem{
        public:
        fixedPoint **x0;
        int sizex0[2];
        fixedPoint **xr;
        int sizexr[2];
        fixedPoint **ur;
        int sizeur[2];
        fixedPoint **xk;
        int sizexk[2];
        fixedPoint **zk;
        int sizezk[2];
        fixedPoint **A;
        int sizeA[2];
        fixedPoint **B;
        int sizeB[2];
        fixedPoint **C;
        int sizeC[2];
        fixedPoint **xtemp;
        int sizextemp[2];
        subSystem(){}

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

        //No noise
        void measureState(fixedPoint **uk){
                fixedPoint *AxkIn = new fixedPoint[this->sizeA[0]];
                fixedPoint *Buk = new fixedPoint[this->sizeB[0]];
                matrixVecMul(this->A, this->xk, AxkIn, this->sizeA);
                matrixVecMul(this->B, uk, Buk, this->sizeB);
                for(int i = 0; i < this->sizeA[0]; i++){
                        this->xk[i][0] = AxkIn[i] + Buk[i] /*+ wk[i][0]*/;
                }
                delete[] AxkIn;
                delete[] Buk;
        }

        //No noise
        void computezk(){
                matrixMul(this->C,this->xk, this->zk, this->sizeC, this->sizexk);
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

                double **dataur;
                getFileSize("Data/ur.txt", this->sizeur);
                dataur = new double*[this->sizeur[0]];
                for(int i = 0; i < this->sizeur[0]; i++){
                        dataur[i] = new double[this->sizeur[1]];
                }

                readFile(dataur, "Data/ur.txt", this->sizeur);
                this->ur = new fixedPoint*[this->sizeur[0]];
                for(int i = 0; i < this->sizeur[0]; i++){
                        this->ur[i] = new fixedPoint[this->sizeur[1]];
                }


                double **dataxr;
                getFileSize("Data/xr.txt", this->sizexr);
                dataxr = new double*[this->sizexr[0]];
                for(int i = 0; i < this->sizexr[0]; i++){
                        dataxr[i] = new double[this->sizexr[1]];
                }

                readFile(dataxr, "Data/xr.txt", this->sizexr);
                this->xr = new fixedPoint*[this->sizexr[0]];
                for(int i = 0; i < this->sizexr[0]; i++){
                        this->xr[i] = new fixedPoint[this->sizexr[1]];
                }


                double **datax0;
                getFileSize("Data/x0.txt", this->sizex0);
                datax0 = new double*[this->sizex0[0]];
                for(int i = 0; i < this->sizex0[0]; i++){
                        datax0[i] = new double[this->sizex0[1]];
                }

                readFile(datax0, "Data/x0.txt", this->sizex0);
                this->x0 = new fixedPoint*[this->sizex0[0]];
                for(int i = 0; i < this->sizex0[0]; i++){
                        this->x0[i] = new fixedPoint[this->sizex0[1]];
                }

                this->sizezk[0] = this->sizeC[0];
                this->sizezk[1] = this->sizex0[1];
                this->zk = new fixedPoint*[this->sizezk[0]];
                for(int i = 0; i < this->sizezk[0]; i++){
                        this->zk[i] = new fixedPoint[this->sizezk[1]];
                }

                this->sizexk[0] = this->sizex0[0];
                this->sizexk[1] = this->sizex0[1];
                this->xk = new fixedPoint*[this->sizexk[0]];
                for(int i = 0; i < this->sizexk[0]; i++){
                        this->xk[i] = new fixedPoint[this->sizexk[1]];
                }

                this->sizextemp[0] = this->sizex0[0];
                this->sizextemp[1] = this->sizex0[1];
                this->xtemp = new fixedPoint*[this->sizextemp[0]];
                for(int i = 0; i < this->sizextemp[0]; i++){
                        this->xtemp[i] = new fixedPoint[this->sizextemp[1]];
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

                for(int i = 0; i < this->sizeur[0]; i++){
                        for(int j = 0; j < this->sizeur[1]; j++){
                                this->ur[i][j] = fixedPoint(dataur[i][j], 24, 24, ALICE);
                        }
                }

                for(int i = 0; i < this->sizexr[0]; i++){
                        for(int j = 0; j < this->sizexr[1]; j++){
                                this->xr[i][j] = fixedPoint(dataxr[i][j], 24, 24, ALICE);
                        }
                }

                for(int i = 0; i < this->sizex0[0]; i++){
                        for(int j = 0; j < this->sizex0[1]; j++){
                                this->x0[i][j] = fixedPoint(datax0[i][j], 24, 24, ALICE);
                        }
                }

                for(int i = 0; i < this->sizexk[0]; i++){
                        for(int j = 0; j < this->sizexk[1]; j++){
                                this->xk[i][j] = this->x0[i][j];
                        }
                }

                for(int i = 0; i < this->sizextemp[0]; i++){
                        for(int j = 0; j < this->sizextemp[1]; j++){
                                this->xtemp[i][j] = this->xk[i][j];
                        }
                }


                delete[] dataA;
                delete[] dataB;
                delete[] dataC;
                delete[] datax0;
                delete[] dataxr;
                delete[] dataur;

        }

};

#endif
