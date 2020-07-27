#include "emp-tool/emp-tool.h"
#include "emp-sh2pc/emp-sh2pc.h"
#include <iomanip>
using namespace std;
using namespace emp;




/*
        24-47 is Integer bits
        0-23 is fractional bits

*/

class fixedPoint {
        public:
        Bit *bits;
        Integer val;
        int length = 0;
        int party;


        fixedPoint(double input, int fractLen, int intLen, int party = PUBLIC){
                length = fractLen + intLen;
                bits = new Bit[length];
                to_binary(input, bits, length, party);
                val = Integer(length, bits);
                party = party;
        }

        fixedPoint(){}
        fixedPoint(const fixedPoint& in){
                length = in.length;
                party=in.party;
                bits = new Bit[length];
                memcpy(bits, in.bits, sizeof(Bit)*length);
        }


        template<typename O>
        O reveal(int party=PUBLIC);

        template<typename... Args>
        size_t bool_s(size_t size, Args... args) {
                return size;
        }

        fixedPoint operator*(fixedPoint & rhs);

        void bool_dat(bool* data, size_t len, long long num, size_t startLoc) {
                string str, bin;
                int temp = 0;
                str = to_string(num);
                bin = dec_to_bin(str);
                reverse(bin.begin(), bin.end());
                int l = (bin.size() > (size_t)len ? len : bin.size());
                for(int i = temp; i < l+temp; ++i)
                        data[startLoc +i-temp] = (bin[i] == '1');
                for (size_t i = l; i < len; ++i){
                        data[startLoc+i] = data[startLoc+l-1];
                }
        }


        void in_it(Bit * bits, const bool* b, int length, int party) {
                block * bbits = (block *) bits;
                if (party == PUBLIC) {
                        block one = CircuitExecution::circ_exec->public_label(true);
                        block zero = CircuitExecution::circ_exec->public_label(false);
                        for(int i = 0; i < length; ++i)
                                bbits[i] = b[i] ? one : zero;
                }
                else {
                        ProtocolExecution::prot_exec->feed((block *)bits, party, b, length);
                }
        }

        void to_binary(double x, Bit *dest, int length, int party){
                double intPart, fractPart;
                bool *b = new bool[48];
                fractPart = modf(x, &intPart);
                for(int i = 0; i < 24; i++){
                        b[i] = 0;
                        b[i+24] = (intPart<0? 1: 0);
                }
                if(x < 0){
                        intPart = floor(x);
                        fractPart = 1 + fractPart;
                }
                bool_dat(b, 24, (long long)intPart, 24);
                bool_dat(b, 24, (long long)(fractPart * pow(2,24)), 0);
                in_it(dest, b, 48, party);
        }

        void full_adder(Bit *a, Bit *b, Bit *cin, Bit *cout, Bit *res, int size){
                Bit carry, bxc, axc, t, carryTemp;
                for(int i = 0; i < size-1; i++){
                        axc = a[i] ^ carry;
                        bxc = b[i] ^ carry;
                        res[i] = a[i] ^ bxc;
                        t = axc & bxc;
                        carry =carry^t;
                }
                carryTemp = carry;
                axc = a[size-1] ^ carry;
                bxc = b[size-1] ^ carry;
                res[size-1] = a[size-1] ^ bxc;
                t = axc & bxc;
                carry = carry^t;
                //Will truncate if there is overflow, but is not needed on the current inputs
//              if(size == 48){
//                      if((carry^carryTemp).reveal<bool>()){
//                              for(int i = 24; i < 47; i++){
//                                      res[i] = 1;
//                              }
//                              res[47] = 0;
//                              std::cout<<"OVERFLOW"<<endl;
//                      }
//              }

        }

        fixedPoint operator+(const fixedPoint & rhs) {
                fixedPoint ret(*this);
                full_adder(this->bits, rhs.bits, nullptr, nullptr,ret.bits, 48);
                return ret;

        }

        void full_sub(Bit *a, Bit *b, Bit *borrowOut, Bit *borrowIn, Bit *ret, int size){
                Bit borrow,bxc,bxa,t;
                int skipLast; int i = 0;
                if(size==0) {
                        if(borrowIn && borrowOut)
                                *borrowOut = *borrowIn;
                        return;
                }
                if(borrowIn)
                        borrow = *borrowIn;
                else
                        borrow = false;
                // skip AND on last bit if borrowOut==NULL
                skipLast = (borrowOut == nullptr);
                while(size-- > skipLast) {
                        bxa = a[i] ^ b[i];
                        bxc = borrow ^ b[i];
                        ret[i] = bxa ^ borrow;
                        t = bxa & bxc;
                        borrow = borrow ^ t;
                        ++i;
                }
                if(borrowOut != nullptr) {
                        *borrowOut = borrow;
                }
                else
                        ret[i] = a[i] ^ b[i] ^ borrow;
        }



        fixedPoint operator-(const fixedPoint& rhs) {
                fixedPoint ret(*this);
                full_sub(this->bits, rhs.bits, nullptr, nullptr, ret.bits, 48);
                return ret;
        }


        fixedPoint operator-() const {
                return (fixedPoint(0,24,24,PUBLIC) - (*this));

        }


        void full_mul(Bit *a, Bit *b, Bit *ret){
                Bit one(1,PUBLIC);
                Bit zero(0, PUBLIC);
                Bit * sum = new Bit[96];
                Bit * temp = new Bit[49];
                for(int i = 0; i < 96; i++){
                        sum[i] = false;
                }
                for(int i = 0; i < 48; i++){
                        for(int k = 0; k < 48; k++){
                                temp[k] = a[i] & b[k];
                        }
                        if ( i == 0){
                                temp[48] = one;
                                temp[47] = !temp[47];
                        }
                        else if (i < 47){
                                temp[47] = !temp[47];
                                temp[48] = zero;
                        }
                        else{
                                for(int j = 0; j < 47; j++){
                                        temp[j] = !temp[j];
                                }
                                temp[48] = one;
                        }
                        full_adder(sum+i, temp, nullptr, nullptr, sum+i, 49);
                }
                memcpy(ret, sum+24, sizeof(Bit)*48);
                delete[] sum;
                delete[] temp;
        }






};

        template<>
        double fixedPoint::reveal<double>(int party) {
                bool *integers = new bool[24];
                bool *decimals = new bool[24];
                string wholenum="", fraction="";
                double result;
                ProtocolExecution::prot_exec->reveal(integers, party, (block *)bits+24,  24);
                ProtocolExecution::prot_exec->reveal(decimals, party, (block *)bits,  24);

                fraction +='0';
                for(int i = 24-1; i >=0; --i){
                        wholenum += (integers[i]? '1' : '0');
                        fraction += (decimals[i]? '1' : '0');
                }
                wholenum = bin_to_dec(wholenum);
                fraction = bin_to_dec(fraction);
                result = stod(wholenum) + stod(fraction)/pow(2,24);
                delete[] integers;
                delete[] decimals;
                return result;
        }

        template<>
        string fixedPoint::reveal<string>(int party)  {
                return to_string(this->reveal<double>(party));
        }

        fixedPoint fixedPoint::operator*( fixedPoint & rhs){
                fixedPoint ret(*this);
                full_mul(this->bits,rhs.bits, ret.bits);
                return ret;
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

//Compute estimate xHat
void xHat(fixedPoint **gamma, fixedPoint **xHatIn, fixedPoint **L, fixedPoint **zk, fixedPoint **xGamma, fixedPoint **ret, int *size){
        fixedPoint *xHatGamma = new fixedPoint[size[0]];
        fixedPoint *Lzk = new fixedPoint[size[0]];
        matrixVecMul(gamma, xHatIn, xHatGamma, size);
        matrixVecMul(L, zk, Lzk, size);
        for(int i = 0; i < size[0]; i++){
                ret[i][0] = xHatGamma[i] + Lzk[i] +xGamma[i][0];
        }
        delete[] xHatGamma;
        delete[] Lzk;
}

//Computes control input uk
void uK(fixedPoint **k, fixedPoint **xHatIn, fixedPoint **ur, fixedPoint **ret, int *size){
        fixedPoint *kxHat = new fixedPoint[size[0]];
        matrixVecMul(k, xHatIn, kxHat, size);
        for(int i = 0; i < size[0]; i++){
                ret[i][0] =  ur[i][0] - kxHat[i];
        }
        delete[] kxHat;
}


//No noise
void measureState(fixedPoint **A, fixedPoint ** xkIn, fixedPoint **B, fixedPoint **uk, /*fixedPoint **wk,*/ fixedPoint **ret, int *sizeA, int *sizeB){
        fixedPoint *AxkIn = new fixedPoint[sizeA[0]];
        fixedPoint *Buk = new fixedPoint[sizeB[0]];
        matrixVecMul(A, xkIn, AxkIn, sizeA);
        matrixVecMul(B, uk, Buk, sizeB);
        for(int i = 0; i < sizeA[0]; i++){
                ret[i][0] = AxkIn[i] + Buk[i] /*+ wk[i][0]*/;
        }
        delete[] AxkIn;
        delete[] Buk;
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


int main(int argc, char** argv){
        int port, party;
        parse_party_and_port(argv, &party, &port);
        NetIO * io = new NetIO(party==ALICE ? nullptr: "127.0.0.1", port);
        setup_semi_honest(io, party);


//Read in .txt files
        double **dataL;
        int sizeL[2];
        getFileSize("L.txt", sizeL);
        dataL = new double*[sizeL[0]];
        for(int i = 0; i < sizeL[0]; i++){
                dataL[i] = new double[sizeL[1]];
        }

        readFile(dataL, "L.txt", sizeL);
        fixedPoint **L = new fixedPoint*[sizeL[0]];
        for(int i = 0; i < sizeL[0]; i++){
                L[i] = new fixedPoint[sizeL[1]];
        }

        for(int i = 0; i <sizeL[0]; i++){
                for(int j = 0; j < sizeL[1]; j++){
                        L[i][j] = fixedPoint(dataL[i][j], 24, 24, ALICE);
                }
        }

        double **dataA;
        int sizeA[2];
        getFileSize("A.txt", sizeA);
        dataA = new double*[sizeA[0]];
        for(int i = 0; i < sizeA[0]; i++){
                dataA[i] = new double[sizeA[1]];
        }

        readFile(dataA, "A.txt", sizeA);
        fixedPoint **A = new fixedPoint*[sizeA[0]];
        for(int i = 0; i < sizeA[0]; i++){
                A[i] = new fixedPoint[sizeA[1]];
        }

        for(int i = 0; i <sizeA[0]; i++){
                for(int j = 0; j < sizeA[1]; j++){
                        A[i][j] = fixedPoint(dataA[i][j], 24, 24, ALICE);
                }
        }

        double **dataB;
        int sizeB[2];
        getFileSize("B.txt", sizeB);
        dataB = new double*[sizeB[0]];
        for(int i = 0; i < sizeB[0]; i++){
                dataB[i] = new double[sizeB[1]];
        }

        readFile(dataB, "B.txt", sizeB);
        fixedPoint **B = new fixedPoint*[sizeB[0]];
        for(int i = 0; i < sizeB[0]; i++){
                B[i] = new fixedPoint[sizeB[1]];
        }

        for(int i = 0; i <sizeB[0]; i++){
                for(int j = 0; j < sizeB[1]; j++){
                        B[i][j] = fixedPoint(dataB[i][j], 24, 24, ALICE);
                }
        }

        double **dataC;
        int sizeC[2];
        getFileSize("C.txt", sizeC);
        dataC = new double*[sizeC[0]];
        for(int i = 0; i < sizeC[0]; i++){
                dataC[i] = new double[sizeC[1]];
        }

        readFile(dataC, "C.txt", sizeC);
        fixedPoint **C = new fixedPoint*[sizeC[0]];
        for(int i = 0; i < sizeC[0]; i++){
                C[i] = new fixedPoint[sizeC[1]];
        }

        for(int i = 0; i <sizeC[0]; i++){
                for(int j = 0; j < sizeC[1]; j++){
                        C[i][j] = fixedPoint(dataC[i][j], 24, 24, ALICE);
                }
        }

        double **dataK;
        int sizeK[2];
        getFileSize("K.txt", sizeK);
        dataK = new double*[sizeK[0]];
        for(int i = 0; i < sizeK[0]; i++){
                dataK[i] = new double[sizeK[1]];
        }

        readFile(dataK, "K.txt", sizeK);
        fixedPoint **K = new fixedPoint*[sizeK[0]];
        for(int i = 0; i < sizeK[0]; i++){
                K[i] = new fixedPoint[sizeK[1]];
        }
        //negative because K.txt holds the nagation of the control gain K
        for(int i = 0; i <sizeK[0]; i++){
                for(int j = 0; j < sizeK[1]; j++){
                        K[i][j] = fixedPoint(-(dataK[i][j]), 24, 24, ALICE);
                }
        }

        double **dataur;
        int sizeur[2];
        getFileSize("ur.txt", sizeur);
        dataur = new double*[sizeur[0]];
        for(int i = 0; i < sizeur[0]; i++){
                dataur[i] = new double[sizeur[1]];
        }

        readFile(dataur, "ur.txt", sizeur);
        fixedPoint **ur = new fixedPoint*[sizeur[0]];
        for(int i = 0; i < sizeur[0]; i++){
                ur[i] = new fixedPoint[sizeur[1]];
        }

        for(int i = 0; i <sizeur[0]; i++){
                for(int j = 0; j < sizeur[1]; j++){
                        ur[i][j] = fixedPoint(dataur[i][j], 24, 24, ALICE);
                }
        }

        double **dataxr;
        int sizexr[2];
        getFileSize("xr.txt", sizexr);
        dataxr = new double*[sizexr[0]];
        for(int i = 0; i < sizexr[0]; i++){
                dataxr[i] = new double[sizexr[1]];
        }

        readFile(dataxr, "xr.txt", sizexr);
        fixedPoint **xr = new fixedPoint*[sizexr[0]];
        for(int i = 0; i < sizexr[0]; i++){
                xr[i] = new fixedPoint[sizexr[1]];
        }

        for(int i = 0; i <sizexr[0]; i++){
                for(int j = 0; j < sizexr[1]; j++){
                        xr[i][j] = fixedPoint(dataxr[i][j], 24, 24, ALICE);
                }
        }

        double **datax0;
        int sizex0[2];
        getFileSize("x0.txt", sizex0);
        datax0 = new double*[sizex0[0]];
        for(int i = 0; i < sizex0[0]; i++){
                datax0[i] = new double[sizex0[1]];
        }

        readFile(datax0, "x0.txt", sizex0);
        fixedPoint **x0 = new fixedPoint*[sizex0[0]];
        for(int i = 0; i < sizex0[0]; i++){
                x0[i] = new fixedPoint[sizex0[1]];
        }

        for(int i = 0; i <sizex0[0]; i++){
                for(int j = 0; j < sizex0[1]; j++){
                        x0[i][j] = fixedPoint(datax0[i][j], 24, 24, ALICE);
                }
        }


//Initialize needed variables
        int *sizegamma3 = new int[2];
        sizegamma3[0] = sizeB[0];
        sizegamma3[1] = sizeB[1];
        fixedPoint **gamma3 = new fixedPoint*[sizegamma3[0]];
        for(int i = 0; i < sizegamma3[0]; i++){
                gamma3[i] = new fixedPoint[sizegamma3[1]];
        }

        int sizeLC[2];
        sizeLC[0] = sizeL[0];
        sizeLC[1] = sizeC[1];
        fixedPoint **LC = new fixedPoint*[sizeLC[0]];
        for(int i = 0; i < sizeL[1]; i++){
                LC[i] = new fixedPoint[sizeLC[1]];
        }

        fixedPoint **LCB = new fixedPoint*[sizeLC[0]];
        for(int i = 0; i < sizeLC[0]; i++){
                LCB[i] = new fixedPoint[sizeB[1]];
        }

        int sizegamma2[2];
        sizegamma2[0] = sizegamma3[0];
        sizegamma2[1] = sizeK[1];
        fixedPoint **gamma2 = new fixedPoint*[sizegamma2[0]];
        for(int i = 0; i < sizegamma2[0]; i++){
                gamma2[i] = new fixedPoint[sizegamma2[1]];
        }


        int sizegamma1[2];
        sizegamma1[0] = sizeA[0];
        sizegamma1[1] = sizeA[1];
        fixedPoint **gamma1 = new fixedPoint*[sizegamma1[0]];
        for(int i = 0; i < sizegamma1[0]; i++){
                gamma1[i] = new fixedPoint[sizegamma1[1]];
        }

        fixedPoint **LCA = new fixedPoint*[sizeA[0]];
        for(int i = 0; i < sizeA[0]; i++){
                LCA[i] = new fixedPoint[sizeA[1]];
        }

        int sizeuTilder[2];
        sizeuTilder[0] = sizeur[0];
        sizeuTilder[1] = sizeur[1];
        fixedPoint **uTilder = new fixedPoint*[sizeuTilder[0]];
        for(int i = 0; i < sizeuTilder[0]; i++){
                uTilder[i] = new fixedPoint[sizeuTilder[1]];
        }
        fixedPoint **Kxr = new fixedPoint*[sizeur[0]];
        for(int i = 0; i < sizeur[0]; i++){
                Kxr[i] = new fixedPoint[sizeur[1]];
        }

        int sizexGamma[2];
        sizexGamma[0] = sizegamma2[0];
        sizexGamma[1] = sizexr[1];
        fixedPoint **xGamma = new fixedPoint*[sizexGamma[0]];
        for(int i = 0; i < sizexGamma[0]; i++){
                xGamma[i] = new fixedPoint[sizexGamma[1]];
        }
        fixedPoint **gamma2xr = new fixedPoint*[sizegamma2[0]];
        fixedPoint **gamma3ur = new fixedPoint*[sizegamma3[0]];
        for(int i = 0; i < sizegamma2[0]; i++){
                gamma2xr[i] = new fixedPoint[sizexr[1]];
                gamma3ur[i] = new fixedPoint[sizeur[1]];
        }

        int sizezk[2];
        sizezk[0] = sizeC[0];
        sizezk[1] = sizex0[1];
        fixedPoint **zk = new fixedPoint*[sizezk[0]];
        for(int i = 0; i < sizezk[0]; i++){
                zk[i] = new fixedPoint[sizezk[1]];
        }



        int sizexHatk[2];
        sizexHatk[0] = sizexGamma[0];
        sizexHatk[1] = sizexGamma[1];
        fixedPoint **xHatk = new fixedPoint*[sizexHatk[0]];
        for(int i = 0; i < sizexHatk[0]; i++){
                xHatk[i] = new fixedPoint[sizexHatk[1]];
        }
        for(int i = 0; i < sizexHatk[0]; i++){
                for(int j = 0; j < sizexHatk[1]; j++){
                        xHatk[i][j] = x0[i][j];
                }
        }


        int sizeuk[2];
        sizeuk[0] = sizeuTilder[0];
        sizeuk[1] = sizeuTilder[1];
        fixedPoint **uk = new fixedPoint*[sizeuk[0]];
        for(int i = 0; i < sizeuk[0]; i++){
                uk[i] = new fixedPoint[sizeuk[1]];
        }
        for(int i = 0; i < sizeuk[0]; i++){
                for(int j = 0; j < sizeuk[1]; j++){
                        uk[i][j] = fixedPoint(0,24,24,ALICE);
                }
        }

        int sizexk[2];
        sizexk[0] = sizex0[0];
        sizexk[1] = sizex0[1];
        fixedPoint **xk = new fixedPoint*[sizexk[0]];
        for(int i = 0; i < sizexk[0]; i++){
                xk[i] = new fixedPoint[sizexk[1]];
        }
        for(int i = 0; i < sizexk[0]; i++){
                for(int j = 0; j < sizexk[1]; j++){
                        xk[i][j] = x0[i][j];
                }
        }
        int sizextemp[2];
        sizextemp[0] = sizex0[0];
        sizextemp[1] = sizex0[1];
        fixedPoint **xtemp = new fixedPoint*[sizextemp[0]];
        for(int i = 0; i < sizextemp[0]; i++){
                xtemp[i] = new fixedPoint[sizextemp[1]];
        }
        for(int i = 0; i < sizextemp[0]; i++){
                for(int j = 0; j < sizextemp[1]; j++){
                        xtemp[i][j] = xk[i][j];
                }
        }
        int sizexHattemp[2];
        sizexHattemp[0] = sizexHatk[0];
        sizexHattemp[1] = sizexHatk[1];
        fixedPoint **xHattemp = new fixedPoint*[sizexHattemp[0]];
        for(int i = 0; i < sizexHattemp[0]; i++){
                xHattemp[i] = new fixedPoint[sizexHattemp[1]];
        }
        for(int i = 0; i < sizexHattemp[0]; i++){
                for(int j = 0; j < sizexHattemp[1]; j++){
                        xHattemp[i][j] = xHatk[i][j];
                }
        }




//Start computations

//Setup computing constants
        //auto setupConStart = clock_start();

        //compute gamma3
        matrixMul(L,C,LC, sizeL, sizeC);
        matrixMul( LC, B, LCB, sizeLC, sizeB);
        for(int i = 0; i < sizegamma3[0]; i++){
                for(int j = 0; j < sizegamma3[1]; j++){
                        gamma3[i][j] = B[i][j] - LCB[i][j];
                }
        }

        //compute gamma2
        matrixMul(gamma3, K, gamma2, sizegamma3, sizeK);

        //compute gamma1
        matrixMul(LC, A, LCA, sizeLC, sizeA);
        for(int i = 0; i < sizegamma1[0]; i++){
                for(int j = 0; j < sizegamma1[0]; j++){
                        gamma1[i][j] = (A[i][j] - LCA[i][j]) - gamma2[i][j];
                }
        }
        //cout<<fixed<<setprecision(7)<<"Time for setup to compute constants: "<<time_from(setupConStart) *.000001 <<" seconds"<<endl;

//Cloud computing constants
        //auto cloudConStart = clock_start();

        //compute uTiler
        matrixMul(K, xr, Kxr, sizeK, sizexr);
        for(int i = 0; i < sizeuTilder[0]; i++){
                for(int j = 0; j < sizeuTilder[1]; j++){
                        uTilder[i][j] =  ur[i][j] + Kxr[i][j];
                }
        }

        //compute xGamma
        matrixMul(gamma2, xr, gamma2xr, sizegamma2, sizexr);
        matrixMul(gamma3, ur, gamma3ur, sizegamma3, sizeur);
        for(int i = 0; i < sizexGamma[0]; i++){
                for(int j = 0; j < sizexGamma[1]; j++){
                        xGamma[i][j] = gamma3ur[i][j] + gamma2xr[i][j];
                }
        }
        //cout<<"Time for cloud to compute constants: "<<time_from(cloudConStart)*.000001<<" senconds"<<endl;


//Start computations
        int k = 0;
                uK(K, x0, uTilder, uk, sizeK);
                cout<<"u"<<0<<":  "<<endl;
                for(int i = 0; i < sizeuk[0]; i++){
                        for(int j = 0; j < sizeuk[1]; j++){
                                cout<<uk[i][j].reveal<double>(BOB)<<", ";
                        }
                        cout<<endl;
                }
                cout<<endl<<endl;

                measureState(A, xtemp, B, uk, xk, sizeA, sizeB);

                //compute zk with no noise
                matrixMul(C, xk, zk, sizeC, sizexk);
                cout<<"z"<<k<<":  "<<endl;
                for(int i = 0; i < sizezk[0]; i++){
                        for(int j = 0; j < sizezk[1]; j++){
                                cout<<zk[i][j].reveal<double>(BOB)<<", ";
                        }
                        cout<<endl;
                }
                cout<<endl<<endl;

//long long time = 0;
//long long hold = 0;

        //Start loop
        for(k=1;k<5;k++){
                for(int i = 0; i < sizexHattemp[0]; i++){
                        for(int j = 0; j < sizexHattemp[1]; j++){
                                xHattemp[i][j] = xHatk[i][j];
                        }
                }
//auto startCloudOnline = clock_start();

                //Compute state estimate xHat
                xHat(gamma1, xHattemp, L, zk, xGamma, xHatk, sizegamma1);
                cout<<"xHat"<<k<<":  "<<endl;
                for(int i = 0; i < sizexHatk[0]; i++){
                        for(int j = 0; j < sizexHatk[1]; j++){
                                cout<<xHatk[i][j].reveal<double>(BOB)<<", ";
                        }
                        cout<<endl;
                }
                cout<<endl<<endl;

                //compute uk
                uK(K, xHatk, uTilder, uk, sizeK);
//hold=time_from(startCloudOnline);
//cout<<"Time for cloud computations: "<<hold<<endl;
//time+=hold;
                cout<<"u"<<k<<":  "<<endl;
                for(int i = 0; i < sizeuk[0]; i++){
                        for(int j = 0; j < sizeuk[1]; j++){
                                cout<<uk[i][j].reveal<double>(BOB)<<", ";
                        }
                        cout<<endl;
                }
                cout<<endl<<endl;

                for(int i = 0; i < sizextemp[0]; i++){
                        for(int j = 0; j < sizextemp[1]; j++){
                                xtemp[i][j] = xk[i][j];
                        }
                }

                //measure state ie. compute xk
                measureState(A, xtemp, B, uk, xk, sizeA, sizeB);

                //compute zk with no noise
                matrixMul(C, xk, zk, sizeC, sizexk);
                cout<<"z"<<k<<":  "<<endl;
                for(int i = 0; i < sizezk[0]; i++){
                        for(int j = 0; j < sizezk[1]; j++){
                                cout<<zk[i][j].reveal<double>(BOB)<<", ";
                        }
                        cout<<endl;
                }
                cout<<endl<<endl;

        }


//cout<<"Average time: "<<time/100<<endl;

        delete io;
        return 0;


}
