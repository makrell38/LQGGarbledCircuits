#ifndef FIXEDPOINT_H
#define FIXEDPOINT_H

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
                bool_dat(b, 24, (long long)round(fractPart * pow(2,24)), 0);
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
                Bit over = zero;
                Bit s;
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

                s = a[47]^b[47];
                for(int i = 72; i < 96; i++){
                        over = over | (sum[i]^sum[71]);
                }
               for(int i = 48; i < 71; i++){
                        ret[i-24] = (over&(!s)) | ((!over)&sum[i]);
                        ret[i-48] = sum[i-24];
                }
                ret[47] = (over&s) | ((!over)&sum[71]);
                ret[23] = sum[47];


//              memcpy(ret, sum+24, sizeof(Bit)*48);
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

#endif
