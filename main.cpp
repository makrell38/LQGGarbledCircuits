#include "subsystem.h"
#include "Setup.h"
#include "Cloud.h"
#include "Actuator.h"



int main(int argc, char** argv){
        int port, party;
        parse_party_and_port(argv, &party, &port);
        NetIO * io = new NetIO(party==ALICE ? nullptr: "127.0.0.1", port);
        setup_semi_honest(io, party);
        bool print = 0;


        Actuator *actuator = new Actuator();
        subSystem *subsystem = new subSystem();
        Setup *setup = new Setup();
        Cloud *cloud = new Cloud();
        subsystem->inputData();
        setup->inputData();
auto startSetupConst = clock_start();
        setup->computeConstants();
cout<<time_from(startSetupConst)<<endl;

auto startCloudConst = clock_start();
        cloud->getInputs(setup->L, setup->sizeL, setup->K, setup->sizeK, setup->gamma3, setup->sizegamma3, setup->gamma2, setup->sizegamma2, setup->gamma1, setup->sizegamma1, subsystem->xr, subsystem->sizexr, subsystem->ur, subsystem->sizeur, subsystem->x0, subsystem->sizex0);
        cloud->computeConstants();
cout<<time_from(startCloudConst)<<endl;

        int k = 0;
auto start = clock_start();
        cloud->computeuk();
        subsystem->measureState(cloud->uk);
        subsystem->computezk();

auto startAct = clock_start();
        for(int i = 0; i <

        if(print){
                cout<<"u"<<0<<":  "<<endl;
                for(int i = 0; i < cloud->sizeuk[0]; i++){
                        for(int j = 0; j < cloud->sizeuk[1]; j++){
                                cout<<fixed<<setprecision(5)<<cloud->uk[i][j].reveal<double>(BOB)<<", ";
                        }
                        cout<<endl;
                }
                cout<<endl<<endl;
                cout<<"z"<<k<<":  "<<endl;
                for(int i = 0; i < subsystem->sizezk[0]; i++){
                        for(int j = 0; j < subsystem->sizezk[1]; j++){
                                cout<<subsystem->zk[i][j].reveal<double>(BOB)<<", ";
                        }
                        cout<<endl;
                }
                cout<<endl<<endl;
        }
//long long holdCloud = 0;
//long long timeCloud = 0;
//long long holdSub = 0;
//long long timeSub = 0;
        for(k = 1; k < 5; k++){
//auto startCloudComp = clock_start();
                cloud->computexHat(subsystem->zk);
                cloud->computeuk();
//holdCloud = time_from(startCloudComp);
//timeCloud+= holdCloud;

//auto startSubComp = clock_start();
                subsystem->measureState(cloud->uk);
                subsystem->computezk();
//holdSub = time_from(startSubComp);
//timeSub+= holdSub;
                if(print){
                        cout<<"xHat"<<k<<":  "<<endl;
                        for(int i = 0; i < cloud->sizexHatk[0]; i++){
                                for(int j = 0; j < cloud->sizexHatk[1]; j++){
                                        cout<<cloud->xHatk[i][j].reveal<double>(BOB)<<", ";
                                }
                                cout<<endl;
                        }
                        cout<<endl<<endl;
                        cout<<"u"<<k<<":  "<<endl;
                        for(int i = 0; i < cloud->sizeuk[0]; i++){
                                for(int j = 0; j < cloud->sizeuk[1]; j++){
                                        cout<<cloud->uk[i][j].reveal<double>(BOB)<<", ";
                                }
                                cout<<endl;
                        }
                        cout<<endl<<endl;
                        cout<<"z"<<k<<":  "<<endl;
                        for(int i = 0; i < subsystem->sizezk[0]; i++){
                                for(int j = 0; j < subsystem->sizezk[1]; j++){
                                        cout<<subsystem->zk[i][j].reveal<double>(BOB)<<", ";
                                }
                                cout<<endl;
                        }
                        cout<<endl<<endl;

                }

        }
cout<<time_from(start)<<endl;
//cout<<"Average Cloud comp time: "<<timeCloud/100<<endl;
//cout<<"Average Sub comp time: "<<timeSub/100<<endl;

        delete io;
        return 0;


}
