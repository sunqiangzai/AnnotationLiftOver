/*
 * =====================================================================================
 *
 *       Filename:  testMultipleThreads.cpp
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  06/21/2017 18:17:56
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  YOUR NAME (), 
 *   Organization:  
 *
 * =====================================================================================
 */
#include <stdlib.h>
/*************************************************************************




 ************************************************************************/

#include <iostream>
#include <pthread.h>
#include "model.h"
#include <unistd.h>

pthread_mutex_t mutex;

struct myPara{
    int i;
    MyThreadCount* myThreadCount;
    myPara(int _i, MyThreadCount& _myThreadCount){
        i=_i;
        myThreadCount=&_myThreadCount;
    }
};

void *doSomething( void *arg){
    struct myPara *pstru = (struct myPara *) arg;
    for( double j=-9999999; j<99999999; j+=0.1){

    }

    pthread_mutex_lock(&mutex);
    (*pstru->myThreadCount).countDown();
    pthread_mutex_unlock(&mutex);
    std::cout << pstru->i << " hh " << (*pstru->myThreadCount).getCount() << std::endl;
    return NULL;
}
int main(){
    int maxThread = 4;
    MyThreadCount myThreadCount;
    pthread_attr_t attr;
    pthread_attr_init( &attr );
    pthread_attr_setscope( &attr, PTHREAD_SCOPE_SYSTEM );

    for(int i=-9999999; i<99999999; i+=1){
        std::cout << "57" << std::endl;
        bool isThisThreadUnrun=true;
        while(isThisThreadUnrun){
            if(myThreadCount.getCount() < maxThread){
                std::cout << 62 << std::endl;
                pthread_t thread;
                struct myPara para(i, myThreadCount);
                pthread_mutex_lock(&mutex);
                myThreadCount.plusOne();
                pthread_mutex_unlock(&mutex);
                pthread_create(&thread, &attr, doSomething, &para);
                //pthread_join(&thread, NULL);
                std::cout << "70" << std::endl;
                isThisThreadUnrun=false;
                break;
            }else{
                //std::cout << "77" << std::endl;
                usleep(10);
                //std::cout << "79" << std::endl;
            }
        }
    }

    while(myThreadCount.hasNext()){// wait for all the thread
        usleep(50);
    }
    pthread_mutex_destroy(&mutex);


    return 0;
}
