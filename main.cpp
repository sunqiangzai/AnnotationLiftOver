#include <iostream>
#include <exception>
#include <string>
#include "trycatchfinally.h"
using namespace std;

int main()
{
    // No exception
    TryFinally(
        []{ cout<<"Try Body"<<endl; },
        [](exception){ cout<<"Catch Body"<<endl; },
        []{ cout<<"Finally Body"<<endl; }
    );
    
    // Exception thrown
    TryFinally(
        []{ cout<<"Try Body"<<endl; throw exception(); },
        [](exception){ cout<<"Catch Body"<<endl; },
        []{ cout<<"Finally Body"<<endl; }
    );
    
    // You can even change values in this scope (notice the [&]?)
    string what;
    TryFinally(
        []{ throw exception(); },
        [&](const exception &ex){ what = ex.what(); },
        []{}
    );
    cout<<what<<endl;
    
    // What about custom exception types?
    class test_exception_t : public std::exception{ public:
        virtual const char *what() const noexcept{ return "My Test Exception"; }
    };
    
    TryFinally(
        []{ throw test_exception_t(); },
        [](const exception &ex){ cout<<ex.what()<<endl; },
        []{}
    );
    
    // Can you tell the difference between two types of exceptions?
    class second_exception_t : public std::exception{ public:
        virtual const char *what() const noexcept{ return "My Second Test Exception"; }
        string name;
    };
    
    TryFinally(
        []{ throw second_exception_t(); },
        [](const exception &ex){ 
            if(NULL == dynamic_cast<test_exception_t const *>(&ex))
                cout<<"This is not the first kind of exception"<<endl;
            if(NULL != dynamic_cast<second_exception_t const *>(&ex))
                cout<<"This is the second kind of exception"<<endl;
        },
        []{}
    );
    
    // What about exception types that aren't derived from std::exception?
    TryFinally<const char *>(
        []{ throw "This is not even an exception!"; },
        [](const char *ex){ cout<<ex<<endl; },
        []{}
    );
    
    // If you fail to catch the exception, the finally should still execute
    try{
        TryFinally(
            []{ throw "This is not even an exception!"; },
            [](exception){ cout<<"This should not execute because we didn't catch the right type"<<endl; },
            []{ cout<<"Finally body still executes"<<endl; }
        );
    } catch(...){
        cout<<"Exception escaped"<<endl;
    }
    
    // Let's do something more complex. I will modify and rethrow the exception
    try{
        TryFinally(
            []{ throw second_exception_t(); },
            [](exception &ex){
                // Let's give the exception a name and rethrow it
                second_exception_t *s = dynamic_cast<second_exception_t *>(&ex);
                s->name = "George";
                
                // I don't even need to tell it what to throw
                throw;
            },
            []{ cout<<"Finally Body"<<endl; }
        );
    }
    catch(const second_exception_t &ex){
        cout<<"Exception named "<<ex.name<<endl;
    }
    
    return 0;
}
