#include<iostream>
using namespace std;

#include<Eigen/Core>
#include<Eigen/Dense>
using namespace Eigen;

#define PI 3.14159265

int main(){
    cout<<acos(-1)<<endl;
    float x=2,y=1;
    Vector3f v(x,y,1);
    cout<<endl<<v<<endl;

    Matrix3f trans = Matrix3f::Identity();
    float axis=-45.0/180.0*PI,tranx=1,transy=2;
    float cosx=cos(axis),sinx=sin(axis);
    cout<<cosx<<endl<<sinx<<endl;
    trans.block<2,2>(0,0)<<cosx,-sinx,sinx,cosx;
    cout<<trans<<endl;

    // Vector3f res=Vector3f::Zero();
    // res=trans*v;
    // cout<<res<<endl;
}
