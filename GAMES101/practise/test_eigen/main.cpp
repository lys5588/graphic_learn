#include<iostream>
using namespace std;

#include<Eigen/Core>
#include<Eigen/Dense>
using namespace Eigen;

#define PI 3.14159265

static bool insideTriangle(int x, int y, const Eigen::Vector3f* _v)
{      
    bool inside=false;
    Eigen::Vector3f v0=_v[0],v1=_v[1],v2=_v[2];
    Eigen::Vector3f point(float(x),float(y),0);
    float res1 = point.transpose() * v0,
        res2 = point.transpose() * v1, 
        res3 = point.transpose() * v2;
    std::cout<<point<<endl<<v0<<endl<<v1<<endl<<v2<<endl;
    std::cout<<res1<<endl<<res2<<endl<<res3<<endl;
    if(res1 * res2 > 0 && res2 * res3 > 0 && res3 * res1 > 0 ){
        inside=true;
    }
    else if(res1 * res2 < 0 && res2 * res3 < 0 && res3 * res1 < 0 ){
        inside=true;
    }
    // TODO : Implement this function to check if the point (x, y) is inside the triangle represented by _v[0], _v[1], _v[2]
    return inside;
}

int main(){
    // float rotation_angle=45.0;
    // Eigen::Matrix4f model = Eigen::Matrix4f::Identity();

    // // TODO: Implement this function
    // // Create the model matrix for rotating the triangle around the Z axis.
    // // Then return it.
    // Eigen::Matrix4f translate;
    // translate << cos(rotation_angle/180.0*acos(-1)), -sin(rotation_angle/180.0*acos(-1)), 0.0, 0.0, sin(rotation_angle/180.0*acos(-1)),
    //     cos(rotation_angle/180.0*acos(-1)), 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0;

    // model = translate * model;
    // cout<<"model"<<model<<endl;

    // cout<<acos(-1)<<endl;
    // float x=2,y=1;
    // Vector3f v(x,y,1);
    // cout<<endl<<v<<endl;

    // Matrix3f trans = Matrix3f::Identity();
    // float axis=-45.0/180.0*PI,tranx=1,transy=2;
    // float cosx=cos(axis),sinx=sin(axis);
    // cout<<cosx<<endl<<sinx<<endl;
    // trans.block<2,2>(0,0)<<cosx,-sinx,sinx,cosx;
    // cout<<trans<<endl;

    // Vector3f res=Vector3f::Zero();
    // res=trans*v;
    // cout<<res<<endl;


    int x=2 ,y=4;
    Eigen::Vector3f p0(3,1,0),p1(4,3,2),p2(5,1,9);
    Eigen::Vector3f _v[3]={p0,p1,p2};
    std::cout<<insideTriangle(x,y,_v);
}
