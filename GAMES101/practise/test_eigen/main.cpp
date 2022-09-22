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
    //triangle edge and point-to-point vector
    Eigen::Vector3f e1,e2,e3,vec1,vec2,vec3;

    e1<<(v1(0)-v0(0)),(v1(1)-v0(1)),0;
    e2<<(v2(0)-v1(0)),(v2(1)-v1(1)),0;
    e3<<(v0(0)-v2(0)),(v0(1)-v2(1)),0;
    
    vec1<<(point(0)-v0(0)),(point(1)-v0(1)),0;
    vec1<<(point(0)-v1(0)),(point(1)-v1(1)),0;
    vec1<<(point(0)-v2(0)),(point(1)-v2(1)),0;

    float res1 = e1.transpose() * vec1,
        res2 = e2.transpose() * vec2, 
        res3 = e3.transpose() * vec3;
    std::cout<<e1<<endl<<e2<<endl<<e3<<endl<<endl<<vec1<<endl<<vec2<<endl<<vec3<<endl;
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
