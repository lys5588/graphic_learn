#include<iostream>
using namespace std;

#include<Eigen/Core>
#include<Eigen/Dense>
using namespace Eigen;

#define PI 3.14159265

static bool insideTriangle_my(float x, float y, const Eigen::Vector3f* _v)
{      
    bool inside=false;
    Eigen::Vector3f v0=_v[0],v1=_v[1],v2=_v[2];
    Eigen::Vector3f point(x,y,0);
    //triangle edge and point-to-point vector
    Eigen::Vector3f e1,e2,e3,vec1,vec2,vec3;

    e1=v1-v0;
    e2=v2-v1;
    e3=v0-v2;
    
    vec1=point-v0;
    vec1=point-v1;
    vec1=point-v2;

    float res1 = e1.cross(vec1).z(),
        res2 = e2.cross(vec1).z(), 
        res3 = e3.cross(vec1).z();
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

static bool insideTriangle(float x, float y, const Vector3f* _v)
{   
    // 0 means negative, 1 means positive
    int flag = -1;

    for(int i = 0; i < 3; i++) {
        // the current point
        Eigen::Vector3f p0 = {x, y, 0};
        // the 1st vertex
        Eigen::Vector3f p1 = _v[i];
        // the 2nd vertex
        Eigen::Vector3f p2 = _v[(i+1)%3];
        
        // the 1st vector (p1-p0)
        Eigen::Vector3f v1 = p0-p1;
        // the 2nd vector (p1-p2)
        Eigen::Vector3f v2 = p2-p1;

        // get the cross product
        float cp = v1.cross(v2).z();
        if(cp == 0) continue;

        cout<<v1<<" "<<v2<<" "<<cp<<endl;
        int sign = cp < 0 ? 0: 1;
        if(flag == -1) flag = sign;
        if(flag != sign) return false;
    }

    return true;
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


    int x=433 ,y=307;
    Eigen::Vector3f p0(542,434,0),p1(138,223,2),p2(434,307,9);
    Eigen::Vector3f _v[3]={p0,p1,p2};
    std::cout<<insideTriangle(x,y,_v);
}
