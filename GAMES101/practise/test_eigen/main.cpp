#include<iostream>
using namespace std;

#include<Eigen/Core>
#include<Eigen/Dense>
using namespace Eigen;

#define PI 3.14159265

static bool insideTriangle_my(float x, float y, const Eigen::Vector3f* _v)
{      
    bool inside=-1;
    for(int i=0;i<3;i++){
        Eigen::Vector3f p(x,y,0);
        Eigen::Vector3f p1=_v[i];
        Eigen::Vector3f p2=_v[(i+1)%3];

        Eigen::Vector3f v1,v2;
        v1=p-p1;
        v2=p-p2;

        float res=v1.cross(v2).z();
        int pos=res>0 ? 1:0;
        if(inside==-1) inside=pos;
        if(pos!=inside)return false;


    }

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
        cout<<p0<<" "<<p1<<" "<<p2<<" "<<cp<<endl;
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


    // int x=433 ,y=307;
    // Eigen::Vector3f p0(542,434,0),p1(138,223,2),p2(434,307,9);
    // Eigen::Vector3f _v[3]={p0,p1,p2};
    // std::cout<<insideTriangle_my(x,y,_v);

    Eigen::Vector3f ka = Eigen::Vector3f(0.005, 0.005, 0.005);
    Eigen::Vector3f amb_light_intensity{10, 10, 10};

    Eigen::Vector3f result=ka*amb_light_intensity;
    cout<<result;
}
