// clang-format off
#include <iostream>
#include <opencv2/opencv.hpp>
#include "rasterizer.hpp"
#include "global.hpp"
#include "Triangle.hpp"

constexpr double MY_PI = 3.1415926;

Eigen::Matrix4f get_view_matrix(Eigen::Vector3f eye_pos)
{
    Eigen::Matrix4f view = Eigen::Matrix4f::Identity();

    Eigen::Matrix4f translate;
    translate << 1,0,0,-eye_pos[0],
                 0,1,0,-eye_pos[1],
                 0,0,1,-eye_pos[2],
                 0,0,0,1;

    view = translate*view;

    return view;
}

Eigen::Matrix4f get_model_matrix(float rotation_angle)
{
    float cosx=cos(rotation_angle[0]/180.0*MY_PI),sinx=sin(rotation_angle[0]/180.0*MY_PI);
    Eigen::Matrix4f m_x = Eigen::Matrix4f::Identity();
    m_x.block<2,2>(1,1)<<cosx,-sinx,sinx,cosx;

    float cosy=cos(rotation_angle[1]/180.0*MY_PI),siny=sin(rotation_angle[1]/180.0*MY_PI);
    Eigen::Matrix4f m_y = Eigen::Matrix4f::Identity();
    m_y.block<3,3>(0,0)<<cosy,0,siny,0,1,0,-siny,0,cosy;

    float cosz=cos(rotation_angle[2]/180.0*MY_PI),sinz=sin(rotation_angle[2]/180.0*MY_PI);
    Eigen::Matrix4f m_z = Eigen::Matrix4f::Identity();
    m_z.block<2,2>(0,0)<<cosz,-sinz,sinz,cosz;

    Eigen::Matrix4f model = Eigen::Matrix4f::Identity();
    model = m_x * m_y * m_z;
    // TODO: Implement this function
    // Create the model matrix for rotating the triangle around the Z axis.
    // Then return it.
    std::cout<<"x"<<std::endl<<m_x<<std::endl;
    std::cout<<"y"<<std::endl<<m_y<<std::endl;
    std::cout<<"z"<<std::endl<<m_z<<std::endl;
    std::cout<<"model"<<std::endl<<model<<std::endl;
    return model;
}

Eigen::Matrix4f get_projection_matrix(float eye_fov, float aspect_ratio, float zNear, float zFar)
{
    //eye_fov  视角 aspect_ratio width/height
    // Students will implement this function
    //step:
    //1. using perspective projection to squish
    //2. using orthographic projection to transform into [-1,1]

    Eigen::Matrix4f projection = Eigen::Matrix4f::Identity();
    //1.
    Eigen::Matrix4f m_p_ersp_to_ortho = Eigen::Matrix4f::Identity();
    m_p_ersp_to_ortho(0,0)=zNear;
    m_p_ersp_to_ortho(1,1)=zNear;
    m_p_ersp_to_ortho(2,2)=zNear+zFar;
    m_p_ersp_to_ortho(2,3)=-zNear*zFar;
    m_p_ersp_to_ortho(3,2)=1;
    m_p_ersp_to_ortho(3,3)=0;


    //2.计算挤压后的l,r,b,t
    float half_fov = eye_fov/2.0*MY_PI/180.0; //radium
    float top = zNear * tan(half_fov);
    float bottum = -top;
    float right = top * aspect_ratio;
    float left = -right;
    
    Eigen::Matrix4f m_ortho_trans = Eigen::Matrix4f::Identity();
    m_ortho_trans(0,3) = -(right + left)/2.0;
    m_ortho_trans(1,3) = -(top + bottum)/2.0;
    m_ortho_trans(2,3) = -(zNear + zFar)/2.0;

    Eigen::Matrix4f m_ortho_scale = Eigen::Matrix4f::Identity();
    m_ortho_trans(0,0) = 2.0/(right - left);
    m_ortho_trans(1,1) = 2.0/(top - bottum);
    m_ortho_trans(2,2) = 2.0/(zNear - zFar);
    
    projection = m_ortho_scale * m_ortho_trans * m_p_ersp_to_ortho;


    // TODO: Implement this function
    // Create the projection matrix for the given parameters.
    // Then return it.

    return projection;
}

int main(int argc, const char** argv)
{
    float angle = 0;
    bool command_line = false;
    std::string filename = "output.png";

    if (argc == 2)
    {
        command_line = true;
        filename = std::string(argv[1]);
    }

    rst::rasterizer r(700, 700);

    Eigen::Vector3f eye_pos = {0,0,5};


    std::vector<Eigen::Vector3f> pos
            {
                    {2, 0, -2},
                    {0, 2, -2},
                    {-2, 0, -2},
                    {3.5, -1, -5},
                    {2.5, 1.5, -5},
                    {-1, 0.5, -5}
            };

    std::vector<Eigen::Vector3i> ind
            {
                    {0, 1, 2},
                    {3, 4, 5}
            };

    std::vector<Eigen::Vector3f> cols
            {
                    {217.0, 238.0, 185.0},
                    {217.0, 238.0, 185.0},
                    {217.0, 238.0, 185.0},
                    {185.0, 217.0, 238.0},
                    {185.0, 217.0, 238.0},
                    {185.0, 217.0, 238.0}
            };

    auto pos_id = r.load_positions(pos);
    auto ind_id = r.load_indices(ind);
    auto col_id = r.load_colors(cols);

    int key = 0;
    int frame_count = 0;

    if (command_line)
    {
        r.clear(rst::Buffers::Color | rst::Buffers::Depth);

        r.set_model(get_model_matrix(angle));
        r.set_view(get_view_matrix(eye_pos));
        r.set_projection(get_projection_matrix(45, 1, 0.1, 50));

        r.draw(pos_id, ind_id, col_id, rst::Primitive::Triangle);
        cv::Mat image(700, 700, CV_32FC3, r.frame_buffer().data());
        image.convertTo(image, CV_8UC3, 1.0f);
        cv::cvtColor(image, image, cv::COLOR_RGB2BGR);

        cv::imwrite(filename, image);

        return 0;
    }

    while(key != 27)
    {
        r.clear(rst::Buffers::Color | rst::Buffers::Depth);

        r.set_model(get_model_matrix(angle));
        r.set_view(get_view_matrix(eye_pos));
        r.set_projection(get_projection_matrix(45, 1, 0.1, 50));

        r.draw(pos_id, ind_id, col_id, rst::Primitive::Triangle);

        cv::Mat image(700, 700, CV_32FC3, r.frame_buffer().data());
        image.convertTo(image, CV_8UC3, 1.0f);
        cv::cvtColor(image, image, cv::COLOR_RGB2BGR);
        cv::imshow("image", image);
        key = cv::waitKey(10);

        std::cout << "frame count: " << frame_count++ << '\n';

        float flex_angle=10;
        switch(key){
            case 'q': 
            angle[0]+=flex_angle;
            break;
            case 'e':
            angle[0]-=flex_angle;
            break;
            case 'a': 
            angle[1]+=flex_angle;
            break;
            case 'd':
            angle[1]-=flex_angle;
            break;
            case 'z': 
            angle[2]+=flex_angle;
            break;
            case 'c':
            angle[2]-=flex_angle;
            break;
            
        }
    }

    return 0;
}
// clang-format on