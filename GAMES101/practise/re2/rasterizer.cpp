// clang-format off
//
// Created by goksu on 4/6/19.
//

#include <algorithm>
#include <vector>
#include "rasterizer.hpp"
#include <opencv2/opencv.hpp>
#include <math.h>


rst::pos_buf_id rst::rasterizer::load_positions(const std::vector<Eigen::Vector3f> &positions)
{
    auto id = get_next_id();
    pos_buf.emplace(id, positions);

    return {id};
}

rst::ind_buf_id rst::rasterizer::load_indices(const std::vector<Eigen::Vector3i> &indices)
{
    auto id = get_next_id();
    ind_buf.emplace(id, indices);

    return {id};
}

rst::col_buf_id rst::rasterizer::load_colors(const std::vector<Eigen::Vector3f> &cols)
{
    auto id = get_next_id();
    col_buf.emplace(id, cols);

    return {id};
}

auto to_vec4(const Eigen::Vector3f& v3, float w = 1.0f)
{
    return Vector4f(v3.x(), v3.y(), v3.z(), w);
}

//to fill
static bool insideTriangle(float x, float y, const Eigen::Vector3f* _v)
{      
    bool inside=false;
    Eigen::Vector3f v0=_v[0],v1=_v[1],v2=_v[2];
    Eigen::Vector3f point(x,y,0);
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
    // std::cout<<e1<<endl<<e2<<endl<<e3<<endl<<endl<<vec1<<endl<<vec2<<endl<<vec3<<endl;
    // std::cout<<res1<<endl<<res2<<endl<<res3<<endl;
    if(res1 * res2 > 0 && res2 * res3 > 0 && res3 * res1 > 0 ){
        inside=true;
    }
    else if(res1 * res2 < 0 && res2 * res3 < 0 && res3 * res1 < 0 ){
        inside=true;
    }
    // TODO : Implement this function to check if the point (x, y) is inside the triangle represented by _v[0], _v[1], _v[2]
    return inside;
}

static std::tuple<float, float, float> computeBarycentric2D(float x, float y, const Vector3f* v)
{
    float c1 = (x*(v[1].y() - v[2].y()) + (v[2].x() - v[1].x())*y + v[1].x()*v[2].y() - v[2].x()*v[1].y()) / (v[0].x()*(v[1].y() - v[2].y()) + (v[2].x() - v[1].x())*v[0].y() + v[1].x()*v[2].y() - v[2].x()*v[1].y());
    float c2 = (x*(v[2].y() - v[0].y()) + (v[0].x() - v[2].x())*y + v[2].x()*v[0].y() - v[0].x()*v[2].y()) / (v[1].x()*(v[2].y() - v[0].y()) + (v[0].x() - v[2].x())*v[1].y() + v[2].x()*v[0].y() - v[0].x()*v[2].y());
    float c3 = (x*(v[0].y() - v[1].y()) + (v[1].x() - v[0].x())*y + v[0].x()*v[1].y() - v[1].x()*v[0].y()) / (v[2].x()*(v[0].y() - v[1].y()) + (v[1].x() - v[0].x())*v[2].y() + v[0].x()*v[1].y() - v[1].x()*v[0].y());
    return {c1,c2,c3};
}

void rst::rasterizer::draw(pos_buf_id pos_buffer, ind_buf_id ind_buffer, col_buf_id col_buffer, Primitive type)
{
    auto& buf = pos_buf[pos_buffer.pos_id];
    auto& ind = ind_buf[ind_buffer.ind_id];
    auto& col = col_buf[col_buffer.col_id];

    float f1 = (50 - 0.1) / 2.0;
    float f2 = (50 + 0.1) / 2.0;

    Eigen::Matrix4f mvp = projection * view * model;
    for (auto& i : ind)
    {
        Triangle t;
        Eigen::Vector4f v[] = {
                mvp * to_vec4(buf[i[0]], 1.0f),
                mvp * to_vec4(buf[i[1]], 1.0f),
                mvp * to_vec4(buf[i[2]], 1.0f)
        };
        //Homogeneous division
        for (auto& vec : v) {
            vec /= vec.w();
        }
        //Viewport transformation
        for (auto & vert : v)
        {
            vert.x() = 0.5*width*(vert.x()+1.0);
            vert.y() = 0.5*height*(vert.y()+1.0);
            vert.z() = vert.z() * f1 + f2;
        }

        for (int i = 0; i < 3; ++i)
        {
            t.setVertex(i, v[i].head<3>());
            t.setVertex(i, v[i].head<3>());
            t.setVertex(i, v[i].head<3>());
        }

        auto col_x = col[i[0]];
        auto col_y = col[i[1]];
        auto col_z = col[i[2]];

        t.setColor(0, col_x[0], col_x[1], col_x[2]);
        t.setColor(1, col_y[0], col_y[1], col_y[2]);
        t.setColor(2, col_z[0], col_z[1], col_z[2]);

        rasterize_triangle(t);
    }
}
//to fill up
//Screen space rasterization
void rst::rasterizer::rasterize_triangle(const Triangle& t) {
    auto v = t.toVector4();
    
    //using INT_MIN as boingbox border
    float alpha,beta,gamma;
    float l_sc=int(std::min(std::min(v[0](0),v[1](0)),v[2](0))),
        r_sc=int(std::max(std::max(v[0](0),v[1](0)),v[2](0)
        ))+1,
        b_sc=int(std::min(std::min(v[0](1),v[1](1)),v[2](1))),
        t_sc=int(std::max(std::max(v[0](1),v[1](1)),v[2](1)))+1;
    // std::cout<<v[0][0]<<" "<<v[0][1]<<std::endl;
    // std::cout<<v[1][0]<<" "<<v[1][1]<<std::endl;
    // std::cout<<v[2][0]<<" "<<v[2][1]<<std::endl;
    
    // std::cout<<l_sc<<" "<<r_sc<<" "<<b_sc<<" "<<t_sc<<std::endl;


    for(float i = l_sc;i<r_sc;i++){
        for(float j=b_sc;j<t_sc;j++){
            if(insideTriangle(i+0.5,j+0.5,t.v)){
                //compute the interpolation result of z
                auto[alpha, beta, gamma] = computeBarycentric2D(i+0.5, j+0.5, t.v);
                float w_reciprocal = 1.0/(alpha / v[0].w() + beta / v[1].w() + gamma / v[2].w());
                float z_interpolated = alpha * v[0].z() / v[0].w() + beta * v[1].z() / v[1].w() + gamma * v[2].z() / v[2].w();
                z_interpolated *= w_reciprocal;

                //setcolor
                //遮挡判断
                int index=get_index(i,j);
                if(z_interpolated < depth_buf[index]){//如果当前z值比像素z值小（这里是把z值换成正数比较的）
                    // TODO : set the current pixel (use the set_pixel function) to the color of the triangle (use getColor function) if it should be painted.
                    Eigen::Vector3f pixel;
                    pixel<< i,j,z_interpolated;
                    set_pixel(pixel,t.getColor());
                    depth_buf[index] = z_interpolated;//设置像素颜色，修改像素当前深度   
                }
            }
        }
    }


    // TODO : Find out the bounding box of current triangle.
    // iterate through the pixel and find if the current pixel is inside the triangle

    // If so, use the following code to get the interpolated z value.
    //auto[alpha, beta, gamma] = computeBarycentric2D(x, y, t.v);
    //float w_reciprocal = 1.0/(alpha / v[0].w() + beta / v[1].w() + gamma / v[2].w());
    //float z_interpolated = alpha * v[0].z() / v[0].w() + beta * v[1].z() / v[1].w() + gamma * v[2].z() / v[2].w();
    //z_interpolated *= w_reciprocal;

    // TODO : set the current pixel (use the set_pixel function) to the color of the triangle (use getColor function) if it should be painted.
}

void rst::rasterizer::set_model(const Eigen::Matrix4f& m)
{
    model = m;
}

void rst::rasterizer::set_view(const Eigen::Matrix4f& v)
{
    view = v;
}

void rst::rasterizer::set_projection(const Eigen::Matrix4f& p)
{
    projection = p;
}

void rst::rasterizer::clear(rst::Buffers buff)
{
    if ((buff & rst::Buffers::Color) == rst::Buffers::Color)
    {
        std::fill(frame_buf.begin(), frame_buf.end(), Eigen::Vector3f{0, 0, 0});
    }
    if ((buff & rst::Buffers::Depth) == rst::Buffers::Depth)
    {
        std::fill(depth_buf.begin(), depth_buf.end(), std::numeric_limits<float>::infinity());
    }
}

rst::rasterizer::rasterizer(int w, int h) : width(w), height(h)
{
    frame_buf.resize(w * h);
    depth_buf.resize(w * h);
}

int rst::rasterizer::get_index(int x, int y)
{
    return (height-1-y)*width + x;
}

void rst::rasterizer::set_pixel(const Eigen::Vector3f& point, const Eigen::Vector3f& color)
{
    //old index: auto ind = point.y() + point.x() * width;
    auto ind = (height-1-point.y())*width + point.x();
    frame_buf[ind] = color;

}

// clang-format on