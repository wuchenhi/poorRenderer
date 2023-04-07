#include <vector>
#include <cstdlib>
#include <limits>
#include <iostream>
#include "tgaimage.h"
#include "model.h"
#include "geometry.h"
#include "our_gl.h"

Model *model = NULL;
Model *modelShadow = NULL;
float *shadowbuffer = NULL; //用于shadow map

const int width  = 800;
const int height = 800;

Vec3f light_dir(1,1,0);     //用于shadow map
Vec3f   eye(1.2,-.8,3);
Vec3f    center(0,0,0);
Vec3f        up(0,1,0);

struct ZShader : public IShader {
    mat<4,3,float> varying_tri;

    virtual Vec4f vertex(int iface, int nthvert) {
        Vec4f gl_Vertex = Projection*ModelView*embed<4>(model->vert(iface, nthvert));
        varying_tri.set_col(nthvert, gl_Vertex);
        return gl_Vertex;
    }

    virtual bool fragment(Vec3f gl_FragCoord, Vec3f bar, TGAColor &color) {
        color = TGAColor(0, 0, 0);
        return false;
    }
};

// 最大斜率方向
float max_elevation_angle(float *zbuffer, Vec2f p, Vec2f dir) {
    float maxangle = 0;
    for (float t=0.; t<1000.; t+=1.) {
        Vec2f cur = p + dir*t;
        if (cur.x>=width || cur.y>=height || cur.x<0 || cur.y<0) return maxangle;

        float distance = (p-cur).norm();
        if (distance < 1.f) continue;
        float elevation = zbuffer[int(cur.x)+int(cur.y)*width]-zbuffer[int(p.x)+int(p.y)*width];
        maxangle = std::max(maxangle, atanf(elevation/distance));
    }
    return maxangle;
}

/******************************************阴影渲染*******************************************/

// z 缓冲区复制到帧缓冲区
struct DepthShader : public Shader2 {
    mat<3,3,float> varying_tri;

    DepthShader() : varying_tri() {}

    virtual Vec4f vertex(int iface, int nthvert) {
        Vec4f gl_Vertex = embed<4>(model->vert(iface, nthvert));      // read the vertex from .obj file
        gl_Vertex = Viewport*Projection*ModelView*gl_Vertex;          // transform it to screen coordinates
        varying_tri.set_col(nthvert, proj<3>(gl_Vertex/gl_Vertex[3]));
        return gl_Vertex;
    }

    virtual bool fragment(Vec3f bar, TGAColor &color) {
        Vec3f p = varying_tri*bar;
        color = TGAColor(255, 255, 255)*(p.z/depth);
        return false;
    }
};

struct Shader : public Shader2 {
    mat<4,4,float> uniform_M;   //  Projection*ModelView
    mat<4,4,float> uniform_MIT; // (Projection*ModelView).invert_transpose()转置+取逆
    mat<4,4,float> uniform_Mshadow; // 将帧缓冲区的屏幕坐标转换为阴影缓冲区的屏幕坐标
    mat<2,3,float> varying_uv;  // 三角形uv坐标, written by the vertex shader, read by the fragment shader
    mat<3,3,float> varying_tri; // 视口变换前的三角坐标，由VS写入，由FS读取

    Shader(Matrix M, Matrix MIT, Matrix MS) : uniform_M(M), uniform_MIT(MIT), uniform_Mshadow(MS), varying_uv(), varying_tri() {}
    
    virtual Vec4f vertex(int iface, int nthvert) {
        varying_uv.set_col(nthvert, model->uv(iface, nthvert));
        Vec4f gl_Vertex = Viewport*Projection*ModelView*embed<4>(model->vert(iface, nthvert));
        varying_tri.set_col(nthvert, proj<3>(gl_Vertex/gl_Vertex[3]));
        return gl_Vertex;
    }

    virtual bool fragment(Vec3f bar, TGAColor &color) {
        Vec4f sb_p = uniform_Mshadow*embed<4>(varying_tri*bar);   // 阴影缓冲区中的对应点 补1为齐次矩阵
        sb_p = sb_p/sb_p[3];
        int idx = int(sb_p[0]) + int(sb_p[1])*width;              // 阴影缓冲器阵列中的索引
        float shadow = .3+.7*(shadowbuffer[idx]<sb_p[2] + 43.34); // 现在要确定当前像素是否被点亮,将z坐标与存储在阴影缓冲区中的值进行比较
        
        Vec2f uv = varying_uv*bar;
        Vec3f n = proj<3>(uniform_MIT*embed<4>(model->normal(uv))).normalize(); // 法线
        Vec3f l = proj<3>(uniform_M  *embed<4>(light_dir        )).normalize(); // 光照方向
        Vec3f r = (n*(n*l*2.f) - l).normalize();   // reflected light  r = 2 n<n，l> - l
        float spec = pow(std::max(r.z, 0.0f), model->specular(uv)); // 镜面
        float diff = std::max(0.f, n*l); //漫反射
        TGAColor c = model->diffuse(uv);
        //环境分量取了 20，为漫反射分量取了 1.2 分，为镜面反射分量取了 .6 分
        for (int i=0; i<3; i++) color[i] = std::min<float>(20 + c[i]*shadow*(1.2*diff + .6*spec), 255);
        return false;
    }
};
/******************************************阴影渲染*******************************************/

struct GouraudShader : public Shader2 {
    Vec3f varying_intensity; // 由顶点着色器写入，由片段着色器读取
    //  iface  第几个三角形  nthvert 三角形的第几个点
    virtual Vec4f vertex(int iface, int nthvert) {
        Vec4f gl_Vertex = embed<4>(model->vert(iface, nthvert)); // 从.obj文件中读取顶点
        gl_Vertex = Viewport*Projection*ModelView*gl_Vertex;     // 将其转换为屏幕坐标
        varying_intensity[nthvert] = std::max(0.f, model->normal(iface, nthvert)*light_dir); // 获得漫反射照明强度
        return gl_Vertex;
    }
    // 接收重心坐标，用于varying_数据的插值
    virtual bool fragment(Vec3f bar, TGAColor &color) {
        float intensity = varying_intensity*bar;   // 对当前像素进行强度插值 简单地计算为两个向量之间的点积
        // if (intensity>.85) intensity = 1;
        // else if (intensity>.60) intensity = .80;
        // else if (intensity>.45) intensity = .60;
        // else if (intensity>.30) intensity = .45;
        // else if (intensity>.15) intensity = .30;
        // else intensity = 0;
        // color = TGAColor(255, 155, 0)*intensity;
        color = TGAColor(255, 255, 255)*intensity; 
        return false;                            
    }
};

// 纹理渲染
struct ShaderWenli : public Shader2 {
    Vec3f          varying_intensity; 
    mat<2,3,float> varying_uv;        

    virtual Vec4f vertex(int iface, int nthvert) {
        varying_uv.set_col(nthvert, model->uv(iface, nthvert));//?
        varying_intensity[nthvert] = std::max(0.f, model->normal(iface, nthvert)*light_dir); // get diffuse lighting intensity
        Vec4f gl_Vertex = embed<4>(model->vert(iface, nthvert)); // read the vertex from .obj file
        return Viewport*Projection*ModelView*gl_Vertex; // transform it to screen coordinates
    }
    
    virtual bool fragment(Vec3f bar, TGAColor &color) {
        float intensity = varying_intensity*bar;   // interpolate intensity for the current pixel
        Vec2f uv = varying_uv*bar;                 // interpolate uv for the current pixel
        color = model->diffuse(uv)*intensity;      
        return false;                             
    }
};

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << "obj/model.obj" << std::endl;
        return 1;
    }
    float *zbuffer = new float[width*height];
    for (int i=width*height; i--; zbuffer[i] = -std::numeric_limits<float>::max());
    model = new Model(argv[1]);

    TGAImage frame(width, height, TGAImage::RGB);
    lookat(eye, center, up);
    viewport(width/8, height/8, width*3/4, height*3/4);
    projection(-1.f/(eye-center).norm());

    // 环境光遮蔽效果
    ZShader zshader;
    for (int i=0; i<model->nfaces(); i++) {
        for (int j=0; j<3; j++) {
            zshader.vertex(i, j);
        }
        triangle(zshader.varying_tri, zshader, frame, zbuffer);
    }

    for (int x=0; x<width; x++) {
        for (int y=0; y<height; y++) {
            if (zbuffer[x+y*width] < -1e5) continue;
            float total = 0;
            //立体角近似为（90°-max_elevation_angle） / 8 的总和
            for (float a=0; a<M_PI*2-1e-4; a += M_PI/4) {
                total += M_PI/2 - max_elevation_angle(zbuffer, Vec2f(x, y), Vec2f(cos(a), sin(a)));
            }
            total /= (M_PI/2)*8;
            total = pow(total, 100.f);//增加对比度
            frame.set(x, y, TGAColor(total*255, total*255, total*255));
        }
    }

    frame.flip_vertically();
    frame.write_tga_file("framebuffer.tga");

    // 阴影渲染
    float *zbufferShadow = new float[width*height];
    shadowbuffer   = new float[width*height];
    for (int i=width*height; --i; ) {
        zbufferShadow[i] = shadowbuffer[i] = -std::numeric_limits<float>::max();
    }   
    modelShadow = new Model(argv[1]);
    light_dir.normalize();
    { // rendering the shadow buffer
        TGAImage depth(width, height, TGAImage::RGB);
        //我们将渲染将相机放置在光源位置的图像。它将允许确定哪些部分被点亮，哪些部分被隐藏在光线之外
        lookat(light_dir, center, up);
        //viewport(width/8, height/8, width*3/4, height*3/4);
        viewportShadow(width/8, height/8, width*3/4, height*3/4);

        projection(0);

        DepthShader depthshader;
        Vec4f screen_coords[3];
        for (int i=0; i<model->nfaces(); i++) {
            for (int j=0; j<3; j++) {
                screen_coords[j] = depthshader.vertex(i, j);
            }
            triangleShadow(screen_coords, depthshader, depth, shadowbuffer);
        }
        depth.flip_vertically(); // 将原点放在图像的左下角
        depth.write_tga_file("depth.tga");
    }
    //保留了对象到屏幕的转换矩阵
    Matrix M = Viewport*Projection*ModelView;

    { // rendering the frame buffer
        TGAImage frameShadow(width, height, TGAImage::RGB);
        lookat(eye, center, up);
        //viewport(width/8, height/8, width*3/4, height*3/4);
        viewportShadow(width/8, height/8, width*3/4, height*3/4);

        projection(-1.f/(eye-center).norm());

        Shader shader(ModelView, (Projection*ModelView).invert_transpose(), M*(Viewport*Projection*ModelView).invert());
        Vec4f screen_coords[3];
        // 迭代模型的所有三角形并光栅化每个三角形
        // 外循环遍历所有三角形
        // 内部循环循环遍历当前三角形的所有顶点，并为每个顶点调用顶点着色器
        for (int i=0; i<model->nfaces(); i++) {
            for (int j=0; j<3; j++) {
                screen_coords[j] = shader.vertex(i, j);
            }
            triangleShadow(screen_coords, shader, frameShadow, zbufferShadow);
        }
        frameShadow.flip_vertically(); 
        frameShadow.write_tga_file("shadow.tga");
    }
    delete [] zbuffer;
    delete [] zbufferShadow;
    delete model;
    delete modelShadow;
    return 0;
}