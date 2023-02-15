#ifndef __OUR_GL_H__
#define __OUR_GL_H__
#include "tgaimage.h"
#include "geometry.h"

extern Matrix ModelView;
extern Matrix Projection;
extern Matrix Viewport;
const float depth = 2000.f;

// 视口矩阵  双单元立方体 [-1，1]*[-1，1]*[-1，1] 映射到屏幕立方体 [x，x+w]*[y，y+h]*[0，d] 上
void viewport(int x, int y, int w, int h);
void viewportShadow(int x, int y, int w, int h);
void projection(float coeff=0.f); // coeff = -1/c
// eye 眼睛（摄像机） 返回view
void lookat(Vec3f eye, Vec3f center, Vec3f up);

void line(int x0, int y0, int x1, int y1, TGAImage &image, TGAColor color);

struct IShader {
    virtual ~IShader();
    // 顶点着色器的主要目标:转换顶点的坐标。次要目标:为片段着色器准备数据
    virtual Vec4f vertex(int iface, int nthvert) = 0;
    // 片段着色器的主要目标:确定当前像素的颜色。次要目标:返回 true 来丢弃当前像素
    virtual bool fragment(Vec3f gl_FragCoord, Vec3f bar, TGAColor &color) = 0;
};
// void triangle(Vec4f *pts, IShader &shader, TGAImage &image, TGAImage &zbuffer);
void triangle(mat<4,3,float> &pts, IShader &shader, TGAImage &image, float *zbuffer);

struct Shader2 {
    virtual ~Shader2();
    // 顶点着色器的主要目标:转换顶点的坐标。次要目标:为片段着色器准备数据
    virtual Vec4f vertex(int iface, int nthvert) = 0;
    // 片段着色器的主要目标:确定当前像素的颜色。次要目标:返回 true 来丢弃当前像素
    virtual bool fragment(Vec3f bar, TGAColor &color) = 0;
};
// void triangle(Vec4f *pts, IShader &shader, TGAImage &image, TGAImage &zbuffer);
void triangleShadow(Vec4f *pts, Shader2 &shader, TGAImage &image, float *zbuffer);

void triangleSao(Vec2i t0, Vec2i t1, Vec2i t2, TGAImage &image, TGAColor color);
#endif //__OUR_GL_H__

