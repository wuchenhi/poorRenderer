#ifndef __OUR_GL_H__
#define __OUR_GL_H__
#include "tgaimage.h"
#include "geometry.h"

extern Matrix ModelView;
extern Matrix Projection;
extern Matrix Viewport;
const float depth = 2000.f;

// �ӿھ���  ˫��Ԫ������ [-1��1]*[-1��1]*[-1��1] ӳ�䵽��Ļ������ [x��x+w]*[y��y+h]*[0��d] ��
void viewport(int x, int y, int w, int h);
void viewportShadow(int x, int y, int w, int h);
void projection(float coeff=0.f); // coeff = -1/c
// eye �۾���������� ����view
void lookat(Vec3f eye, Vec3f center, Vec3f up);

struct IShader {
    virtual ~IShader();
    // ������ɫ������ҪĿ��:ת����������ꡣ��ҪĿ��:ΪƬ����ɫ��׼������
    virtual Vec4f vertex(int iface, int nthvert) = 0;
    // Ƭ����ɫ������ҪĿ��:ȷ����ǰ���ص���ɫ����ҪĿ��:���� true ��������ǰ����
    virtual bool fragment(Vec3f gl_FragCoord, Vec3f bar, TGAColor &color) = 0;
};
// void triangle(Vec4f *pts, IShader &shader, TGAImage &image, TGAImage &zbuffer);
void triangle(mat<4,3,float> &pts, IShader &shader, TGAImage &image, float *zbuffer);

struct Shader2 {
    virtual ~Shader2();
    // ������ɫ������ҪĿ��:ת����������ꡣ��ҪĿ��:ΪƬ����ɫ��׼������
    virtual Vec4f vertex(int iface, int nthvert) = 0;
    // Ƭ����ɫ������ҪĿ��:ȷ����ǰ���ص���ɫ����ҪĿ��:���� true ��������ǰ����
    virtual bool fragment(Vec3f bar, TGAColor &color) = 0;
};
// void triangle(Vec4f *pts, IShader &shader, TGAImage &image, TGAImage &zbuffer);
void triangleShadow(Vec4f *pts, Shader2 &shader, TGAImage &image, float *zbuffer);
#endif //__OUR_GL_H__

