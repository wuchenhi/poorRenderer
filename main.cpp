#include <vector>
#include <cstdlib>
#include <limits>
#include <iostream>
#include "tgaimage.h"
#include "model.h"
#include "geometry.h"
#include "our_gl.h"

Model *model = NULL;

const int width  = 800;
const int height = 800;

Vec3f       eye(1.2,-.8,3);
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

// ���б�ʷ���
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

// struct Shader : public IShader {
//     mat<4,4,float> uniform_M;   //  Projection*ModelView
//     mat<4,4,float> uniform_MIT; // (Projection*ModelView).invert_transpose()ת��+ȡ��
//     mat<4,4,float> uniform_Mshadow; // ��֡����������Ļ����ת��Ϊ��Ӱ����������Ļ����
//     mat<2,3,float> varying_uv;  // ������uv����, written by the vertex shader, read by the fragment shader
//     mat<3,3,float> varying_tri; // �ӿڱ任ǰ���������꣬��VSд�룬��FS��ȡ

//     Shader(Matrix M, Matrix MIT, Matrix MS) : uniform_M(M), uniform_MIT(MIT), uniform_Mshadow(MS), varying_uv(), varying_tri() {}
    
//     virtual Vec4f vertex(int iface, int nthvert) {
//         varying_uv.set_col(nthvert, model->uv(iface, nthvert));
//         Vec4f gl_Vertex = Viewport*Projection*ModelView*embed<4>(model->vert(iface, nthvert));
//         varying_tri.set_col(nthvert, proj<3>(gl_Vertex/gl_Vertex[3]));
//         return gl_Vertex;
//     }

//     virtual bool fragment(Vec3f bar, TGAColor &color) {
//         Vec4f sb_p = uniform_Mshadow*embed<4>(varying_tri*bar); // ��Ӱ�������еĶ�Ӧ�� ��1Ϊ��ξ���
//         sb_p = sb_p/sb_p[3];
//         int idx = int(sb_p[0]) + int(sb_p[1])*width; // ��Ӱ�����������е�����
//         float shadow = .3+.7*(shadowbuffer[idx]<sb_p[2] + 43.34); //����Ҫȷ����ǰ�����Ƿ񱻵���,��z������洢����Ӱ�������е�ֵ���бȽ�
        
//         Vec2f uv = varying_uv*bar;
//         Vec3f n = proj<3>(uniform_MIT*embed<4>(model->normal(uv))).normalize(); // ����
//         Vec3f l = proj<3>(uniform_M  *embed<4>(light_dir        )).normalize(); // ���շ���
//         Vec3f r = (n*(n*l*2.f) - l).normalize();   // reflected light  r = 2 n<n��l> - l
//         float spec = pow(std::max(r.z, 0.0f), model->specular(uv)); // ����
//         float diff = std::max(0.f, n*l); //������
//         TGAColor c = model->diffuse(uv);
//         //��������ȡ�� 20��Ϊ���������ȡ�� 1.2 �֣�Ϊ���淴�����ȡ�� .6 ��
//         for (int i=0; i<3; i++) color[i] = std::min<float>(20 + c[i]*shadow*(1.2*diff + .6*spec), 255);
//         return false;
//     }
// };

// // z ���������Ƶ�֡������
// struct DepthShader : public IShader {
//     mat<3,3,float> varying_tri;

//     DepthShader() : varying_tri() {}

//     virtual Vec4f vertex(int iface, int nthvert) {
//         Vec4f gl_Vertex = embed<4>(model->vert(iface, nthvert)); // read the vertex from .obj file
//         gl_Vertex = Viewport*Projection*ModelView*gl_Vertex;          // transform it to screen coordinates
//         varying_tri.set_col(nthvert, proj<3>(gl_Vertex/gl_Vertex[3]));
//         return gl_Vertex;
//     }

//     virtual bool fragment(Vec3f bar, TGAColor &color) {
//         Vec3f p = varying_tri*bar;
//         color = TGAColor(255, 255, 255)*(p.z/depth);
//         return false;
//     }
// };

int main(int argc, char** argv) {
    if (2>argc) {
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
            //����ǽ���Ϊ��90��-max_elevation_angle�� / 8 ���ܺ�
            for (float a=0; a<M_PI*2-1e-4; a += M_PI/4) {
                total += M_PI/2 - max_elevation_angle(zbuffer, Vec2f(x, y), Vec2f(cos(a), sin(a)));
            }
            total /= (M_PI/2)*8;
            total = pow(total, 100.f);//���ӶԱȶ�
            frame.set(x, y, TGAColor(total*255, total*255, total*255));
        }
    }

    frame.flip_vertically();
    frame.write_tga_file("framebuffer.tga");
    delete [] zbuffer;
    delete model;
    return 0;
}

//     { // rendering the shadow buffer
//         TGAImage depth(width, height, TGAImage::RGB);
//         //���ǽ���Ⱦ����������ڹ�Դλ�õ�ͼ����������ȷ����Щ���ֱ���������Щ���ֱ������ڹ���֮��
//         lookat(light_dir, center, up);
//         viewport(width/8, height/8, width*3/4, height*3/4);
//         projection(0);

//         DepthShader depthshader;
//         Vec4f screen_coords[3];
//         for (int i=0; i<model->nfaces(); i++) {
//             for (int j=0; j<3; j++) {
//                 screen_coords[j] = depthshader.vertex(i, j);
//             }
//             triangle(screen_coords, depthshader, depth, shadowbuffer);
//         }
//         depth.flip_vertically(); // ��ԭ�����ͼ������½�
//         depth.write_tga_file("depth.tga");
//     }
//     //�����˶�����Ļ��ת������
//     Matrix M = Viewport*Projection*ModelView;

//     { // rendering the frame buffer
//         TGAImage frame(width, height, TGAImage::RGB);
//         lookat(eye, center, up);
//         viewport(width/8, height/8, width*3/4, height*3/4);
//         projection(-1.f/(eye-center).norm());

//         Shader shader(ModelView, (Projection*ModelView).invert_transpose(), M*(Viewport*Projection*ModelView).invert());
//         Vec4f screen_coords[3];
//         // ����ģ�͵����������β���դ��ÿ��������
//         // ��ѭ����������������
//         // �ڲ�ѭ��ѭ��������ǰ�����ε����ж��㣬��Ϊÿ��������ö�����ɫ��
//         for (int i=0; i<model->nfaces(); i++) {
//             for (int j=0; j<3; j++) {
//                 screen_coords[j] = shader.vertex(i, j);
//             }
//             triangle(screen_coords, shader, frame, zbuffer);
//         }
//         frame.flip_vertically(); 
//         frame.write_tga_file("framebuffer.tga");
//     }


// struct GouraudShader : public IShader {
//     Vec3f varying_intensity; // �ɶ�����ɫ��д�룬��Ƭ����ɫ����ȡ
//     //  iface  �ڼ���������  nthvert �����εĵڼ�����
//     virtual Vec4f vertex(int iface, int nthvert) {
//         Vec4f gl_Vertex = embed<4>(model->vert(iface, nthvert)); // ��.obj�ļ��ж�ȡ����
//         gl_Vertex = Viewport*Projection*ModelView*gl_Vertex;     // ����ת��Ϊ��Ļ����
//         varying_intensity[nthvert] = std::max(0.f, model->normal(iface, nthvert)*light_dir); // �������������ǿ��
//         return gl_Vertex;
//     }
//     // �����������꣬����varying_���ݵĲ�ֵ
//     virtual bool fragment(Vec3f bar, TGAColor &color) {
//         float intensity = varying_intensity*bar;   // �Ե�ǰ���ؽ���ǿ�Ȳ�ֵ �򵥵ؼ���Ϊ��������֮��ĵ��
//         // if (intensity>.85) intensity = 1;
//         // else if (intensity>.60) intensity = .80;
//         // else if (intensity>.45) intensity = .60;
//         // else if (intensity>.30) intensity = .45;
//         // else if (intensity>.15) intensity = .30;
//         // else intensity = 0;
//         // color = TGAColor(255, 155, 0)*intensity;
//         color = TGAColor(255, 255, 255)*intensity; 
//         return false;                            
//     }
// };

// // ������Ⱦ
// struct ShaderWenli : public IShader {
//     Vec3f          varying_intensity; 
//     mat<2,3,float> varying_uv;        

//     virtual Vec4f vertex(int iface, int nthvert) {
//         varying_uv.set_col(nthvert, model->uv(iface, nthvert));//?
//         varying_intensity[nthvert] = std::max(0.f, model->normal(iface, nthvert)*light_dir); // get diffuse lighting intensity
//         Vec4f gl_Vertex = embed<4>(model->vert(iface, nthvert)); // read the vertex from .obj file
//         return Viewport*Projection*ModelView*gl_Vertex; // transform it to screen coordinates
//     }
    
//     virtual bool fragment(Vec3f bar, TGAColor &color) {
//         float intensity = varying_intensity*bar;   // interpolate intensity for the current pixel
//         Vec2f uv = varying_uv*bar;                 // interpolate uv for the current pixel
//         color = model->diffuse(uv)*intensity;      
//         return false;                             
//     }
// };


/* zbuffer use

// 4ά����תΪ��
Vec3f m2v(Matrix m) {
    return Vec3f(m[0][0]/m[3][0], m[1][0]/m[3][0], m[2][0]/m[3][0]);
}

//��תΪ4ά����
Matrix v2m(Vec3f v) {
    Matrix m(4, 1);
    m[0][0] = v.x;
    m[1][0] = v.y;
    m[2][0] = v.z;
    m[3][0] = 1.f;
    return m;
}

Vec3f world2screen(Vec3f v) {
    return Vec3f(int((v.x+1.)*width/2.+.5), int((v.y+1.)*height/2.+.5), v.z);
}

void line(int x0, int y0, int x1, int y1, TGAImage &image, TGAColor color) { 
    // ȷ���߶�С�ڿ��  ��Ȼ�п�
    bool steep = false;
    if(std::abs(x1-x0) < std::abs(y1-y0)) {
        std::swap(x0,y0);
        std::swap(x1,y1);
        steep = true;
    } 
    // ȷ���������һ�
    if(x0>x1) {
        std::swap(x0,x1);
        std::swap(y0,y1);
    }
    int dx = x1-x0;
    int dy = y1-y0;
    //float derror = std::abs(dy/float(dx));//����ֱ��б��
    int derror2 = std::abs(dy*2); //б�� * dx *2 ��ΪҪ��.5�Ƚϣ���2��Ϊ�˿�����int
    int error2 = 0;
    int y = y0;
    for (int x=x0; x<=x1; x++) { 
        if(steep) {
            image.set(y, x, color);    
        } else {
            image.set(x, y, color); 
        } 
        error2 += derror2;
        if(error2>dx) {
            y += (y1>y0?1:-1);
            error2 -= dx*2;
        }
    } 
}

//ɨ�߷�  1��Ϊ2
void triangle1(Vec2i t0, Vec2i t1, Vec2i t2, TGAImage &image, TGAColor color) {
    if (t0.y==t1.y && t0.y==t2.y) return;
    // ��y����
    if(t0.y > t1.y) std::swap(t0, t1);
    if(t0.y > t2.y) std::swap(t0, t2);
    if(t1.y > t2.y) std::swap(t1, t2);
    int total_height = t2.y - t0.y;
    for(int i=0; i<total_height; ++i) {
        bool second_half = i>t1.y-t0.y || t1.y==t0.y; //���ڶ���
        int segment_height = second_half ? t2.y-t1.y : t1.y-t0.y;
        float alpha = (float)i/total_height;
        float beta  = (float)(i-(second_half? t1.y-t0.y : 0))/segment_height;
        Vec2i A = t0 + (t2-t0)*alpha;
        Vec2i B = second_half ? t1 + (t2-t1)*beta : t0 + (t1-t0)*beta; 
        if(A.x > B.x) std::swap(A,B);
        for (int j=A.x; j<=B.x; j++) { 
            image.set(j, t0.y+i, color); 
        } 
    }
}

void triangle2(Vec3f *pts, float *zbuffer, TGAImage &image, TGAColor color) {
    Vec2f bboxmin( std::numeric_limits<float>::max(),  std::numeric_limits<float>::max());
    Vec2f bboxmax(-std::numeric_limits<float>::max(), -std::numeric_limits<float>::max());
    Vec2f clamp(image.get_width()-1, image.get_height()-1);
    for (int i=0; i<3; i++) {
        for (int j=0; j<2; j++) {
            bboxmin[j] = std::max(0.f,      std::min(bboxmin[j], pts[i][j]));
            bboxmax[j] = std::min(clamp[j], std::max(bboxmax[j], pts[i][j]));
        }
    }
    Vec3f P;
    for (P.x=bboxmin.x; P.x<=bboxmax.x; P.x++) {
        for (P.y=bboxmin.y; P.y<=bboxmax.y; P.y++) {
            Vec3f bc_screen  = barycentric(pts[0], pts[1], pts[2], P);
            if (bc_screen.x<0 || bc_screen.y<0 || bc_screen.z<0) continue;
            P.z = 0;
            for (int i=0; i<3; i++) P.z += pts[i][2]*bc_screen[i];//����������Ҫ���Ƶ�ÿ�����أ�ֻ�轫�����������������դ�񻯵������ζ���� z ֵ
            if (zbuffer[int(P.x+P.y*width)]<P.z) {
                zbuffer[int(P.x+P.y*width)] = P.z;
                image.set(P.x, P.y, color);
            }
        }
    }
}
*/