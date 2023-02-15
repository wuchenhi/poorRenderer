
#include <cmath>
#include <limits>
#include <cstdlib>
#include "our_gl.h"

Matrix ModelView;
Matrix Viewport;
Matrix Projection;

IShader::~IShader() {}
Shader2::~Shader2() {}

void viewport(int x, int y, int w, int h) {
    Viewport = Matrix::identity();
    Viewport[0][3] = x+w/2.f;
    Viewport[1][3] = y+h/2.f;
    Viewport[2][3] = 1.f;
//  Viewport[2][3] = 255.f/2.f;
    Viewport[0][0] = w/2.f;
    Viewport[1][1] = h/2.f;
    Viewport[2][2] = 0;
//  Viewport[2][2] = 255.f/2.f;
}

void viewportShadow(int x, int y, int w, int h) {
    Viewport = Matrix::identity();
    Viewport[0][3] = x+w/2.f;
    Viewport[1][3] = y+h/2.f;
    Viewport[2][3] = depth/2.f;
    Viewport[0][0] = w/2.f;
    Viewport[1][1] = h/2.f;
    Viewport[2][2] = depth/2.f;
}

void projection(float coeff) {
    Projection = Matrix::identity();
    Projection[3][2] = coeff;
}

// 更改摄像机视角=更改物体位置和角度，操作为互逆矩阵
// 摄像机变换是先旋转再平移，所以物体需要先平移后旋转，且都是逆矩阵
void lookat(Vec3f eye, Vec3f center, Vec3f up) {
    Vec3f z = (eye-center).normalize();// z 为眼睛指向中心
    Vec3f x = cross(up,z).normalize();
    Vec3f y = cross(z,x).normalize();
    Matrix Minv = Matrix::identity();
    Matrix Tr   = Matrix::identity();
    for (int i=0; i<3; i++) {
        Minv[0][i] = x[i];
        Minv[1][i] = y[i];
        Minv[2][i] = z[i];
        Tr[i][3] = -center[i];//矩阵第四列用于平移,因为观察位置从原点变为了center，所以需要将物体平移-cente
    }
    ModelView = Minv*Tr;
}

// 重心坐标
Vec3f barycentric(Vec2f A, Vec2f B, Vec2f C, Vec2f P) {
    Vec3f s[2];
    for (int i=2; i--; ) {
        s[i][0] = C[i]-A[i];
        s[i][1] = B[i]-A[i];
        s[i][2] = A[i]-P[i];
    }
    Vec3f u = cross(s[0], s[1]);// 通过叉积寻找一个同时与 （ABx，ACx，PAx）和 （ABy，ACy，PAy）正交的向量（u，v，1）
    if (std::abs(u[2])>1e-2)    // 三点共线时，会导致u[2]为0
        return Vec3f(1.f-(u.x+u.y)/u.z, u.y/u.z, u.x/u.z);
    return Vec3f(-1,1,1);       // 在这种情况下产生负坐标，它将被光栅化器丢弃
}

void line(int x0, int y0, int x1, int y1, TGAImage &image, TGAColor color) { 
    // 确保高度小于宽度  不然有孔
    bool steep = false;
    if(std::abs(x1-x0) < std::abs(y1-y0)) {
        std::swap(x0,y0);
        std::swap(x1,y1);
        steep = true;
    } 
    // 确保从左往右画
    if(x0>x1) {
        std::swap(x0,x1);
        std::swap(y0,y1);
    }
    int dx = x1-x0;
    int dy = y1-y0;
    //float derror = std::abs(dy/float(dx));//所画直线斜率
    int derror2 = std::abs(dy*2); //斜率 * dx *2 因为要和.5比较，乘2是为了可以用int
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

void triangleShadow(Vec4f *pts, Shader2 &shader, TGAImage &image, float *zbuffer) {
    Vec2f bboxmin( std::numeric_limits<float>::max(),  std::numeric_limits<float>::max());
    Vec2f bboxmax(-std::numeric_limits<float>::max(), -std::numeric_limits<float>::max());
    for (int i=0; i<3; i++) {
        for (int j=0; j<2; j++) {
            bboxmin[j] = std::min(bboxmin[j], pts[i][j]/pts[i][3]);
            bboxmax[j] = std::max(bboxmax[j], pts[i][j]/pts[i][3]);
        }
    }
    Vec2i P;
    TGAColor color;
    for (P.x=bboxmin.x; P.x<=bboxmax.x; P.x++) {
        for (P.y=bboxmin.y; P.y<=bboxmax.y; P.y++) {
            Vec3f c = barycentric(proj<2>(pts[0]/pts[0][3]), proj<2>(pts[1]/pts[1][3]), proj<2>(pts[2]/pts[2][3]), proj<2>(P));
            float z = pts[0][2]*c.x + pts[1][2]*c.y + pts[2][2]*c.z;
            float w = pts[0][3]*c.x + pts[1][3]*c.y + pts[2][3]*c.z;
            int frag_depth = z/w;
            if (c.x<0 || c.y<0 || c.z<0 || zbuffer[P.x+P.y*image.get_width()]>frag_depth) continue;
            bool discard = shader.fragment(c, color);
            if (!discard) {
                zbuffer[P.x+P.y*image.get_width()] = frag_depth;
                image.set(P.x, P.y, color);
            }
        }
    }
}

void triangle(mat<4,3,float> &clipc, IShader &shader, TGAImage &image, float *zbuffer) {
    mat<3,4,float> pts  = (Viewport*clipc).transpose(); // transposed 以方便访问每个点s
    mat<3,2,float> pts2;
    for (int i=0; i<3; i++) pts2[i] = proj<2>(pts[i]/pts[i][3]);

    Vec2f bboxmin( std::numeric_limits<float>::max(),  std::numeric_limits<float>::max());
    Vec2f bboxmax(-std::numeric_limits<float>::max(), -std::numeric_limits<float>::max());
    Vec2f clamp(image.get_width()-1, image.get_height()-1);
    for (int i=0; i<3; i++) {
        for (int j=0; j<2; j++) {
            bboxmin[j] = std::max(0.f,      std::min(bboxmin[j], pts2[i][j]));
            bboxmax[j] = std::min(clamp[j], std::max(bboxmax[j], pts2[i][j]));
        }
    }
    Vec2i P;
    TGAColor color;
    for (P.x=bboxmin.x; P.x<=bboxmax.x; P.x++) {
        for (P.y=bboxmin.y; P.y<=bboxmax.y; P.y++) {
            Vec3f bc_screen  = barycentric(pts2[0], pts2[1], pts2[2], P);
            Vec3f bc_clip    = Vec3f(bc_screen.x/pts[0][3], bc_screen.y/pts[1][3], bc_screen.z/pts[2][3]);
            bc_clip = bc_clip/(bc_clip.x+bc_clip.y+bc_clip.z);
            float frag_depth = clipc[2]*bc_clip;
            if (bc_screen.x<0 || bc_screen.y<0 || bc_screen.z<0 || zbuffer[P.x+P.y*image.get_width()]>frag_depth) continue;
            bool discard = shader.fragment(Vec3f(P.x, P.y, frag_depth), bc_clip, color);
            if (!discard) {
                zbuffer[P.x+P.y*image.get_width()] = frag_depth;
                image.set(P.x, P.y, color);
            }
        }
    }
}

//扫线法  1分为2
void triangleSao(Vec2i t0, Vec2i t1, Vec2i t2, TGAImage &image, TGAColor color) {
    if (t0.y==t1.y && t0.y==t2.y) return;
    // 按y升序
    if(t0.y > t1.y) std::swap(t0, t1);
    if(t0.y > t2.y) std::swap(t0, t2);
    if(t1.y > t2.y) std::swap(t1, t2);
    int total_height = t2.y - t0.y;
    for(int i=0; i<total_height; ++i) {
        bool second_half = i>t1.y-t0.y || t1.y==t0.y; //画第二半
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

// void triangleZbuffer(Vec3f *pts, float *zbuffer, TGAImage &image, TGAColor color) {
//     Vec2f bboxmin( std::numeric_limits<float>::max(),  std::numeric_limits<float>::max());
//     Vec2f bboxmax(-std::numeric_limits<float>::max(), -std::numeric_limits<float>::max());
//     Vec2f clamp(image.get_width()-1, image.get_height()-1);
//     for (int i=0; i<3; i++) {
//         for (int j=0; j<2; j++) {
//             bboxmin[j] = std::max(0.f,      std::min(bboxmin[j], pts[i][j]));
//             bboxmax[j] = std::min(clamp[j], std::max(bboxmax[j], pts[i][j]));
//         }
//     }
//     Vec3f P;
//     for (P.x=bboxmin.x; P.x<=bboxmax.x; P.x++) {
//         for (P.y=bboxmin.y; P.y<=bboxmax.y; P.y++) {
//             Vec3f bc_screen  = barycentric(pts[0], pts[1], pts[2], P);
//             if (bc_screen.x<0 || bc_screen.y<0 || bc_screen.z<0) continue;
//             P.z = 0;
//             for (int i=0; i<3; i++) P.z += pts[i][2]*bc_screen[i];//对于我们想要绘制的每个像素，只需将其重心坐标乘以我们栅格化的三角形顶点的 z 值
//             if (zbuffer[int(P.x+P.y*width)]<P.z) {
//                 zbuffer[int(P.x+P.y*width)] = P.z;
//                 image.set(P.x, P.y, color);
//             }
//         }
//     }
// }

// 光栅器，对于三角形中的每个点调用片段着色器，然后执行深度检查（Z-buffer）
// void triangle1(Vec4f *pts, IShader &shader, TGAImage &image, TGAImage &zbuffer) {
//     Vec2f bboxmin( std::numeric_limits<float>::max(),  std::numeric_limits<float>::max());
//     Vec2f bboxmax(-std::numeric_limits<float>::max(), -std::numeric_limits<float>::max());
//     for (int i=0; i<3; i++) {
//         for (int j=0; j<2; j++) {
//             bboxmin[j] = std::min(bboxmin[j], pts[i][j]/pts[i][3]);
//             bboxmax[j] = std::max(bboxmax[j], pts[i][j]/pts[i][3]);
//         }
//     }
//     Vec2i P;
//     TGAColor color;
//     for (P.x=bboxmin.x; P.x<=bboxmax.x; P.x++) {
//         for (P.y=bboxmin.y; P.y<=bboxmax.y; P.y++) {
//             Vec3f c = barycentric(proj<2>(pts[0]/pts[0][3]), proj<2>(pts[1]/pts[1][3]), proj<2>(pts[2]/pts[2][3]), proj<2>(P));
//             float z = pts[0][2]*c.x + pts[1][2]*c.y + pts[2][2]*c.z;
//             float w = pts[0][3]*c.x + pts[1][3]*c.y + pts[2][3]*c.z;
//             int frag_depth = std::max(0, std::min(255, int(z/w+.5)));
//             if (c.x<0 || c.y<0 || c.z<0 || zbuffer.get(P.x, P.y)[0]>frag_depth) continue;
//             bool discard = shader.fragment(c, color);
//             if (!discard) {
//                 zbuffer.set(P.x, P.y, TGAColor(frag_depth));
//                 image.set(P.x, P.y, color);
//             }
//         }
//     }
// }
