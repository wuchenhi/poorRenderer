
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

// ����������ӽ�=��������λ�úͽǶȣ�����Ϊ�������
// ������任������ת��ƽ�ƣ�����������Ҫ��ƽ�ƺ���ת���Ҷ��������
void lookat(Vec3f eye, Vec3f center, Vec3f up) {
    Vec3f z = (eye-center).normalize();// z Ϊ�۾�ָ������
    Vec3f x = cross(up,z).normalize();
    Vec3f y = cross(z,x).normalize();
    Matrix Minv = Matrix::identity();
    Matrix Tr   = Matrix::identity();
    for (int i=0; i<3; i++) {
        Minv[0][i] = x[i];
        Minv[1][i] = y[i];
        Minv[2][i] = z[i];
        Tr[i][3] = -center[i];//�������������ƽ��,��Ϊ�۲�λ�ô�ԭ���Ϊ��center��������Ҫ������ƽ��-cente
    }
    ModelView = Minv*Tr;
}

// ��������
Vec3f barycentric(Vec2f A, Vec2f B, Vec2f C, Vec2f P) {
    Vec3f s[2];
    for (int i=2; i--; ) {
        s[i][0] = C[i]-A[i];
        s[i][1] = B[i]-A[i];
        s[i][2] = A[i]-P[i];
    }
    Vec3f u = cross(s[0], s[1]);// ͨ�����Ѱ��һ��ͬʱ�� ��ABx��ACx��PAx���� ��ABy��ACy��PAy��������������u��v��1��
    if (std::abs(u[2])>1e-2)    // ���㹲��ʱ���ᵼ��u[2]Ϊ0
        return Vec3f(1.f-(u.x+u.y)/u.z, u.y/u.z, u.x/u.z);
    return Vec3f(-1,1,1);       // ����������²��������꣬��������դ��������
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
    mat<3,4,float> pts  = (Viewport*clipc).transpose(); // transposed �Է������ÿ����s
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

//ɨ�߷�  1��Ϊ2
void triangleSao(Vec2i t0, Vec2i t1, Vec2i t2, TGAImage &image, TGAColor color) {
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
//             for (int i=0; i<3; i++) P.z += pts[i][2]*bc_screen[i];//����������Ҫ���Ƶ�ÿ�����أ�ֻ�轫�����������������դ�񻯵������ζ���� z ֵ
//             if (zbuffer[int(P.x+P.y*width)]<P.z) {
//                 zbuffer[int(P.x+P.y*width)] = P.z;
//                 image.set(P.x, P.y, color);
//             }
//         }
//     }
// }

// ��դ���������������е�ÿ�������Ƭ����ɫ����Ȼ��ִ����ȼ�飨Z-buffer��
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
