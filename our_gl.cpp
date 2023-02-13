
#include <cmath>
#include <limits>
#include <cstdlib>
#include "our_gl.h"

Matrix ModelView;
Matrix Viewport;
Matrix Projection;

IShader::~IShader() {}

void viewport(int x, int y, int w, int h) {
    Viewport = Matrix::identity();
    Viewport[0][3] = x+w/2.f;
    Viewport[1][3] = y+h/2.f;
    Viewport[2][3] = 1.f;
    Viewport[0][0] = w/2.f;
    Viewport[1][1] = h/2.f;
    Viewport[2][2] = 0;
}

void projection(float coeff) {
    Projection = Matrix::identity();
    Projection[3][2] = coeff;
}

void lookat(Vec3f eye, Vec3f center, Vec3f up) {
    Vec3f z = (eye-center).normalize();
    Vec3f x = cross(up,z).normalize();
    Vec3f y = cross(z,x).normalize();
    Matrix Minv = Matrix::identity();
    Matrix Tr   = Matrix::identity();
    for (int i=0; i<3; i++) {
        Minv[0][i] = x[i];
        Minv[1][i] = y[i];
        Minv[2][i] = z[i];
        Tr[i][3] = -center[i];
    }
    ModelView = Minv*Tr;
}

Vec3f barycentric(Vec2f A, Vec2f B, Vec2f C, Vec2f P) {
    Vec3f s[2];
    for (int i=2; i--; ) {
        s[i][0] = C[i]-A[i];
        s[i][1] = B[i]-A[i];
        s[i][2] = A[i]-P[i];
    }
    Vec3f u = cross(s[0], s[1]);
    if (std::abs(u[2])>1e-2) // dont forget that u[2] is integer. If it is zero then triangle ABC is degenerate
        return Vec3f(1.f-(u.x+u.y)/u.z, u.y/u.z, u.x/u.z);
    return Vec3f(-1,1,1); // in this case generate negative coordinates, it will be thrown away by the rasterizator
}

void triangle(mat<4,3,float> &clipc, IShader &shader, TGAImage &image, float *zbuffer) {
    mat<3,4,float> pts  = (Viewport*clipc).transpose(); // transposed to ease access to each of the points
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


// #include <cmath>
// #include <limits>
// #include <cstdlib>
// #include "our_gl.h"

// Matrix ModelView;
// Matrix Viewport;
// Matrix Projection;

// IShader::~IShader() {}

// void viewport(int x, int y, int w, int h) {
//     Viewport = Matrix::identity();
//     Viewport[0][3] = x+w/2.f;
//     Viewport[1][3] = y+h/2.f;
//     Viewport[2][3] = 255.f/2.f;
//     Viewport[0][0] = w/2.f;
//     Viewport[1][1] = h/2.f;
//     Viewport[2][2] = 255.f/2.f;
// }

// void projection(float coeff) {
//     Projection = Matrix::identity();
//     Projection[3][2] = coeff;
// }

// void lookat(Vec3f eye, Vec3f center, Vec3f up) {
//     Vec3f z = (eye-center).normalize();//z为眼睛指向中心
//     Vec3f x = cross(up,z).normalize();
//     Vec3f y = cross(z,x).normalize();
//     ModelView = Matrix::identity();
//     for (int i=0; i<3; i++) {
//         ModelView[0][i] = x[i];
//         ModelView[1][i] = y[i];
//         ModelView[2][i] = z[i];
//         ModelView[i][3] = -center[i];
//     }
// }

// Vec3f barycentric(Vec2f A, Vec2f B, Vec2f C, Vec2f P) {
//     Vec3f s[2];
//     for (int i=2; i--; ) {
//         s[i][0] = C[i]-A[i];
//         s[i][1] = B[i]-A[i];
//         s[i][2] = A[i]-P[i];
//     }
//     Vec3f u = cross(s[0], s[1]);// 寻找一个同时与 （ABx，ACx，PAx） 和 （ABy，ACy，PAy） 正交的向量 （u，v，1）  叉积满足
//     if (std::abs(u[2])>1e-2) // u[2]是整数。如果它是零，那么三角形ABC是退化的
//         return Vec3f(1.f-(u.x+u.y)/u.z, u.y/u.z, u.x/u.z); //（1 ? u ? v，u，v）
//     return Vec3f(-1,1,1); // 在这种情况下产生负坐标，它将被光栅化器丢弃
// }

// void triangle(mat<4,3,float> &clipc, IShader &shader, TGAImage &image, float *zbuffer) {
//     mat<3,4,float> pts  = (Viewport*clipc).transpose(); // transposed to ease access to each of the points
//     mat<3,2,float> pts2;
//     for (int i=0; i<3; i++) pts2[i] = proj<2>(pts[i]/pts[i][3]);

//     Vec2f bboxmin( std::numeric_limits<float>::max(),  std::numeric_limits<float>::max());
//     Vec2f bboxmax(-std::numeric_limits<float>::max(), -std::numeric_limits<float>::max());
//     Vec2f clamp(image.get_width()-1, image.get_height()-1);
//     for (int i=0; i<3; i++) {
//         for (int j=0; j<2; j++) {
//             bboxmin[j] = std::max(0.f,      std::min(bboxmin[j], pts2[i][j]));
//             bboxmax[j] = std::min(clamp[j], std::max(bboxmax[j], pts2[i][j]));
//         }
//     }
//     Vec2i P;
//     TGAColor color;
//     for (P.x=bboxmin.x; P.x<=bboxmax.x; P.x++) {
//         for (P.y=bboxmin.y; P.y<=bboxmax.y; P.y++) {
//             Vec3f bc_screen  = barycentric(pts2[0], pts2[1], pts2[2], P);
//             Vec3f bc_clip    = Vec3f(bc_screen.x/pts[0][3], bc_screen.y/pts[1][3], bc_screen.z/pts[2][3]);
//             bc_clip = bc_clip/(bc_clip.x+bc_clip.y+bc_clip.z);
//             float frag_depth = clipc[2]*bc_clip;
//             if (bc_screen.x<0 || bc_screen.y<0 || bc_screen.z<0 || zbuffer[P.x+P.y*image.get_width()]>frag_depth) continue;
//             bool discard = shader.fragment(Vec3f(P.x, P.y, frag_depth), bc_clip, color);
//             if (!discard) {
//                 zbuffer[P.x+P.y*image.get_width()] = frag_depth;
//                 image.set(P.x, P.y, color);
//             }
//         }
//     }
// }

// void triangle2(Vec4f *pts, IShader &shader, TGAImage &image, float *zbuffer) {
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
//             int frag_depth = z/w;
//             if (c.x<0 || c.y<0 || c.z<0 || zbuffer[P.x+P.y*image.get_width()]>frag_depth) continue;
//             bool discard = shader.fragment(c, color);
//             if (!discard) {
//                 zbuffer[P.x+P.y*image.get_width()] = frag_depth;
//                 image.set(P.x, P.y, color);
//             }
//         }
//     }
// }

// //光栅器，对于三角形中的每个点调用片段着色器，然后执行深度检查（Z-buffer）
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
