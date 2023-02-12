#include <vector>
#include <cmath>
#include "tgaimage.h"
#include "model.h"
#include "geometry.h"

const TGAColor white = TGAColor(255, 255, 255, 255);
const TGAColor red   = TGAColor(255, 0,   0,   255);
const TGAColor green = TGAColor(0,   255, 0,   255);
Model *model = NULL;
const int width  = 800;
const int height = 800;

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
void triangle(Vec2i t0, Vec2i t1, Vec2i t2, TGAImage &image, TGAColor color) {
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

int main(int argc, char** argv) {
    
    if(argc == 2) {
        model = new Model(argv[1]);
    }else {
        model = new Model("obj/african_head/african_head.obj");
    }

    TGAImage image(width, height, TGAImage::RGB);

    Vec3f light_dir(0, 0, -1);

    for(int i=0; i<model->nfaces(); ++i) {
        std::vector<int> face = model->face(i);
        Vec2i screen_coords[3]; 
        Vec3f world_coords[3];
        for(int j=0; j<3; ++j) {
            Vec3f v = model->vert(face[j]);
            screen_coords[j] = Vec2i((v.x+1.)*width/2., (v.y+1.)*height/2.);//��Ļ����
            world_coords[j]  = v; //��������
        } 
        Vec3f n = (world_coords[2]-world_coords[0])^(world_coords[1]-world_coords[0]);//�����εķ��߿��Լ򵥵ؼ���Ϊ�������ߵĲ��
        n.normalize();
        float intensity = n*light_dir;//����ǿ�ȵ��ڹ�ʸ���ı����˻��͸��������εķ���
        if (intensity>0) {
            triangle(screen_coords[0], screen_coords[1], screen_coords[2], image, TGAColor(intensity*255, intensity*255, intensity*255, 255));
        }
    }
    image.flip_vertically(); // i want to have the origin at the left bottom corner of the image
    image.write_tga_file("output.tga");
    delete model;
    return 0;
   
}