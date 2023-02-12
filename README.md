复制过来查看 sudo cp framebuffer.tga /mnt/e

g++ -ggdb -g -pg -O0

```
int main(int argc, char** argv) {
    TGAImage image(100, 100, TGAImage::RGB);
    for(int i= 0; i<1000000; ++i) {
        //line(13, 20, 80, 40, image, white);
        line1(20, 13, 40, 80, image, red); 
        line1(80, 40, 13, 20, image, red);
        image.flip_vertically();
    }
    image.write_tga_file("output.tga");
    return 0;
}
```

./main 运行后 gprof main gmon.out

```
time   seconds   seconds    calls  Ts/call  Ts/call  name  
 56.56      1.09     1.09                             line1(int, int, int, int, TGAImage&, TGAColor)
 18.42      1.45     0.36                             TGAImage::set(int, int, TGA
```

画线处需要优化 优化后

```
time   seconds   seconds    calls  Ts/call  Ts/call  name  
 32.91      0.46     0.46                             line1(int, int, int, int, TGAImage&, TGAColor)
 32.55      0.92     0.46                             TGAImage::set(i
```

float 换int   seconds   seconds    calls  Ts/call  Ts/call  name
