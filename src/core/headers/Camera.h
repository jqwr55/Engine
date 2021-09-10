#pragma once
#include <Common.h>

template<typename T> struct v2 {
    union {
        struct {
            T x;
            T y;
        };
        T arr[2];
    };
    T& operator [](u32 i) {
        return arr[i];
    }
    v2 operator + (T scalar) {
        return {x + scalar , y + scalar};
    }
    v2 operator - (T scalar) {
        return {x - scalar , y - scalar};
    }
    v2 operator * (T scalar) {
        return {x * scalar , y * scalar};
    }
    v2 operator / (T scalar) {
        return {x / scalar , y / scalar};
    }
    v2 operator + (v2<T> other) {
        return {x + other.x , y + other.y};
    }
    v2 operator - (v2<T> other) {
        return {x - other.x , y - other.y};
    }
    v2 operator * (v2<T> other) {
        return {x * other.x , y * other.y};
    }
    v2 operator / (v2<T> other) {
        return {x / other.x , y / other.y};
    }
    v2 operator == (v2<T> other) {
        return {x == other.x && y == other.y};
    }
    v2 operator != (v2<T> other) {
        return {x != other.x || y != other.y};
    }
    T len() {
        return sqrt(x * x + y * y);
    }
    v2 normalize() {
        T l = len();
        return {x/l, y/l};
    }
};
template<typename T> struct v3 {
    union {
        struct {
            T x;
            T y;
            T z;
        };
        T arr[3];
    };
    T& operator [](u32 i) {
        return arr[i];
    }
    v3 operator + (T scalar) {
        return {x + scalar , y + scalar , z + scalar};
    }
    v3 operator - (T scalar) {
        return {x - scalar , y - scalar , z - scalar};
    }
    v3 operator * (T scalar) {
        return {x * scalar , y * scalar , z * scalar};
    }
    v3 operator / (T scalar) {
        return {x / scalar , y / scalar , z / scalar};
    }
    
    v3 operator + (v3<T> other) {
        return {x + other.x , y + other.y , z + other.z};
    }
    v3 operator - (v3<T> other) {
        return {x - other.x , y - other.y , z - other.z};
    }
    v3 operator * (v3<T> other) {
        return {x * other.x , y * other.y , z * other.z};
    }
    v3 operator / (v3<T> other) {
        return {x / other.x , y / other.y , z / other.z};
    }
    v3 operator == (v3<T> other) {
        return {x == other.x && y == other.y && z == other.z};
    }
    v3 operator != (v3<T> other) {
        return {x != other.x || y != other.y || z != other.z};
    }
    T len() {
        return sqrt(x * x + y * y + z * z);
    }
    v3 normalize() {
        T l = len();
        return {x/l, y/l, z/l};
    }
};

template<typename T> struct v4 {
    union {
        struct {
            T x;
            T y;
            T z;
            T w;
        };
        T arr[4];
    };
    T& operator [](u32 i) {
        return arr[i];
    }
    v4 operator + (T scalar) {
        return {x + scalar , y + scalar , z + scalar, w + scalar};
    }
    v4 operator - (T scalar) {
        return {x - scalar , y - scalar , z - scalar, w - scalar};
    }
    v4 operator * (T scalar) {
        return {x * scalar , y * scalar , z * scalar , w * scalar};
    }
    v4 operator / (T scalar) {
        return {x / scalar , y / scalar , z / scalar, w / scalar};
    }
    
    v4 operator + (v4<T> other) {
        return {x + other.x , y + other.y , z + other.z, w + other.w};
    }
    v4 operator - (v4<T> other) {
        return {x - other.x , y - other.y , z - other.z, w - other.w};
    }
    v4 operator * (v4<T> other) {
        return {x * other.x , y * other.y , z * other.z, w * other.w};
    }
    v4 operator / (v4<T> other) {
        return {x / other.x , y / other.y , z / other.z, w / other.w};
    }
    v4 operator == (v4<T> other) {
        return {x == other.x && y == other.y && z == other.z && w == other.w};
    }
    v4 operator != (v4<T> other) {
        return {x != other.x || y != other.y || z != other.z || w != other.w};
    }
    T len() {
        return sqrt(x * x + y * y + z * z + w * w);
    }
    v4 normalize() {
        T l = len();
        return {x/l, y/l, z/l, w/l};
    }
};
template<typename T> T dot(v3<T> l, v3<T> r) {
    return (l.x * r.x) + (l.y * r.y) + (l.z * r.z);
}
template<typename T> T dot(v4<T> l, v4<T> r) {
    return (l.x * r.x) + (l.y * r.y) + (l.z * r.z) + (l.w * r.w);
}
template<typename T> v3<T> cross(v3<T> l , v3<T> r) {
    return {
        l.y * r.z - l.z * r.y,
        l.z * r.x - l.x * r.z,
        l.x * r.y - l.y * r.x,
    };
}
template<typename T> struct Mat4 {
    v4<T> arr[4];

    v4<T>& operator [](u32 i) {
        return arr[i];
    }
    Mat4<T> operator * (Mat4<T> other) {
        Mat4<T> ret;

        ret.arr[0] = (arr[0] * other.arr[0].arr[0]) + (arr[1] * other.arr[0].arr[1]) + (arr[2] * other.arr[0].arr[2]) + (arr[3] * other.arr[0].arr[3]);
        ret.arr[1] = (arr[0] * other.arr[1].arr[0]) + (arr[1] * other.arr[1].arr[1]) + (arr[2] * other.arr[1].arr[2]) + (arr[3] * other.arr[1].arr[3]);
        ret.arr[2] = (arr[0] * other.arr[2].arr[0]) + (arr[1] * other.arr[2].arr[1]) + (arr[2] * other.arr[2].arr[2]) + (arr[3] * other.arr[2].arr[3]);
        ret.arr[3] = (arr[0] * other.arr[3].arr[0]) + (arr[1] * other.arr[3].arr[1]) + (arr[2] * other.arr[3].arr[2]) + (arr[3] * other.arr[3].arr[3]);
        return ret;
    }
};

struct Camera {
    v3<f32> position;
    v3<f32> direction;
    v3<f32> vel;
    u8 keys;
};

void MoveCameraAlong(Camera& cam) {

    u8 keys = cam.keys;

    u8 w = ((keys >> 0) & 1);
    u8 a = ((keys >> 1) & 1);
    u8 s = ((keys >> 2) & 1);
    u8 d = ((keys >> 3) & 1);
    u8 space = ((keys >> 4) & 1);
    u8 shift = ((keys >> 5) & 1);

    glm::vec2 horizontalForwardDir{cam.direction.x ,cam.direction.z};
    horizontalForwardDir = glm::normalize( horizontalForwardDir );
    glm::vec2 horizontalOrtoDir{horizontalForwardDir.y , -horizontalForwardDir.x};

    int8_t forward = w-s;
    int8_t ortogonal = a-d;

    const float speed = 0.00003f;

    cam.vel.x += (horizontalForwardDir.x * forward * speed) + (horizontalOrtoDir.x * ortogonal * speed);
    cam.vel.z += (horizontalForwardDir.y * forward * speed) + (horizontalOrtoDir.y * ortogonal * speed);

    cam.keys = 0;
}


void RotateCamera(Camera* cam , float vertRotAngle , float horizRotAngle) {

    f32 cosHoriz = cos(horizRotAngle);
    f32 sinHoriz = sin(horizRotAngle);

    f32 cosVert = cos(vertRotAngle);
    f32 sinVert = sin(vertRotAngle);

    cam->direction.x = cam->direction.x * cosHoriz - cam->direction.z * sinHoriz;
    cam->direction.z = cam->direction.x * sinHoriz + cam->direction.z * cosHoriz;

    cam->direction = cam->direction.normalize();
    v3<f32> right = cross(cam->direction, {0,1,0}).normalize();
    v3<f32> w = cross(right, cam->direction).normalize();

    cam->direction = cam->direction * cosVert + w * sinVert;
    cam->direction = cam->direction.normalize();
}

Mat4<f32> LookAt(v3<f32> from, v3<f32> to, v3<f32> worldUp = {0,1,0}) {

    v3<f32> forward{ (to - from).normalize() };
    v3<f32> right{ cross(forward, worldUp).normalize() };
    v3<f32> up{ cross(right, forward).normalize() };

    return {
        right.x , up.x , -forward.x , 0,
        right.y , up.y , -forward.y , 0,
        right.z , up.z , -forward.z , 0,
        -dot(right, from), -dot(up, from), dot(forward, from) , 1
    };
}

Mat4<f32> PerspectiveMatrix_(f32 fov , f32 aspect , f32 near , f32 far) {

    f32 tanFov = tan( fov * 0.5 );
    f32 x = 1 / ( aspect * tanFov );
    f32 y = 1 / ( tanFov );

    f32 z = -(far + near) / (far - near);
    f32 w = (-2 * far * near) / (far - near);

    return Mat4<f32> {
        x,0,0,0,
        0,y,0,0,
        0,0,z,-1,
        0,0,w,0
    };
}
void PrintglmMat(glm::mat4 mat) {
    for(u32 i = 0; i < 16; i++) {
        std::cout << ((f32*)&mat)[i] << " ";
        if((i+1) % 4 == 0) std::cout << std::endl;
    }
}
void PrintMat4(Mat4<f32> mat) {
    for(u32 i = 0; i < 16; i++) {
        std::cout << ((f32*)&mat)[i] << " ";
        if((i+1) % 4 == 0) std::cout << std::endl;
    }
}

/*
glm::mat4 LookAt(glm::vec3 from , glm::vec3 to , glm::vec3 worldUp = {0,1,0} ) {
    glm::vec3 forward{ glm::normalize(to - from) };
    glm::vec3 right{glm::normalize(glm::cross( forward , worldUp ))};
    glm::vec3 up{glm::normalize(glm::cross(right , forward))};

    return glm::mat4 {
       right.x , up.x , -forward.x , 0,
       right.y , up.y , -forward.y , 0,
       right.z , up.z , -forward.z , 0,
       -glm::dot(right , from),-glm::dot(up , from),glm::dot(forward , from),1
    };
}

glm::mat4 PerspectiveMatrix(f32 fov , f32 aspect , f32 near , f32 far) {

    f32 tanFov = tan( fov * 0.5 );
    f32 x = 1 / ( aspect * tanFov );
    f32 y = 1 / ( tanFov );
    f32 z = -(far + near) / (far - near);
    f32 w = (-2 * far * near) / (far - near);

    return glm::mat4{
        x,0,0,0,
        0,y,0,0,
        0,0,z,-1,
        0,0,w,0
    };
}
*/