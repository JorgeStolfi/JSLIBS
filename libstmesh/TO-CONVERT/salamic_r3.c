
/* See {salamic_r3.h} */

#include <math.h>
#include <stdint.h>
/*
#include <glm/glm.h>
#include <glm/gtc/matrix_transform.h>
*/
#include <salamic_r3.h>
#include <salamic_utils.h>

int32_t salamic_r3_hash (const salamic_stl_r3_t *v) {
  int32_t h = 0;
  salamic_utils_hash_combine(&h, salamic_utils_float_hash(v->x));
  salamic_utils_hash_combine(&h, salamic_utils_float_hash(v->y));
  salamic_utils_hash_combine(&h, salamic_utils_float_hash(v->z));
  return h;
};


/*
salamic_r3_read(string label) {
  
		istringstream ss(label);
		string coord;
		getline(ss, coord, '|');
		x = atof(coord.c_str());
		getline(ss, coord, '|');
		y = atof(coord.c_str());
		getline(ss, coord, '|');
		z = atof(coord.c_str());
	}

    bool seen;
    

	void setCoords(array<float,3> a) {
		x = a[0]; y = a[1]; z = a[2];
	}
	void setCoords(string label) {
		istringstream ss(label);
		string coord;
		getline(ss, coord, '|');
		x = atof(coord.c_str());
		getline(ss, coord, '|');
		y = atof(coord.c_str());
		getline(ss, coord, '|');
		z = atof(coord.c_str());
	}
        void set_seen (bool value) {
           seen = value;
        }

	array<float,3> getCoords() { array<float,3> c = {{x, y, z}}; return c; }
    float dotproduct(const salamic_stl_r3_t &v) const { return (x*v.x + y*v.y + z*v.z); }

float salamic_r3_distTo(const salamic_stl_r3_t *a, const salamic_stl_r3_t *b) {
  float dx = a->x - b->x;
  float dy = a->y - b->y;
  float dz = a->z - b->z;
  return sqrt(dx*dx + dy*dy + dz*dz); 
}

void transform(salamic_stl_r3_t *v, const glm::mat4 *mat) { 
  glm::vec4 v4 = glm::vec4(v->x, v->y, v->z, 1.0);
  glm::vec4 vt = (*mat)*v4; 
  v->x = vt.x; v->y = vt.y; v->z = vt.z;
  }
*/

/*
    salamic_stl_r3_t& operator-=(const salamic_stl_r3_t &pt) { x= (x-pt.x); y=(y-pt.y); z=(z-pt.z); return *this;    }
    salamic_stl_r3_t operator-(const salamic_stl_r3_t &pt) { return salamic_stl_r3_t((x-pt.x), (y-pt.y), (z-pt.z)); }
    salamic_stl_r3_t operator+(const salamic_stl_r3_t &pt) { return salamic_stl_r3_t((x+pt.x), (y+pt.y), (z+pt.z)); }
    salamic_stl_r3_t operator/(float a) { return salamic_stl_r3_t((x/a), (y/a), (z/a)); }
    salamic_stl_r3_t operator*(float a) { return salamic_stl_r3_t((x*a), (y*a), (z*a)); }
	bool operator<(const salamic_stl_r3_t &pt) const { return z < pt.z; }
	bool operator>(const salamic_stl_r3_t &pt) const { return z > pt.z; }
	//bool operator==(const salamic_stl_r3_t &pt) {return compf(x, pt.x) && compf(y, pt.y) && compf(z, pt.z); }
	bool operator==(const salamic_stl_r3_t &pt) const {return distTo(pt) < 0.005; }
	//bool operator!=(const salamic_stl_r3_t &pt) {return !(compf(x, pt.x) && compf(y, pt.y) && compf(z, pt.z)); }
	bool operator!=(const salamic_stl_r3_t &pt) const { return distTo(pt) > 0.005; }
    float normalize() const { return sqrt(x*x+y*y+z*z); }
	string getLabel() const {
		stringstream ss;
		ss << roundIt(x) << "|" << roundIt(y) << "|" << roundIt(z);
		//ss << x << "|" << y << "|" << z;
		return ss.str();
	}
salamic_stl_r3_t operator-(const salamic_stl_r3_t &a, const salamic_stl_r3_t &b) {return salamic_stl_r3_t((a.x-b.x), (a.y-b.y), (a.z-b.z)); }
salamic_stl_r3_t operator+(const salamic_stl_r3_t &a, const salamic_stl_r3_t &b) {return salamic_stl_r3_t((a.x+b.x), (a.y+b.y), (a.z+b.z)); }

ostream& operator<<(ostream& os, const salamic_stl_r3_t& v) {
		os << "x: " << v.x << "; y: " << v.y << "; z: " << v.z;
		return os;
	}

class salamic_r3_Plane_t {
public:
    salamic_r3_Plane_t() : mDistance(0) {}
    float distance() const { return mDistance; }
    float distanceToPoint(const salamic_stl_r3_t &vertex) const {    return vertex.dotproduct(mNormal) - mDistance; }
    void setNormal(salamic_stl_r3_t normal) { mNormal = normal; }
    void setDistance(float distance) { mDistance = distance; }
protected:
    salamic_stl_r3_t        mNormal;    // normalized Normal-Vector of the plane
    float   mDistance;  // shortest distance from plane to Origin
};
*/

/*
public:
    salamic_r3_Segment_t(salamic_stl_r3_t p0=salamic_stl_r3_t(), salamic_stl_r3_t p1=salamic_stl_r3_t(), int32_t l=0) { v[0]=p0; v[1]=p1; label=l; }
    int32_t label;
	bool operator==(const salamic_r3_Segment_t &ls) const { return v[0] == ls.v[0] && v[1] == ls.v[1]; }
	int32_t interceptPlane(salamic_r3_Plane_t &plane, salamic_stl_r3_t &point) {
        salamic_stl_r3_t a = v[0];
		salamic_stl_r3_t b = v[1];
		const float da = plane.distanceToPoint(a);
        const float db = plane.distanceToPoint(b);
        if (da*db<0) {
            const float s = da/(da-db); // intersection factor (between 0 and 1)
            salamic_stl_r3_t bMinusa = b-a;
            point = a+bMinusa*s;
			return 0;
        }
        else if (0==da) { // plane falls exactly on one of the three Triangle vertices
            point = a;
			return 0;
        }
        else if (0==db) { // plane falls exactly on one of the three Triangle vertices
            point = b;
			return 0;
        }
		return 0;
	}
	
	friend ostream& operator<<(ostream& os, const salamic_r3_Segment_t& ls) {
		os << "V0: (" << ls.v[0] << "); V1: (" << ls.v[1] << ")"	;
		return os;
	}
*/


