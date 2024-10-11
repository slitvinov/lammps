#ifndef LMP_RIGID_CONST_H
#define LMP_RIGID_CONST_H 
namespace LAMMPS_NS {
  namespace RigidConst {
    enum{SINGLE, MOLECULE, GROUP};
    enum{NONE, XYZ, XY, YZ, XZ};
    enum{ISO, ANISO, TRICLINIC};
    enum{FULL_BODY, INITIAL, FINAL, FORCE_TORQUE, VCM_ANGMOM, XCM_MASS, ITENSOR, DOF};
    enum {POINT = 1<<0,
          SPHERE = 1<<1,
          ELLIPSOID = 1<<2,
          LINE = 1<<3,
          TRIANGLE = 1<<4,
          DIPOLE = 1<<5,
          OMEGA = 1<<6,
          ANGMOM = 1<<7,
          TORQUE = 1<<8
    };
    static constexpr double TOLERANCE = 1.0e-6;
    static constexpr double EPSILON = 1.0e-7;
    static constexpr double BIG = 1.0e20;
    static constexpr double SINERTIA = 0.4;
    static constexpr double EINERTIA = 0.2;
    static constexpr double LINERTIA = 1.0/12.0;
    static constexpr int MAXLINE = 1024;
    static constexpr int CHUNK = 1024;
    static constexpr int DELTA_BODY = 10000;
    static constexpr int ATTRIBUTE_PERBODY = 20;
  }
}
#endif
