#pragma once

#include "medusa_random.h"

#include <array>
#include <string>
#include <vector>

namespace medusa {

using namespace ::std;

#define MaxFloat 1.0e+38   /* Good enough for most applications */
#define MAXLINT 2147483647 /* Assuming 32 bit 2's complement */
#define MAXSINT 32767      /* Assuming 16 bit 2's complement */
#define GAMMA 0.57721566490153286060
#define DEG 57.29577951308232087680
#define PHI 1.61803398874989484820
#define ENAT 2.71828182845904523536
#define LOG2E 1.4426950408889634074
#define LOG10E 0.43429448190325182765
#define LN2 0.69314718055994530942
#define LN10 2.30258509299404568402
#define SQRT2 1.41421356237309504880
#define SQRT1_2 0.70710678118654752440

const static double INF = 1.0e+38;
const static double PI = 4.0 * atan(1.0);
const static double rpi = PI / 180.0;
const static double DEG2GRAD = PI / 180.0;
const static double GRAD2DEG = 1 / DEG2GRAD;
const static double Diheral_Angle = 2.0 * atan(sqrt(2.0));

// -- PHYSICAL CONSTANTS --

const double e_MKS = 1.60219e-19; // electron charge (coulomb)
const double e_CGS = 4.80324e-10; // electron charge (esu)

const double eV_MKS = 1.60219e-19; // Electron volt (J)
const double eV_CGS = 1.60219e-12; // Electron volt (erg)

const double m_MKS = 9.1095e-31; // Electron rest mass (kg)
const double m_CGS = 9.1095e-28; // Electron rest mass (gr)

const double h_MKS = 6.6262e-34; // Planck's constant (J*s)
const double h_CGS = 6.6262e-27; // Planck's constant (erg*s)
const double h_eV = 4.1357e-15;  // Planck's constant (eV*s)

const double hbar_MKS = 1.05459e-34; // Planck's constant (J*s)
const double hbar_CGS = 1.05459e-27; // Planck's constant (erg*s)
const double hbar_eV = 6.5822e-16;   // Planck's constant (eV*s)

const double k_MKS = 1.380658e-23; // Boltzmann constant ( J/K)
const double k_eV = 8.617318e-05;  // Boltzmann constant (eV/K)

const double DBL2 = DBL_MAX * 10e-10;
const double DBL1 = DBL_MAX * 0.5e-10;

constexpr double INF = 1.0e+38;

class hbond_acceptor;

class vec {
public:
  double x;
  double y;
  double z;
  // default constructor
  vec(double ix = 0, double iy = 0, double iz = 0) {
    x = ix;
    y = iy;
    z = iz;
  }
  // clone
  vec(const vec &old) {
    x = old.x;
    y = old.y;
    z = old.z;
  }

  ~vec() {
    // cout << "delete vector" << endl;
  }

  // operator overlading
  friend const void operator<<(ostream &ostr, vec &a) {
    ostr << a.getX() << " " << a.getY() << " " << a.getZ() << endl;
  }

  const double &operator[](int i) const {
    switch (i) {
    case 0:
      return x;
    case 1:
      return y;
    case 2:
      return z;
    default:
      std::cerr << "out of range" << endl;
      exit(1);
    }
  }

  double &operator[](int i) {
    switch (i) {
    case 0:
      return x;
    case 1:
      return y;
    case 2:
      return z;
    default:
      std::cerr << "out of range" << endl;
      exit(1);
    }
  }

  friend const vec operator+(const vec &left, const vec &right);

  friend const vec operator-(const vec &left, const vec &right);

  friend const vec operator*(const vec &left, const double &right);

  friend const vec operator*(const double &left, const vec &right);

  friend const double operator*(const vec &left,
                                const vec &right); // dot product

  friend const vec operator^(const vec &left,
                             const vec &right); // cross product

  friend const vec operator/(const vec &left, const double &right);

  vec &operator=(const vec &right) {
    if (this == &right)
      return *this;
    x = right.x;
    y = right.y;
    z = right.z;
    return *this;
  }

  vec &operator+=(const vec &right) {
    x += right.x;
    y += right.y;
    z += right.z;
    return *this;
  }

  vec &operator-=(const vec &right) {
    x -= right.x;
    y -= right.y;
    z -= right.z;
    return *this;
  }

  vec &operator*=(double right) {
    x *= right;
    y *= right;
    z *= right;
    return *this;
  }

  vec &operator/=(double right) {
    x /= right;
    y /= right;
    z /= right;
    return *this;
  }

  int operator==(vec &right) const {
    if (x == right.x && y == right.y && z == right.z)
      return 1;
    else
      return 0;
  }

  void print(ostream &out) const {
    out << "(" << x << ", " << y << ", " << z << ")";
  }

  double getX() const { return x; }

  double getY() const { return y; }

  double getZ() const { return z; }

  void setX(double ix) { x = ix; }

  void setY(double iy) { y = iy; }

  void setZ(double iz) { z = iz; }

  void addX(double ix) { x += ix; }

  void minusX(double ix) { x -= ix; }

  void addY(double iy) { y += iy; }

  void minusY(double iy) { y -= iy; }

  void addZ(double iz) { z += iz; }

  void minusZ(double iz) { z -= iz; }

  double getDist(const vec &a) const {
    double tmp = 0;
    tmp += (x - a.x) * (x - a.x);
    tmp += (y - a.y) * (y - a.y);
    tmp += (z - a.z) * (z - a.z);
    return sqrt(tmp);
  }

  // get the Magnitude
  double getMag() const { return sqrt(x * x + y * y + z * z); }

  // to normalize the vector
  void Normalize() {
    double tmp = sqrt(x * x + y * y + z * z);
    if (tmp > 0) {
      x /= tmp;
      y /= tmp;
      z /= tmp;
    }
  }

  // return a normalized vector without changed the original vector
  vec norm() const {
    double tmp = sqrt(x * x + y * y + z * z);
    if (tmp > 0) {
      return vec(x / tmp, y / tmp, z / tmp);
    } else {
      return vec(INF, INF, INF);
    }
  }

  // get magnitude square
  double mod2() const { return x * x + y * y + z * z; }

  void Rotate(const vec &p1, const vec &p2, double theta);

  void EularRotate(double **m);

  // do DimondRotate
  void DimondRotate(const vec &p1, const vec &p2, double theta);
};

// operator overlading
inline const vec operator+(const vec &left, const vec &right) {
  return vec(left.x + right.x, left.y + right.y, left.z + right.z);
}

inline const vec operator-(const vec &left, const vec &right) {
  return vec(left.x - right.x, left.y - right.y, left.z - right.z);
}

inline const vec operator*(const vec &left, const double &right) {
  return vec(left.x * right, left.y * right, left.z * right);
}

inline const vec operator*(const double &left, const vec &right) {
  return vec(left * right.x, left * right.y, left * right.z);
}

inline const double operator*(const vec &left,
                              const vec &right) { // dot product
  return (left.x * right.x + left.y * right.y + left.z * right.z);
}

inline const vec operator^(const vec &left,
                           const vec &right) { // cross product
  return vec(left.y * right.z - left.z * right.y,
             left.z * right.x - left.x * right.z,
             left.x * right.y - left.y * right.x);
}

inline const vec operator/(const vec &left, const double &right) {
  return vec(left.x / right, left.y / right, left.z / right);
}

const static int NPOLY = 16;
const static double internal_weight = 0.5E0;

class poly {
private:
  int degree;
  double offset;
  double scale;
  double normalize;
  double coeff[8];

public:
  poly(int d, double o, double s, double c1, double c2, double c3, double c4,
       double c5, double c6, double c7, double c8)
      : degree(d), offset(o), scale(s) {
    coeff[0] = c1;
    coeff[1] = c2;
    coeff[2] = c3;
    coeff[3] = c4;
    coeff[4] = c5;
    coeff[5] = c6;
    coeff[6] = c7;
    coeff[7] = c8;
    normalize = internal_weight / scale;
  }

  ~poly() {}

  double getValue(double x) const {
    double v = 0;
    for (int i = degree - 1; i >= 0; i--) {
      v = coeff[i] + v * x;
    }
    v -= offset;
    v *= normalize;
    return v;
  }

  double getDerivative(double x) const {
    double v = 0;
    for (int i = degree - 2; i >= 0; i--) {
      v = (i + 1) * coeff[i + 1] + v * x;
    }
    v *= normalize;
    return v;
  }
};

const static poly POLYS[] = {
    poly(8,                    // degree
         9.939481771384827E0,  // offset
         1.949735682035433E0,  // scale
         -20897.2963242561E0,  // 0
         64247.6438427233E0,   // 1
         -83351.2088271679E0,  // 2
         59275.3348525644E0,   // 3
         -24993.0213101615E0,  // 4
         6256.27794663695E0,   // 5
         -861.843894662537E0,  // 6
         50.4501160458263E0),  // 7
    poly(8,                    // degree
         10.19618640787044E0,  // offset
         2.356555305029287E0,  // scale
         -20548.3707305135E0,  // 0
         65573.1420852652E0,   // 1
         -88281.0899077548E0,  // 2
         65135.8946527841E0,   // 3
         -28487.2557887243E0,  // 4
         7394.5552104088E0,    // 5
         -1055.95001101041E0,  // 6
         64.0502420663362E0),  // 7
    poly(8,                    // degree
         3.95898E0,            // offset
         1.55E0,               // scale
         887.47797126217E0,    // 0
         -1585.99448214285E0,  // 1
         1057.0173554807E0,    // 2
         -310.744439916284E0,  // 3
         34.0605598761848E0,   // 4
         0E0,                  // 5
         0E0,                  // 6
         0E0),                 // 7
    poly(8,                    // degree
         0.17E0,               // offset
         1E0,                  // scale
         477.645791505638E0,   // 0
         -806.71862668605E0,   // 1
         503.64186195361E0,    // 2
         -137.999376572433E0,  // 3
         14.0234710763302E0,   // 4
         0E0,                  // 5
         0E0,                  // 6
         0E0),                 // 7
    poly(8,                    // degree
         11.0144095070536E0,   // offset
         3.88539106080337E0,   // scale
         13.1761011064371E0,   // 0
         6.67268543875679E0,   // 1
         -294.055688032283E0,  // 2
         1866.10364337877E0,   // 3
         -5355.8036552271E0,   // 4
         7923.03899964551E0,   // 5
         -5888.11906511642E0,  // 6
         1736.97818014834E0),  // 7
    poly(8,                    // degree
         9.62421538318229E0,   // offset
         2.75547929573618E0,   // scale
         13.127927303063E0,    // 0
         0.822868665987038E0,  // 1
         -109.910665370086E0,  // 2
         686.786517638211E0,   // 3
         -1942.11269671157E0,  // 4
         2836.84871202931E0,   // 5
         -2095.85402633764E0,  // 6
         617.160098870171E0),  // 7
    poly(3,                    // degree
         7.15765555313026E0,   // offset
         3.1773184908633E0,    // scale
         12.3002419354918E0,   // 0
         -4.7448972963389E0,   // 1
         -3.57500757688594E0,  // 2
         0E0,                  // 3
         0E0,                  // 4
         0E0,                  // 5
         0E0,                  // 6
         0E0),                 // 7
    poly(3,                    // degree
         6.573635101318359E0,  // offset
         8.319905042648315E0,  // scale
         6.57363492607565E0,   // 0
         -4.23172288336822E0,  // 1
         1.58016437009784E0,   // 2
         0E0,                  // 3
         0E0,                  // 4
         0E0,                  // 5
         0E0,                  // 6
         0E0),                 // 7
    poly(3,                    // degree
         7.56221962058537E0,   // offset
         3.1445927360857E0,    // scale
         12.2116835353074E0,   // 0
         -7.16425671041666E0,  // 1
         -0.629799940391145E0, // 2
         0E0,                  // 3
         0E0,                  // 4
         0E0,                  // 5
         0E0,                  // 6
         0E0),                 // 7
    poly(3,                    // degree
         6.787431716918945E0,  // offset
         7.794056534767151E0,  // scale
         6.78743191739409E0,   // 0
         -3.2345923209685E0,   // 1
         1.39326898971812E0,   // 2
         0E0,                  // 3
         0E0,                  // 4
         0E0,                  // 5
         0E0,                  // 6
         0E0),                 // 7
    poly(8,                    // degree
         11.926299544787E0,    // offset
         5.85119894086317E0,   // scale
         11.9350231894526E0,   // 0
         18.409859118441E0,    // 1
         -241.176866650372E0,  // 2
         1135.03014758504E0,   // 3
         -2716.19210529742E0,  // 4
         3463.22973054283E0,   // 5
         -2303.26203622031E0,  // 6
         641.36914663319E0),   // 7
    poly(8,                    // degree
         8.09262412710125E0,   // offset
         2.16379301251132E0,   // scale
         11.5196524630225E0,   // 0
         -0.982514077302885E0, // 1
         -38.7140999410775E0,  // 2
         204.61420850299E0,    // 3
         -573.17550726292E0,   // 4
         837.762308008926E0,   // 5
         -624.832463761649E0,  // 6
         190.175068834576E0),  // 7
    poly(3,                    // degree
         7.482917785644531E0,  // offset
         3.084799265279046E0,  // scale
         7.48291779029918E0,   // 0
         -11.4461758882338E0,  // 1
         10.6177843320035E0,   // 2
         0E0,                  // 3
         0E0,                  // 4
         0E0,                  // 5
         0E0,                  // 6
         0E0),                 // 7
    poly(3,                    // degree
         4.760147571563721E0,  // offset
         3.084799265279046E0,  // scale
         4.76014765119106E0,   // 0
         -2.74345435413267E0,  // 1
         2.71214471299479E0,   // 2
         0E0,                  // 3
         0E0,                  // 4
         0E0,                  // 5
         0E0,                  // 6
         0E0),                 // 7
    poly(3,                    // degree
         5.085277557373047E0,  // offset
         1.831269767104503E0,  // scale
         5.08527774066672E0,   // 0
         -6.87443803002661E0,  // 1
         6.45151885957664E0,   // 2
         0E0,                  // 3
         0E0,                  // 4
         0E0,                  // 5
         0E0,                  // 6
         0E0),                 // 7
    poly(8,                    // degree
         -0.5667E0,            // offset
         1.000000000000000E0,  // scale
         2.582609E0,           // 0
         -27.898499E0,         // 1
         184.031158E0,         // 2
         -497.581482E0,        // 3
         490.462250E0,         // 4
         74.764015E0,          // 5
         -412.507416E0,        // 6
         184.580704E0)         // 7
};

#define ABS(x) (((x) < 0) ? (-(x)) : (x))

class DCLPDB {

  char pfile[100]; // PDB file
  ifstream in;     // stream

  class atom {
  public:
    vec r;        // position
    char name[6]; // actual name
    int CA;       // CA flag, default = 0 (not CA)
    int CB;       // CB flag, default = 0 (not CB)
    int heavy;    // heavy atom (0 if H, 1 otherwise)
    double T;     // temperature factor

    atom() {}
    ~atom() {}
  };

  class residue {
  public:
    char name[100]; // name of an amino acid
    int na;         // number of atoms
    atom *a;        // pointer to an array of atoms belonging to aa
    int pCA;        // position of CA within the aa
    int pCB;        // position of CB within the aa

    residue(){};
    ~residue() { delete[] a; };
  };

  class helix { // structure defines the alpha-helices
  public:
    int start;
    int end;
    int len;

    helix() {}
    ~helix() {}
  };

  class sheet { // structure defines the beta-sheets
  public:
    int start;
    int end;
    int len;
    int sence; // sence =1,-1,0 if parallel,antiparallel,1st strand

    sheet() {}
    ~sheet() {}
  };

  class turn { // structure defines the turns
  public:
    int start;
    int end;
    int len;

    turn() {}
    ~turn() {}
  };

  class ssbond { // structure defines the -S-S- bonds
  public:
    int cys1;
    int cys2; // positions of Cysteins

    ssbond() {}
    ~ssbond() {}
  };

public:
  int N;        // Number of residues
  int Na;       // Number of atoms
  int ipos;     // initial number of the first residue
  residue *res; // residues
  int **CM;     // Contact matrix
  vec COM;      // Center of Mass vector
  char **Seq3;  // Sequence in a 3-letter representation
  char *Seq;    // Sequence in a 1-letter representation
  short *SeqID; // integer identifier of a given sequence

  // Secondary strucutre elements
  int Nh;      // Number of helices
  int Ns;      // Number of sheets
  int Nt;      // Number of turns
  int Nss;     // Number of -S-S- bonds
  helix *ahx;  // allocate helix structure
  sheet *bsh;  // allocate sheet structure
  turn *trn;   // allocate turn  structure
  ssbond *ssb; // allocate -S-S- bond  structure

  // chain
  char CHAIN; // chain ID

  DCLPDB();
  DCLPDB(char *);
  DCLPDB(char *, char);
  ~DCLPDB();

  void Init(char *afname) { strcpy(pfile, afname); };
  void InitRes(int, int *);         // initializes the residues and atoms
  int GetNR(void);                  // counts number of residues
  void Read(void);                  // reads the PDB file
  void ReadChain(void);             // reads the PDB file but only first chain
  int ReadStructure2(void);         // reads secondary structure from PDB file
  double GetHelixFraction(void);    // computes the fraction of helices
  double GetSheetFraction(void);    // computes the fraction of beta-sheets
  void Write(ostream &ostr = cout); // writes to a stream, default=cout
  void
  WriteCA(ostream &ostr = cout); // writes only CA to a stream, default=cout
  void
  WriteCB(ostream &ostr = cout);   // writes only CB to a stream, default=cout
  int GetAtomNC(int, int, double); // returns number of contacts for atoms
  int GetCANC(double);             // returns number of contacts for CAs
  int GetCANC(int, double);        // returns number of contacts for CAs
  int GetCBNC(int, double);        // returns number of contacts for CBs
  int GetCBNC(double);             // returns number of contacts for CBs
  int GetNC(int, double);          // returns number of contacts for all atoms
  int GetNCncn(int, double);       // returns number of contacts for all
                                   // atoms, excluding chain i,i+1 nighbours
  int GetHeavyNC(int, double);     // returns number of contacts for all
                                   // heavy atoms

  void GetCAContMap(double);           // returns contact map based on CAs
  void GetCBContMap(double);           // returns contact map based on CBs
  double GetHeavyContactOrder(double); // returns contact order \cite{Plaxco98}
                                       // if any of heavy atoms of aa in the
                                       // specific range
  void GetCACenterMass(void); // returns xyz of the center of mass based on CA
  void GetCBCenterMass(void); // returns xyz of the center of mass based on CB
  void Shift2CACenterMass(void); // shifts to the center of mass based on CA
  void Shift2CBCenterMass(void); // shifts to the center of mass based on CB
  double D(int, int, int, int);  // returns the distance between two atoms
  double D2(int, int, int,
            int);               // returns the distance square between two atoms
  double DCA(int, int);         // returns the distance between two CAs
  double DCA2(int, int);        // returns the distance square between two CAs
  double DCB(int, int);         // returns the distance between two CBs
  double DCB2(int, int);        // returns the distance square between two CBs
  double End2EndCA(void);       // returns end-to-end distance between CAs
  double End2EndCB(void);       // returns end-to-end distance between CAs
  double RgCA(void);            // returns Rg based on CAs
  double RgCB(void);            // returns Rg based on CBs
  double Dist2Center(int, int); // returns distance to the center of an atom
  double aveCADist(int);        // returns the ave length between a given aa
  double aveCBDist(int);        // and the rest of aa in CA and CB represent.
  double aveCADist(void);       // returns the ave length between all aa
  double aveCBDist(void);       // in CA and CB representations
  double CaCb_CaCbAngle(int, int); // returns cos(CaCb_i, CaCb_j)

  void Reset(void);          // resets all the values
  void SetN(int);            // sets N
  void SetNa(int);           // sets Na
  void MemAllocateRes(void); // Allocates memory for residues N has been set
  void MemAllocateRes(int);  // Allocates memory for residues with argv=N
  void MemAllocateRes(int, int *); // Allocates memory for residues with
                                   // argv[1]=N and argv[2]=*na
  void MemAllocateAtoms(int *);    // Allocates memory for atoms with argv=*na
  void SetpCA(int *v);             // Sets the index of CA atom
  void SetpCB(int *v);             // Sets the index of CB atom
  void SetDefaultName(void);    // Sets default name for the amino acids (GLY)
  void SetDefaultName(char[]);  // Sets default name for the amino acids to argv
  void SetName(char[], int);    // Sets name for an amino acid
  void SetDefaultAtomName(int); // Sets default name for ith atoms (CA)
  void SetDefaultAtomName(int, char[]); // Sets default name for ith atoms
  void SetAtomName(int, int, char[]);   // Sets default name for ith atoms of aa

  void GetSeq(void);             // initializes a seq string
  double Dist2(vec *);           // returns the square of vector size
  double Dist2(vec *, vec *);    // returns the square of diff between 2 vectors
  double ScalProd(vec *, vec *); // returns scalar product of two vectors
  double VecAngle(vec *, vec *); // returns the angle between two vectors
  char Convert3to1(char[]); // returns the name of an amino acid in 1 letter
                            // representation
  int IsHeavy(char *);      // checks if atom is H (return=0) or not (1)

  void ClearMem(void); // clears memory

private:
  void squeeze(char[], int);
  void clean(char s[], int n) {
    for (int i = n; i < strlen(s); i++)
      s[i] = '\0';
  };
  void clear(char s[]) {
    for (int i = 0; i < strlen(s); i++)
      s[i] = '\0';
  };
  void clear(char s[], int n) {
    for (int i = 0; i < n; i++)
      s[i] = '\0';
  };
  void Error(char *msg) {
    cerr << msg << endl << flush;
    exit(1);
  };
  char tmp[100];
  char aa1[20], aa3[20][4];
};

// define atom types (see Tsai et al, JMB 290, 253--266 (1999))
typedef enum {
  C3H0 = 0,
  C3H1,
  C4H1,
  C4H2,
  C4H3,
  N3H0,
  N3H1,
  N3H2,
  N4H3,
  O1H0,
  O1H1,
  S2H0,
  S2H1,
  MOL2, // added for all small molecule atoms
  H,
} ATMtype;

// Van der Waals radii for various atom types (in A)
const static double VDWR[] = {
    1.61, //_C3H0
    1.76, //_C3H1
    1.88, //_C4H1
    1.88, //_C4H2
    1.88, //_C4H3
    1.64, //_N3H0
    1.64, //_N3H1
    1.64, //_N3H2
    1.64, //_N4H3
    1.42, //_O1H0
    1.46, //_O1H1
    1.77, //_S2H0
    1.77, //_S2H1
    0.00, // added for all small molecule atoms
    1.00, //_H
};

typedef enum {
  _N_ = 0, // 1 backbone
  _CA_,    // 2 backbone
  _C_,     // 3 backbone
  _O_,     // 4 backbone
  _CB_,    // 5 CB
  _CG_,    // 6
  _CG1_,   // 7
  _CG2_,   // 8
  _CD_,    // 9
  _CD1_,   // 10
  _CD2_,   // 11
  _CE_,    // 12
  _CE1_,   // 13
  _CE2_,   // 14
  _CE3_,   // 15
  _CH2_,   // 16
  _CZ_,    // 17
  _CZ2_,   // 18
  _CZ3_,   // 19
  _SD_,    // 20
  _SG_,    // 21
  _ND1_,   // 22
  _ND2_,   // 23
  _NE_,    // 24
  _NE1_,   // 25
  _NE2_,   // 26
  _NH1_,   // 27
  _NH2_,   // 28
  _NZ_,    // 29
  _OD1_,   // 30
  _OD2_,   // 31
  _OE1_,   // 32
  _OE2_,   // 33
  _OH_,    // 34
  _OG_,    // 35
  _OG1_,   // 36

  // Hydrogens
  // CYS
  _HG1_,
  // HIS
  _HD1_,
  // TRP
  _HE1_,
  // lys
  _HZ1_,
  _HZ2_,
  _HZ3_,
  // arg
  _HE_,
  _HH11_,
  _HH12_,
  _HH21_,
  _HH22_,
  // backbone
  _HN_,
  // ASN
  _HND1_,
  _HND2_,
  // GLN
  _HNE1_,
  _HNE2_,
  // SER,THR,TYR
  _HO_,
  // added for small molecule atoms:
  _H_,
  _S_,
  _P_,

  _F_,
  _Cl_,
  _Br_,
  _I_,

  _Li_,
  _Na_,
  _K_,
  _Mg_,
  _Ca_,
  _Al_,
  _Si_,
  _Sn_,
  _Se_,

  _Cr_,
  _Mo_,
  _Mn_,
  _Fe_,
  _Co_,
  _Cu_,
  _Zn_,
  _Eh_,
  _Aum_,
  _Ce_,

  _LP_,  // lone pair            SPECIAL SYMBOLS
  _Du_,  // dummy
  _DuC_, // dummy with the weight of carbon
  _ANY_, // any atom
  _HEV_, // heavy (non H) atom
  _HET_, // heteroatom (N, O, S, P)
  _HAL_, // halogen

  // Backbone N terminal hydrogens [ Added by Jian 11/2/18 ]
  _HN1_,
  _HN2_,
  _HN3_,
  _HD2_,
  _HE2_
} PDBAtype;

const static char *PDBAname[] = {
    "N", "CA", "C", "O", "CB", "CG", "CG1", "CG2", "CD", "CD1", "CD2", "CE",
    "CE1", "CE2", "CE3", "CH2", "CZ", "CZ2", "CZ3", "SD", "SG", "ND1", "ND2",
    "NE", "NE1", "NE2", "NH1", "NH2", "NZ", "OD1", "OD2", "OE1", "OE2", "OH",
    "OG", "OG1",

    // Hydrogens
    // CYS
    "HG1",
    // HIS
    "HD1",
    // TRP
    "HE1",
    // lys
    "HZ1", "HZ2", "HZ3",
    // arg
    "HE", "HH11", "HH12", "HH21", "HH22",
    // backbone
    "HN",
    // ASN
    "HND1", "HND2",
    // GLN
    "HNE1", "HNE2",
    // SER,THR,TYR
    "HO",
    // added for small molecule atoms
    "H", "S", "P",

    "F", "Cl", "Br", "I",

    "Li", "Na", "K", "Mg", "Ca", "Al", "Si", "Sn", "Se",

    "Cr", "Mo", "Mn", "Fe", "Co", "Cu", "Zn", "Eh", "Aum", "Ce",

    "LP",  // lone pair            SPECIAL SYMBOLS
    "Du",  // dummy
    "DuC", // dummy with the weight of carbon
    "ANY", // any atom
    "HEV", // heavy (non H) atom
    "HET", // heteroatom (N, O, S, P)
    "HAL", // halogen

    // Backbone N terminal Hydrogens
    "HN1", "HN2", "HN3", "HD2", "HE2"};

// returns the atom id based on its PDA name
int atomid(const std::string &name);

const char *PDBAtype2Name(PDBAtype t);

double radius(ATMtype type);

//--------------------------------------------------------------------
//
//  A T O M   C L A S S
//
class atom {
public:
  vec r; // position
  vec p; // old position
  vec gme_r;
  vec origin_r; // original pos

  ATMtype ATOMtype; // atom type (see above)
  PDBAtype PDBtype; // PDBtype

  short topAtomIndex;

  int FF_t;      // ForceField type
  double mass;   // MASS
  double charge; // partial charge

  bool _isCharged;
  int icharge; // the integer charge, simplied to account for the
  // salt-bridge interaction between charged residues,
  // such as ARG(CZ+),LYS(NZ1) and ASP(CG-),GLU(CD-)

  atom *attachedHA;      // attached heavy atom; NULL for non-proton atom
  int n_attachedProtons; // number of attached protons;
  atom *attachedProtons; // attached protons;

  // used for cell indexing
  int cell_index;
  atom *cell_next;
  atom *cell_prev;
  int resID;
  int subResID; // residue id counted within each chain
  int chainID;  // which chain the atom belongs to in the whole complex
  int index;    // atom index, numbered within the complex

  // hydrogen bonding parameter
  short isSideChain;
  hbond_acceptor *acceptor;         // for acceptors; NULL for the non-acceptor;
  hbond_acceptor *hbonded_acceptor; // for protons; NULL for no-hbond partner
  int fine_cell_index;
  atom *fine_cell_next;
  atom *fine_cell_prev;

  // assign the pdb name
  string _name;
  void name(string n) { _name = n; }
  string &name() { return _name; }

  // functions or METHODS:
  atom(double x = 0, double y = 0, double z = 0) {
    r = vec(0, 0, 0);
    topAtomIndex = -1;
    attachedHA = NULL;
    n_attachedProtons = 0;
    attachedProtons = NULL;
    cell_next = NULL;
    cell_prev = NULL;
    fine_cell_next = NULL;
    fine_cell_prev = NULL;
    acceptor = NULL;
    hbonded_acceptor = NULL;
    isSideChain = -1;
    icharge = 0;
    _isCharged = false;
  }

  atom(const char *PDBname) {
    r = vec(0, 0, 0);
    topAtomIndex = -1;
    PDBtype = static_cast<PDBAtype>(atomid(PDBname));
    attachedHA = NULL;
    n_attachedProtons = 0;
    attachedProtons = NULL;
    cell_next = NULL;
    cell_prev = NULL;
    fine_cell_next = NULL;
    fine_cell_prev = NULL;
    acceptor = NULL;
    hbonded_acceptor = NULL;
    isSideChain = -1;
    icharge = 0;
    _isCharged = false;
  }

  atom(const ATMtype type) {
    r = vec(0, 0, 0);
    topAtomIndex = -1;
    ATOMtype = type;
    attachedHA = NULL;
    n_attachedProtons = 0;
    attachedProtons = NULL;
    cell_next = NULL;
    cell_prev = NULL;
    fine_cell_next = NULL;
    fine_cell_prev = NULL;
    acceptor = NULL;
    hbonded_acceptor = NULL;
    isSideChain = -1;
    icharge = 0;
    _isCharged = false;
  }

  atom(const char *PDBname, vec &p) {
    r = p;
    topAtomIndex = -1;
    PDBtype = static_cast<PDBAtype>(atomid(PDBname));
    attachedHA = NULL;
    n_attachedProtons = 0;
    attachedProtons = NULL;
    cell_next = NULL;
    cell_prev = NULL;
    fine_cell_next = NULL;
    fine_cell_prev = NULL;
    acceptor = NULL;
    hbonded_acceptor = NULL;
    isSideChain = -1;
    icharge = 0;
    _isCharged = false;
  }

  atom(const atom &right) { *this = right; }

  /*
     USE with extreme care!!!
     this equal operator is not a compelte imprelemtnation!
     */
  atom &operator=(const atom &right) {
    if (this == &right)
      return *this;
    //    if (attachedProtons) delete[] attachedProtons;
    if (n_attachedProtons > 1)
      delete[] attachedProtons;
    else if (attachedProtons)
      delete attachedProtons;
    r = right.r;
    p = right.p;
    topAtomIndex = right.topAtomIndex;
    ATOMtype = right.ATOMtype;
    PDBtype = right.PDBtype;
    attachedHA = NULL;
    n_attachedProtons = 0;
    attachedProtons = NULL;
    return *this;
  }

  void init(const ATMtype type) { ATOMtype = type; }
  void init(const ATMtype atm, const PDBAtype pdb, const vec &a) {
    ATOMtype = atm;
    PDBtype = pdb;
    r = a;
  }

  void setATMtype(ATMtype type) { ATOMtype = type; }

  void setPDBAtype(PDBAtype type) { PDBtype = type; }

  int getIndex() { return index; }

  /***************************
    Here the destructor assumes that where only one proton attached,
    the attachedProtons is crated by:
    attachedProtons = new atom;
    instead of
    attachedProtons = new atom[1];
   ***************************/
  ~atom() { removeHydrogens(); }

  void removeHydrogens() {
    if (n_attachedProtons > 1) {
      delete[] attachedProtons;
    } else if (attachedProtons) {
      delete attachedProtons;
    }
  }

  vec *getR() { return &r; }

  void setR(const vec &a) { r = a; }

  void setR(double x, double y, double z) { r = vec(x, y, z); }

  void setMass(double im) { mass = im; }

  double getMass() { return mass; }

  void setCharge(double iq) { charge = iq; }

  double getCharge() { return charge; }

  void setICharge(int iq) { icharge = iq; }
  int getICharge() { return icharge; }
  bool isCharged() { return _isCharged; }
  void charged() { _isCharged = true; }

  void setFFtype(int type) { FF_t = type; }

  int getFFtype() { return FF_t; }

  void savePosition() { p = r; }

  void restorPosition() { r = p; }

  void deprotonate() { n_attachedProtons--; }

  const char *getName() { return PDBAname[PDBtype]; }

  void updateSC_TYPE() {
    if (PDBtype == _N_ || PDBtype == _CA_ || PDBtype == _C_ || PDBtype == _O_ ||
        PDBtype == _HN_ || PDBtype == _HN1_ || PDBtype == _HN2_ ||
        PDBtype == _CB_) {
      isSideChain = 0;
    } else
      isSideChain = 1;
  }

  short isSC() { return isSideChain; }

  int clearHBAcceptor() {
    if (attachedHA) { // a proton
      if (hbonded_acceptor) {
        hbonded_acceptor = NULL;
        return 1;
      } else
        return 0;
    }
    return 0;
  }

  int isHBonded() {
    if (hbonded_acceptor)
      return 1;
    return 0;
  }

  int isH() {
    // if(PDBtype>=_HG1_ && PDBtype<_S_) return 1;
    if (attachedHA)
      return 1;
    return 0;
  }
};

typedef enum {
  NATURAL_AA,
  GEN_NATURAL_AA,
  OTHER_AA,
  GEN_OTHER_AA,
  SMALL_MOLECULE,
  NATURAL_RNA
} BaseResidueType;

// for all types

class BaseResidue {
private:
  /// residue id given in the pdb
  int _resid;
  /// residue insertion flag in pdb
  char _resins;

public:
  /// Save the current position of each atom.
  /// Recover the previously save positions for each atom.

  int getResid() { return _resid; };
  char getResins() { return _resins; };
  void setResid(int resid) { _resid = resid; };
  void setResins(char resins) { _resins = resins; };
  // interface methods
  virtual int getNA() = 0;
  virtual int
  getNH() = 0; // count number of attached protons, should be same as getNPH();
  virtual int getNChi() = 0;   // number of SC chi angles
  virtual int getNBBChi() = 0; // number of BB dihedral angles
  virtual atom *getAtom(int i) = 0;
  // get atoms is not required now to allows vector-based contaners
  // virtual atom* getAtoms() = 0;
  virtual const char *getName() = 0;
  virtual atom *getAtom(const string &pdbname, const topology &top) = 0;

  // set the array of polar hydrogens
  virtual int getNPH() = 0;
  virtual atom **getPH_ARRAY() = 0;

  virtual int getNHBAcceptor() = 0;
  virtual atom **getHBAcceptor_ARRAY() = 0;
  virtual void getBoundary(vec &min, vec &max) = 0;
  virtual void print(ostream &out) = 0;
  virtual BaseResidueType getBaseResidueType() = 0;
  // by doing updateTopology, every residue should have the fftype, mass and
  // index assigned and ready for energy evaluation
  virtual void updateTopology(const topology &top) = 0;
  virtual void saveGME_Conf() = 0;
  virtual void recoverGME_Conf() = 0;

  void savePos() {
    for (int ia = 0; ia < this->getNA(); ia++)
      this->getAtom(ia)->savePosition();
    for (int ia = 0; ia < this->getNPH(); ia++)
      this->getPH_ARRAY()[ia]->savePosition();
  }

  void recoverPos() {
    for (int ia = 0; ia < this->getNA(); ia++)
      this->getAtom(ia)->restorPosition();
    for (int ia = 0; ia < this->getNPH(); ia++)
      this->getPH_ARRAY()[ia]->restorPosition();
  }
};

typedef enum { SP2 = 0, SP3, RING } ACCEPTOR_HYBRID_t;
typedef enum { BB_HELIX = 0, BB_SHEET, SC_SP2, SC_SP3, SC_RING } HB_STAT_t;
typedef enum { allRHA = 0, shortAHD, shortXAH, longAHD, longXAH } HB_TERM_t;
// variable used to seperate angular dependence for sc hbond
const static double RHA_cutoff = 2.10;
// for SC only, if the HA distance are locate close the cut-off range
// so to say, plus/minus a small range, the interpolation will happen
const static double RHA_INTERP_RANGE = 0.20;
const static double RHA_INTERP_MAX = RHA_cutoff + RHA_INTERP_RANGE;
const static double RHA_INTERP_MIN = RHA_cutoff - RHA_INTERP_RANGE;
// the distance, cos(angles) interpolation edge
const static double RHA_INTERP_EDGE = 2.10;
const static double ANGLE_INTERP_EDGE = 0.05;
// max and min distance, cos(angels) to define hbond
const static double MAX_R = 3.0;
const static double MAX_R2 = MAX_R * MAX_R;
const static double MAX_ROTH_R = MAX_R + 1;
const static double MAX_ROTH_R2 = MAX_ROTH_R * MAX_ROTH_R;
const static double MIN_R = 1.4;
const static double MAX_CAHD = 1.0;
const static double MIN_CAHD = 0.0;
const static double MAX_CXAH = 1.0;
const static double MIN_CXAH = 0.0;
// maxr to defined E(r) as zero
const static double maxr[] = {2.800000E0, 2.745000E0, 2.5000E0, 2.500000E0,
                              2.500000E0};

const static int POLY_TABLE[][5] = {
    // helix: bb-bb
    0,
    4,
    10,
    -1,
    -1,
    // sheet: bb-bb
    1,
    5,
    11,
    -1,
    -1,
    // sp2:  sc-sc/sc-bb
    2,
    6,
    12,
    7,
    13,
    // sp3:  sc-sc/sc-bb
    3,
    8,
    14,
    9,
    14,
    // ring N:sc-sc/sc-bb
    2,
    6,
    15,
    7,
    15,
};

class hbond_acceptor {
private:
  atom *A; // acceptor atom
  atom *X; // the atom covalently attached to A; non-hydorgen
  atom *Y;
  ACCEPTOR_HYBRID_t hybrid;
  unsigned char MAX_HBONDS;
  atom **Hs;             // linked hboned protons
  double *EHBs;          // the bonding energy
  unsigned char nhbonds; // number of hbonded protons
public:
  hbond_acceptor() {
    A = NULL;
    X = NULL;
    Y = NULL;
    hybrid = static_cast<ACCEPTOR_HYBRID_t>(-1);
    MAX_HBONDS = 0;
    Hs = NULL;
    EHBs = NULL;
    nhbonds = 0;
  }

  ~hbond_acceptor() {
    if (Hs)
      delete[] Hs;
    if (EHBs)
      delete[] EHBs;
  }

  void init(atom *iA, atom *iX, ACCEPTOR_HYBRID_t ihybrid, int max_hbonds) {
    A = iA;
    X = iX;
    hybrid = ihybrid;
    MAX_HBONDS = max_hbonds;
    Hs = new atom *[MAX_HBONDS];
    EHBs = new double[MAX_HBONDS];
    nhbonds = 0;
    for (int i = 0; i < MAX_HBONDS; i++) {
      Hs[i] = NULL;
      EHBs[i] = INF;
    }
  }

  atom **getHs() { return Hs; }
  double *getEHBs() { return EHBs; }

  unsigned char getMAX_HBONDS() { return MAX_HBONDS; }

  unsigned char getNBonds() { return nhbonds; }

  atom *getA() { return A; }

  atom *getX() { return X; }

  void setY(atom *y) { Y = y; }
  atom *getY() { return Y; }

  ACCEPTOR_HYBRID_t getHybrid() { return hybrid; }

  int addH(atom *h, const double &e) {
    if (nhbonds < MAX_HBONDS) {
      for (int i = 0; i < MAX_HBONDS; i++) {
        if (!Hs[i]) {
          nhbonds++;
          Hs[i] = h;
          EHBs[i] = e;
          return 1;
        }
      }
    }
    return 0;
  }

  double getE(atom *h) {
    for (int i = 0; i < MAX_HBONDS; i++) {
      if (Hs[i] == h) {
        return EHBs[i];
      }
    }
    return INF;
  }

  int clearH(atom *h, double &e) {
    for (int i = 0; i < MAX_HBONDS; i++) {
      if (Hs[i] == h) {
        nhbonds--;
        Hs[i]->clearHBAcceptor();
        Hs[i] = NULL;
        e = EHBs[i];
        EHBs[i] = INF;
        return 1;
      }
    }
    return 0;
  }

  int isSaturated() {
    if (nhbonds < MAX_HBONDS)
      return 0;
    return 1;
  }

  bool hasProtonAt(int iRes) {
    for (int i = 0; i < MAX_HBONDS; i++) {
      if (Hs[i]) {
        if (Hs[i]->resID == iRes) {
          return true;
        }
      }
    }
    return false;
  }

  int allowedHBonds() { return MAX_HBONDS - nhbonds; }

  double clearHs(double &E_bb_bb, double &E_bb_sc, double &E_sc_sc) {
    double Ehb = 0;
    for (int i = 0; i < MAX_HBONDS; i++) {
      if (Hs[i]) {
        if (Hs[i]->isSC() && A->isSC())
          E_sc_sc += EHBs[i];
        else if (!Hs[i]->isSC() && !A->isSC())
          E_bb_bb += EHBs[i];
        else
          E_bb_sc += EHBs[i];
        Ehb += EHBs[i];
        Hs[i]->clearHBAcceptor();
        Hs[i] = NULL;
      }
    }
    nhbonds = 0;
    return Ehb;
  }

  void currentHBE(double &hbe_bb_bb, double &hbe_bb_sc, double &hbe_sc_sc) {
    hbe_bb_bb = hbe_bb_sc = hbe_sc_sc = 0;
    int n = 0;
    for (int i = 0; i < MAX_HBONDS; i++) {
      if (Hs[i]) {
        if (Hs[i]->isSC() && A->isSC())
          hbe_sc_sc += EHBs[i];
        else if (!Hs[i]->isSC() && !A->isSC()) {
          // cout << Hs[i]->resID << " " << A->resID << endl;
          hbe_bb_bb += EHBs[i];
        } else
          hbe_bb_sc += EHBs[i];
        n++;
      }
    }
    if (nhbonds != n) {
      cerr << "fatal error: mismatch of hbonds" << endl;
      exit(1);
    }
  }
};

class rotRecord;
class resRotLib;
class BBDep_rotRecord;
class BBDep_resRotLib;

// define the residue name, type and interconversion
// the NATURAL amino acid
typedef enum {
  CYS = 0,
  MET,
  PHE,
  ILE,
  LEU,
  VAL,
  TRP,
  TYR,
  ALA,
  GLY,
  THR,
  SER,
  GLN,
  ASN,
  GLU,
  ASP,
  HIS,
  ARG,
  LYS,
  PRO
} AAtype;

const static char *AAname3[] = {"CYS", "MET", "PHE", "ILE", "LEU", "VAL", "TRP",
                                "TYR", "ALA", "GLY", "THR", "SER", "GLN", "ASN",
                                "GLU", "ASP", "HIS", "ARG", "LYS", "PRO"};

const static char *AAname1[] = {"C", "M", "F", "I", "L", "V", "W",
                                "Y", "A", "G", "T", "S", "Q", "N",
                                "E", "D", "H", "R", "K", "P"};

typedef enum { _N_CA_C_ = 0, _N_CA_CB_, _CB_CA_C_ } CEN_ANGLE_t;
const static double CEN_ANGLES[][6] = {
    1.933, 0.0368, 1.929, 0.0235, 1.926, 0.0286, // CYS=0,
    1.933, 0.0368, 1.929, 0.0235, 1.926, 0.0286, // MET,
    1.933, 0.0368, 1.929, 0.0235, 1.926, 0.0286, // PHE,
    1.933, 0.0368, 1.929, 0.0235, 1.926, 0.0286, // ILE,
    1.933, 0.0368, 1.929, 0.0235, 1.926, 0.0286, // LEU,
    1.933, 0.0368, 1.929, 0.0235, 1.926, 0.0286, // VAL,
    1.933, 0.0368, 1.929, 0.0235, 1.926, 0.0286, // TRP,
    1.933, 0.0368, 1.929, 0.0235, 1.926, 0.0286, // TYR,
    1.933, 0.0368, 1.929, 0.0235, 1.926, 0.0286, // ALA,
    1.972, 0.0439, 0.000, 0.0000, 0.000, 0.0000, // GLY,
    1.933, 0.0368, 1.929, 0.0235, 1.926, 0.0286, // THR,
    1.933, 0.0368, 1.929, 0.0235, 1.926, 0.0286, // SER,
    1.933, 0.0368, 1.929, 0.0235, 1.926, 0.0286, // GLN,
    1.933, 0.0368, 1.929, 0.0235, 1.926, 0.0286, // ASN,
    1.933, 0.0368, 1.929, 0.0235, 1.926, 0.0286, // GLU,
    1.933, 0.0368, 1.929, 0.0235, 1.926, 0.0286, // ASP,
    1.933, 0.0368, 1.929, 0.0235, 1.926, 0.0286, // HIS,
    1.933, 0.0368, 1.929, 0.0235, 1.926, 0.0286, // ARG,
    1.933, 0.0368, 1.929, 0.0235, 1.926, 0.0286, // LYS,
    1.964, 0.0394, 1.800, 0.0132, 1.939, 0.0256, // PRO
};
const static double KCEN_CONST = 1.00;
static double KCEN = KCEN_CONST;
void setKCEN(const double k);
void resetKCEN();
double getKCEN();

const static double KBBDEV = 1.00;

int residueid(const char *aa_name);
const char *residuename(int t);

const static double N_CA = 1.458;
const static double CA_C = 1.530;
const static double C_N = 1.329;
const static double C_O = 1.230;

const static double N_CA_C = 111.59;
const static double CA_C_N = 116.20;
const static double C_N_CA = 121.70;

//
// R E S I D U E    C L A S S
//

class residue : public BaseResidue {
private:
  AAtype type;

  short nHA; // number of heavy atoms
  atom *HA;  // heavy atoms

  short nPH;
  atom **PH;

  short nHAcceptor;
  atom **HAcceptor;
  hbond_acceptor *Acceptors;

  int nCHI;
  double *CHI;
  double *CHI_bak;
  double *CHI_gme;
  unsigned char *iCHI;
  vec *chi_uv; // unit vector of rotamer angles;
  bool is_Aromatic;

  char chainName; // the actual name in the PDB file
  atom ***scAtomDiv;
  int *scAtomDivN;
  /*scAtomDiv[i] is the array of the atoms affective by dihedral angle Chi[i].
  To be used in side chain atom division with Go trick.
  scAtomDivN[i] stores number of atoms in the array.
  Abe, Braun, Noguti and Go, Computers & Chemistry 1984*/

  double cen_angles[3]; // N_CA_C, N_CA_CB, CB_CA_C central angles
  double cen_angles_back[3];
  double Ecen;
  double Ecen_back;

  double E_bbdev;

  double _CEN_ANGLES[3];

public:
  residue();

  residue(const residue &right);

  residue &operator=(const residue &right) {
    if (this == &right)
      return *this;
    clear();
    init(residuename(right.getType()));
    for (int i = 0; i < nHA; i++) {
      HA[i] = right.HA[i];
    }
    chainName = right.chainName;
    return *this;
  }

  ~residue();
  /*init function initialized/allocate the atoms.*/
  void init(const char *name);
  /*initTemplate intialized the xyz positions for all the atoms,
   the residue is initialized in the origin*/
  void initTemplate_old();
  void initTemplate();
  void clear();
  void setXYZ(int at, const vec &r);
  void setChainName(char name) { chainName = name; }
  char getChainName() { return chainName; }

  bool checkSize(const topology &top);

  atom *getN() { return HA; }
  atom *getCA() { return HA + 1; }
  atom *getC() { return HA + 2; }
  atom *getO() { return HA + 3; }
  atom *getCB() { return nHA > 4 ? HA + 4 : NULL; }
  atom *getHN() {
    if (type != PRO)
      return HA->attachedProtons;
    else
      return NULL;
  }

  int getType() const { return type; }
  int getNA() { return nHA; }
  int getNH();
  int getNChi() { return nCHI; }
  int getNBBChi() { return 3; } // number of backbone dihydrals
  atom *getAtom(int i) { return HA + i; }
  atom *getAtoms() { return HA; }

  /**
   * Get Heavy atom based on PDBAtype.
   */
  atom *getHAtom(int pdb_type) {
    for (int i = 0; i < nHA; i++) {
      if (HA[i].PDBtype == pdb_type) {
        return HA + i;
      }
      for (int j = 0; j < HA[i].n_attachedProtons; j++) {
        if (HA[i].attachedProtons[j].PDBtype == pdb_type) {
          return HA[i].attachedProtons + j;
        }
      }
    }
    return NULL;
  }
  atom *getAtom(const string &pdbname) {
    for (int i = 0; i < nHA; i++) {
      for (int k = -1; k < HA[i].n_attachedProtons; k++) {
        atom *theAtom = (k < 0) ? (HA + i) : (HA[i].attachedProtons + k);
        if (string(PDBAtype2Name(theAtom->PDBtype)) == pdbname)
          return theAtom;
      }
    }
    return NULL;
  }

  atom *getAtom(const string &pdbname, const topology &top) {
    for (int i = 0; i < nHA; i++) {
      for (int k = -1; k < HA[i].n_attachedProtons; k++) {
        atom *theAtom = (k < 0) ? (HA + i) : (HA[i].attachedProtons + k);
        if (string(PDBAtype2Name(theAtom->PDBtype)) == pdbname) {
          return theAtom;
        }
      }
    }
    return NULL;
  }

  // In this method, we remove one of the attached protons
  void deprotonate(const string &pdbname) {
    // cout << "hi, there" << AAname3[type] << " " << pdbname << endl;
    atom *theAtom = getAtom(pdbname);
    if (theAtom != NULL && theAtom->n_attachedProtons > 0) {
      theAtom->deprotonate();
      setPH_ARRAY();
    }
  }

  // In this method, we add one proton to a heavyatom
  void protonate(const string &pdbname, bool is_n_terminal = false);

  // in this method, we turn an non-hbond acceptor into an hbond-acceptor
  void addToHBondAcceptor(const string &pdbname, const topology &top);

  double *getChis() { return CHI; }
  unsigned char *getIChis() { return iCHI; }
  double getChi(int i) { return CHI[i]; }
  unsigned char getIChi(int i) { return iCHI[i]; }
  const char *getName() { return residuename(type); }

  // return the largest chi angles affecting the postion of iAtom's
  // atom. Input iAtom, output the number of Chi angles that will
  // affect the position of it.
  int getNChi_AFFECTED(int iAtom);
  // calculate the unit vectors along each rotamers
  void updateChiUV();
  vec *getChiUV(int i);
  vec *getChiP2(int i);

  void updateRotamers();
  int rotateChi(int i, double theta, double **m);
  void updateICHIs(resRotLib &rotlib);
  void updateICHIs(BBDep_resRotLib &rotlib);
  /**************************************************************************
   stachasticly rotate the rotamers:
   Each rotamer angle i is determined by the average rotamer angle aveChi[i]
   and perturbed by a gaussian distribution with a standard derivation of
   the sigChi[i].
   For the rotatable protons: THR,TYR,SER, and LYS
   THR,SER, has three staggered positions(-180, -60, 60)
   TYR has only 0,180(in the ring plane)
   LYS protons have the 3-fold symmetry, therefore, only in 0.
   But for each one, an extra 20 degree is allowed with a gaussian
  distributions.
  ***************************************************************************/
  // also return the relative energy from the deviation (E/kT)
  double stachastic_rotate(BBDep_rotRecord &rot, randomGenerator &ran,
                           double **m, const double &z = 1.0);
  void stachastic_rotate_nonGaussian(BBDep_rotRecord &rot, randomGenerator &ran,
                                     double **m, const double &z = 1.0);
  /****************************************************************************
    Within a given rotamers, make random rotation around the previous chi
  angles. The variations is gaussion. The STD is 2 degrees. The Allowed range
  for a given chi angles is around z.
    ---Added by Yin---
    Also return the off-rotamer energy after the random rotation
    ---End addition---
  *****************************************************************************/
  double stachastic_rotate_Gaussian_var(BBDep_rotRecord &rot,
                                        randomGenerator &ran, double **m,
                                        const double &z = 1.0);
  /// Return the off rotamer energy.
  double getEOffRot(BBDep_rotRecord &rot);
  void fixed_rotate(BBDep_rotRecord &rot, double **m);
  void fixed_rotateH(double angle, double **m);

  void getIAxis(vec &xhat, vec &yhat, vec &zhat);
  vec nextN(vec &xbar, vec &ybar, vec &zbar);
  void rotate(const vec &xbar, const vec &ybar, const vec &zbar);
  void shift(const vec &d);

  void saveConfiguration() {
    for (int i = 0; i < nHA; i++) {
      HA[i].savePosition();
      for (int j = 0; j < HA[i].n_attachedProtons; j++)
        HA[i].attachedProtons[j].savePosition();
    }
    for (int i = 0; i < nCHI; i++)
      CHI_bak[i] = CHI[i];
    if (type == SER || type == THR || type == TYR || type == LYS) {
      CHI_bak[nCHI] = CHI[nCHI];
    }
  }
  void recoverConfiguration() {
    for (int i = 0; i < nHA; i++) {
      if (!(HA[i].r == HA[i].p))
        HA[i].r = HA[i].p;
      for (int j = 0; j < HA[i].n_attachedProtons; j++)
        if (!(HA[i].attachedProtons[j].r == HA[i].attachedProtons[j].p))
          HA[i].attachedProtons[j].r = HA[i].attachedProtons[j].p;
    }
    for (int i = 0; i < nCHI; i++)
      CHI[i] = CHI_bak[i];
    if (type == SER || type == THR || type == TYR || type == LYS) {
      CHI[nCHI] = CHI_bak[nCHI];
    }
  }
  void saveGME_Conf() {
    for (int i = 0; i < nHA; i++) {
      HA[i].gme_r = HA[i].r;
      for (int j = 0; j < HA[i].n_attachedProtons; j++)
        HA[i].attachedProtons[j].gme_r = HA[i].attachedProtons[j].r;
    }
    for (int i = 0; i < nCHI; i++)
      CHI_gme[i] = CHI[i];
    if (type == SER || type == THR || type == TYR || type == LYS) {
      CHI_gme[nCHI] = CHI[nCHI];
    }
  }
  void recoverGME_Conf() {
    for (int i = 0; i < nHA; i++) {
      if (!(HA[i].r == HA[i].gme_r))
        HA[i].r = HA[i].gme_r;
      for (int j = 0; j < HA[i].n_attachedProtons; j++)
        if (!(HA[i].attachedProtons[j].r == HA[i].attachedProtons[j].gme_r))
          HA[i].attachedProtons[j].r = HA[i].attachedProtons[j].gme_r;
    }
    for (int i = 0; i < nCHI; i++)
      CHI[i] = CHI_gme[i];
    if (type == SER || type == THR || type == TYR || type == LYS) {
      CHI[nCHI] = CHI_gme[nCHI];
    }
  }
  void saveOrigin_Conf() {
    for (int i = 0; i < nHA; i++) {
      HA[i].origin_r = HA[i].r;
      for (int j = 0; j < HA[i].n_attachedProtons; j++)
        HA[i].attachedProtons[j].origin_r = HA[i].attachedProtons[j].r;
    }
  }
  void recoverOrigin_Conf() {
    for (int i = 0; i < nHA; i++) {
      if (!(HA[i].r == HA[i].origin_r))
        HA[i].r = HA[i].origin_r;
      for (int j = 0; j < HA[i].n_attachedProtons; j++)
        if (!(HA[i].attachedProtons[j].r == HA[i].attachedProtons[j].origin_r))
          HA[i].attachedProtons[j].r = HA[i].attachedProtons[j].origin_r;
    }
  }

  // call with extreme care!
  // the two residues should have the same atoms and
  // belong to the same type
  void copyConfiguration(residue *theRes) {
    for (int i = 0; i < nHA; i++) {
      HA[i].r = theRes->HA[i].r;
      for (int j = 0; j < HA[i].n_attachedProtons; j++)
        HA[i].attachedProtons[j].r = theRes->HA[i].attachedProtons[j].r;
    }
  }

  void addPolarSH();
  void addBBH(residue &prev);
  void addBBH(bool is_n_terminal = false);

  void updatePolarSH();
  void updateBBH(residue &prev);
  void updateBBH(bool is_n_terminal = false);
  void updateBBH(const vec &r_h);

  // set the array of polar hydrogens
  void setPH_ARRAY();
  int getNPH() { return nPH; }
  atom **getPH_ARRAY() { return PH; }

  // set and get the array of hydrogen acceptors
  void initHBAcceptor_ARRAY(int n) {
    if (HAcceptor)
      delete[] HAcceptor;
    if (Acceptors) {
      for (int i = 0; i < nHAcceptor; i++) {
        Acceptors[i].getA()->acceptor = NULL;
      }
      delete[] Acceptors;
    }
    nHAcceptor = n;
    HAcceptor = new atom *[n];
    Acceptors = new hbond_acceptor[n];
  }

  void setHBAcceptor(int i, atom *iA, atom *iX, ACCEPTOR_HYBRID_t it, int mhb) {
    if (i < nHAcceptor) {
      HAcceptor[i] = iA;
      Acceptors[i].init(iA, iX, it, mhb);
      if (iA->acceptor) {
        cout << "acceptor pre-assigned" << endl;
        exit(1);
      }
      iA->acceptor = Acceptors + i;
    } else {
      cout << "error in assigning acceptors" << endl;
      exit(1);
    }
  }

  void initHBAcceptors(const topology &top);

  int getNHBAcceptor() { return nHAcceptor; }

  atom **getHBAcceptor_ARRAY() { return HAcceptor; }

  void print(ostream &out) {
    for (int i = 0; i < nHA; i++) {
      getAtom(i)->r.print(out);
      out << endl;
    }
  }

  int isNAN() {
    for (int i = 0; i < nHA; i++) {
      if (std::isnan(getAtom(i)->r.getX()) ||
          std::isnan(getAtom(i)->r.getY()) ||
          std::isnan(getAtom(i)->r.getZ())) {
        // cout << "atom " << i << endl;
        return 1;
      }
      for (int j = 0; j < HA[i].n_attachedProtons; j++) {
        if (std::isnan(HA[i].attachedProtons[j].r.getX()) ||
            std::isnan(HA[i].attachedProtons[j].r.getY()) ||
            std::isnan(HA[i].attachedProtons[j].r.getZ())) {
          // cout << "atom: " << i << " proton: " << j << endl;
          return 1;
        }
      }
    }
    return 0;
  }

  bool isAROMATIC() { return is_Aromatic; }

  int *getscAtomDivN() { return scAtomDivN; }
  atom ***getscAtomDiv() { return scAtomDiv; }

  /// Rebuild the whole side chain coordinates from the template coordinates.
  /// Align the template backbone with the current backbone.
  void rebuildSCFromTemplate();

  /// For aromatic+ ARG, to replace the planar atoms with ideal coordinate
  void planarize();
  void writePDB(ostream &out, int startRes, char chainID) {
    char line[81];
    char chain = chainID;
    int iatom = 1;
    // print the heavy atoms first
    for (int j = 0; j < this->getNA(); j++) {
      sprintf(&line[0],
              "ATOM  %5d  %-4.4s%3.3s %1c%4d    %8.3lf%8.3lf%8.3lf  1.00  0.00 "
              "          %c  ",
              iatom++, this->getAtom(j)->getName(), this->getName(), chain,
              startRes, this->getAtom(j)->getR()->getX(),
              this->getAtom(j)->getR()->getY(),
              this->getAtom(j)->getR()->getZ(), this->getAtom(j)->getName()[0]);
      out << line << endl;
    }
    // print the polar protons
    for (int j = 0; j < this->getNA(); j++) {
      for (int k = 0; k < this->getAtom(j)->n_attachedProtons; k++) {
        char atomname[10];
        int len;
        if (len = strlen(this->getAtom(j)->attachedProtons[k].getName()) == 4) {
          atomname[0] = this->getAtom(j)->attachedProtons[k].getName()[3];
          for (int c = 0; c < 3; c++)
            atomname[c + 1] = this->getAtom(j)->attachedProtons[k].getName()[c];
          atomname[4] = '\0';
          sprintf(&line[0],
                  "ATOM  %5d %-5.5s%3.3s %1c%4d    %8.3lf%8.3lf%8.3lf  1.00  "
                  "0.00           %c  ",
                  iatom++, atomname, this->getName(), chain, startRes,
                  this->getAtom(j)->attachedProtons[k].getR()->getX(),
                  this->getAtom(j)->attachedProtons[k].getR()->getY(),
                  this->getAtom(j)->attachedProtons[k].getR()->getZ(),
                  this->getAtom(j)->attachedProtons[k].getName()[0]);
        } else {
          sprintf(&line[0],
                  "ATOM  %5d  %-4.4s%3.3s %1c%4d    %8.3lf%8.3lf%8.3lf  1.00  "
                  "0.00           %c  ",
                  iatom++, this->getAtom(j)->attachedProtons[k].getName(),
                  this->getName(), chain, startRes,
                  this->getAtom(j)->attachedProtons[k].getR()->getX(),
                  this->getAtom(j)->attachedProtons[k].getR()->getY(),
                  this->getAtom(j)->attachedProtons[k].getR()->getZ(),
                  this->getAtom(j)->attachedProtons[k].getName()[0]);
        }
        out << line << endl;
      }
    }
  }

  // for compatiblity with the abstract layer
  void getBoundary(vec &min, vec &max) {
    min = vec(INF, INF, INF);
    max = vec(-INF, -INF, -INF);
    for (int j = 0; j < getNA(); j++) {
      if (min.x > getAtom(j)->r.x)
        min.x = getAtom(j)->r.x;
      if (min.y > getAtom(j)->r.y)
        min.y = getAtom(j)->r.y;
      if (min.z > getAtom(j)->r.z)
        min.z = getAtom(j)->r.z;
      if (max.x < getAtom(j)->r.x)
        max.x = getAtom(j)->r.x;
      if (max.y < getAtom(j)->r.y)
        max.y = getAtom(j)->r.y;
      if (max.z < getAtom(j)->r.z)
        max.z = getAtom(j)->r.z;
      for (int k = 0; k < getAtom(j)->n_attachedProtons; k++) {
        vec &tmp = getAtom(j)->attachedProtons[k].r;
        if (min.x > tmp.x)
          min.x = tmp.x;
        if (min.y > tmp.y)
          min.y = tmp.y;
        if (min.z > tmp.z)
          min.z = tmp.z;
        if (max.x < tmp.x)
          max.x = tmp.x;
        if (max.y < tmp.y)
          max.y = tmp.y;
        if (max.z < tmp.z)
          max.z = tmp.z;
      }
    }
  };
  BaseResidueType getBaseResidueType() { return NATURAL_AA; };
  void updateTopology(const topology &top);

  void updateCEN_ANGLES();
  void updateCEN_ANGLES0(); // use the native angles as references
  void initCEN_ANGLES();
  double getCEN_ANGLE(int i) { return _CEN_ANGLES[i]; }

  void saveCEN_ANGLES() {
    cen_angles_back[0] = cen_angles[0];
    cen_angles_back[1] = cen_angles[1];
    cen_angles_back[2] = cen_angles[2];
    Ecen_back = Ecen;
  }

  void recoverCEN_ANGLES() {
    cen_angles[0] = cen_angles_back[0];
    cen_angles[1] = cen_angles_back[1];
    cen_angles[2] = cen_angles_back[2];
    Ecen = Ecen_back;
  }

  double getEcen() { return Ecen; }

  void initE_BBDEV() {
    for (int i = 0; i < getNA(); i++) {
      if (!getAtom(i)->isSC()) {
        getAtom(i)->origin_r = getAtom(i)->r;
      }
    }
    E_bbdev = 0;
  }

  void updateE_BBDEV() {
    E_bbdev = 0;
    for (int i = 0; i < getNA(); i++) {
      if (!getAtom(i)->isSC()) {
        E_bbdev += (getAtom(i)->origin_r - getAtom(i)->r).mod2();
      }
    }
    E_bbdev *= KBBDEV / 2.0;
  }

  double getE_BBDEV() { return E_bbdev; }
};

void getEulerRotateMatrix(double alpha, double beta, double gamma,
                          double **Rot) {
  double sa = sin(alpha);
  double ca = cos(alpha);
  double sb = sin(beta);
  double cb = cos(beta);
  double sg = sin(gamma);
  double cg = cos(gamma);

  Rot[0][0] = ca * cg - sa * cb * sg;
  Rot[0][1] = -ca * sg - sa * cb * cg;
  Rot[0][2] = sb * sa;

  Rot[1][0] = sa * cg + ca * cb * sg;
  Rot[1][1] = -sa * sg + ca * cb * cg;
  Rot[1][2] = -sb * ca;

  Rot[2][0] = sb * sg;
  Rot[2][1] = sb * cg;
  Rot[2][2] = cb;
}

/**
 * A chain is defined as a collection of one or more BaseResidue connected.
 * It is an independent dockable object in a complex, it can be a protein (means
 * peptide) chain, a RNA chain , or a small molecule chain
 */
class BaseChain {
protected:
  /// center of heavy atoms
  vec center;
  vec center_bak; // for backup
private:
  /// chain Name given in the PDB
  /// record the chain name given in the pdb
  /// if the chain name is not specified it will be "_"
  char chainName;
  /// segmentation in the chain (same chain name, seperated by TER record)
  /// if one chain (e.g. "A") is artificially seperated by several TER record,
  /// the names will be "A","A1","A2",... for the artificially seperated chains
  short int chainSeg;

public:
  /// Save current atom position.
  void savePos() {
    for (int i = 0; i < size(); i++)
      getResidue(i)->savePos();
    center_bak = center;
  }
  /// Recover previously save atom position.
  void recoverPos() {
    for (int i = 0; i < size(); i++)
      getResidue(i)->recoverPos();
    center = center_bak;
  }
  /// get the chain name
  char getChainName() { return chainName; };
  /// set the chain name
  void setChainName(char _chainName) { chainName = _chainName; };
  /// get the chain segmentation id
  short int getChainSeg() { return chainSeg; };
  /// set the chain segmentation id
  void setChainSeg(short int _id) { chainSeg = _id; };
  virtual void updateTopology(const topology &top) = 0;
  virtual BaseResidue *getResidue(int ires) = 0;
  virtual int size() = 0; // the number of residues
  virtual void getBoundary(vec &min, vec &max) = 0;
  virtual double **getM() = 0; // every chain has an M matrix to facilitate
                               // external call for rotations
  virtual void setAtomIndex(int &index) = 0;
  const vec &getCenter() { return center; };
  // roate around an axis through the center

  void updateCenter() {
    center.x = center.y = center.z = 0;
    int NHA = 0; // number of heavy atoms
    for (int ires = 0; ires < this->size(); ires++) {
      BaseResidue *theRes = this->getResidue(ires);
      for (int ia = 0; ia < theRes->getNA(); ia++) {
        center += theRes->getAtom(ia)->r;
        NHA++;
      }
    }
    center /= NHA;
  }

  void shift(const vec &dist) {
    for (int ires = 0; ires < this->size(); ires++) {
      BaseResidue *theRes = this->getResidue(ires);
      for (int ia = 0; ia < theRes->getNA(); ia++)
        theRes->getAtom(ia)->r += dist;
      for (int ia = 0; ia < theRes->getNPH(); ia++)
        theRes->getPH_ARRAY()[ia]->r += dist;
    }
    center += dist;
  }

  void rotateX(double theta) {
    vec axis(1, 0, 0);
    this->rotate(axis, theta);
  }

  void rotateY(double theta) {
    vec axis(0, 1, 0);
    this->rotate(axis, theta);
  }
  void rotateZ(double theta) {
    vec axis(0, 0, 1);
    this->rotate(axis, theta);
  }

  void rotate(vec &axis, double theta) {
    vec center2 = center + axis;
    getDimondRotateMatrix(center, center2, theta, this->getM());
    for (int ires = 0; ires < this->size(); ires++) {
      BaseResidue *theRes = this->getResidue(ires);
      for (int ia = 0; ia < theRes->getNA(); ia++) {
        vec rRel = (theRes->getAtom(ia)->r - center);
        rRel.EularRotate(this->getM());
        theRes->getAtom(ia)->r = rRel + center;
      }
      for (int ia = 0; ia < theRes->getNPH(); ia++) {
        vec rRel = theRes->getPH_ARRAY()[ia]->r - center;
        rRel.EularRotate(this->getM());
        theRes->getPH_ARRAY()[ia]->r = rRel + center;
      }
    }
  }

  void rotate(double alpha, double beta, double gamma) {
    getEulerRotateMatrix(alpha, beta, gamma, this->getM());
    for (int ires = 0; ires < this->size(); ires++) {
      BaseResidue *theRes = this->getResidue(ires);
      for (int ia = 0; ia < theRes->getNA(); ia++) {
        vec rRel = (theRes->getAtom(ia)->r - center);
        rRel.EularRotate(this->getM());
        theRes->getAtom(ia)->r = rRel + center;
      }
      for (int ia = 0; ia < theRes->getNPH(); ia++) {
        vec rRel = theRes->getPH_ARRAY()[ia]->r - center;
        rRel.EularRotate(this->getM());
        theRes->getPH_ARRAY()[ia]->r = rRel + center;
      }
    }
  }
};

inline void getDimondRotateMatrix(vec &O, vec &v, double theta, double **Rot) {
  vec b = v - O;
  b.Normalize();

  double angle = theta - static_cast<int>(theta / (2 * jnc::pi)) * 2 * jnc::pi;
  if (angle < 0)
    angle += 2 * jnc::pi;

  if (angle != jnc::pi) {
    double t = tan(0.5 * angle);
    double uu = t * b.x;
    double vv = t * b.y;
    double ww = t * b.z;
    double f = 2 / (1 + t * t);

    // computing the elements of the rotational matrix
    Rot[0][0] = 0.5 * f * (1 + uu * uu - vv * vv - ww * ww);
    Rot[0][1] = f * (uu * vv - ww);
    Rot[0][2] = f * (uu * ww + vv);
    Rot[1][0] = f * (uu * vv + ww);
    Rot[1][1] = 0.5 * f * (1 - uu * uu + vv * vv - ww * ww);
    Rot[1][2] = f * (vv * ww - uu);
    Rot[2][0] = f * (uu * ww - vv);
    Rot[2][1] = f * (vv * ww + uu);
    Rot[2][2] = 0.5 * f * (1 - uu * uu - vv * vv + ww * ww);
  } else {
    double uu = b.x;
    double vv = b.y;
    double ww = b.z;
    Rot[0][0] = uu * uu - vv * vv - ww * ww;
    Rot[0][1] = 2 * uu * vv;
    Rot[0][2] = 2 * uu * ww;
    Rot[1][0] = 2 * uu * vv;
    Rot[1][1] = -uu * uu + vv * vv - ww * ww;
    Rot[1][2] = 2 * vv * ww;
    Rot[2][0] = 2 * uu * ww;
    Rot[2][1] = 2 * vv * ww;
    Rot[2][2] = -uu * uu - vv * vv + ww * ww;
  }
}

const static int MAX_BONDS = 12;

class bond_table {
private:
  int na;
  int **table;
  int nbonds;

public:
  bond_table() {
    na = 0;
    nbonds = 0;
    table = NULL;
  }

  ~bond_table() {
    if (table && table[0]) {
      delete[] table[0];
      delete[] table;
    }
  }

  void clear() {
    if (table && table[0]) {
      delete[] table[0];
      delete[] table;
    }
    na = 0;
    nbonds = 0;
    table = NULL;
  }

  void init(int n) {
    na = n;
    table = new int *[na];
    table[0] = new int[na * (MAX_BONDS + 1)];
    for (int i = 1; i < na; i++) {
      table[i] = table[i - 1] + MAX_BONDS + 1;
    }
    for (int i = 0; i < na; i++) {
      table[i][0] = 0;
    }
  }

  void setBond(int i, int j) {
    int prev_n = table[i][0];
    table[i][prev_n + 1] = j;
    table[i][0]++;
    prev_n = table[j][0];
    table[j][prev_n + 1] = i;
    table[j][0]++;
    nbonds++;
  }

  int isBonded(int i, int j) {
    for (int a = 1; a < table[i][0] + 1; a++) {
      if (table[i][a] == j)
        return 1;
    }
    return 0;
  }

  int getNBonds() { return nbonds; }

  int getNBonds(int i) { return table[i][0]; }

  int *getBondList(int i) { return table[i] + 1; }

  int getNA() { return na; }
};
//------------------------------------------------------
// The table of angles are derived from the bond table
// i<->j<->k; assuming i<k
//------------------------------------------------------
class angle_table {
private:
  int nangles;
  int **table;
  int *search;
  class triple {
  public:
    int i;
    int j;
    int k;
    triple(int a, int b, int c) : i(a), j(b), k(c){};
    ~triple(){};
  };

public:
  angle_table() {
    table = NULL;
    search = NULL;
  }

  ~angle_table() {
    if (table && search && table[0]) {
      delete[] table[0];
      delete[] table;
      delete[] search;
    }
  }

  void clear() {
    if (table && search && table[0]) {
      delete[] table[0];
      delete[] table;
      delete[] search;
    }
    table = NULL;
    search = NULL;
  }

  void init(bond_table &bonds) {
    vector<triple> tmp_pool;
    search = new int[bonds.getNA()];
    int index = 0;
    for (int ii = 0; ii < bonds.getNA(); ii++) {
      search[ii] = index;
      // cout << ii << endl;
      for (int jj = ii + 1; jj < bonds.getNA(); jj++) {
        if (!bonds.isBonded(ii, jj)) { // ii and jj should be bonded
          int ii_nb = bonds.getNBonds(ii);
          int jj_nb = bonds.getNBonds(jj);
          for (int iip = 0; iip < ii_nb; iip++) {
            int kk = bonds.getBondList(ii)[iip];
            for (int jjp = 0; jjp < jj_nb; jjp++) {
              if (kk == bonds.getBondList(jj)[jjp]) {
                tmp_pool.push_back(triple(ii, kk, jj));
                index++;
              }
            }
          }
        }
      } // end of looping jj
    }   // end of looping ii
    nangles = tmp_pool.size();
    table = new int *[nangles];
    table[0] = new int[nangles * 3];
    for (int i = 1; i < nangles; i++)
      table[i] = table[i - 1] + 3;
    for (int i = 0; i < tmp_pool.size(); i++) {
      table[i][0] = tmp_pool[i].i;
      table[i][1] = tmp_pool[i].j;
      table[i][2] = tmp_pool[i].k;
    }
    tmp_pool.clear();
  }

  int getNAngles() { return nangles; }

  int **getTable() { return table; }

  int isAngled(int i, int k) {
    int pt;
    if (i < k) {
      for (pt = search[i]; pt < search[i + 1]; pt++) {
        if (table[pt][2] == k)
          return 1;
      }
    } else if (i > k) {
      for (pt = search[k]; pt < search[k + 1]; pt++) {
        if (table[pt][2] == i)
          return 1;
      }
    }
    return 0;
  }
};
//-------------------------------------------------------
// 1-4 lists, also derived from bond tables and also its derivative angle table
// which have different VDW and/or EEL
// i<->?<->?<->j; assuming i<j
//-------------------------------------------------------
class link14_table {
private:
  int nlinks;
  int **table;
  int *search;
  class dimer {
  public:
    int i;
    int j;
    dimer(int a, int b) : i(a), j(b){};
    ~dimer() {}
  };

public:
  link14_table() {
    table = NULL;
    search = NULL;
  }

  ~link14_table() {
    if (table && search && table[0]) {
      delete[] table[0];
      delete[] table;
      delete[] search;
    }
  }

  void clear() {
    if (table && search && table[0]) {
      delete[] table[0];
      delete[] table;
      delete[] search;
    }
    table = NULL;
    search = NULL;
  }

  void init(bond_table &bonds, angle_table &angles) {
    vector<dimer> tmp_pool;
    search = new int[bonds.getNA()];
    int index = 0;
    for (int ii = 0; ii < bonds.getNA(); ii++) {
      search[ii] = index;
      for (int ll = ii + 1; ll < bonds.getNA(); ll++) {
        if (!bonds.isBonded(ii, ll) &&  // not bonded
            !angles.isAngled(ii, ll)) { // not angled
          for (int iip = 0; iip < bonds.getNBonds(ii); iip++) {
            int jj = bonds.getBondList(ii)[iip];
            for (int llp = 0; llp < bonds.getNBonds(ll); llp++) {
              int kk = bonds.getBondList(ll)[llp];
              if (bonds.isBonded(jj, kk)) { // ii's bp and jj'bp are also bonded
                tmp_pool.push_back(dimer(ii, ll));
                index++;
              }
            } // end of loop jj's bond partner
          }   // end of loop ii's bond partner
        }
      } // end of loop jj
    }   // end of loop ii
    // cout << tmp_pool.size() << endl;
    nlinks = tmp_pool.size();
    table = new int *[nlinks];
    table[0] = new int[nlinks * 2];
    for (int i = 1; i < nlinks; i++)
      table[i] = table[i - 1] + 2;
    for (int i = 0; i < tmp_pool.size(); i++) {
      table[i][0] = tmp_pool[i].i;
      table[i][1] = tmp_pool[i].j;
    }
    tmp_pool.clear();
  }

  int getNLinks() { return nlinks; }

  int **getTable() { return table; }

  int isLink14(int i, int j) {
    if (i < j) {
      for (int a = search[i]; a < search[i + 1]; a++)
        if (table[a][1] == j)
          return 1;
    } else if (j < i) {
      for (int a = search[j]; a < search[j + 1]; a++)
        if (table[a][1] == i)
          return 1;
    }
    return 0;
  }
};

class tetra {
public:
  int i;
  int j;
  int k;
  int l;
  tetra(int a, int b, int c, int d) : i(a), j(b), k(c), l(d) {}
  ~tetra(){};
};
//-----------------------------------------------------
// Dihedral terms should be arranged as i<->j<->k<->l
// assuming i<l
//-----------------------------------------------------
class dihedral_table {
private:
  vector<tetra> vtable;
  int ndihedral;
  int **table;

public:
  dihedral_table() {
    ndihedral = 0;
    table = NULL;
  }

  ~dihedral_table() {
    vtable.clear();
    if (table && table[0]) {
      delete[] table[0];
      delete[] table;
    }
  }

  void clear() {
    vtable.clear();
    if (table && table[0]) {
      delete[] table[0];
      delete[] table;
    }
    ndihedral = 0;
    table = NULL;
  }

  void addTetra(int i, int j, int k, int l) {
    if (i < l)
      vtable.push_back(tetra(i, j, k, l));
    else
      vtable.push_back(tetra(l, k, j, i));
    ndihedral++;
  }

  int checkConnectivity(bond_table &bonds) {
    for (int n = 0; n < vtable.size(); n++) {
      tetra &theTetra = vtable[n];
      if (bonds.isBonded(theTetra.i, theTetra.j) &&
          bonds.isBonded(theTetra.j, theTetra.k) &&
          bonds.isBonded(theTetra.k, theTetra.l))
        ;
      else
        return 0;
    }
    return 1;
  }

  void constructTable() {
    // ndihedral = vtable.size();
    table = new int *[ndihedral];
    table[0] = new int[ndihedral * 4];
    for (int i = 1; i < ndihedral; i++)
      table[i] = table[i - 1] + 4;
    for (int i = 0; i < ndihedral; i++) {
      table[i][0] = vtable[i].i;
      table[i][1] = vtable[i].j;
      table[i][2] = vtable[i].k;
      table[i][3] = vtable[i].l;
    }
    vtable.clear();
  }

  int getNDihedral() { return ndihedral; }

  int **getTable() { return table; }
};
//---------------------------------------------------
// Improper dihedrals i, j, k, l
// i<->j, k<->j, l<->j; j is the centera atom
//--------------------------------------------------
class improper_table {
private:
  vector<tetra> vtable;
  int nimproper;
  int **table;

public:
  improper_table() {
    nimproper = 0;
    table = NULL;
  }

  ~improper_table() {
    vtable.clear();
    if (table && table[0]) {
      delete[] table[0];
      delete[] table;
    }
  }

  void clear() {
    vtable.clear();
    if (table && table[0]) {
      delete[] table[0];
      delete[] table;
    }
    nimproper = 0;
    table = NULL;
  }

  void addTetra(int i, int j, int k, int l) {
    nimproper++;
    vtable.push_back(tetra(i, j, k, l));
  }

  int checkConnectivity(bond_table &bonds) {
    for (int n = 0; n < vtable.size(); n++) {
      tetra &theTetra = vtable[n];
      if (bonds.isBonded(theTetra.i, theTetra.j) &&
          bonds.isBonded(theTetra.k, theTetra.j) &&
          bonds.isBonded(theTetra.l, theTetra.j))
        ;
      return 0;
    }
    return 1;
  }

  void constructTable() {
    // nimproper = vtable.size();
    table = new int *[nimproper];
    table[0] = new int[nimproper * 4];
    for (int i = 1; i < nimproper; i++)
      table[i] = table[i - 1] + 4;
    for (int i = 0; i < nimproper; i++) {
      table[i][0] = vtable[i].i;
      table[i][1] = vtable[i].j;
      table[i][2] = vtable[i].k;
      table[i][3] = vtable[i].l;
    }
    vtable.clear();
  }

  int getNImproper() { return nimproper; }

  int **getTable() { return table; }
};

class connect {
private:
  int na;
  bond_table bonds;
  angle_table angles;
  link14_table links1_4;
  dihedral_table dihedrals;
  improper_table impropers;

public:
  connect() { na = 0; }
  ~connect() {}

  bond_table &getBonds() { return bonds; }

  angle_table &getAngles() { return angles; }

  dihedral_table &getDihedrals() { return dihedrals; }

  improper_table &getImpropers() { return impropers; }

  // init process before all the addin gof BOND, DIHE, IMPR
  void init(int na) { bonds.init(na); }
  void clear() {
    bonds.clear();
    angles.clear();
    links1_4.clear();
    dihedrals.clear();
    impropers.clear();
  }

  // bonds
  void initBonds(int n) { bonds.init(n); }
  void setBond(int i, int j) { bonds.setBond(i, j); }
  int isBonded(int i, int j) { return bonds.isBonded(i, j); }
  int getNBonds() { return bonds.getNBonds(); }
  int getNBondsOf(int i) { return bonds.getNBonds(i); }
  int *getBondListOf(int i) { return bonds.getBondList(i); }
  // angles
  void initAngles() { angles.init(bonds); }
  int getNAngles() { return angles.getNAngles(); }
  int **getAngleTable() { return angles.getTable(); }
  int isAngled(int i, int j) { return angles.isAngled(i, j); }
  // 1-4links
  void init1_4Links() { links1_4.init(bonds, angles); }
  int getN1_4Links() { return links1_4.getNLinks(); }
  int **get1_4LinkTable() { return links1_4.getTable(); }
  int is1_4Linked(int i, int j) { return links1_4.isLink14(i, j); }
  // dihedrals
  void addDihedral(int i, int j, int k, int l) {
    dihedrals.addTetra(i, j, k, l);
  }
  int finalizeDihedrals() {
    int tmp = dihedrals.checkConnectivity(bonds);
    if (tmp)
      dihedrals.constructTable();
    return tmp;
  }
  int getNDihedrals() { return dihedrals.getNDihedral(); }
  int **getDihedralTable() { return dihedrals.getTable(); }
  // impropers
  void addImproper(int i, int j, int k, int l) {
    impropers.addTetra(i, j, k, l);
  }
  int finalizeImpropers() {
    int tmp = impropers.checkConnectivity(bonds);
    if (tmp)
      impropers.constructTable();
    return tmp;
  }
  int getNImpropers() { return impropers.getNImproper(); }
  int **getImproperTable() { return impropers.getTable(); }

  // post process after all the input of BOND, DIHE, IMPR
  int postInitialize() {
    angles.init(bonds);
    links1_4.init(bonds, angles);
    int isDIHE = finalizeDihedrals();
    int isIMPR = finalizeImpropers();
    if (isDIHE && isIMPR)
      return 1;
    else
      return 0;
  }
};

typedef enum { PROTEIN = 0, RNA, DNA, MOLECULE, OTHER } chain_t;

class atm {
public:
  std::string _name;
  vec _r;
  int _index;
  double _bfactor;

  atm() {}

  atm(const std::string &n, double x, double y, double z, double b = 0)
      : _r(x, y, z) {
    _name = n;
    _bfactor = b;
    //        std::cout << _name << std::endl;
  }

  atm(const atm &atm2) {
    //        std::cout << "atm copy construct" << std::endl;
    _name = atm2._name;
    _r = atm2._r;
    _index = atm2._index;
    _bfactor = atm2._bfactor;
  }

  atm &operator=(const atm &atm2) {
    //        std::cout << "atm copy" << std::endl;
    _name = atm2._name;
    _r = atm2._r;
    _index = atm2._index;
    _bfactor = atm2._bfactor;
    return *this;
  }

  const string name() const { return _name; }
  void name(const string &n) { _name = n; }
  const vec pos() const { return _r; }
  void pos(const vec &r) { _r = r; }

  void index(const int i) { _index = i; }
  const int index() const { return _index; }

  void bfactor(const double b) { _bfactor = b; }
  const double bfactor() const { return _bfactor; }
};

class res {
public:
  std::string _name;
  std::deque<atm> _res = {};
  int _index;
  short _isHETATM;
  int _resID;
  //    friend class pdb;

  res() { _isHETATM = 0; }

  res(const string &n, int index) : _name(n), _index(index) { _isHETATM = 0; }

  const deque<atm> &atoms() const { return _res; }
  const atm &atom(int i) const {
    if (i >= 0 && i < _res.size())
      return _res[i];
    else {
      cerr << "wrong index" << endl;
      exit(1);
    }
  }
  const atm &atom(string name) const {
    for (int i = 0; i < _res.size(); i++) {
      if (_res[i].name() == name) {
        return _res[i];
      }
    }
    cerr << "can not find the specified atom:" << name << endl;
    exit(1);
  }
  const atm &gamma() const {
    for (int i = 0; i < _res.size(); i++) {
      if (_res[i].name().substr(0, 2) == "CG" ||
          _res[i].name().substr(0, 2) == "OG" ||
          _res[i].name().substr(0, 2) == "SG") {
        return _res[i];
      }
    }
    cerr << "can not find the specified gamma atom" << endl;
    exit(1);
  }
  void addAtom(const atm &theAtom) { _res.push_back(theAtom); }
  void name(const string &n) { _name = n; }
  const string &name() const { return _name; }
  int nAtoms() const { return _res.size(); }
  int index() const { return _index; }
  void index(int i) { _index = i; }
  short isHETATM() const { return _isHETATM; }
  void isHETATM(int i) { _isHETATM = i; }
  int resID() const { return _resID; }
  void resID(int i) { _resID = i; }
};

class mol {
public:
  string _name;
  std::deque<res> _mol = {};
  char _ID;
  short _ct;

  mol() {}
  mol(char ID) : _ID(ID) { _ct = -1; }
  ~mol() {}

  const deque<res> &residues() const { return _mol; }
  const res &residue(int i) const { return _mol[i]; }
  void addResidue(const res &theRes) { _mol.push_back(theRes); }
  const string &name() const { return _name; }
  void name(const string &n) { _name = n; }
  const char ID() const { return _ID; }
  void ID(char a) { _ID = a; }
  int nResidues() const { return _mol.size(); }
  void write(ostream &out) const;
  void write(ostream &out, const vec &trans) const;

  chain_t chain_type() { return static_cast<chain_t>(_ct); }
};

class pdb {
public:
  string _name;
  std::deque<mol> _molecules = {};

  pdb(string n) : _name(n) {}

  const deque<mol> &molecules() { return _molecules; }
  mol &molecule(int i) { return _molecules[i]; }
  void addMolecules(const mol &theMol) { _molecules.push_back(theMol); }
  const string &name() { return _name; }
  void name(const string &n) { _name = n; }
  int nMolecules() { return _molecules.size(); }

  void init(const char *file);
  void write(ostream &out);
  void assignBFactor(const map<string, double> &bf_array);

  void setGoIndex();
  int setDopeIndex();
  int setDope2BIndex();
  int setDopeHIndex();
  void addBBH();
  bool hasBBH();
};

typedef enum { _PDB_ = 0, _PDB_BB_, _PDB_MIX_, _SEQ_ } generate_t;

class protein : public BaseChain {
private:
  char name[100];
  //  char chainName;
  int length;
  residue *res;
  double *phi;
  double *psi;

  double **m;
  connect con;
  friend class Complex;

public:
  protein() {
    res = NULL;
    phi = NULL;
    psi = NULL;
    m = NULL;
  }

  void init(const mol &theMol, const bool hasBBH = false);
  void init(const string &name, const vector<string> &seq);

  protein(const char *thename, char *file, generate_t t);
  protein(const char *thename, char *file, vector<int> &missingSCs,
          generate_t t);
  ~protein() {
    if (res)
      delete[] res;
    if (phi)
      delete[] phi;
    if (psi)
      delete[] psi;
    if (m && m[0]) {
      delete[] m[0];
      delete[] m;
    }
  }
  int pointMutation(int i, const char *newRes);
  double **getM() { return m; }

  // update the force field type, mass and charge of each atom
  // and initialize the hydrogen bond related poloar protons and acceptors
  void updateTopology(const topology &top);
  // in the complex, also assign the res ids
  void updateTopology(const topology &top, const int iP,
                      const int startIndexOfChain);
  // update the FF_t, mass and charge for a given residue
  // and initialize the hydrogen bond related poloar protons and acceptors
  void updateTopology(const topology &top, const int ires);

  residue *getResidue(int i) {
    if (i >= 0 && i < length)
      return res + i;
    return NULL;
  }

  residue *getResidues() { return res; }

  // rotate the residue to a new rotatmer state
  int rotateResidue(int resi_index, const rotRecord &newRotamer);
  void addPH();
  void updatePH();

  void writePDB(ostream &out);
  void
  writePDB(ostream &out, int startRes,
           char chainID); // write PDB file with given start resid and chain ID
  void writeSTATE(ostream &out);
  void getBoundary(vec &min, vec &max);

  int size() { return length; }
  void setResidueIndex();
  int setAtomIndex();
  int setHAtomIndex();
  void setAtomIndex(int &index);

  void
  initConnect(const topology &top); // after point mutation;
                                    // 1) updatetopology of the residue
                                    // 2) indexing 3) initialize connectivity
  void clearConnect() { con.clear(); }
  connect &getConnnect() { return con; }

  int isBonded(atom *ai, atom *aj) {
    return con.isBonded(ai->index, aj->index);
  }

  int isAngled(atom *ai, atom *aj) {
    return con.isAngled(ai->index, aj->index);
  }

  int is1_4Linked(atom *ai, atom *aj) {
    return con.is1_4Linked(ai->index, aj->index);
  }

  void translateTo(const vec &cm) {
    vec old_cm(0, 0, 0);
    int iatom = 0;
    for (int i = 0; i < length; i++) {
      for (int j = 0; j < res[i].getNA(); j++) {
        old_cm += res[i].getAtom(j)->r;
        iatom++;
        for (int k = 0; k < res[i].getAtom(j)->n_attachedProtons; k++) {
          old_cm += res[i].getAtom(j)->attachedProtons[k].r;
          iatom++;
        }
      }
    }
    old_cm /= static_cast<double>(iatom);
    old_cm = cm - old_cm;
    for (int i = 0; i < length; i++) {
      for (int j = 0; j < res[i].getNA(); j++) {
        res[i].getAtom(j)->r += old_cm;
        for (int k = 0; k < res[i].getAtom(j)->n_attachedProtons; k++) {
          res[i].getAtom(j)->attachedProtons[k].r += old_cm;
        }
      }
    }
  }

  void rotatePhi(int ires, const double &phi_i) {
    double phi_old = phi[ires];
    double dth = phi_i - phi_old;
    getDimondRotateMatrix(getResidue(ires)->getN()->r,
                          getResidue(ires)->getCA()->r, dth, getM());
    for (int i = 0; i < getResidue(ires)->getNA(); i++) {
      if (i != 0 && i != 1 &&
          !(i == 6 && getResidue(ires)->getType() == PRO)) { // not N, CA
        atom *theAtom = getResidue(ires)->getAtom(i);
        vec tmp = theAtom->r - getResidue(ires)->getCA()->r;
        tmp.EularRotate(getM());
        theAtom->r = tmp + getResidue(ires)->getCA()->r;
        for (int k = 0; k < theAtom->n_attachedProtons; k++) {
          tmp = theAtom->attachedProtons[k].r - getResidue(ires)->getCA()->r;
          tmp.EularRotate(getM());
          theAtom->attachedProtons[k].r = tmp + getResidue(ires)->getCA()->r;
        }
      }
    }
    if (getResidue(ires)->getType() == PRO)
      getResidue(ires)->updateRotamers();
    for (int i = ires + 1; i < size(); i++) {
      for (int j = 0; j < getResidue(i)->getNA(); j++) {
        atom *theAtom = getResidue(i)->getAtom(j);
        vec tmp = theAtom->r - getResidue(ires)->getCA()->r;
        tmp.EularRotate(getM());
        theAtom->r = tmp + getResidue(ires)->getCA()->r;
        for (int k = 0; k < theAtom->n_attachedProtons; k++) {
          tmp = theAtom->attachedProtons[k].r - getResidue(ires)->getCA()->r;
          tmp.EularRotate(getM());
          theAtom->attachedProtons[k].r = tmp + getResidue(ires)->getCA()->r;
        }
      }
    }

    phi[ires] = phi_i;
  }

  void rotatePsi(int ires, const double &psi_i) {
    double psi_old = psi[ires];
    double dth = psi_i - psi_old;
    getDimondRotateMatrix(getResidue(ires)->getCA()->r,
                          getResidue(ires)->getC()->r, dth, getM());

    vec tmp = getResidue(ires)->getO()->r - getResidue(ires)->getC()->r;
    tmp.EularRotate(getM());
    getResidue(ires)->getO()->r = tmp + getResidue(ires)->getC()->r;

    for (int i = ires + 1; i < size(); i++) {
      for (int j = 0; j < getResidue(i)->getNA(); j++) {
        atom *theAtom = getResidue(i)->getAtom(j);
        vec tmp = theAtom->r - getResidue(ires)->getC()->r;
        tmp.EularRotate(getM());
        theAtom->r = tmp + getResidue(ires)->getC()->r;
        for (int k = 0; k < theAtom->n_attachedProtons; k++) {
          tmp = theAtom->attachedProtons[k].r - getResidue(ires)->getC()->r;
          tmp.EularRotate(getM());
          theAtom->attachedProtons[k].r = tmp + getResidue(ires)->getC()->r;
        }
      }
    }

    psi[ires] = psi_i;
  }

  //  char getChainName(){ return chainName;}

  /************************
    BackRub motion.
    Ref: The Backrub Mition: How Protein Backbone Shrugs when a sidechain
    Dances. --Davis, IW; Arendal II, WE; Richardson&Richardson; Structure 2006.
  ********************************/
  /* Rotate around the axis of CAs of residues i-1 and i+1.
     dtheta is in RADIAN. */
  void rotateTheta13(int i, const double &dtheta);

  /*Rotate around the axis of CAs of residues i-1 and i.
    dtheta is in RADIAN. */
  void rotateTheta12(int i, const double &dtheta);

  /*Rotate around the axis of CAs of residues i+1 and i*/
  void rotateTheta32(int i, const double &dtheta);

  /*update the CEN-ANGLES*/
  double updateCEN_ANGLES() {
    for (int i = 0; i < size(); i++) {
      getResidue(i)->updateCEN_ANGLES();
    }
  }

  bool updatePhiPsi();
  void getPhiPsi(int ires, double &f, double &y);

  void setPhiPsi(int i, const double &f, const double &y) {
    if (i >= 0 && i < length) {
      phi[i] = f;
      psi[i] = y;
    }
  }
};

class gresidue : public BaseResidue {
private:
  AAtype type;
  residue *instances;
  AAtype _GME_type;

public:
  gresidue();
  ~gresidue();

  /*return the current amino acid type*/
  AAtype getType() { return type; }
  AAtype getGMEType() { return _GME_type; }
  /*return the current residue*/
  residue *getResidue();
  void mutate(AAtype newType) { type = newType; }

  residue *get_residue(AAtype theType) { return &instances[theType]; }
  /*create the accosicate*/
  /*init function intialized all the possible amino acids in the position.
    the initial type corresponds to the input template.*/
  void init(residue &temp);
  /*recalculate the rotamer angles for each amino acid*/
  void updateRotamers_ALL();
  /*for each of candidate residues, update the iCHIs*/
  void updateICHIs_ALL(BBDep_RotLib &rotlib, int &iF, int &iY);
  /*for each of candidate residues, adding the polar sidechain protons*/
  void addPolarSH_ALL();
  /*for each of candidate residues, add the backbone protons*/
  void addBBH_ALL(residue &prev);
  void addBBH_ALL();
  void updateBBH_ALL(const vec &r);

  /*for each of the candidate residues, initialize the polar proton array*/
  void setPH_ARRAY_ALL();
  /*for each of the candidate residues, initialize the acceptor array*/
  void initHBAcceptors_ALL(const topology &top);

  /*after adding the pHs for all, the function assign the charge, mass and
    force field type and sidechain bit. In addition, it will initialize the
    proton and  acceotpr arrays*/
  void updateTopology(const topology &top); // after add the pHs
  void updateTopology_ALL(const topology &top) { return updateTopology(top); };
  /*initialize the residue index*/
  void setResidueIndex_ALL(int idx) {
    for (int i = 0; i < 20; i++) {
      for (int j = 0; j < instances[i].getNA(); j++) {
        instances[i].getAtom(j)->resID = idx;
        for (int k = 0; k < instances[i].getAtom(j)->n_attachedProtons; k++) {
          instances[i].getAtom(j)->attachedProtons[k].resID = idx;
        }
      }
    }
  }

  /* initialize the sub-residue index
   * In complex, atom::resID includes all chains, atom::resSubID are labeled
   * within each chain*/
  void setSubResidueIndex_ALL(int idx) {
    for (int i = 0; i < 20; i++) {
      for (int j = 0; j < instances[i].getNA(); j++) {
        instances[i].getAtom(j)->subResID = idx;
        for (int k = 0; k < instances[i].getAtom(j)->n_attachedProtons; k++) {
          instances[i].getAtom(j)->attachedProtons[k].subResID = idx;
        }
      }
    }
  }
  void setChainIndex_ALL(int idx) {
    for (int i = 0; i < 20; i++) {
      for (int j = 0; j < instances[i].getNA(); j++) {
        instances[i].getAtom(j)->chainID = idx;
        for (int k = 0; k < instances[i].getAtom(j)->n_attachedProtons; k++) {
          instances[i].getAtom(j)->attachedProtons[k].chainID = idx;
        }
      }
    }
  }

  void updateICHIs_ALL(BBDep_RotLib &rotlib, double &phi, double &psi);

  // make the backbone of other residues aligned to the active one
  void updateGResidue() {
    vec new_xhat, new_yhat, new_zhat;
    vec old_xhat, old_yhat, old_zhat;
    getResidue()->getIAxis(new_xhat, new_yhat, new_zhat);

    double **m = new double *[3];
    m[0] = new double[9];
    m[1] = m[0] + 3;
    m[2] = m[1] + 3;

    double **rr = new double *[3];
    rr[0] = new double[9];
    rr[1] = rr[0] + 3;
    rr[2] = rr[1] + 3;

    vec new_ca_r = getResidue()->getCA()->r;
    // new_ca_r.print(cout);
    // cout << endl;
    for (int i = 0; i < 20; i++) {
      if (i != type) { // non active residue
        // instances[i].getCA()->r.print(cout);
        instances[i].getIAxis(old_xhat, old_yhat, old_zhat);
        // calculate the revers matrox of (old_xhat, old_yhat, old_zhat) !!!
        // should be just the transpose
        m[0][0] = old_yhat.getY() * old_zhat.getZ() -
                  old_yhat.getZ() * old_zhat.getY();
        m[1][0] = -(old_xhat.getY() * old_zhat.getZ() -
                    old_xhat.getZ() * old_zhat.getY());
        m[2][0] = old_xhat.getY() * old_yhat.getZ() -
                  old_xhat.getZ() * old_yhat.getY();

        m[0][1] = -(old_yhat.getX() * old_zhat.getZ() -
                    old_yhat.getZ() * old_zhat.getX());
        m[1][1] = old_xhat.getX() * old_zhat.getZ() -
                  old_xhat.getZ() * old_zhat.getX();
        m[2][1] = -(old_xhat.getX() * old_yhat.getZ() -
                    old_xhat.getZ() * old_yhat.getX());

        m[0][2] = old_yhat.getX() * old_zhat.getY() -
                  old_yhat.getY() * old_zhat.getX();
        m[1][2] = -(old_xhat.getX() * old_zhat.getY() -
                    old_xhat.getY() * old_zhat.getX());
        m[2][2] = old_xhat.getX() * old_yhat.getY() -
                  old_xhat.getY() * old_yhat.getX();

        // to calculate the rotation matrix
        // (new_xhat, new_yhat, new_zhat)*rr
        rr[0][0] = new_xhat.getX() * m[0][0] + new_yhat.getX() * m[1][0] +
                   new_zhat.getX() * m[2][0];
        rr[0][1] = new_xhat.getX() * m[0][1] + new_yhat.getX() * m[1][1] +
                   new_zhat.getX() * m[2][1];
        rr[0][2] = new_xhat.getX() * m[0][2] + new_yhat.getX() * m[1][2] +
                   new_zhat.getX() * m[2][2];

        rr[1][0] = new_xhat.getY() * m[0][0] + new_yhat.getY() * m[1][0] +
                   new_zhat.getY() * m[2][0];
        rr[1][1] = new_xhat.getY() * m[0][1] + new_yhat.getY() * m[1][1] +
                   new_zhat.getY() * m[2][1];
        rr[1][2] = new_xhat.getY() * m[0][2] + new_yhat.getY() * m[1][2] +
                   new_zhat.getY() * m[2][2];

        rr[2][0] = new_xhat.getZ() * m[0][0] + new_yhat.getZ() * m[1][0] +
                   new_zhat.getZ() * m[2][0];
        rr[2][1] = new_xhat.getZ() * m[0][1] + new_yhat.getZ() * m[1][1] +
                   new_zhat.getZ() * m[2][1];
        rr[2][2] = new_xhat.getZ() * m[0][2] + new_yhat.getZ() * m[1][2] +
                   new_zhat.getZ() * m[2][2];

        // perform the rotations and translations.
        vec ca_r = instances[i].getCA()->r;
        for (int j = 0; j < instances[i].getNA(); j++) {
          atom *theAtom = instances[i].getAtom(j);
          vec tmp = theAtom->r - ca_r;
          tmp.EularRotate(rr);
          theAtom->r = tmp + new_ca_r;
          for (int k = 0; k < theAtom->n_attachedProtons; k++) {
            tmp = theAtom->attachedProtons[k].r - ca_r;
            tmp.EularRotate(rr);
            theAtom->attachedProtons[k].r = tmp + new_ca_r;
          }
        } // end the rotatiom and translation;
      }

      // checking
      {
        double EPS = 1.0e-6;
        if ((getResidue()->getN()->r.getDist(instances[i].getN()->r) > EPS) ||
            (getResidue()->getCA()->r.getDist(instances[i].getCA()->r) > EPS) ||
            (getResidue()->getC()->r.getDist(instances[i].getC()->r) > EPS)) {
          cout << "error " << endl;
        }
      }
      // ensumre the backbone H & O are the same
      if (instances[i].getHN() && getResidue()->getHN())
        instances[i].getHN()->r = getResidue()->getHN()->r;
      instances[i].getO()->r = getResidue()->getO()->r;
    }
    // cout << endl;
    delete[] m[0];
    delete[] m;
    delete[] rr[0];
    delete[] rr;
  }

  int isNAN() {
    for (int i = 0; i < 20; i++) {
      if (instances[i].isNAN()) {
        // cout << " error in type: " << i+1 << " " << type << endl;
        return 1;
      }
    }
    return 0;
  }

  /** to comply with BaseResidue abstract class interface **/
  int getNA() { return getResidue()->getNA(); };
  int getNH() { return getResidue()->getNH(); };
  atom *getAtom(int i) { return getResidue()->getAtom(i); };
  atom *getAtoms() { return getResidue()->getAtoms(); };

  atom *getAtom(const string &pdbname, const topology &top) {
    return getResidue()->getAtom(pdbname, top);
  }

  const char *getName() { return getResidue()->getName(); }
  int getNPH() { return getResidue()->getNPH(); };
  atom **getPH_ARRAY() { return getResidue()->getPH_ARRAY(); };
  int getNHBAcceptor() { return getResidue()->getNHBAcceptor(); };
  atom **getHBAcceptor_ARRAY() { return getResidue()->getHBAcceptor_ARRAY(); };
  void print(ostream &out) { return getResidue()->print(out); };
  BaseResidueType getBaseResidueType() { return GEN_NATURAL_AA; };
  int getNChi() { return getResidue()->getNChi(); }
  int getNBBChi() { return getResidue()->getNBBChi(); }
  void getBoundary(vec &min, vec &max) {
    return getResidue()->getBoundary(min, max);
  };
  void saveGME_Conf() {
    _GME_type = type;
    getResidue()->saveGME_Conf();
  }
  void recoverGME_Conf() {
    type = _GME_type;
    getResidue()->recoverGME_Conf();
  }
};

class aa_rotamer {
private:
  int type;
  BBDep_rotRecord *record;
  double min_bb_vdwr;
  bool available;
  int nHit;

public:
  aa_rotamer() {
    type = -1;
    record = NULL;
    available = true;
    min_bb_vdwr = -1;
    nHit = 0;
  }
  aa_rotamer(int itype, BBDep_rotRecord *irecord) {
    type = itype;
    record = irecord;
    available = true;
    min_bb_vdwr = -1;
    nHit = 0;
  }
  aa_rotamer(const aa_rotamer &old) {
    type = old.type;
    record = old.record;
    min_bb_vdwr = old.min_bb_vdwr;
    available = old.available;
    nHit = old.nHit;
  }

  ~aa_rotamer() {}

  void init(int itype, BBDep_rotRecord *irecord) {
    type = itype;
    record = irecord;
    available = true;
    nHit = 0;
  }

  void init(aa_rotamer &iaa_rot) {
    type = iaa_rot.type;
    record = iaa_rot.record;
    min_bb_vdwr = iaa_rot.min_bb_vdwr;
    available = true;
    nHit = iaa_rot.nHit;
  }

  void getDATA(int &otype, BBDep_rotRecord *&orecord) {
    // if(available){
    otype = type;
    orecord = record;
    //}
    // else{
    //  otype = -1;
    //  orecord = NULL;
    //}
  }

  void disable() { available = false; }

  bool isAvailable() { return available; }

  void setMIN_BB_VDWR(double &vdwr) {
    if (min_bb_vdwr < 0) {
      min_bb_vdwr = vdwr;
    } else if (vdwr < min_bb_vdwr) {
      min_bb_vdwr = vdwr;
    }
  }

  double getMIN_BB_VDWR() { return min_bb_vdwr; }

  int getNHit() { return nHit; }

  void increaseHit() { nHit++; }

  int resetHit() { nHit = 0; }

  int setHit(const int &in) { nHit = in; }

  BBDep_rotRecord *getRotRecord() { return record; }
  int getType() { return type; }

  void updateRotamer(BBDep_RotLib &rotLib, const int &newF, const int &newY) {
    if (type != ALA && type != GLY) {
      record = rotLib.newRotamer(record, type, newF, newY);
    }
  }
};

class con_res {
private:
  char *type2index;
  char na;
  char **bond_table;
  char **angle_table;

public:
  con_res() {
    type2index = NULL;
    na = 0;
    bond_table = NULL;
    angle_table = NULL;
  }

  ~con_res() {
    if (type2index)
      delete[] type2index;
    if (bond_table)
      delete[] bond_table[0];
    if (bond_table)
      delete[] bond_table;
    if (angle_table)
      delete[] angle_table[0];
    if (angle_table)
      delete[] angle_table;
  }

  con_res &operator=(const con_res &right) {
    if (this == &right)
      return *this;

    if (type2index)
      delete[] type2index;
    if (bond_table)
      delete[] bond_table[0];
    if (bond_table)
      delete[] bond_table;
    if (angle_table)
      delete[] angle_table[0];
    if (angle_table)
      delete[] angle_table;

    type2index = new char[_HAL_ + 1];
    for (int i = 0; i < _HAL_ + 1; i++)
      type2index[i] = right.type2index[i];
    bond_table = new char *[na];
    bond_table[0] = new char[na * na];
    for (int i = 0; i < na * na; i++)
      bond_table[0][i] = right.bond_table[0][i];
    for (int i = 1; i < na; i++)
      bond_table[i] = right.bond_table[i];

    angle_table = new char *[na];
    angle_table[0] = new char[na * na];
    for (int i = 0; i < na * na; i++)
      angle_table[0][i] = right.angle_table[0][i];
    for (int i = 1; i < na; i++)
      angle_table[i] = right.angle_table[i];
    return *this;
  }

  // copy constructor
  con_res(const con_res &right) : na(right.na) {
    type2index = new char[_HAL_ + 1];
    for (int i = 0; i < _HAL_ + 1; i++)
      type2index[i] = right.type2index[i];
    bond_table = new char *[na];
    bond_table[0] = new char[na * na];
    for (int i = 0; i < na * na; i++)
      bond_table[0][i] = right.bond_table[0][i];
    for (int i = 1; i < na; i++)
      bond_table[i] = right.bond_table[i];

    angle_table = new char *[na];
    angle_table[0] = new char[na * na];
    for (int i = 0; i < na * na; i++)
      angle_table[0][i] = right.angle_table[0][i];
    for (int i = 1; i < na; i++)
      angle_table[i] = right.angle_table[i];
  }

  void init(top_residue &ires, int aa_t) {
    type2index = new char[_HAL_ + 1];
    for (int i = 0; i < _HAL_ + 1; i++)
      type2index[i] = -1;

    na = static_cast<char>(ires.size());

    bond_table = new char *[na];
    bond_table[0] = new char[na * na];
    for (int i = 0; i < na * na; i++)
      bond_table[0][i] = 0;
    for (int i = 1; i < na; i++)
      bond_table[i] = bond_table[i - 1] + na;

    angle_table = new char *[na];
    angle_table[0] = new char[na * na];
    for (int i = 0; i < na * na; i++)
      angle_table[0][i] = 0;
    for (int i = 1; i < na; i++)
      angle_table[i] = angle_table[i - 1] + na;

    // initialize the indexing array
    {
      residue tempRes;
      tempRes.init(residuename(aa_t));
      tempRes.addPolarSH();
      tempRes.addBBH();

      int theNAs = tempRes.getNA() + tempRes.getNH();
      if (theNAs != na) {
        cerr << "error in atom numbers" << endl;
        exit(1);
      }

      for (int i = 0; i < na; i++) {
        int found = 0;
        for (int j = 0; j < tempRes.getNA(); j++) {
          atom *theAtom = tempRes.getAtom(j);
          if (ires.getAtomIndex(theAtom->getName()) == i) {
            found = 1;
            type2index[theAtom->PDBtype] = i;
            break;
          }
          for (int k = 0; k < theAtom->n_attachedProtons; k++) {
            atom *theH = &theAtom->attachedProtons[k];
            if (ires.getAtomIndex(theH->getName()) == i) {
              found = 1;
              type2index[theH->PDBtype] = i;
              break;
            }
          }
        }

        if (!found) {
          cout << ires.PDBname << " " << ires.getAtom(i)->PDBname << endl;
          cout << atomid(ires.getAtom(i)->PDBname) << " " << i << endl;
          exit(1);
        }
      }
    }

    // initialize the bonds
    for (int i = 0; i < ires.nbond; i++) {
      int i1 = ires.bond[i][0];
      int i2 = ires.bond[i][1];
      // cout << i1 << " " << i2 << endl;
      if (i1 == i2) {
        cerr << i1 << " bond with itself " << endl;
        exit(1);
      }
      bond_table[i1][i2] = 1;
      bond_table[i2][i1] = 1;
    }

    // initialize the angles
    for (int i = 0; i < na; i++) {
      for (int j = i + 1; j < na; j++) {
        if (!bond_table[i][j]) { // not a bonding pair
          for (int k = 0; k < na; k++) {
            if (bond_table[i][k] && bond_table[j][k]) {
              angle_table[i][j] = 1;
              angle_table[j][i] = 1;
              break;
            }
          }
        }
      }
    }
  }

  int isBonded(PDBAtype ai_t, PDBAtype aj_t) {
    int i1 = type2index[ai_t];
    int i2 = type2index[aj_t];
    if (i1 < 0 || i2 < 0) {
      cerr << "fatal error: mistake in the atom typing" << endl;
      exit(1);
    }
    return bond_table[i1][i2];
  }

  int isAngled(PDBAtype ai_t, PDBAtype aj_t) {
    int i1 = type2index[ai_t];
    int i2 = type2index[aj_t];
    if (i1 < 0 || i2 < 0) {
      cerr << "fatal error: mistake in the atom typing" << endl;
      exit(1);
    }
    return angle_table[i1][i2];
  }
};

class RES_Connect {
private:
  con_res *con;

public:
  RES_Connect() { con = new con_res[20]; }
  ~RES_Connect() { delete[] con; }

  RES_Connect &operator=(const RES_Connect &right) {
    if (this == &right)
      return *this;
    con = new con_res[20];
    for (int i = 0; i < 20; i++) {
      con[i] = right.con[i];
    }
    return *this;
  }

  RES_Connect(const RES_Connect &right) {
    con = new con_res[20];
    for (int i = 0; i < 20; i++) {
      con[i] = right.con[i];
    }
  }

  void init(const topology &top) {
    for (int i = 0; i < 20; i++) {
      con[i].init(*top.getTopResi(residuename(i)), i);
    }
  }

  int isBonded(AAtype res_t, PDBAtype ai_t, PDBAtype aj_t) {
    if (res_t >= 0 and res_t < 19)
      return con[res_t].isBonded(ai_t, aj_t);
  }

  int isAngled(AAtype res_t, PDBAtype ai_t, PDBAtype aj_t) {
    if (res_t >= 0 and res_t < 19)
      return con[res_t].isAngled(ai_t, aj_t);
  }
};

typedef enum {
  NODES = -1,
  ALLAA = 0,
  NATAA,
  NAROT,
  FIXNR,
  POLAR,
  HYDPH,
  AROMA,
  PIKAA,
} DesignTable_type;

const static char *DesignTable_keywords[] = {
    "ALLAA", "NATAA", "NAROT", "FIXNR", "POLAR", "HYDPH", "AROMA", "PIKAA"};

int DT_Keyword2Int(const char *keyword);
const char *DT_Int2Keyword(int t);

class TableEntry {
private:
  DesignTable_type type;
  vector<AAtype> available_aa;

public:
  TableEntry() { type = NODES; }

  TableEntry(const DesignTable_type &itype) { type = itype; }

  TableEntry(const TableEntry &entry) {
    type = entry.type;
    for (int i = 0; i < entry.available_aa.size(); i++) {
      available_aa.push_back(entry.available_aa[i]);
    }
  }

  ~TableEntry() {}

  void setDT_type(const DesignTable_type &itype) { type = itype; }

  int getDT_type() { return type; }

  void addAA_type(const AAtype &aa) { available_aa.push_back(aa); }

  vector<AAtype> &getAvailableAA() { return available_aa; }

  AAtype getAvailableAA(int ia) { return available_aa[ia]; }

  int getNAvailableAA() { return available_aa.size(); }

  void clearAA_Array() { available_aa.clear(); }

  TableEntry &operator=(const TableEntry &right) {
    type = right.type;
    available_aa.clear();
    for (int i = 0; i < right.available_aa.size(); i++) {
      available_aa.push_back(right.available_aa[i]);
    }
  }
};

class DesignTable {
private:
  vector<TableEntry> table;
  string DT_file;

public:
  DesignTable(int nRES) { table.assign(nRES, TableEntry(NODES)); }

  ~DesignTable() {}

  void readFile(const string &file);

  void print(ostream &out) {
    for (int i = 0; i < table.size(); i++) {
      out << i << " " << DT_Int2Keyword(table[i].getDT_type()) << " ";
      switch (table[i].getDT_type()) {
      case PIKAA:
      case ALLAA:
      case POLAR:
      case HYDPH:
      case AROMA:
        for (int j = 0; j < table[i].getNAvailableAA(); j++) {
          out << residuename(table[i].getAvailableAA()[j]) << " ";
        }
        break;
      } // break;
      out << endl;
    }
  }

  TableEntry &getDT_Entry(int iRes) {
    if (iRes >= 0 && iRes < table.size()) {
      return table[iRes];
    } else {
      cerr << "error in DT: out of range" << endl;
      exit(1);
    }
  }

  int getDT_type(int iRes) {
    if (iRes >= 0 && iRes < table.size()) {
      return table[iRes].getDT_type();
    } else {
      cerr << "error in DT: out of range" << endl;
      exit(1);
    }
  }
};

// to calculate the Dihedral angle
double getDihedralAngle(const vec &p1, const vec &p2, const vec &p3,
                        const vec &p4) {
  double tmp = 0;
  vec v12 = p2 - p1;
  vec v23 = p3 - p2;
  vec v34 = p4 - p3;

  vec t1 = v12 ^ v23;
  vec t2 = v23 ^ v34;

  tmp = (t1 * t2) / sqrt(t1.mod2() * t2.mod2());
  if (tmp >= 1)
    tmp = 1.0;
  if (tmp <= -1)
    tmp = -1.0;

  if ((t1 ^ t2) * v23 > 0)
    return acos(tmp);
  else
    return -acos(tmp);
}

class gprotein : public BaseChain {
private:
  protein &initProtein;
  int length;
  gresidue *gres;
  double *phi;
  double *psi;
  double *omiga;
  aa_rotamer **list;
  int *naa_rot;
  RES_Connect res_con;
  vec *fy_uv;
  vec *omiga_uv;

  friend class gcomplex;

public:
  gprotein(protein &p);
  ~gprotein();

  gprotein &operator=(const gprotein &right) {
    if (this == &right)
      return *this;
    delete[] gres;
    delete[] phi;
    delete[] psi;
    delete[] omiga;
    delete[] list;
    delete[] naa_rot;
    delete[] fy_uv;
    delete[] omiga_uv;
    length = right.length;
    initProtein = right.initProtein;
    phi = new double[length];
    psi = new double[length];
    omiga = new double[length];
    list = new aa_rotamer *[length];
    gres = new gresidue[length];
    naa_rot = new int[length];
    res_con = right.res_con;
    fy_uv = new vec[length * 2];
    omiga_uv = new vec[length];
    for (int i = 0; i < length; i++) {
      gres[i] = right.gres[i];
      phi[i] = right.phi[i];
      psi[i] = right.psi[i];
      omiga[i] = right.omiga[i];
      list[i] = right.list[i];
      naa_rot[i] = right.naa_rot[i];
      omiga_uv[i] = right.omiga_uv[i];
    }
    for (int i = 0; i < length * 2; i++) {
      fy_uv[i] = right.fy_uv[i];
    }
    return *this;
  }

  gprotein(const gprotein &right)
      : length(right.length), initProtein(right.initProtein),
        res_con(right.res_con) {
    phi = new double[length];
    psi = new double[length];
    omiga = new double[length];
    list = new aa_rotamer *[length];
    gres = new gresidue[length];
    naa_rot = new int[length];
    fy_uv = new vec[length * 2];
    omiga_uv = new vec[length];
    for (int i = 0; i < length; i++) {
      gres[i] = right.gres[i];
      phi[i] = right.phi[i];
      psi[i] = right.psi[i];
      omiga[i] = right.omiga[i];
      list[i] = right.list[i];
      naa_rot[i] = right.naa_rot[i];
      omiga_uv[i] = right.omiga_uv[i];
    }
    for (int i = 0; i < length * 2; i++) {
      fy_uv[i] = right.fy_uv[i];
    }
  }
  /*set the force field, mass, charge and sidechain type bit for each atom,
    also, the residue index is assigned. The proton, acceptor array is assigned
    for the calculation of hydrogen bonds interaction.*/
  void updateTopology(const topology &top);
  // the residue index is assigned with a shift in order to index all residues
  // in a complex subResId in an atom is counted for each chain
  void updateTopology(const topology &top, int chainID, int shift);

  bool updateFY(int ires, int &newF, int &newY);

  residue *getResidue(int ires) { return gres[ires].getResidue(); }

  gresidue *getGResidue(int ires) { return gres + ires; }

  /*F->phi, Y->psi*/
  void getFY(int ires, double &iphi, double &ipsi) {
    iphi = phi[ires];
    ipsi = psi[ires];
  }

  // between residue ires, ires+1;
  double &getOmiga(int ires) { return omiga[ires]; }

  void updateFY_UV() {
    for (int i = 0; i < length; i++) {
      residue *theRes = getResidue(i);
      fy_uv[2 * i] = theRes->getCA()->r - theRes->getN()->r;
      fy_uv[2 * i].Normalize();
      fy_uv[2 * i + 1] = theRes->getC()->r - theRes->getCA()->r;
      fy_uv[2 * i + 1].Normalize();
      if (i < length - 1) {
        omiga_uv[i] = getResidue(i + 1)->getN()->r - theRes->getC()->r;
        omiga_uv[i].Normalize();
      }
    }
  }

  vec &getPhi_UV(int ires) { return fy_uv[2 * ires]; }

  vec &getPsi_UV(int ires) { return fy_uv[2 * ires + 1]; }

  vec &getOmiga_UV(int ires) { return omiga_uv[ires]; }

  protein &getProtein() { return initProtein; }
  int size() { return length; }

  double **getM() { return initProtein.getM(); }

  /************************************************************************
  For each residue position, create a list of available amino acid/side-chain
  rotamer. Thus, the search in the space will be facilitated
  *************************************************************************/
  int build_AA_Rotamer_List(BBDep_RotLib &rotlib, BBDep_AA &aalib,
                            const double p_AA = 0.001,
                            const double p_ROT = 0.001);
  void build_AA_Rotamer_List(BBDep_RotLib &rotlib, vector<int> &missing_idx);
  int build_AA_Rotamer_List(BBDep_RotLib &rotlib, BBDep_AA &aalib,
                            DesignTable &dt, BBDep_rotRecord **nativeRotamers,
                            const double p_AA = 0.001,
                            const double p_ROT = 0.001);
  int build_AA_Rotamer_List(BBDep_RotLib &rotlib, BBDep_AA &aalib,
                            DesignTable &dt, const double p_AA = 0.001,
                            const double p_ROT = 0.001);

  int rebuild_AA_Rotamer_List();
  int rebuild_AA_Rotamer_List(const double &vdwr_cutoff);
  void clear_Rotamer_List();

  /*******************
   due to change in phi/psi, the rotamer list needs to be updated
  ****/
  void updateList(BBDep_RotLib &rotlib, const int &iRes, const int &newF,
                  const int &newY);
  void updateList(BBDep_RotLib &rotlib, const int &iRes, const int &newF,
                  const int &newY, BBDep_rotRecord *nativeRotamer);

  int getNList(int iaa) { return naa_rot[iaa]; }
  aa_rotamer *getList(int iaa) { return list[iaa]; }
  void setList(int iaa, vector<aa_rotamer> &pool) {
    naa_rot[iaa] = pool.size();
    delete[] list[iaa];
    list[iaa] = new aa_rotamer[pool.size()];
    for (int j = 0; j < pool.size(); j++) {
      list[iaa][j].init(pool[j]);
    }
  }

  void writePDB(ostream &out);
  void
  writePDB(ostream &out, int startRes,
           char chainID); // write PDB file with given start resid and chain ID
                          /*weather two atoms within a residue is bonded*/
  int isBonded(AAtype res_t, PDBAtype ai_t, PDBAtype aj_t) {
    return res_con.isBonded(res_t, ai_t, aj_t);
  }
  int isAngled(AAtype res_t, PDBAtype ai_t, PDBAtype aj_t) {
    return res_con.isAngled(res_t, ai_t, aj_t);
  }

  int isBonded(int iRes, PDBAtype ai_t, PDBAtype aj_t) {
    // cout << getResidue(iRes)->getType() << endl;
    return res_con.isBonded(static_cast<AAtype>(getResidue(iRes)->getType()),
                            ai_t, aj_t);
  }
  int isAngled(int iRes, PDBAtype ai_t, PDBAtype aj_t) {
    // cout << getResidue(iRes)->getType() << endl;
    return res_con.isAngled(static_cast<AAtype>(getResidue(iRes)->getType()),
                            ai_t, aj_t);
  }

  void rotatePhi(int ires, const double &phi_i) {
    double phi_old = phi[ires];
    double dth = phi_i - phi_old;
    getDimondRotateMatrix(getResidue(ires)->getN()->r,
                          getResidue(ires)->getCA()->r, dth, getM());
    for (int i = 0; i < getResidue(ires)->getNA(); i++) {
      if (i != 0 && i != 1 &&
          !(i == 6 && getResidue(ires)->getType() == PRO)) { // not N, CA
        atom *theAtom = getResidue(ires)->getAtom(i);
        vec tmp = theAtom->r - getResidue(ires)->getCA()->r;
        tmp.EularRotate(getM());
        theAtom->r = tmp + getResidue(ires)->getCA()->r;
        for (int k = 0; k < theAtom->n_attachedProtons; k++) {
          tmp = theAtom->attachedProtons[k].r - getResidue(ires)->getCA()->r;
          tmp.EularRotate(getM());
          theAtom->attachedProtons[k].r = tmp + getResidue(ires)->getCA()->r;
        }
      }
    }
    if (getResidue(ires)->getType() == PRO)
      getResidue(ires)->updateRotamers();
    for (int i = ires + 1; i < size(); i++) {
      for (int j = 0; j < getResidue(i)->getNA(); j++) {
        atom *theAtom = getResidue(i)->getAtom(j);
        vec tmp = theAtom->r - getResidue(ires)->getCA()->r;
        tmp.EularRotate(getM());
        theAtom->r = tmp + getResidue(ires)->getCA()->r;
        for (int k = 0; k < theAtom->n_attachedProtons; k++) {
          tmp = theAtom->attachedProtons[k].r - getResidue(ires)->getCA()->r;
          tmp.EularRotate(getM());
          theAtom->attachedProtons[k].r = tmp + getResidue(ires)->getCA()->r;
        }
      }
    }

    phi[ires] = phi_i;
  }

  void rotatePsi(int ires, const double &psi_i) {
    double psi_old = psi[ires];
    double dth = psi_i - psi_old;
    getDimondRotateMatrix(getResidue(ires)->getCA()->r,
                          getResidue(ires)->getC()->r, dth, getM());

    vec tmp = getResidue(ires)->getO()->r - getResidue(ires)->getC()->r;
    tmp.EularRotate(getM());
    getResidue(ires)->getO()->r = tmp + getResidue(ires)->getC()->r;

    for (int i = ires + 1; i < size(); i++) {
      for (int j = 0; j < getResidue(i)->getNA(); j++) {
        atom *theAtom = getResidue(i)->getAtom(j);
        vec tmp = theAtom->r - getResidue(ires)->getC()->r;
        tmp.EularRotate(getM());
        theAtom->r = tmp + getResidue(ires)->getC()->r;
        for (int k = 0; k < theAtom->n_attachedProtons; k++) {
          tmp = theAtom->attachedProtons[k].r - getResidue(ires)->getC()->r;
          tmp.EularRotate(getM());
          theAtom->attachedProtons[k].r = tmp + getResidue(ires)->getC()->r;
        }
      }
    }

    psi[ires] = psi_i;
  }

  void rotateOmiga(int ires, double &omiga_i) {
    if (ires < size() - 1) {
      double omiga_old = omiga[ires];
      double dth = omiga_i - omiga_old;
      getDimondRotateMatrix(getResidue(ires)->getC()->r,
                            getResidue(ires + 1)->getN()->r, dth, getM());
      for (int i = ires + 1; i < size(); i++) {
        for (int j = 0; j < getResidue(i)->getNA(); j++) {
          atom *theAtom = getResidue(i)->getAtom(j);
          vec tmp = theAtom->r - getResidue(ires + 1)->getN()->r;
          tmp.EularRotate(getM());
          theAtom->r = tmp + getResidue(ires + 1)->getN()->r;
          for (int k = 0; k < theAtom->n_attachedProtons; k++) {
            tmp =
                theAtom->attachedProtons[k].r - getResidue(ires + 1)->getN()->r;
            tmp.EularRotate(getM());
            theAtom->attachedProtons[k].r =
                tmp + getResidue(ires + 1)->getN()->r;
          }
        }
      }
      omiga[ires] = omiga_i;
    }
  }

  /************************
  BackRub motion.
  Ref: The Backrub Mition: How Protein Backbone Shrugs when a sidechain
   Dances. --Davis, IW; Arendal II, WE; Richardson&Richardson; Structure 2006.
   Rotate around the axis of CAs of residues ires-1 and ires+1.
   dtheta is in RADIAN.

  Theta13 is defined as: Ni-1, CAi-1, CAi+1, CAi
  Theta12 is defined as: CAi+1, CAi-1, CAi, Ni
  Theta32 is defined as: CAi-1, CAi+1, CAi, Ci
  */
  void rotateTheta13(int ires, const double &dtheta);

  /*Rotate around the axis of CAs of residues ires-1 and ires.
    dtheta is in RADIAN. */
  void rotateTheta12(int ires, const double &dtheta);

  /*Rotate around the axis of CAs of residues ires+1 and ires.
    dtheta is in RADIAN. */
  void rotateTheta32(int ires, const double &dtheta);

  /***************************************************
   during the backbone minimization, the backbone is moved around. However,
  since the protein(generalized protein) has copies of other residues. It is
  necessary to move other residues around with respect to the active residues.
  ****************************************************/
  void updateGResidues() {
    for (int i = 0; i < size(); i++) {
      gres[i].updateGResidue();
    }
  }

  int isNAN() {
    int nan_res;
    for (int i = 0; i < size(); i++) {
      if (gres[i].isNAN()) {
        // cout << "error in residue: " << i << endl;
        return 1;
      }
    }
    return 0;
  }

  void debugOut() {
    for (int ires = 0; ires < size(); ires++) {
      residue *theRes = getResidue(ires);
      cout << "ires = " << ires << "  " << theRes->getName() << endl;
      for (int iHA = 0; iHA < theRes->getNA(); iHA++) {
        atom *theAtom = theRes->getAtom(iHA);
        printf("%f %f %f \n", theAtom->r.x, theAtom->r.y, theAtom->r.z);
        for (int iP = 0; iP < theAtom->n_attachedProtons; iP++) {
          atom theH = theAtom->attachedProtons[iP];
          printf("\t %f %f %f\n", theH.r.x, theH.r.y, theH.r.z);
        }
      }
    }
  }

  void getBoundary(vec &min, vec &max) {
    min = vec(INF, INF, INF);
    max = vec(-INF, -INF, -INF);
    vec mintmp, maxtmp;
    for (int i = 0; i < size(); i++) {
      gres->getBoundary(mintmp, maxtmp);
      if (mintmp.x < min.x)
        min.x = mintmp.x;
      if (mintmp.y < min.y)
        min.y = mintmp.y;
      if (mintmp.z < min.z)
        min.z = mintmp.z;
      if (maxtmp.x > max.x)
        max.x = maxtmp.x;
      if (maxtmp.y > max.y)
        max.y = maxtmp.y;
      if (maxtmp.z > max.z)
        max.z = maxtmp.z;
    }
  }

  /// Idealize the bond and angle. All values are copied from template.
  void idealize();
  void setAtomIndex(int &index) {}

  bool updatePhiPsi() {
    for (int i = 0; i < size(); i++) {
      if (i > 0) {
        phi[i] = getDihedralAngle(gres[i - 1].getResidue()->getC()->r,
                                  gres[i].getResidue()->getN()->r,
                                  gres[i].getResidue()->getCA()->r,
                                  gres[i].getResidue()->getC()->r);
      }
      if (i < size() - 1) {
        psi[i] = getDihedralAngle(gres[i].getResidue()->getN()->r,
                                  gres[i].getResidue()->getCA()->r,
                                  gres[i].getResidue()->getC()->r,
                                  gres[i + 1].getResidue()->getN()->r);

        omiga[i] = getDihedralAngle(gres[i].getResidue()->getCA()->r,
                                    gres[i].getResidue()->getC()->r,
                                    gres[i + 1].getResidue()->getN()->r,
                                    gres[i + 1].getResidue()->getCA()->r);
      }
    }
  }

  void getPhiPsi(int i, double &f, double &y) {
    if (i >= 0 && i < length) {
      f = phi[i];
      y = psi[i];
    }
  }

  void setPhiPsi(int i, double &f, double &y) {
    if (i >= 0 && i < length) {
      phi[i] = f;
      psi[i] = y;
    }
  }

  void setOmiga(int i, double &o) {
    if (i >= 0 && i < length)
      omiga[i] = o;
  }

}; // end gresidue class

class grid {
protected:
  vec max;                   // maximum x,y,z
  vec min;                   // minimum x,y,z
  vec cell_box;              // dimension of unit cell along different axis
  int ncell[3];              // number of cells along each axises
  protein &p;                // the protein
  double max_ir;             // the maximum interactino range
  const static int step = 2; // each max_ir is roughly devied by 2;
  double max_ir2;            // the squre of the max interaction range

  int ncells;        // number of cells
  atom **cell_atoms; // the linked list of atoms in each cell
  VDW &vdw;          // the VDW intertation tables
  // double** EEF1_SOLV;
  const vector<vector<double>> &EEF1_SOLV;

public:
  grid(protein &theProtein, VDW &theVDW)
      : p(theProtein), vdw(theVDW), EEF1_SOLV(theVDW.getEEF1_SOLV()) {
    // EEF1_SOLV = vdw.getEEF1_SOLV();
  }
  grid(protein &theProtein, VDW &theVDW, double mir);
  ~grid() { delete[] cell_atoms; }

  void constructCells();
  void dropCellResidue(int iRes);
  void insertCellResidue(int iRes);
  int updateCellResidue(int iRes);

  int getNCell() { return ncells; }
  atom *getCellAtomsOf(int i) { return cell_atoms[i]; }

  /**************************************************************
   Get the energy related to ith atom. For the SOLVation energy calculation,
   only the dGi = dGi0 - sum_j{f(ij)Vj} is calculated.
   For VDW energy calculation, all other terms are taken considered.
   Thus, the total vdw will be doubled if all of the i is summed while
   solvation energy is not!!!!!!!!!!!!!! Use with care!
  ***************************************************************/
  double getEVS_iatom(atom *ai, double &E_VDW_attr, double &E_VDW_rep,
                      double &E_SOLV);

  /****************************************************************
   Get the energy related to ith atom with j>i.
   For solvation energy, dGi' = dGi0-sum_j>i{f(ij)Vj+f(ji)Vi}.
   Thus, the total energy can be simply summed for all terms
   using this fuction!
  ****************************************************************/
  double getEVS_iatom_half(atom *ai, double &E_VDW_attr, double &E_VDW_rep,
                           double &E_SOLV);

  /***************************************************************
   Get the free energy with respect to the sidechain of a specific residue.
   The internal energy of the residue is not included!!!! The internal energy
   will be taken into account for the backbone dependent rotermer
   energies.
  ***************************************************************/
  double getEVS_RES_SC(int iRes, double &E_VDW_attr, double &E_VDW_rep,
                       double &E_SOLV);
  /*The function to indentify the internally connected atoms, the first argument
    is the current sidechain atoms, the second argement is the atom to be
    determined. */
  int isIncluded(atom *theAtom, atom *pt) {
    if (abs(theAtom->resID - pt->resID) >= 2)
      return 1;
    else {
      if (theAtom->resID == pt->resID)
        return 0;
      else if (theAtom->resID == pt->resID + 1 &&
               (pt->PDBtype == _CA_ || pt->PDBtype == _C_ ||
                pt->PDBtype == _O_)) // previous residue
        return 0;
      else if (theAtom->resID == pt->resID - 1 &&
               (pt->PDBtype == _CA_ || pt->PDBtype == _N_ ||
                pt->PDBtype == _HN_ || pt->PDBtype == _HN1_ ||
                pt->PDBtype == _HN2_)) // next residue
        return 0;
      else
        return 1;
    }
  }
  /*Entened, return
    1, if the two atoms are in the following
    CAi-CO-NH-CAi+1-CO-NH-CAi+2
              |
              CRi+1
    0, otherwise
    do not depend on the order in the input atom pair */
  int isInternal(atom *iatom, atom *jatom);

  double getEVS() {
    double VDW_rep = 0, VDW_attr = 0, SOLV = 0;
    double vdw_rep, vdw_attr, solv;
    for (int i = 0; i < p.size(); i++) {
      residue *theRes = p.getResidue(i);
      atom *theAtom = NULL;
      for (int j = 0; j < theRes->getNA(); j++) {
        theAtom = theRes->getAtom(j);
        getEVS_iatom_half(theAtom, vdw_attr, vdw_rep, solv);
        VDW_rep += vdw_rep;
        VDW_attr += vdw_attr;
        SOLV += solv;
      }
    }
    cout << VDW_rep << " " << VDW_attr << " " << SOLV << endl;
    // VDW_rep=VDW_attr=SOLV=0;
    // for(int i=0; i<p.size(); i++){
    //   residue* theRes = p.getResidue(i);
    //   atom* theAtom = NULL;
    //   for(int j=0; j<theRes->getNA(); j++){
    //	theAtom = theRes->getAtom(j);
    //	getEVS_iatom(theAtom, vdw_attr, vdw_rep, solv);
    //	VDW_rep += vdw_rep;
    //	VDW_attr+= vdw_attr;
    //	SOLV    += solv;
    //   }
    // }
    // cout << VDW_rep/2 << " " << VDW_attr/2 << " " << SOLV << endl;
  }

  // debug
  void printCells(int index) {
    atom *pt = cell_atoms[index];
    while (pt) {
      cout << pt->getName();
      pt->r.print(cout);
      cout << endl;
      pt = pt->cell_next;
    }
  }
};

class gen_grid : public grid {
private:
  gprotein &gp;

public:
  gen_grid(protein &theProtein, gprotein &theGProtein, VDW &theVDW, double mir);
  ~gen_grid() {}

  void dropCellResidue(int iRes);
  void insertCellResidue(int iRes);
  int updateCellResidue(int iRes);
  void updateCellAtom(atom *theAtom);

  double getEVS_Atom(atom *theAtom, double &E_VDW_attr, double &E_VDW_rep,
                     double &E_SOLV);

  /***************************************************************************
    The VDW,SOLV are arrays to stored the reside pairwise interactions.
   ***************************************************************************/

  double getEVS_RES_SC(int iRes, double &E_VDW_attr, double &E_VDW_rep,
                       double &E_SOLV, double *VDWA_SS = NULL,
                       double *VDWR_SS = NULL, double *SOLV_SS = NULL,
                       double *VDWA_SB = NULL, double *VDWR_SB = NULL,
                       double *SOLV_SB = NULL);

  /*include the SB interaction*/
  double getEVS_RES_SC(int iRes, double &E_VDW_attr, double &E_VDW_rep,
                       double &E_SOLV, double &E_SB, double *VDWA_SS = NULL,
                       double *VDWR_SS = NULL, double *SOLV_SS = NULL,
                       double *VDWA_SB = NULL, double *VDWR_SB = NULL,
                       double *SOLV_SB = NULL, double *SB_SS = NULL);

  /*Get the VDW, SOLV energy and it derivative to the chi angles with respect to
   * the sidechain*/
  double getE_DEV_RES_SC(int iRes, double &E_VDW_attr, double &E_VDW_rep,
                         double &E_SOLV, double *DEV_VDW_attr,
                         double *DEV_VDW_rep, double *DEV_SOLV);

  /*include the SB interaction*/
  double getE_DEV_RES_SC(int iRes, double &E_VDW_attr, double &E_VDW_rep,
                         double &E_SOLV, double &E_SB, double *DEV_VDW_attr,
                         double *DEV_VDW_rep, double *DEV_SOLV, double *DEV_SB);

  /***************************************************************
    obtain the VDW and SOLV for GLY and PRO with respect to a
    reference aa,

    N(NZ with on proton) --- CA(CA1 with one proton) --- CO

    for GLY, CA1--> CA2;
    for PRO, NZ --> NZNP;
   ****************************************************************/
  void getEVS_SPECIAL_RES(int iRes, double &E_VDW_attr, double &E_VDW_rep,
                          double &E_SOLV, double *VDW_A = NULL,
                          double *VDW_R = NULL, double *SOLV = NULL);

  double getE_DEV(double &E_VDWA, double &E_VDWR, double &E_SOLV, double *DVDWA,
                  double *DVDWR, double *DSOLV);

  /*calculate HBOND DERIVATIVE using Abe&Go trick*/
  double getE_DEV_GOTRICK(double &E_VDWA, double &E_VDWR, double &E_SOLV,
                          double *DVDWA, double *DVDWR, double *DSOLV);

  /*calculate HBOND DERIVATIVE using Abe&Go trick, using only one loop*/
  double getE_DEV_GOTRICK1(double &E_VDWA, double &E_VDWR, double &E_SOLV,
                           double *DVDWA, double *DVDWR, double *DSOLV);

  /*extend the FY to FYO, including Omiga angle*/
  double getE_DEV_FYO_GOTRICK(double &E_VDWA, double &E_VDWR, double &E_SOLV,
                              double *DVDWA, double *DVDWR, double *DSOLV);

  void getGF_GOTRICK(atom *theAtom, vec &F_VDWR, vec &F_VDWA, vec &F_SOLV,
                     vec &G_VDWR, vec &G_VDWA, vec &G_SOLV, double &SOLVi,
                     double &E_VDWR, double &E_VDWA);
  /* extended to include all sidechain angles*/
  double getE_DEV_FYOChi_GOTRICK(double &E_VDWA, double &E_VDWR, double &E_SOLV,
                                 double *DVDWA, double *DVDWR, double *DSOLV,
                                 int nDihedral, int *INDEX_MAPPING);

  double getE(double &E_VDWA, double &E_VDWR, double &E_SOLV);

  /*compute the derivate for backrub motion*/
  void getDE_BACKRUB(const int ires, double &E_VDWA, double &E_VDWR,
                     double &E_SOLV, double *DVDWA, double *DVDWR,
                     double *DSOLV, vec *backrub_uv);
};

class fine_grid {
protected:
  vec max;
  vec min;
  vec cell_box;
  int ncell[3];
  protein &p;

  double max_r;
  double max_r2;

  int ncells;
  atom **cell_acceptors;
  atom **cell_protons;

  void constructCells();

public:
  fine_grid(protein &theProtein, int DUMMY) : p(theProtein) {}
  fine_grid(protein &theProtein);
  ~fine_grid() {
    delete[] cell_protons;
    delete[] cell_acceptors;
  }

  void updateFineCellResidue(int iRes);
  void dropFineCellResidue(int iRes);
  void insertFineCellResidue(int iRes);

  void getHBE_RES(int iRes, double &EHB_bb_bb, double &EHB_bb_sc,
                  double &EHB_sc_sc);
  void clearHBE_RES(int iRes, double &EHB_bb_bb, double &EHB_bb_sc,
                    double &EHB_sc_sc);
  void clearSC_HBE_RES(int iRes, double &EHB_bb_sc, double &EHB_sc_sc);
};

class gen_fine_grid : public fine_grid {
private:
  gprotein &gp;

public:
  gen_fine_grid(protein &theProtein, gprotein &theGProtein);
  ~gen_fine_grid() {}

  void updateFineCellProton(atom *);
  void updateFineCellAcceptor(atom *);

  void updateFineCellResidue(int iRes);
  void dropFineCellResidue(int iRes);
  void insertFineCellResidue(int iRes);

  void getHBE_Proton(atom *theH, double &EHB_bb_bb, double &EHB_bb_sc,
                     double &EHB_sc_sc);
  void getHBE_Acceptor(atom *theA, double &EHB_bb_bb, double &EHB_bb_sc,
                       double &EHB_sc_sc);
  void getHBE_RES(int iRes, double &EHB_bb_bb, double &EHB_bb_sc,
                  double &EHB_sc_sc);

  void clearHBE_Proton(atom *theH, double &EHB_bb_bb, double &EHB_bb_sc,
                       double &EHB_sc_sc);
  void clearHBE_Acceptor(atom *theA, double &EHB_bb_bb, double &EHB_bb_sc,
                         double &EHB_sc_sc);
  void clearHBE_RES(int iRes, double &EHB_bb_bb, double &EHB_bb_sc,
                    double &EHB_sc_sc);

  void clearSC_HBE_RES(int iRes, double &EHB_bb_sc, double &EHB_sc_sc);
  /*obtain the derivative also!*/
  void getHBE_RES_DEV(int iRes, double &EHB_bb_bb, double &EHB_bb_sc,
                      double &EHB_sc_sc, double *DEV_bb_sc, double *DEV_sc_sc);
  /*clear the hydrogen bond status table*/
  void clearHBE();
  /*obtain the derivative for whole structure*/
  void getEHB_DEV(double &EHB_bb_bb, double &EHB_bb_sc, double &EHB_sc_sc,
                  double *DEV_bb_bb, double *DEV_bb_sc, double *DEV_sc_sc);
  /*extend the derive with respect to FYO instead of FY*/
  void getEHB_FYO_DEV(double &EHB_bb_bb, double &EHB_bb_sc, double &EHB_sc_sc,
                      double *DEV_bb_bb, double *DEV_bb_sc, double *DEV_sc_sc);
  /*obtain the energy*/
  void getEHB(double &EHB_bb_bb, double &EHB_bb_sc, double &EHB_sc_sc);

  /*calculate the repulsion between two-close protons*/
  void getHH_EVR_RES_SC(int iRes, double &VDW_REP, double *VDWR_SS_1D = NULL,
                        double *VDWR_SB_1D = NULL);

  /*calculate the derative of the repulsion between protons*/
  void getHH_EVR_DEV_SC(int iRes, double &VDWR, double *DEV_VDWR);

  /*compute the DERIVATIVE of hydrogenbonds for BACKRUB motion*/
  void getEHB_DEV_BACKRUB(const int ires, double &EHB_bb_bb, double &EHB_bb_sc,
                          double &EHB_sc_sc, double *DEV_bb_bb,
                          double *DEV_bb_sc, double *DEV_sc_sc,
                          vec *backrub_uv);
};

int isComment(const char *line) {
  int pt = 0;
  while ((line[pt] == ' ' || line[pt] == '\t') && line[pt] != '\0') {
    pt++;
  }

  if (line[pt] == '#') {
    return 1;
  }

  if (line[pt] == '/' && line[pt + 1] == '/') {
    return 1;
  }

  return 0;
}

typedef enum {
  VDW_ATTR_W = 0,
  VDW_REP_W,
  SOLV_EEF1_W,
  SB_W,
  HB_BB_BB_W,
  HB_BB_SC_W,
  HB_SC_SC_W,
  FY_AA_W,
  FY_AA_CHI_W,
  CYS_W,
  MET_W,
  PHE_W,
  ILE_W,
  LEU_W,
  VAL_W,
  TRP_W,
  TYR_W,
  ALA_W,
  GLY_W,
  THR_W,
  SER_W,
  GLN_W,
  ASN_W,
  GLU_W,
  ASP_W,
  HIS_W,
  ARG_W,
  LYS_W,
  PRO_W,
  WATER_HB_W,
  CONS_W
} WEIGHT_t;

static const int nCOEF1 = 31;

class WEIGHT_COEF {
private:
  double *coef;

public:
  WEIGHT_COEF() {
    coef = new double[nCOEF1];
    for (int i = 0; i < nCOEF1; i++)
      coef[i] = 0;
  }

  ~WEIGHT_COEF() { delete[] coef; }

  int readFile(const char *file) {
    char buf[1024];
    int i = 0;
    ifstream in(file, ios::in);
    while (in.getline(buf, 1024)) {
      if (i >= nCOEF1)
        return 0;
      sscanf(buf, "%lf", &coef[i]);
      i++;
    }
    return 1;
  }

  double &operator()(int t) { return coef[t]; }

  const double &operator()(int t) const { return coef[t]; }

  double en_base(double vdwa, double vdwr, double solv, double esb,
                 double cons) {
    return vdwa + vdwr * coef[VDW_REP_W] + solv * coef[SOLV_EEF1_W] +
           esb * coef[SB_W] + cons * coef[CONS_W];
  }

  double en_vdw(double vdwa, double vdwr) {
    return vdwa + vdwr * coef[VDW_REP_W];
  }

  double en_vdw(double vdwa, double vdwr, double solv) {
    return vdwa + vdwr * coef[VDW_REP_W] + solv * coef[SOLV_EEF1_W];
  }

  double en_hb(double bb_bb, double bb_sc, double sc_sc) {
    return bb_bb * coef[HB_BB_BB_W] + bb_sc * coef[HB_BB_SC_W] +
           sc_sc * coef[HB_SC_SC_W];
  }

  double en_hb(double bb_sc, double sc_sc) {
    return bb_sc * coef[HB_BB_SC_W] + sc_sc * coef[HB_SC_SC_W];
  }

  double en_fy(double aa, double chi) {
    return aa * coef[FY_AA_W] + chi * coef[FY_AA_CHI_W];
  }

  double *getWeight() { return coef; }

  int nCoefs() { return nCOEF1; }

  void correctREF(const char *file) {
    ifstream in(file, ios::in);
    string line;
    char res[10];
    double ref;
    int rid;
    while (getline(in, line)) {
      if (!isComment(line.c_str())) {
        sscanf(line.c_str(), "%s %lf", res, &ref);
        if ((rid = residueid(res)) < 0) {
          cerr << "REFERENCE ENERGY FILE: wrong residue name" << endl;
          exit(1);
        }
        coef[CYS_W + rid] += ref;
      }
    }
    in.close();
  }
};

struct atom_mass_t {
  int FFTindex = -1;
  double atm_mass = INF;
  string FFname;

  void init(const string &iFFname, double imass) {
    FFname = iFFname;
    atm_mass = imass;
  }

  //   void setFFtype(int i) { FFTindex = i; }
  //   int getFFtype() const { return FFTindex; }

  //   void setMass(double imass) { atm_mass = imass; }
  //   double getMass() const { return atm_mass; }
};

struct top_atom {
  int FFtype = -1;
  double charge = INF;
  vec ic;
  string PDBname;

  void init(const string &name, int type, double q) {
    PDBname = name;
    FFtype = type;
    charge = q;
  }

  //   void setType(int c) { FFtype = c; }
  //   int getType() { return FFtype; }
  //   void setCharge(double q) { charge = q; }
  //   double getCharge() { return charge; }
  //   char *getName() { return PDBname; }
  //   void setIC(double x, double y, double z) { ic = vec(x, y, z); }
  //   const vec &getIC() { return ic; }
};

struct top_residue {
  string PDBname;
  top_atom *a = NULL;
  int natom = 0;

  int nbond = 0;
  int **bond; // bond[i][0]--bond[i][1]
  int ndihedral = 0;
  int **dihedral; // diheral[i][0]-dih[i][1]-dih[i][2]-dih[i][3]
  int nimproper = 0;
  int **improper;

  top_residue(const string &name, int n) { init(name, n); }

  ~top_residue() {
    if (a != NULL) {
      delete[] a;
    }
    if (nbond) {
      delete[] bond[0];
      delete[] bond;
    }
    if (ndihedral) {
      delete[] dihedral[0];
      delete[] dihedral;
    }
    if (nimproper) {
      delete[] improper[0];
      delete[] improper;
    }
  }

  void init(const string &name, int n) {
    PDBname = name;
    natom = n;
    a = new top_atom[n];
  }

  int getNAtom() { return natom; }

  top_atom *getAtom(int i) { return a + i; }

  int size() { return natom; }

  //   char *getName() { return PDBname; }

  int getIndex(const char *pdbname) {
    for (int i = 0; i < natom; i++) {
      if (a[i].PDBname == pdbname)
        //   if (!strcmp(a[i].PDBname, pdbname))
        return i;
    }
    return -1;
  }

  int top_residue::getAtomIndex(const string &atm_name) {
    for (int j = 0; j < size(); j++) {
      if (getAtom(j)->PDBname == atm_name)
        //   if (!strcmp(getAtom(j)->getName(), atm_name)) // corr atm name
        return j;
    }
    string aname;
    aname = atm_name;
    uint len = aname.size();
    // int len = strlen(aname);
    if (aname[len - 1] >= '0' && aname[len - 1] <= '9') { // the last digit
      //   aname[len - 1] = '\0';
      aname.pop_back();
      int nmatch = 0;
      int atom_index = -1;
      for (int j = 0; j < size(); j++) {
        if (getAtom(j)->PDBname == aname) { // corr atm name
          nmatch++;
          atom_index = j;
        }
      }
      if (nmatch == 0) {
        cerr << "No match of atom '" << atm_name << "' in '" << PDBname << "'"
             << endl;
        return -1;
      } else if (nmatch > 1) {
        cerr << "Multiple match of atom '" << atm_name << "' in '" << PDBname
             << "'" << endl;
        return -1;
      } else {
        return atom_index;
      }
    }
  }
  //   int getNBonds() { return nbond; }
  //   int **getBonds() { return bond; }

  //   int getNDihedrals() { return ndihedral; }
  //   int **getDihedrals() { return dihedral; }

  //   int getNImpropers() { return nimproper; }
  //   int **getImpropers() { return improper; }
};

struct topology {
  char file[100];
  // atom_mass_t* mass_table;
  vector<atom_mass_t> mass_table;
  int ntype = 0;
  vector<top_residue *> resi_table;

  topology(const string &name) { addTopology(name); }

  ~topology() {
    // delete [] mass_table;
    for (int i = 0; i < resi_table.size(); i++)
      delete resi_table[i];
    resi_table.clear();
  }

  void addTopology(const string &name);

  int fftid(const string &fftype) const;
  double getMass(int fftid) const { return mass_table[fftid].atm_mass; }

  int getNTopRecords() { return resi_table.size(); }

  void addRes(top_residue *newRes) { resi_table.push_back(newRes); }

  top_residue *getTopResi(int i) const {
    if (i >= 0 && i < resi_table.size())
      return resi_table[i];
    return NULL;
  }

  top_residue *getTopResi(const char *res_name) const {
    int res_index = -1;
    for (int i = 0; i < resi_table.size(); i++) {
      if (resi_table[i]->PDBname == res_name) {
        res_index = i;
        break;
      }
    }
    if (res_index == -1) {
      cerr << "No match of residue name:" << res_name << endl;
      return NULL;
    }
    return resi_table[res_index];
  }

  int getTopResIndex(const char *res_name) const {
    int res_index = -1;
    for (int i = 0; i < resi_table.size(); i++) {
      if (resi_table[i]->PDBname == res_name) {
        res_index = i;
        break;
      }
    }
    if (res_index == -1) {
      cerr << "No match of residue name:" << res_name << endl;
      return -1;
    }
    return res_index;
  }

  // atom_mass_t* getMassTable()const{return mass_table;}
  //   const vector<atom_mass_t> &getMassTable() const { return mass_table; }

  //   int getNType() const { return ntype; }
  // return the atom type, mass and chage for a given atom of a given residue
  // return 0,if no match
  int getAtomInfor(const string &res_name, const string &atm_name, int &type,
                   double &mass, double &charge) const;
};

class rotRecord {
private:
  int iChi[4];
  double dChi[4];
  double sChi[4];
  double p, sp;
  double cp1, scp1;

public:
  rotRecord() {
    p = INF;
    cp1 = INF;
    for (int i = 0; i < 4; i++) {
      iChi[i] = 0;
      dChi[i] = INF;
      sChi[i] = INF;
    }
  }

  ~rotRecord(){};

  void setIChi(int i1, int i2, int i3, int i4) {
    iChi[0] = i1;
    iChi[1] = i2;
    iChi[2] = i3;
    iChi[3] = i4;
  }

  void setChi(double &chi1, double sig1, double &chi2, double sig2,
              double &chi3, double sig3, double &chi4, double sig4) {
    dChi[0] = chi1;
    sChi[0] = sig1;
    dChi[1] = chi2;
    sChi[1] = sig2;
    dChi[2] = chi3;
    sChi[2] = sig3;
    dChi[3] = chi4;
    sChi[3] = sig4;
  }

  void setP(double ip, double isp = 0) {
    p = ip;
    sp = isp;
  }
  void setCP1(double icp1, double iscp1 = 0) {
    cp1 = icp1;
    scp1 = iscp1;
  }
  double getP() { return p; }
  double getCP1() { return cp1; }
  const int *getIChi() const { return iChi; }
  const double *getDChi() const { return dChi; }
  const double *getSChi() const { return sChi; }

  void print(ostream &out) {
    char buf[1000];
    int pt = sprintf(buf,
                     "%1d %1d %1d %1d                %6.2lf %6.2lf  %6.2lf "
                     "%6.2lf %6.1lf %4.1lf",
                     iChi[0], iChi[1], iChi[2], iChi[3], p, sp, cp1, scp1,
                     dChi[0], sChi[0]);
    for (int i = 0; i < 4; i++) {
      if (iChi[i]) {
        pt += sprintf(buf + pt, " %8.1lf %4.1lf", dChi[i], sChi[i]);
      }
    }
    out << buf << endl;
  }
};

class BBDep_AA {
private:
  double ***fy_aa;
  double ***e_stat;
  void updateCAP();

public:
  BBDep_AA();
  ~BBDep_AA();

  void readFile(const char *file);
  void updateStatisticalEnergy();

  double getE_STAT(int iF, int iY, int ia) { return e_stat[iF][iY][ia]; }
};

const static double FY_MIN = -180.0;
const static double FY_MAX = 180.0;
const static int FY_NBIN = 36;
const static double FY_BIN = 10.0;

class BBDep_rotRecord {
private:
  unsigned char iChi[4];
  double aveChi[4];
  double sigChi[4];
  double p;
  double E_STAT;

public:
  BBDep_rotRecord() {
    p = INF;
    E_STAT = INF;
    for (int i = 0; i < 4; i++) {
      iChi[i] = 0;
      aveChi[i] = INF;
      sigChi[i] = INF;
    }
  }

  ~BBDep_rotRecord() {}

  void setIChi(int i1, int i2, int i3, int i4) {
    iChi[0] = i1;
    iChi[1] = i2;
    iChi[2] = i3;
    iChi[3] = i4;
  }

  void setChi(const double &chi1, const double &sig1, const double &chi2,
              const double &sig2, const double &chi3, const double &sig3,
              const double &chi4, const double &sig4) {
    aveChi[0] = chi1;
    sigChi[0] = sig1;
    aveChi[1] = chi2;
    sigChi[1] = sig2;
    aveChi[2] = chi3;
    sigChi[2] = sig3;
    aveChi[3] = chi4;
    sigChi[3] = sig4;
  }

  void setP(double ip) { p = ip; }
  void setMinimumP() {
    p = 1E-6;
    E_STAT = -log(p);
  }
  double getP() { return p; }
  const unsigned char *getIChi() const { return iChi; }
  const double *getAveChi() const { return aveChi; }
  const double *getSigChi() const { return sigChi; }

  void setE_STAT(double e) { E_STAT = e; }
  void setE_STAT() {
    if (p != INF && p != 0) {
      E_STAT = -log(p);
      // cout << E_STAT << endl;
    }
  }
  double get_ESTAT() { return E_STAT; }
  friend class BBDep_RotLib;
};

class BBDep_resRotLib {
private:
  int id;
  int nRecords;
  int nRotamers;
  BBDep_rotRecord *rotamers;
  int base[4];

public:
  BBDep_resRotLib() {
    id = -1;
    nRecords = 0;
    nRotamers = 0;
    rotamers = NULL;
  }

  BBDep_resRotLib(int type, int nRot) {
    id = type;
    nRecords = 0;
    nRotamers = nRot;
    rotamers = new BBDep_rotRecord[nRot];
  }

  ~BBDep_resRotLib() { delete[] rotamers; }

  void init(int type, int nRot) {
    id = type;
    nRecords = 0;
    nRotamers = nRot;
    rotamers = new BBDep_rotRecord[nRot];
  }

  void setNRecords(int iRecords) { nRecords = iRecords; }
  int getNRecords() { return nRecords; }

  int getNRotamers() { return nRotamers; }
  BBDep_rotRecord *getRotamers() { return rotamers; }

  BBDep_rotRecord *getRotamer(int i) {
    if (i < 0 && i >= nRotamers)
      return NULL;
    return rotamers + i;
  }

  void setE_STAT() {
    for (int i = 0; i < nRotamers; i++) {
      rotamers[i].setE_STAT();
    }
  }

  void ordering() {
    /*find the base for each integer chi*/
    for (int i = 0; i < 4; i++)
      base[i] = 0;
    for (int irot = 0; irot < nRotamers; irot++) {
      const unsigned char *pt = rotamers[irot].getIChi();
      for (int i = 0; i < 4; i++) {
        if (pt[i] > base[i])
          base[i] = pt[i];
      }
    }
    // cout << residuename(id) << endl;
    // for(int i=0; i<4; i++) cout << base[i] << " " ;
    // cout << endl;
    for (int i = 0; i < nRotamers; i++) {
      const unsigned char *pt = rotamers[i].getIChi();
      int add = 0;
      for (int k = 0; k < 4; k++) {
        if (pt[k] > 0) {
          add = add * base[k] + pt[k] - 1;
        }
      }
      if (i != add) {
        cout << "fatal error in BBDep RotLib" << endl;
        exit(1);
      }
    }
  }

  BBDep_rotRecord *getRotamer(const int &nCHI, const unsigned char *iCHI) {
    int add = 0;
    for (int i = 0; i < nCHI; i++) {
      add = add * base[i] + (iCHI[i] - 1);
    }
    if (add > nRotamers) {
      cerr << "fatal error" << endl;
      exit(1);
    }
    return rotamers + add;
  }

  BBDep_rotRecord *getRotamer(const unsigned char *iCHI) {
    int add = 0;
    for (int i = 0; i < 4; i++) {
      if (iCHI[i] <= 0)
        break;
      add = add * base[i] + (iCHI[i] - 1);
    }
    if (add > nRotamers) {
      cerr << "fatal error" << endl;
      exit(1);
    }
    return rotamers + add;
  }
};

class BBDep_RotLib {
private:
  BBDep_resRotLib ***records;
  void set_record(int phi, int psi, int id, vector<string> &pool);

public:
  BBDep_RotLib();
  ~BBDep_RotLib();

  void readFile(const std::string &rotFile);
  void updateStatisticalEnergy();
  void updateCAP();
  BBDep_resRotLib *get_resRotLib(const int &iF, const int &iY,
                                 const int &type) {
    return &records[iF][iY][type];
  }
  BBDep_rotRecord *newRotamer(BBDep_rotRecord *oldRotamer, const int &type,
                              const int &newF, const int &newY);
};

const static double cutoff_ratio = 0.93;
const static double RMIN_HBD_HBA = 2.95;
// for EEF1 solvation energy parameters, now in VDW class
const static int N_EEF1_t = 5;
typedef enum {
  EEF1_VOL = 0,
  EEF1_DG_REF,
  EEF1_DG_FREE,
  EEF1_LAMBDA,
  CHARMM_RADIUS
} EEF1_t;

class VDWpair {
private:
  double B;
  double A;
  double sigma;
  double epsilon;
  double rmin;

  double d_cutoff;
  double slope;
  double y_intercept;

  double B1_4;
  double A1_4;
  double sigma1_4;
  double epsilon1_4;
  double rmin1_4;

  double d_cutoff1_4;
  double slope1_4;
  double y_intercept1_4;

public:
  VDWpair(double iB = 0, double iA = 0, double iB1_4 = 0, double iA1_4 = 0) {
    B = iB;
    A = iA;
    B1_4 = iB1_4;
    A1_4 = iA1_4;

    if (A && B) {
      epsilon = (A * A) / (4.0 * B);
      sigma = pow(B / A, 1.0 / 6.0);
      rmin = sigma * pow(2.0, 1.0 / 6.0);

      d_cutoff = sigma * cutoff_ratio;
      slope = -24 * epsilon *
              (2 * pow(cutoff_ratio, -13.0) - pow(cutoff_ratio, -7.0)) / sigma;
      y_intercept =
          4.0 * epsilon * (pow(cutoff_ratio, -12.0) - pow(cutoff_ratio, -6.0)) -
          slope * d_cutoff;
    } else {
      sigma = INF;
      epsilon = 0;
      rmin = INF;

      d_cutoff = INF;
      slope = INF;
      y_intercept = INF;
    }

    if (A1_4 && B1_4) {
      epsilon1_4 = (A1_4 * A1_4) / (4.0 * B1_4);
      sigma1_4 = pow(B1_4 / A1_4, 1.0 / 6.0);
      rmin1_4 = sigma1_4 * pow(2.0, 1.0 / 6.0);

      d_cutoff1_4 = sigma1_4 * cutoff_ratio;
      slope1_4 = -24 * epsilon1_4 *
                 (2 * pow(cutoff_ratio, -13.0) - pow(cutoff_ratio, -7.0)) /
                 sigma1_4;
      y_intercept1_4 =
          4.0 * epsilon1_4 *
              (pow(cutoff_ratio, -12.0) - pow(cutoff_ratio, -6.0)) -
          slope1_4 * d_cutoff1_4;
    } else {
      sigma1_4 = INF;
      epsilon1_4 = 0;
      rmin1_4 = INF;

      d_cutoff1_4 = INF;
      slope1_4 = INF;
      y_intercept1_4 = INF;
    }
  }

  ~VDWpair() {}

  void setBA(double iB, double iA) {
    B = iB;
    A = iA;
    if (A && B) {
      epsilon = (A * A) / (4 * B);
      sigma = pow(B / A, 1.0 / 6.0);
      rmin = sigma * pow(2.0, 1.0 / 6.0);

      d_cutoff = sigma * cutoff_ratio;
      slope = -24 * epsilon *
              (2 * pow(cutoff_ratio, -13.0) - pow(cutoff_ratio, -7.0)) / sigma;
      y_intercept =
          4.0 * epsilon * (pow(cutoff_ratio, -12.0) - pow(cutoff_ratio, -6.0)) -
          slope * d_cutoff;
    } else {
      sigma = INF;
      epsilon = 0;
      rmin = INF;

      d_cutoff = INF;
      slope = INF;
      y_intercept = INF;
    }
  }

  void setBA1_4(double iB, double iA) {
    B1_4 = iB;
    A1_4 = iA;
    if (A1_4 && B1_4) {
      epsilon1_4 = (iA * iA) / (4 * iB);
      sigma1_4 = pow(iB / iA, 1.0 / 6.0);
      rmin1_4 = sigma1_4 * pow(2.0, 1.0 / 6.0);

      d_cutoff1_4 = sigma1_4 * cutoff_ratio;
      slope1_4 = -24 * epsilon1_4 *
                 (2 * pow(cutoff_ratio, -13.0) - pow(cutoff_ratio, -7.0)) /
                 sigma1_4;
      y_intercept1_4 =
          4.0 * epsilon1_4 *
              (pow(cutoff_ratio, -12.0) - pow(cutoff_ratio, -6.0)) -
          slope1_4 * d_cutoff1_4;
    } else {
      sigma1_4 = INF;
      epsilon1_4 = 0;
      rmin1_4 = INF;

      d_cutoff1_4 = INF;
      slope1_4 = INF;
      y_intercept1_4 = INF;
    }
  }

  void getBA(double &oB, double &oA, double &oE, double &oS) {
    oB = B;
    oA = A;
    oE = epsilon;
    oS = sigma;
  }

  void getBA(double &oB, double &oA, double &oE, double &oS, double &oRmin) {
    oB = B;
    oA = A;
    oE = epsilon;
    oS = sigma;
    oRmin = rmin;
  }

  void getREP_PARAM(double &o_cutoff, double &o_slope, double &o_intercept) {
    o_cutoff = d_cutoff;
    o_slope = slope;
    o_intercept = y_intercept;
  }

  void getREP_PARAM1_4(double &o_cutoff, double &o_slope, double &o_intercept) {
    o_cutoff = d_cutoff1_4;
    o_slope = slope1_4;
    o_intercept = y_intercept1_4;
  }

  void getBA1_4(double &oB, double &oA, double &oE, double &oS) {
    oB = B1_4;
    oA = A1_4;
    oE = epsilon1_4;
    oS = sigma1_4;
  }

  void getBA1_4(double &oB, double &oA, double &oE, double &oS, double &oRmin) {
    oB = B1_4;
    oA = A1_4;
    oE = epsilon1_4;
    oS = sigma1_4;
    oRmin = rmin1_4;
  }

  VDWpair &operator=(const VDWpair &right) {
    if (this == &right)
      return *this;
    B = right.B;
    A = right.A;
    B1_4 = right.B1_4;
    A1_4 = right.A1_4;

    d_cutoff = right.d_cutoff;
    slope = right.slope;
    y_intercept = right.y_intercept;

    d_cutoff1_4 = right.d_cutoff1_4;
    slope1_4 = right.slope1_4;
    y_intercept1_4 = right.y_intercept1_4;

    sigma = right.sigma;
    epsilon = right.epsilon;
    rmin = right.rmin;
    sigma1_4 = right.sigma1_4;
    epsilon1_4 = right.epsilon1_4;
    rmin1_4 = right.rmin1_4;
  }

  void shrink(const double ratio) {
    sigma *= ratio;
    rmin *= ratio;
    B *= pow(ratio, 12);
    A *= pow(ratio, 6);
    if (A && B) {
      d_cutoff = sigma * cutoff_ratio;
      slope = -24 * epsilon *
              (2 * pow(cutoff_ratio, -13.0) - pow(cutoff_ratio, -7.0)) / sigma;
      y_intercept =
          4.0 * epsilon * (pow(cutoff_ratio, -12.0) - pow(cutoff_ratio, -6.0)) -
          slope * d_cutoff;
    }

    sigma1_4 *= ratio;
    rmin1_4 *= ratio;
    B1_4 *= pow(ratio, 12);
    A1_4 *= pow(ratio, 6);
    if (A1_4 && B1_4) {
      d_cutoff1_4 = sigma1_4 * cutoff_ratio;
      slope1_4 = -24 * epsilon1_4 *
                 (2 * pow(cutoff_ratio, -13.0) - pow(cutoff_ratio, -7.0)) /
                 sigma1_4;
      y_intercept1_4 =
          4.0 * epsilon1_4 *
              (pow(cutoff_ratio, -12.0) - pow(cutoff_ratio, -6.0)) -
          slope1_4 * d_cutoff1_4;
    }
  }

  void shrink_non14(const double ratio) {
    sigma *= ratio;
    rmin *= ratio;
    B *= pow(ratio, 12);
    A *= pow(ratio, 6);
    if (A && B) {
      d_cutoff = sigma * cutoff_ratio;
      slope = -24 * epsilon *
              (2 * pow(cutoff_ratio, -13.0) - pow(cutoff_ratio, -7.0)) / sigma;
      y_intercept =
          4.0 * epsilon * (pow(cutoff_ratio, -12.0) - pow(cutoff_ratio, -6.0)) -
          slope * d_cutoff;
    }

    // sigma1_4 *= ratio;
    // rmin1_4 *= ratio;
    // B1_4 *= pow(ratio, 12);
    // A1_4 *= pow(ratio, 6);
    // if(A1_4&&B1_4){
    //   d_cutoff1_4 = sigma1_4*cutoff_ratio;
    //   slope1_4 =
    //   -24*epsilon1_4*(2*pow(cutoff_ratio,-13.0)-pow(cutoff_ratio,-7.0))/sigma1_4;
    //   y_intercept1_4 = 4.0*epsilon1_4*(pow(cutoff_ratio,
    //   -12.0)-pow(cutoff_ratio,-6.0))
    //	- slope1_4*d_cutoff1_4;
    // }
  }

  void shrink_14(const double ratio) {
    // sigma *= ratio;
    // rmin *= ratio;
    // B *= pow(ratio, 12);
    // A *= pow(ratio, 6);
    // if(A&&B){
    //   d_cutoff = sigma*cutoff_ratio;
    //   slope =
    //   -24*epsilon*(2*pow(cutoff_ratio,-13.0)-pow(cutoff_ratio,-7.0))/sigma;
    //   y_intercept = 4.0*epsilon*(pow(cutoff_ratio,
    //   -12.0)-pow(cutoff_ratio,-6.0))
    //	- slope*d_cutoff;
    // }

    sigma1_4 *= ratio;
    rmin1_4 *= ratio;
    B1_4 *= pow(ratio, 12);
    A1_4 *= pow(ratio, 6);
    if (A1_4 && B1_4) {
      d_cutoff1_4 = sigma1_4 * cutoff_ratio;
      slope1_4 = -24 * epsilon1_4 *
                 (2 * pow(cutoff_ratio, -13.0) - pow(cutoff_ratio, -7.0)) /
                 sigma1_4;
      y_intercept1_4 =
          4.0 * epsilon1_4 *
              (pow(cutoff_ratio, -12.0) - pow(cutoff_ratio, -6.0)) -
          slope1_4 * d_cutoff1_4;
    }
  }

  void copyTo_14() {
    sigma1_4 = sigma;
    rmin1_4 = rmin;
    B1_4 = B;
    A1_4 = A;
    d_cutoff1_4 = d_cutoff;
    slope1_4 = slope;
    y_intercept1_4 = y_intercept;
    epsilon1_4 = epsilon;
  }

  void scaleEPS_14(double ratio) {
    sigma1_4 *= ratio;
    B1_4 *= ratio;
    A1_4 *= ratio;
  }

  void setCHARMM_PARAM(double Emin_ij, double Rmin_ij, double Emin_14,
                       double Rmin_14) {
    if (A && B) { // do not account PROTON's vdw
      epsilon = Emin_ij;
      sigma = Rmin_ij / pow(2.0, 1.0 / 6.0);
      rmin = Rmin_ij;

      epsilon1_4 = Emin_14;
      sigma1_4 = Rmin_14 / pow(2.0, 1.0 / 6.0);
      rmin1_4 = Rmin_14;

      B = Emin_ij * pow(Rmin_ij, 12);
      A = Emin_ij * pow(Rmin_ij, 6) * 2;

      B1_4 = Emin_14 * pow(Rmin_14, 12);
      A1_4 = Emin_14 * pow(Rmin_14, 6) * 2;

      d_cutoff = sigma * cutoff_ratio;
      slope = -24 * epsilon *
              (2 * pow(cutoff_ratio, -13.0) - pow(cutoff_ratio, -7.0)) / sigma;
      y_intercept =
          4.0 * epsilon * (pow(cutoff_ratio, -12.0) - pow(cutoff_ratio, -6.0)) -
          slope * d_cutoff;

      d_cutoff1_4 = sigma1_4 * cutoff_ratio;
      slope1_4 = -24 * epsilon1_4 *
                 (2 * pow(cutoff_ratio, -13.0) - pow(cutoff_ratio, -7.0)) /
                 sigma1_4;
      y_intercept1_4 =
          4.0 * epsilon1_4 *
              (pow(cutoff_ratio, -12.0) - pow(cutoff_ratio, -6.0)) -
          slope1_4 * d_cutoff1_4;
    }
  }

  void resetRMin(const double &newRMIN) {
    // cout << epsilon << " " << sigma << endl;
    sigma = newRMIN / pow(2.0, 1.0 / 6.0);
    sigma1_4 = sigma;

    rmin = newRMIN;
    rmin1_4 = rmin;

    B = 4.0 * epsilon * pow(sigma, 12);
    A = 4.0 * epsilon * pow(sigma, 6);

    B1_4 = 4.0 * epsilon * pow(sigma1_4, 12);
    A1_4 = 4.0 * epsilon * pow(sigma1_4, 6);
    // cout << epsilon << " " << sigma << endl;
    // cout << endl;

    d_cutoff = sigma * cutoff_ratio;
    slope = -24 * epsilon *
            (2 * pow(cutoff_ratio, -13.0) - pow(cutoff_ratio, -7.0)) / sigma;
    y_intercept =
        4.0 * epsilon * (pow(cutoff_ratio, -12.0) - pow(cutoff_ratio, -6.0)) -
        slope * d_cutoff;

    d_cutoff1_4 = sigma1_4 * cutoff_ratio;
    slope1_4 = -24 * epsilon1_4 *
               (2 * pow(cutoff_ratio, -13.0) - pow(cutoff_ratio, -7.0)) /
               sigma1_4;
    y_intercept1_4 = 4.0 * epsilon1_4 *
                         (pow(cutoff_ratio, -12.0) - pow(cutoff_ratio, -6.0)) -
                     slope1_4 * d_cutoff1_4;
  }
};

const static int N_CEDU_t = 24;
typedef enum {
  CEDU_H = 0,
  CEDU_HC,
  CEDU_CR,   //-> C   (EEF1)
  CEDU_CRNP, //-> CR
  CEDU_CA,   //-> CH1E
  CEDU_CA1,  //-> CH1E
  CEDU_CA1P, //-> CH1E
  CEDU_CA2,  //-> CH2E
  CEDU_CA2P, //-> CH2E
  CEDU_CA3,  //-> CH3E
  CEDU_CM1,  //-> CR1E
  CEDU_CM1P, //-> CR1E
  CEDU_NCK,  //-> NH3
  CEDU_NCR,  //-> NC2
  CEDU_NZ,   //-> NH1
  CEDU_NR,   //-> NR
  CEDU_NZNP, //-> N
  CEDU_NZNQ, //-> NH2
  CEDU_OC,   //-> OC
  CEDU_OZ,   //-> O
  CEDU_OW,   //-> OT
  CEDU_OZH,  //-> OH1
  CEDU_SG,   //-> SH1E
  CEDU_SGNP, //-> S
} CEDU_t;

class VDW {
private:
  int nType;
  // VDWpair** pairs;
  // double** EEF1_SOLV;
  vector<vector<VDWpair>> pairs;
  vector<vector<double>> EEF1_SOLV;
  vector<vector<double>> E_S;

public:
  VDW() {
    nType = 0;
    // pairs = NULL;
    // EEF1_SOLV=NULL;
  }

  ~VDW() {
    // delete [] pairs[0];
    // delete [] pairs;
    // delete [] EEF1_SOLV[0];
    // delete [] EEF1_SOLV;
  }

  void init(const char *file, const topology &top);
  void addParam(const char *file, const topology &top);

  int getVDWdata(int type1, int type2, double &B, double &A, double &eps,
                 double &sigma) {
    if (type1 >= 0 && type1 < nType && type2 >= 0 && type2 < nType) {
      pairs[type1][type2].getBA(B, A, eps, sigma);
      return 1;
    }
    return 0;
  }

  int getVDWdata(int type1, int type2, double &B, double &A, double &eps,
                 double &sigma, double &rep_cutoff, double &rep_slope,
                 double &rep_intercept) {
    if (type1 >= 0 && type1 < nType && type2 >= 0 && type2 < nType) {
      pairs[type1][type2].getBA(B, A, eps, sigma);
      pairs[type1][type2].getREP_PARAM(rep_cutoff, rep_slope, rep_intercept);
      return 1;
    }
    return 0;
  }

  int getVDWdata(int type1, int type2, double &B, double &A, double &eps,
                 double &sigma, double &rmin, double &rep_cutoff,
                 double &rep_slope, double &rep_intercept) {
    if (type1 >= 0 && type1 < nType && type2 >= 0 && type2 < nType) {
      pairs[type1][type2].getBA(B, A, eps, sigma, rmin);
      pairs[type1][type2].getREP_PARAM(rep_cutoff, rep_slope, rep_intercept);
      return 1;
    }
    return 0;
  }

  int getVDWdata1_4(int type1, int type2, double &B, double &A, double &eps,
                    double &sigma) {
    if (type1 >= 0 && type1 < nType && type2 >= 0 && type2 < nType) {
      pairs[type1][type2].getBA1_4(B, A, eps, sigma);
      return 1;
    }
    return 0;
  }

  int getVDWdata1_4(int type1, int type2, double &B, double &A, double &eps,
                    double &sigma, double &rep_cutoff, double &rep_slope,
                    double &rep_intercept) {
    if (type1 >= 0 && type1 < nType && type2 >= 0 && type2 < nType) {
      pairs[type1][type2].getBA1_4(B, A, eps, sigma);
      pairs[type1][type2].getREP_PARAM1_4(rep_cutoff, rep_slope, rep_intercept);
      return 1;
    }
    return 0;
  }
  int getVDWdata1_4(int type1, int type2, double &B, double &A, double &eps,
                    double &sigma, double &rmin, double &rep_cutoff,
                    double &rep_slope, double &rep_intercept) {
    if (type1 >= 0 && type1 < nType && type2 >= 0 && type2 < nType) {
      pairs[type1][type2].getBA1_4(B, A, eps, sigma, rmin);
      pairs[type1][type2].getREP_PARAM1_4(rep_cutoff, rep_slope, rep_intercept);
      return 1;
    }
    return 0;
  }
  int getNType() { return nType; }

  void shrink(const double ratio) {
    for (int i = 0; i < nType; i++) {
      for (int j = 0; j < nType; j++) {
        pairs[i][j].shrink(ratio);
      }
    }
  }

  void shrink_non14(const double ratio) {
    for (int i = 0; i < nType; i++) {
      for (int j = 0; j < nType; j++) {
        pairs[i][j].shrink_non14(ratio);
      }
    }
  }

  void shrink_14(const double ratio) {
    for (int i = 0; i < nType; i++) {
      for (int j = 0; j < nType; j++) {
        pairs[i][j].shrink_14(ratio);
      }
    }
  }

  void copyTo_14() {
    for (int i = 0; i < nType; i++) {
      for (int j = 0; j < nType; j++) {
        pairs[i][j].copyTo_14();
      }
    }
  }

  void scaleEPS_14(double ratio) {
    for (int i = 0; i < nType; i++) {
      for (int j = 0; j < nType; j++) {
        pairs[i][j].scaleEPS_14(ratio);
      }
    }
  }

  /* this method is obsolete, will be replace by HBond_RFix(top) */
  void HBond_RFix() {
    // Fix the RMin between the HBond Donar and the HBond Acceptor
    // to favor the hydrogen bonds
    // Donar list:
    // NCK, NCR, NZ, NZNQ, OZH, OW
    // Acceptor list:
    // NR, OC, OZ, OW, OZH
    vector<int> donors;
    donors.push_back(CEDU_NCK);
    donors.push_back(CEDU_NCR);
    donors.push_back(CEDU_NZ);
    donors.push_back(CEDU_NZNQ);
    donors.push_back(CEDU_OZH);
    donors.push_back(CEDU_OW);

    vector<int> acceptors;
    acceptors.push_back(CEDU_NR);  // HIS, N
    acceptors.push_back(CEDU_OC);  // CHARGED Oxygen
    acceptors.push_back(CEDU_OZ);  // PEPTIDE Oxygen
    acceptors.push_back(CEDU_OW);  // WATER Oxygen
    acceptors.push_back(CEDU_OZH); // SER,THR,TYR Oxygen

    for (int i = 0; i < donors.size(); i++) {
      for (int j = 0; j < acceptors.size(); j++) {
        pairs[donors[i]][acceptors[j]].resetRMin(RMIN_HBD_HBA);
        pairs[acceptors[j]][donors[i]] = pairs[donors[i]][acceptors[j]];
      }
    }
  }
  void HBond_RFix(topology &top) {
    // Fix the RMin between the HBond Donar and the HBond Acceptor
    // to favor the hydrogen bonds
    // Donar list:
    // NCK, NCR, NZ, NZNQ, OZH, OW
    // Acceptor list:
    // NR, OC, OZ, OW, OZH
    int NDonar = 6;
    const char *Donar[] = {"NCK", "NCR", "NZ", "NZNQ", "OZH", "OW"};
    int NAcceptor = 5;
    const char *Acceptor[] = {"NR", "OC", "OZ", "OW", "OZH"};
    vector<int> donors;
    int fft;
    for (int i = 0; i < NDonar; i++) {
      fft = top.fftid(Donar[i]);
      if (fft >= 0)
        donors.push_back(fft);
      else {
        cout << "unkown donar type " << Donar[i]
             << ",VDW radius adjustment fails" << endl;
        exit(1);
      }
    }
    vector<int> acceptors;
    for (int i = 0; i < NAcceptor; i++) {
      fft = top.fftid(Acceptor[i]);
      if (fft >= 0)
        acceptors.push_back(fft);
      else {
        cout << "unkown acceptor type " << Acceptor[i]
             << ", VDW radius adjustment fails" << endl;
        exit(1);
      }
    }
    for (int i = 0; i < donors.size(); i++) {
      for (int j = 0; j < acceptors.size(); j++) {
        pairs[donors[i]][acceptors[j]].resetRMin(RMIN_HBD_HBA);
        pairs[acceptors[j]][donors[i]] = pairs[donors[i]][acceptors[j]];
      }
    }
  }
  void HBond_RFix_Complex(topology &top) {
    // Fix the RMin between the HBond Donar and the HBond Acceptor
    // to favor the hydrogen bonds
    // Donar list:
    // NCK, NCR, NZ, NZNQ, OZH, OW
    // Acceptor list:
    // NR, OC, OZ, OW, OZH
    int NDonar = 18;
    const char *Donar[] = {"NCK", "NCR", "NZ", "NZNQ", "OZH", "OW",
                           // small molecule type
                           "_OH1", "_ONH", "_OS", "_NC2", "_NH2", "_NH3",
                           "_NNH2", "_NH1", "_NCCH2", "_NCNH", "_NCCCH",
                           "_NXX-H"};
    int NAcceptor = 26;
    const char *Acceptor[] = {"NR", "OC", "OZ", "OW", "OZH",
                              // small molecule type
                              "_O", "_OH1", "_OPO", "_ONO", "_ON", "_ONH",
                              "_OS", "_OCC", "_OCP", "_OCN", "_NC", "_NC2",
                              "_NH2", "_NH1", "_NR", "_NCN", "_NCNH", "_NNN",
                              "_NXO", "_NXSO", "_NXPO"};
    vector<int> donors;
    int fft;
    for (int i = 0; i < NDonar; i++) {
      fft = top.fftid(Donar[i]);
      if (fft >= 0)
        donors.push_back(fft);
      else {
        cout << "unkown donar type " << Donar[i]
             << ",VDW radius adjustment fails" << endl;
        exit(1);
      }
    }
    vector<int> acceptors;
    for (int i = 0; i < NAcceptor; i++) {
      fft = top.fftid(Acceptor[i]);
      if (fft >= 0)
        acceptors.push_back(fft);
      else {
        cout << "unkown acceptor type " << Acceptor[i]
             << ", VDW radius adjustment fails" << endl;
        exit(1);
      }
    }
    for (int i = 0; i < donors.size(); i++) {
      for (int j = 0; j < acceptors.size(); j++) {
        pairs[donors[i]][acceptors[j]].resetRMin(RMIN_HBD_HBA);
        pairs[acceptors[j]][donors[i]] = pairs[donors[i]][acceptors[j]];
      }
    }
  }

  // double** getEEF1_SOLV(){return EEF1_SOLV;}
  const vector<vector<double>> &getEEF1_SOLV() const { return EEF1_SOLV; }
};

} // namespace medusa
