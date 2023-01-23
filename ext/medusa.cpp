// #include "jnc/core"
// #include "core/MEDUSAProtein.h"
// #include "core/MEDUSAGrid.h"
// #include "core/MEDUSARotLib.h"
// #include "core/MEDUSATop.h"
// #include "core/MEDUSANonBond.h"
// #include "core/MEDUSAFGrid.h"
// #include "core/MEDUSA_BBDEP_RotLib.h"
// #include "core/MEDUSA_BBDEP_AA.h"
// #include "core/MEDUSAGENProtein.h"
// #include "core/MEDUSAGENGrid.h"
// #include "core/MEDUSAGENFGrid.h"
// #include "core/MEDUSA_COEF.h"
// #include "util/MEDUSAUtility.h"
// #include "util/MEDUSAargument.h"

#include "medusa.h"
#include "jnc.h"

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <string>

namespace medusa {

void protein::init(const string &thename, const vector<string> &seq) {
  strcpy(name, thename.c_str());
  setChainName(' ');

  // the rotation matrix
  m = new double *[3];
  m[0] = new double[9];
  for (int i = 0; i < 9; i++)
    m[0][i] = 0;
  m[1] = m[0] + 3;
  m[2] = m[0] + 6;

  length = seq.size();
  res = new residue[length];
  phi = new double[length];
  psi = new double[length];
  res[0].init(seq[0].c_str());
  res[0].initTemplate();
  vec xbar, ybar, zbar;
  vec shift;
  for (int i = 1; i < length; i++) {
    res[i].init(seq[i].c_str());
    // cout << AAname3[seq[i]] << endl;
    res[i].initTemplate();
    shift = res[i - 1].nextN(xbar, ybar, zbar);
    // xbar.print(cout);
    res[i].rotate(xbar, ybar, zbar);
    shift = shift - res[i].getN()->r;
    res[i].shift(shift);
  }
  // add all the polar hydrogen
  addPH();
  // update the state variables:
  // including chi angles, phi, psis
  for (int i = 0; i < length; i++) {
    phi[i] = INF;
    psi[i] = INF;
    res[i].updateRotamers();
    if (i > 0)
      phi[i] = getDihedralAngle(res[i - 1].getC()->r, res[i].getN()->r,
                                res[i].getCA()->r, res[i].getC()->r);
    if (i < length - 1)
      psi[i] = getDihedralAngle(res[i].getN()->r, res[i].getCA()->r,
                                res[i].getC()->r, res[i + 1].getN()->r);
  }
}

void protein::init(const mol &theMol, const bool hasBBH) {
  strcpy(name, theMol.name().c_str());
  setChainName(theMol.ID());
  // cout << chainName << endl;

  // the rotation matrix
  m = new double *[3];
  m[0] = new double[9];
  for (int i = 0; i < 9; i++)
    m[0][i] = 0;
  m[1] = m[0] + 3;
  m[2] = m[0] + 6;

  length = theMol.nResidues();
  res = new residue[length];
  phi = new double[length];
  psi = new double[length];
  for (int i = 0; i < length; i++) {
    res[i].init(theMol.residue(i).name().c_str());
    res[i].setResid(theMol.residue(i).index());
    res[i].setChainName(theMol.ID());
    int nha = 0;
    for (int j = 0; j < theMol.residue(i).nAtoms(); j++) {
      if (theMol.residue(i).atoms()[j].name().c_str()[0] != 'H' &&
          theMol.residue(i).atoms()[j].name().c_str()[0] !=
              'h') { // heavy atoms
        int aid = atomid(theMol.residue(i).atoms()[j].name().c_str());
        if (aid != -1) {
          nha++;
          res[i].setXYZ(aid, theMol.residue(i).atoms()[j].pos());
        }
      }
    }
    if (nha != res[i].getNA()) {
      cout << theMol.residue(i).nAtoms() << " " << res[i].getNA() << endl;
      cout << " residue:" << i + 1 << theMol.residue(i).name()
           << " contains missing heavy atom" << endl;
    }
  }

  // add all the polar hydrogen
  addPH();
  // update the state variables:
  // including chi angles, phi, psis
  for (int i = 0; i < length; i++) {
    phi[i] = INF;
    psi[i] = INF;
    res[i].updateRotamers();
    if (i > 0)
      phi[i] = getDihedralAngle(res[i - 1].getC()->r, res[i].getN()->r,
                                res[i].getCA()->r, res[i].getC()->r);
    if (i < length - 1)
      psi[i] = getDihedralAngle(res[i].getN()->r, res[i].getCA()->r,
                                res[i].getC()->r, res[i + 1].getN()->r);
  }
  if (hasBBH) { // if input PDB contains backbone hydrogen information
    // cout << " assume the native H?" << endl;
    if (res[0].getType() != PRO) {
      res[0].getHN()->setR(theMol.residue(0).atom("HN").pos());
      // res[0].getHN()->getR()->print(cout); cout << endl;
    }
  }
}

void readSEQ(const char *file, vector<int> &seq) {
  ifstream in(file);
  char tmp[100];
  int id;
  string line;
  while (getline(in, line)) {
    if (!isComment(line.c_str())) {
      // while(in>>tmp){
      sscanf(line.c_str(), "%s", tmp);
      id = residueid(tmp);
      // cout << tmp << " " << id << endl;
      if (id != -1) {
        seq.push_back(id);
      } else {
        cout << "Sequence file:" << file << " has non-standard residue name"
             << endl;
        cout << " omit " << tmp << endl;
      }
    }
  }
  in.close();
}

protein::protein(const char *thename, char *file, generate_t t = _PDB_) {
  strcpy(name, thename);
  setChainName(' ');

  // the rotation matrix
  m = new double *[3];
  m[0] = new double[9];
  for (int i = 0; i < 9; i++)
    m[0][i] = 0;
  m[1] = m[0] + 3;
  m[2] = m[0] + 6;

  if (t == _PDB_) {
    DCLPDB *p = new DCLPDB(file);
    p->Read();
    length = p->N;
    res = new residue[length];
    phi = new double[length];
    psi = new double[length];
    for (int i = 0; i < length; i++) {
      res[i].init(p->res[i].name);
      for (int j = 0; j < p->res[i].na; j++) {
        int aid = atomid(p->res[i].a[j].name);
        if (aid != -1) {
          res[i].setXYZ(aid, p->res[i].a[j].r);
        }
      }
      if (p->res[i].na < res[i].getNA()) {
        cout << p->res[i].na << endl;
        cout << res[i].getNA() << endl;
        cout << "In PDB:" << file << " residue:" << i + 1 << p->res[i].name
             << " contains missing heavy atom" << endl;
      }
    }
    delete p;
  } else if (t == _PDB_BB_) {
    DCLPDB *p = new DCLPDB(file);
    p->Read();
    length = p->N;
    res = new residue[length];
    phi = new double[length];
    psi = new double[length];
    vec N, CA, C, O;
    vec xbar, ybar, zbar;
    for (int i = 0; i < length; i++) {
      res[i].init(p->res[i].name);
      for (int j = 0; j < p->res[i].na; j++) {
        int aid = atomid(p->res[i].a[j].name);
        if (aid == _N_)
          N = p->res[i].a[j].r;
        else if (aid == _CA_)
          CA = p->res[i].a[j].r;
        else if (aid == _C_)
          C = p->res[i].a[j].r;
      }
      xbar = (CA - N).norm();
      zbar = xbar ^ (C - CA);
      zbar.Normalize();
      ybar = zbar ^ xbar;
      res[i].initTemplate();
      res[i].rotate(xbar, ybar, zbar);
      // res[i].shift(N);
      res[i].shift(CA);
      if (i > 0) {
        double p =
            getDihedralAngle(res[i - 1].getAtom(0)->r, res[i - 1].getAtom(1)->r,
                             res[i - 1].getAtom(2)->r, res[i].getAtom(0)->r);
        // cout << psi << endl;
        res[i - 1].getAtom(3)->getR()->Rotate(res[i - 1].getAtom(1)->r,
                                              res[i - 1].getAtom(2)->r, p + PI);
      }
    }
    // overlap the backbone atoms with the original PDB
    for (int i = 0; i < length; i++) {
      for (int j = 0; j < p->res[i].na; j++) {
        int aid = atomid(p->res[i].a[j].name);
        if (aid == _N_)
          *res[i].getAtom(0)->getR() = p->res[i].a[j].r;
        else if (aid == _CA_)
          *res[i].getAtom(1)->getR() = p->res[i].a[j].r;
        else if (aid == _C_)
          *res[i].getAtom(2)->getR() = p->res[i].a[j].r;
        else if (aid == _O_)
          *res[i].getAtom(3)->getR() = p->res[i].a[j].r;
      }
    }
    delete p;
  } else if (t == _PDB_MIX_) {
    DCLPDB *p = new DCLPDB(file);
    p->Read();
    length = p->N;
    res = new residue[length];
    phi = new double[length];
    psi = new double[length];
    vec xbar, ybar, zbar;
    for (int i = 0; i < length; i++) {
      res[i].init(p->res[i].name);
      if (p->res[i].na == res[i].getNA()) {
        res[i].init(p->res[i].name);
        for (int j = 0; j < p->res[i].na; j++) {
          int aid = atomid(p->res[i].a[j].name);
          if (aid != -1) {
            res[i].setXYZ(aid, p->res[i].a[j].r);
          }
        }
      } else { // take n,ca,c,o only
        vec *N = NULL;
        vec *CA = NULL;
        vec *C = NULL;
        vec *O = NULL;
        for (int j = 0; j < p->res[i].na; j++) {
          int aid = atomid(p->res[i].a[j].name);
          switch (aid) {
          case _N_:
            N = new vec(p->res[i].a[j].r);
            N->print(cout);
            cout << endl;
            break;
          case _C_:
            C = new vec(p->res[i].a[j].r);
            break;
          case _CA_:
            CA = new vec(p->res[i].a[j].r);
            break;
          case _O_:
            O = new vec(p->res[i].a[j].r);
            break;
          default:
            break;
          }
        }
        if (N && CA && C) { // initialize
          xbar = (*CA - *N).norm();
          zbar = xbar ^ (*C - *CA);
          zbar.Normalize();
          ybar = zbar ^ xbar;
          res[i].initTemplate();
          res[i].rotate(xbar, ybar, zbar);
          // res[i].shift(N);
          res[i].shift(*CA);
          if (i > 0 && psi[i - 1] == -INF) {
            double p = getDihedralAngle(
                res[i - 1].getAtom(0)->r, res[i - 1].getAtom(1)->r,
                res[i - 1].getAtom(2)->r, res[i].getAtom(0)->r);
            // cout << psi << endl;
            res[i - 1].getAtom(3)->getR()->Rotate(
                res[i - 1].getAtom(1)->r, res[i - 1].getAtom(2)->r, p + PI);
          }
          *res[i].getAtom(0)->getR() = *N;
          *res[i].getAtom(1)->getR() = *CA;
          *res[i].getAtom(2)->getR() = *C;
          if (O) {
            *res[i].getAtom(3)->getR() = *O;
            if (i > 0) {
              psi[i] = getDihedralAngle(
                  res[i - 1].getAtom(0)->r, res[i - 1].getAtom(1)->r,
                  res[i - 1].getAtom(2)->r, res[i].getAtom(0)->r);
            } else {
              psi[i] = -INF;
            }
          } else {
            psi[i] = -INF;
          }
          delete N;
          delete CA;
          delete C;
          delete O;
          cout << p->res[i].na << endl;
          cout << res[i].getNA() << endl;
          cout << "In PDB:" << file << " residue:" << i + 1 << p->res[i].name
               << " contains missing/additional heavy atom" << endl;
        } else {
          cerr << "In PDB:" << file << " residue:" << i + 1 << p->res[i].name
               << " miss backbone essential atoms, N,CA,C" << endl;
          exit(1);
        }
      }
    }
    delete p;

  } else if (t == _SEQ_) {
    vector<int> seq;
    readSEQ(file, seq);
    length = seq.size();
    res = new residue[length];
    phi = new double[length];
    psi = new double[length];
    res[0].init(residuename(seq[0]));
    res[0].initTemplate();
    vec xbar, ybar, zbar;
    vec shift;
    for (int i = 1; i < length; i++) {
      res[i].init(residuename(seq[i]));
      // cout << AAname3[seq[i]] << endl;
      res[i].initTemplate();
      shift = res[i - 1].nextN(xbar, ybar, zbar);
      // xbar.print(cout);
      res[i].rotate(xbar, ybar, zbar);
      shift = shift - res[i].getN()->r;
      res[i].shift(shift);
    }
  }
  // add all the polar hydrogen
  addPH();
  // update the state variables:
  // including chi angles, phi, psis
  for (int i = 0; i < length; i++) {
    phi[i] = INF;
    psi[i] = INF;
    res[i].updateRotamers();
    if (i > 0)
      phi[i] = getDihedralAngle(res[i - 1].getC()->r, res[i].getN()->r,
                                res[i].getCA()->r, res[i].getC()->r);
    if (i < length - 1)
      psi[i] = getDihedralAngle(res[i].getN()->r, res[i].getCA()->r,
                                res[i].getC()->r, res[i + 1].getN()->r);
  }
}

protein::protein(const char *thename, char *file, vector<int> &missingSCs,
                 generate_t t = _PDB_) {
  strcpy(name, thename);

  // the rotation matrix
  m = new double *[3];
  m[0] = new double[9];
  for (int i = 0; i < 9; i++)
    m[0][i] = 0;
  m[1] = m[0] + 3;
  m[2] = m[0] + 6;

  if (t == _PDB_) {
    DCLPDB *p = new DCLPDB(file);
    p->Read();
    length = p->N;
    res = new residue[length];
    for (int i = 0; i < length; i++) {
      res[i].init(p->res[i].name);
      for (int j = 0; j < p->res[i].na; j++) {
        int aid = atomid(p->res[i].a[j].name);
        if (aid != -1) {
          res[i].setXYZ(aid, p->res[i].a[j].r);
        }
      }
      if (p->res[i].na < res[i].getNA()) {
        cout << p->res[i].na << endl;
        cout << res[i].getNA() << endl;
        cout << "In PDB:" << file << " residue:" << i + 1 << p->res[i].name
             << " contains missing heavy atom" << endl;
      }
    }
    delete p;
  } else if (t == _PDB_BB_) {
    DCLPDB *p = new DCLPDB(file);
    p->Read();
    length = p->N;
    res = new residue[length];
    vec N, CA, C, O;
    vec xbar, ybar, zbar;
    double psi = 0;
    for (int i = 0; i < length; i++) {
      res[i].init(p->res[i].name);
      for (int j = 0; j < p->res[i].na; j++) {
        int aid = atomid(p->res[i].a[j].name);
        if (aid == _N_)
          N = p->res[i].a[j].r;
        else if (aid == _CA_)
          CA = p->res[i].a[j].r;
        else if (aid == _C_)
          C = p->res[i].a[j].r;
      }
      xbar = (CA - N).norm();
      zbar = xbar ^ (C - CA);
      zbar.Normalize();
      ybar = zbar ^ xbar;
      res[i].initTemplate();
      res[i].rotate(xbar, ybar, zbar);
      // res[i].shift(N);
      res[i].shift(CA);
      if (i > 0) {
        psi =
            getDihedralAngle(res[i - 1].getAtom(0)->r, res[i - 1].getAtom(1)->r,
                             res[i - 1].getAtom(2)->r, res[i].getAtom(0)->r);
        // cout << psi << endl;
        res[i - 1].getAtom(3)->getR()->Rotate(
            res[i - 1].getAtom(1)->r, res[i - 1].getAtom(2)->r, psi + PI);
      }
    }
    delete p;
  } else if (t == _PDB_MIX_) {
    DCLPDB *p = new DCLPDB(file);
    p->Read();
    length = p->N;
    res = new residue[length];
    vec xbar, ybar, zbar;
    double psi = 0;
    for (int i = 0; i < length; i++) {
      res[i].init(p->res[i].name);
      if (p->res[i].na == res[i].getNA()) {
        res[i].init(p->res[i].name);
        for (int j = 0; j < p->res[i].na; j++) {
          int aid = atomid(p->res[i].a[j].name);
          if (aid != -1) {
            res[i].setXYZ(aid, p->res[i].a[j].r);
          }
        }
      } else { // take n,ca,c,o only
        vec *N = NULL;
        vec *CA = NULL;
        vec *C = NULL;
        vec *O = NULL;
        for (int j = 0; j < p->res[i].na; j++) {
          int aid = atomid(p->res[i].a[j].name);
          switch (aid) {
          case _N_:
            N = new vec(p->res[i].a[j].r);
            break;
          case _C_:
            C = new vec(p->res[i].a[j].r);
            break;
          case _CA_:
            CA = new vec(p->res[i].a[j].r);
            break;
          case _O_:
            O = new vec(p->res[i].a[j].r);
            break;
          default:
            break;
          }
        }
        if (N && O && CA && C) { // initialize
          xbar = (*CA - *N).norm();
          zbar = xbar ^ (*C - *CA);
          zbar.Normalize();
          ybar = zbar ^ xbar;
          res[i].initTemplate();
          res[i].rotate(xbar, ybar, zbar);
          // res[i].shift(N);
          res[i].shift(*CA);
          {
            res[i].getN()->r = *N;
            res[i].getC()->r = *C;
            res[i].getO()->r = *O;
          }

          // if(i>0){
          //  psi = getDihedralAngle(res[i-1].getAtom(0)->r,
          //  res[i-1].getAtom(1)->r,
          //			   res[i-1].getAtom(2)->r,
          // res[i].getAtom(0)->r);
          //  cout << psi << endl;
          //  res[i-1].getO()->getR()->Rotate(res[i-1].getCA()->r,
          //  res[i-1].getC()->r, psi);
          //}
        } else {
          cerr << "missing backbone heavy-atoms" << endl;
          exit(1);
        }
        delete N;
        delete CA;
        delete C;
        delete O;
        missingSCs.push_back(i);
        cout << p->res[i].na << endl;
        cout << res[i].getNA() << endl;
        cout << "In PDB:" << file << " residue:" << i + 1 << p->res[i].name
             << " contains missing/additional heavy atom" << endl;
      }
    }
    delete p;

  } else if (t == _SEQ_) {
    vector<int> seq;
    readSEQ(file, seq);
    length = seq.size();
    res = new residue[length];
    res[0].init(residuename(seq[0]));
    res[0].initTemplate();
    vec xbar, ybar, zbar;
    vec shift;
    for (int i = 1; i < length; i++) {
      res[i].init(residuename(seq[i]));
      // cout << AAname3[seq[i]] << endl;
      res[i].initTemplate();
      shift = res[i - 1].nextN(xbar, ybar, zbar);
      xbar.print(cout);
      res[i].rotate(xbar, ybar, zbar);
      shift = shift - res[i].getN()->r;
      res[i].shift(shift);
    }
  }
  // add all the polar hydrogen
  addPH();
  // update the state variables:
  // including chi angles, phi, psis
  phi = new double[length];
  psi = new double[length];
  for (int i = 0; i < length; i++) {
    phi[i] = INF;
    psi[i] = INF;
    res[i].updateRotamers();
    if (i > 0)
      phi[i] = getDihedralAngle(res[i - 1].getC()->r, res[i].getN()->r,
                                res[i].getCA()->r, res[i].getC()->r);
    if (i < length - 1)
      psi[i] = getDihedralAngle(res[i].getN()->r, res[i].getCA()->r,
                                res[i].getC()->r, res[i + 1].getN()->r);
  }
}

void protein::writePDB(ostream &out) {
  writePDB(out, 1, ' '); // default resid start from 1 and no chainID
  out << "END" << endl;
}

void protein::writePDB(ostream &out, int startRes, char chainID) {
  char line[81];
  // char name = chainName;
  char name = chainID;
  int iatom = 1;
  for (int i = 0; i < length; i++) {
    // print the heavy atoms first
    for (int j = 0; j < res[i].getNA(); j++) {
      sprintf(&line[0],
              "ATOM  %5d  %-4.4s%3.3s %1c%4d    %8.3lf%8.3lf%8.3lf  1.00  0.00 "
              "          %c  ",
              iatom++, res[i].getAtom(j)->getName(), res[i].getName(), name,
              i + startRes, res[i].getAtom(j)->getR()->getX(),
              res[i].getAtom(j)->getR()->getY(),
              res[i].getAtom(j)->getR()->getZ(),
              res[i].getAtom(j)->getName()[0]);
      out << line << endl;
    }
    // print the polar protons
    for (int j = 0; j < res[i].getNA(); j++) {
      for (int k = 0; k < res[i].getAtom(j)->n_attachedProtons; k++) {
        char atomname[10];
        int len;
        if (len =
                strlen(
                    getResidue(i)->getAtom(j)->attachedProtons[k].getName()) ==
                4) {
          atomname[0] =
              getResidue(i)->getAtom(j)->attachedProtons[k].getName()[3];
          for (int c = 0; c < 3; c++)
            atomname[c + 1] =
                getResidue(i)->getAtom(j)->attachedProtons[k].getName()[c];
          atomname[4] = '\0';
          sprintf(&line[0],
                  "ATOM  %5d %-5.5s%3.3s %1c%4d    %8.3lf%8.3lf%8.3lf  1.00  "
                  "0.00           %c  ",
                  iatom++, atomname, getResidue(i)->getName(), name,
                  i + startRes,
                  getResidue(i)->getAtom(j)->attachedProtons[k].getR()->getX(),
                  getResidue(i)->getAtom(j)->attachedProtons[k].getR()->getY(),
                  getResidue(i)->getAtom(j)->attachedProtons[k].getR()->getZ(),
                  getResidue(i)->getAtom(j)->attachedProtons[k].getName()[0]);
        } else {
          sprintf(&line[0],
                  "ATOM  %5d  %-4.4s%3.3s %1c%4d    %8.3lf%8.3lf%8.3lf  1.00  "
                  "0.00           %c  ",
                  iatom++,
                  getResidue(i)->getAtom(j)->attachedProtons[k].getName(),
                  getResidue(i)->getName(), name, i + startRes,
                  getResidue(i)->getAtom(j)->attachedProtons[k].getR()->getX(),
                  getResidue(i)->getAtom(j)->attachedProtons[k].getR()->getY(),
                  getResidue(i)->getAtom(j)->attachedProtons[k].getR()->getZ(),
                  getResidue(i)->getAtom(j)->attachedProtons[k].getName()[0]);
        }
        out << line << endl;
        //              iatom++,
        //              res[i].getAtom(j)->attachedProtons[k].getName(),
        //              res[i].getName(),
        //              i+1,
        //              res[i].getAtom(j)->attachedProtons[k].getR()->getX(),
        //              res[i].getAtom(j)->attachedProtons[k].getR()->getY(),
        //              res[i].getAtom(j)->attachedProtons[k].getR()->getZ(),
        //              res[i].getAtom(j)->attachedProtons[k].getName()[0]);
        //      out << line << endl;
      }
    }
  }
  out << "TER" << endl;
}
void protein::writeSTATE(ostream &out) {
  for (int i = 0; i < length; i++) {
    out << i + 1 << " " << residuename(res[i].getType()) << ":" << endl;
    ;
    if (i > 0)
      out << " phi: " << phi[i];
    else
      out << " phi: 0 ";
    if (i < length - 1)
      out << " psi: " << psi[i];
    else
      out << " psi: 0";
    double *chis = res[i].getChis();
    out << " Chis: " << res[i].getNChi() << " ";
    for (int j = 0; j < res[i].getNChi(); j++) {
      out << chis[j] << " ";
    }
    out << endl;
  }
}

int protein::pointMutation(int i, const char *newRes) {
  int newType = residueid(newRes);
  if (newType != -1) {
    if (newType == res[i].getType()) {
      return 0;
    } else {
      vec N, CA, C, O, H;
      vec xbar, ybar, zbar;
      N = res[i].getN()->r;
      CA = res[i].getCA()->r;
      C = res[i].getC()->r;
      O = res[i].getO()->r;
      // obtain the new internal axis
      xbar = (CA - N).norm();
      zbar = xbar ^ (C - CA);
      zbar.Normalize();
      ybar = zbar ^ xbar;

      // add the protns
      res[i].addPolarSH();
      if (i > 0)
        res[i].addBBH(res[i - 1]);
      else
        res[i].addBBH(true);
      //                res[i].addBBH();

      res[i].clear();                    // remove the old corrdinates
      res[i].init(residuename(newType)); // initialize
      res[i].initTemplate();             // also the cooridnates
      res[i].rotate(xbar, ybar, zbar);   // rotate to the new internal axis
      // res[i].shift(N);               //translate
      res[i].shift(CA);
      res[i].updateRotamers(); // update the rotamer angles
      res[i].getO()->r = O;

      return 1;
    }
  }
  return -1;
}

void protein::updateTopology(const topology &top) {
  for (int i = 0; i < length; i++) {
    res[i].updateTopology(top);
  }
}

void protein::updateTopology(const topology &top, const int iP,
                             const int startIndexOfChain) {
  updateTopology(top);
  for (int i = 0; i < length; i++) {
    residue *theRes = getResidue(i);
    for (int j = 0; j < theRes->getNA(); j++) {
      theRes->getAtom(j)->resID = i + startIndexOfChain;
      theRes->getAtom(j)->subResID = i;
      theRes->getAtom(j)->chainID = iP;
      for (int k = 0; k < theRes->getAtom(j)->n_attachedProtons; k++) {
        theRes->getAtom(j)->attachedProtons[k].resID = i + startIndexOfChain;
        theRes->getAtom(j)->attachedProtons[k].subResID = i;
        theRes->getAtom(j)->attachedProtons[k].chainID = iP;
      }
    }
  }
}

void protein::updateTopology(const topology &top, const int ires) {
  // construct the proton and acceptor arrays
  // used for hydrogen bond interaction
  res[ires].updateTopology(top);
}

int protein::rotateResidue(int index, const rotRecord &newRotamer) {
  const int *ichi = newRotamer.getIChi();
  const double *dchi = newRotamer.getDChi();
  residue *resi = getResidue(index);
  if (!resi) {
    cerr << "Can not find residue " << index + 1 << endl;
    return 0;
  }
  for (int i = 0; i < resi->getNChi(); i++) {
    if (ichi[i]) {
      if (!resi->rotateChi(i, dchi[i] * rpi, m)) {
        cerr << "Rotate resiue" << index + 1 << " rotamer" << i + 1
             << " unsuccessfully" << endl;
        return 0;
      }
    } else {
      cerr << "residue " << index + 1 << " have wrong rotamers" << endl;
      return 0;
    }
  }
}

void protein::getBoundary(vec &min, vec &max) {
  min = vec(INF, INF, INF);
  max = vec(-INF, -INF, -INF);
  for (int i = 0; i < length; i++) {
    for (int j = 0; j < res[i].getNA(); j++) {
      if (min.x > res[i].getAtom(j)->r.x)
        min.x = res[i].getAtom(j)->r.x;
      if (min.y > res[i].getAtom(j)->r.y)
        min.y = res[i].getAtom(j)->r.y;
      if (min.z > res[i].getAtom(j)->r.z)
        min.z = res[i].getAtom(j)->r.z;
      if (max.x < res[i].getAtom(j)->r.x)
        max.x = res[i].getAtom(j)->r.x;
      if (max.y < res[i].getAtom(j)->r.y)
        max.y = res[i].getAtom(j)->r.y;
      if (max.z < res[i].getAtom(j)->r.z)
        max.z = res[i].getAtom(j)->r.z;
      for (int k = 0; k < res[i].getAtom(j)->n_attachedProtons; k++) {
        vec &tmp = res[i].getAtom(j)->attachedProtons[k].r;
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
  }
}

void protein::addPH() {
  for (int i = 0; i < length; i++) {
    res[i].addPolarSH();
    if (i > 0) {
      //            if (res[i].getN()->r.getDist(res[i - 1].getC()->r) < 4) {
      res[i].addBBH(res[i - 1]);
      //            } else {  // treat as first residue if the protein chains is
      //            not connected
      // res[i].addBBH(true);
      //                res[i].addBBH();
      //            }
    }
  }
  res[0].addBBH(true);
  //    res[0].addBBH();
  for (int i = 0; i < length; i++) {
    res[i].setPH_ARRAY();
  }
}

void protein::updatePH() {
  for (int i = 0; i < length; i++) {
    res[i].updatePolarSH();
    if (i > 0) {
      res[i].updateBBH(res[i - 1]);
    }
  }
  res[0].updateBBH(true);
  //    res[0].updateBBH();
}

// each atom has a residue index
void protein::setResidueIndex() {
  for (int i = 0; i < length; i++) {
    for (int j = 0; j < res[i].getNA(); j++) {
      res[i].getAtom(j)->resID = i;
      for (int k = 0; k < res[i].getAtom(j)->n_attachedProtons; k++) {
        res[i].getAtom(j)->attachedProtons[k].resID = i;
      }
    }
  }
}

// each atom also has an atomic index
int protein::setAtomIndex() {
  int index = 0;
  for (int i = 0; i < length; i++) {
    for (int j = 0; j < res[i].getNA(); j++) {
      res[i].getAtom(j)->index = index++;
      for (int k = 0; k < res[i].getAtom(j)->n_attachedProtons; k++) {
        res[i].getAtom(j)->attachedProtons[k].index = index++;
      }
    }
  }
  return index;
}

int protein::setHAtomIndex() {
  int index = 1;
  for (int i = 0; i < length; i++) {
    for (int j = 0; j < res[i].getNA(); j++) {
      res[i].getAtom(j)->index = index++;
      for (int k = 0; k < res[i].getAtom(j)->n_attachedProtons; k++) {
        res[i].getAtom(j)->attachedProtons[k].index = index++;
      }
    }
  }
  return index - 1;
}

void protein::setAtomIndex(int &index) {
  for (int i = 0; i < length; i++) {
    for (int j = 0; j < res[i].getNA(); j++) {
      res[i].getAtom(j)->index = ++index;
      for (int k = 0; k < res[i].getAtom(j)->n_attachedProtons; k++) {
        res[i].getAtom(j)->attachedProtons[k].index = ++index;
      }
    }
  }
}

void protein::initConnect(const topology &top) {
  int natoms = 0;
  for (int i = 0; i < length; i++) {
    natoms += res[i].getNA(); // heavy atom
    natoms += res[i].getNH(); // polar hydrogen atom
  }
  con.init(natoms);

  TopResidue *theTopRes = NULL;
  residue *theRes = NULL;
  atom *theAtom = NULL;
  atom *theH = NULL;

  int prev_N = -1;
  int prev_CA = -1;
  int prev_C = -1;
  int prev_O = -1;

  int curr_N = -1;
  int curr_H = -1;
  int curr_CA = -1;
  int curr_C = -1;
  int curr_O = -1;

  for (int i = 0; i < length; i++) {
    theRes = res + i;
    theTopRes = top.getTopResi(theRes->getName());
    int numAtom = theRes->getNA() + theRes->getNH();

    if (numAtom != theTopRes->size()) {
      cerr << "warning: mismatch in the topology record and input protein "
              "(possibly due to protontation/deprotonation"
           << endl;
      // exit(1);
    }
    // create temperary mapping arrays
    int *top2atm = new int[numAtom];

    for (int j = 0; j < theRes->getNA(); j++) {
      theAtom = theRes->getAtom(j);
      top2atm[theTopRes->getAtomIndex(theAtom->getName())] = theAtom->index;
      for (int k = 0; k < theAtom->n_attachedProtons; k++) {
        theH = theAtom->attachedProtons + k;
        top2atm[theTopRes->getAtomIndex(theH->getName())] = theH->index;
      }
    }

    int **bonds = theTopRes->bond;
    /*set the intra-residue bonds*/
    for (int ibond = 0; ibond < theTopRes->nbond; ibond++) {
      con.setBond(top2atm[bonds[ibond][0]], top2atm[bonds[ibond][1]]);
    }
    /*set the intra-residue dihedrals*/
    int **dihedrals = theTopRes->dihedral;
    for (int idihe = 0; idihe < theTopRes->ndihedral; idihe++) {
      con.addDihedral(
          top2atm[dihedrals[idihe][0]], top2atm[dihedrals[idihe][1]],
          top2atm[dihedrals[idihe][2]], top2atm[dihedrals[idihe][3]]);
    }
    /*set the intra-residue impropers*/
    int **impropers = theTopRes->improper;
    for (int iimpr = 0; iimpr < theTopRes->nimproper; iimpr++) {
      con.addImproper(
          top2atm[impropers[iimpr][0]], top2atm[impropers[iimpr][1]],
          top2atm[impropers[iimpr][2]], top2atm[impropers[iimpr][3]]);
    }
    /*set inter peptide bonds, dihedral, impropers*/
    curr_N = res[i].getN()->index;
    curr_CA = res[i].getCA()->index;
    curr_C = res[i].getC()->index;
    curr_O = res[i].getO()->index;
    if (i > 0) {
      /*add C-N*/
      con.setBond(prev_C, curr_N);
      if (theRes->getType() != PRO) {
        curr_H = res[i].getHN()->index;
        /*dihedrals*/
        con.addDihedral(prev_C, curr_N, curr_CA, curr_C);
        con.addDihedral(prev_N, prev_CA, prev_C, curr_N);
        con.addDihedral(prev_CA, prev_C, curr_N, curr_CA);
        con.addDihedral(prev_CA, prev_C, curr_N, curr_H);
        con.addDihedral(prev_O, prev_C, curr_N, curr_CA);
        con.addDihedral(prev_O, prev_C, curr_N, curr_H);
        /*improper*/
        con.addImproper(curr_N, prev_C, prev_CA, prev_O);
        con.addImproper(curr_CA, curr_N, prev_C, curr_H);
      } else { /*PRO*/
        int curr_CG = curr_CA + 4;
        int curr_CD = curr_CA + 5;
        /*dihedral*/
        con.addDihedral(prev_C, curr_N, curr_CD, curr_CG);
        con.addDihedral(prev_C, curr_N, curr_CA, curr_C);
        con.addDihedral(prev_N, prev_CA, prev_C, curr_N);
        con.addDihedral(prev_CA, prev_C, curr_N, curr_CA);
        con.addDihedral(prev_CA, prev_C, curr_N, curr_CD);
        con.addDihedral(prev_O, prev_C, curr_N, curr_CA);
        con.addDihedral(prev_O, prev_C, curr_N, curr_CD);
        /*improper*/
        con.addImproper(curr_N, prev_C, prev_CA, prev_O);
        con.addImproper(curr_CA, curr_N, prev_C, curr_CD);
      }
    }
    prev_N = curr_N;
    prev_CA = curr_CA;
    prev_C = curr_C;
    prev_O = curr_O;
    delete[] top2atm;
  }
  con.postInitialize();
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
void protein::rotateTheta13(int ires, const double &dtheta) {
  if (ires > 0 && ires < (size() - 1)) {
    getDimondRotateMatrix(getResidue(ires - 1)->getCA()->r,
                          getResidue(ires + 1)->getCA()->r, dtheta, getM());
    // C,O of Residue ires-1
    atom *theAtom = getResidue(ires - 1)->getC();
    vec tmp = theAtom->r - getResidue(ires + 1)->getCA()->r;
    tmp.EularRotate(getM());
    theAtom->r = tmp + getResidue(ires + 1)->getCA()->r;

    theAtom = getResidue(ires - 1)->getO();
    tmp = theAtom->r - getResidue(ires + 1)->getCA()->r;
    tmp.EularRotate(getM());
    theAtom->r = tmp + getResidue(ires + 1)->getCA()->r;

    // N,H of Residue ires+1;
    theAtom = getResidue(ires + 1)->getN();
    tmp = theAtom->r - getResidue(ires + 1)->getCA()->r;
    tmp.EularRotate(getM());
    theAtom->r = tmp + getResidue(ires + 1)->getCA()->r;
    theAtom = theAtom->attachedProtons;
    if (getResidue(ires + 1)->getType() == PRO)
      theAtom = getResidue(ires + 1)->getAtom(6); // PROLINE, the CD Atom
    if (theAtom) { // Proline does not have backbone HN
      tmp = theAtom->r - getResidue(ires + 1)->getCA()->r;
      tmp.EularRotate(getM());
      theAtom->r = tmp + getResidue(ires + 1)->getCA()->r;
    }
    if (getResidue(ires + 1)->getType() == PRO)
      getResidue(ires + 1)->updateRotamers();

    // Residue ires
    for (int i = 0; i < getResidue(ires)->getNA(); i++) {
      theAtom = getResidue(ires)->getAtom(i);
      tmp = theAtom->r - getResidue(ires + 1)->getCA()->r;
      tmp.EularRotate(getM());
      theAtom->r = tmp + getResidue(ires + 1)->getCA()->r;
      for (int j = 0; j < getResidue(ires)->getAtom(i)->n_attachedProtons;
           j++) {
        theAtom = getResidue(ires)->getAtom(i)->attachedProtons + j;
        tmp = theAtom->r - getResidue(ires + 1)->getCA()->r;
        tmp.EularRotate(getM());
        theAtom->r = tmp + getResidue(ires + 1)->getCA()->r;
      }
    }
  }
}

/*Rotate around the axis of CAs of residues ires-1 and ires.
  dtheta is in RADIAN. */
void protein::rotateTheta12(int ires, const double &dtheta) {
  if (ires > 0) {
    getDimondRotateMatrix(getResidue(ires - 1)->getCA()->r,
                          getResidue(ires)->getCA()->r, dtheta, getM());
    // C,O of Residue ires-1
    atom *theAtom = getResidue(ires - 1)->getC();
    vec tmp = theAtom->r - getResidue(ires)->getCA()->r;
    tmp.EularRotate(getM());
    theAtom->r = tmp + getResidue(ires)->getCA()->r;

    theAtom = getResidue(ires - 1)->getO();
    tmp = theAtom->r - getResidue(ires)->getCA()->r;
    tmp.EularRotate(getM());
    theAtom->r = tmp + getResidue(ires)->getCA()->r;

    // N,H of Residue ires
    theAtom = getResidue(ires)->getN();
    tmp = theAtom->r - getResidue(ires)->getCA()->r;
    tmp.EularRotate(getM());
    theAtom->r = tmp + getResidue(ires)->getCA()->r;
    theAtom = theAtom->attachedProtons;
    if (getResidue(ires)->getType() == PRO)
      theAtom = getResidue(ires)->getAtom(6); // PROLINE, the CD Atom
    if (theAtom) {
      tmp = theAtom->r - getResidue(ires)->getCA()->r;
      tmp.EularRotate(getM());
      theAtom->r = tmp + getResidue(ires)->getCA()->r;
    }
    if (getResidue(ires)->getType() == PRO)
      getResidue(ires)->updateRotamers();
  }
}

/*Rotate around the axis of CAs of residues ires+1 and ires.
  dtheta is in RADIAN. */
void protein::rotateTheta32(int ires, const double &dtheta) {
  if (ires > 0) {
    getDimondRotateMatrix(getResidue(ires + 1)->getCA()->r,
                          getResidue(ires)->getCA()->r, dtheta, getM());
    // C,O of Residue ires
    atom *theAtom = getResidue(ires)->getC();
    vec tmp = theAtom->r - getResidue(ires)->getCA()->r;
    tmp.EularRotate(getM());
    theAtom->r = tmp + getResidue(ires)->getCA()->r;

    theAtom = getResidue(ires)->getO();
    tmp = theAtom->r - getResidue(ires)->getCA()->r;
    tmp.EularRotate(getM());
    theAtom->r = tmp + getResidue(ires)->getCA()->r;

    // N,H of Residue ires+1
    theAtom = getResidue(ires + 1)->getN();
    tmp = theAtom->r - getResidue(ires)->getCA()->r;
    tmp.EularRotate(getM());
    theAtom->r = tmp + getResidue(ires)->getCA()->r;
    theAtom = theAtom->attachedProtons;
    if (getResidue(ires + 1)->getType() == PRO)
      theAtom = getResidue(ires + 1)->getAtom(6); // PROLINE, the CD Atom
    if (theAtom) {
      tmp = theAtom->r - getResidue(ires)->getCA()->r;
      tmp.EularRotate(getM());
      theAtom->r = tmp + getResidue(ires)->getCA()->r;
    }
    if (getResidue(ires + 1)->getType() == PRO)
      getResidue(ires + 1)->updateRotamers();
  }
}

bool protein::updatePhiPsi() {
  for (int i = 0; i < length; i++) {
    if (i > 0) {
      phi[i] = getDihedralAngle(res[i - 1].getC()->r, res[i].getN()->r,
                                res[i].getCA()->r, res[i].getC()->r);
    }
    if (i < length - 1) {
      psi[i] = getDihedralAngle(res[i].getN()->r, res[i].getCA()->r,
                                res[i].getC()->r, res[i + 1].getN()->r);
    }
  }
}

void protein::getPhiPsi(int i, double &f, double &y) {
  if (i >= 0 && i < length) {
    f = phi[i];
    y = psi[i];
  }
}

gprotein::gprotein(protein &p) : initProtein(p) {
  length = p.size();
  residue *theResidue;
  gres = new gresidue[length];
  phi = new double[length];
  psi = new double[length];
  omiga = new double[length];
  list = new aa_rotamer *[length];
  naa_rot = new int[length];
  fy_uv = new vec[length * 2];
  omiga_uv = new vec[length];
  //  rotamers = new BBDep_rotRecord*[length];
  //  native_rotamers = new BBDep_rotRecord*[length];
  for (int i = 0; i < length; i++) {
    theResidue = p.getResidue(i);
    gres[i].init(*theResidue);
    gres[i].setResid(theResidue->getResid());
    gres[i].setResins(theResidue->getResins());
    phi[i] = INF;
    psi[i] = INF;
    omiga[i] = INF;
    list[i] = NULL;
    naa_rot[i] = 0;
    fy_uv[2 * i] = fy_uv[2 * i + 1] = vec(0, 0, 0);
    omiga_uv[i] = vec(0, 0, 0);
  }
  // add all the polar hydrogens
  for (int i = 0; i < length; i++) {
    gres[i].addPolarSH_ALL();
    if (i > 0) {
      if (gres[i].getResidue()->getN()->r.getDist(
              gres[i - 1].getResidue()->getC()->r) < 4) {
        gres[i].addBBH_ALL(*gres[i - 1].getResidue());
      } else { // treat as first residue if the protein chains is not connected
        cout << i << " r="
             << gres[i].getResidue()->getN()->r.getDist(
                    gres[i - 1].getResidue()->getC()->r)
             << "broken chain found\n";
        gres[i].addBBH_ALL();
      }
    };
  }
  gres[0].addBBH_ALL();
  // update the hydrogen atoms of the first residue as the input
  if (p.getResidue(0)->getHN()) {
    gres[0].updateBBH_ALL(p.getResidue(0)->getHN()->r);
  }

  // update the rotamers, backhone dihedrals
  for (int i = 0; i < length; i++) {
    gres[i].updateRotamers_ALL();
    if (i > 0 && gres[i].getResidue()->getN()->r.getDist(
                     gres[i - 1].getResidue()->getC()->r) < 4) {
      phi[i] = getDihedralAngle(
          gres[i - 1].getResidue()->getC()->r, gres[i].getResidue()->getN()->r,
          gres[i].getResidue()->getCA()->r, gres[i].getResidue()->getC()->r);
    }
    if (i < length - 1 && gres[i].getResidue()->getC()->r.getDist(
                              gres[i + 1].getResidue()->getN()->r) < 4) {
      psi[i] = getDihedralAngle(
          gres[i].getResidue()->getN()->r, gres[i].getResidue()->getCA()->r,
          gres[i].getResidue()->getC()->r, gres[i + 1].getResidue()->getN()->r);
      omiga[i] = getDihedralAngle(gres[i].getResidue()->getCA()->r,
                                  gres[i].getResidue()->getC()->r,
                                  gres[i + 1].getResidue()->getN()->r,
                                  gres[i + 1].getResidue()->getCA()->r);
    }
  }
}

gprotein::~gprotein() {
  delete[] gres;
  delete[] phi;
  delete[] psi;
  delete[] omiga;
  delete[] fy_uv;
  delete[] omiga_uv;
  for (int i = 0; i < length; i++)
    if (list[i])
      delete[] list[i];
  delete[] list;
  delete[] naa_rot;
  //  delete [] rotamers;
  //  delete [] native_rotamers;
}

void gprotein::updateTopology(const topology &top) {
  for (int i = 0; i < length; i++) {
    gres[i].updateTopology_ALL(top);
    gres[i].setResidueIndex_ALL(i);
  }
  res_con.init(top);
}

void gprotein::updateTopology(const topology &top, int chainID, int shift) {
  for (int i = 0; i < length; i++) {
    gres[i].updateTopology_ALL(top);
    gres[i].setResidueIndex_ALL(i + shift);
    gres[i].setSubResidueIndex_ALL(i);
    gres[i].setChainIndex_ALL(chainID);
  }
  res_con.init(top);
}

int gprotein::build_AA_Rotamer_List(BBDep_RotLib &rotlib, BBDep_AA &aalib,
                                    const double p_AA, const double p_ROT) {
  // use p=1% as the cutoff for amino acid
  double cutoff_E_AA = -log(p_AA);
  double cutoff_E_ROT = -log(p_ROT);
  int iF, iY;
  vector<aa_rotamer> pool;
  int total_rotamers = 0;
  for (int i = 0; i < length; i++) {
    FY_double2int(phi[i], psi[i], iF, iY);
    // if(phi[i]==INF) iF = FY_NBIN;
    // else iF=static_cast<int>((phi[i]/rpi-FY_MIN)/FY_BIN);
    // if(psi[i]==INF) iY = FY_NBIN;
    // else iY=static_cast<int>((psi[i]/rpi-FY_MIN)/FY_BIN);
    if ((iF == FY_NBIN) && (iY == FY_NBIN)) {
      cerr << "fatal error: phi/psi are not realistic at position" << i + 1
           << endl;
      exit(1);
    }

    /*build the list*/
    for (int ia = 0; ia < 20; ia++) {
      if (aalib.getE_STAT(iF, iY, ia) < cutoff_E_AA) {
        if ((ia != GLY) && (ia != ALA)) {
          int nRotamers = rotlib.get_resRotLib(iF, iY, ia)->getNRotamers();
          for (int irot = 0; irot < nRotamers; irot++) {
            BBDep_rotRecord *theRotamer =
                rotlib.get_resRotLib(iF, iY, ia)->getRotamer(irot);
            if (theRotamer->get_ESTAT() < cutoff_E_ROT) {
              pool.push_back(aa_rotamer(ia, theRotamer));
            }
          }
        } else {
          pool.push_back(aa_rotamer(ia, NULL));
        }
      }
    }
    // cout << i+1 << " " << phi[i]/rpi << " " << psi[i]/rpi << " " <<
    // pool.size() << endl;
    naa_rot[i] = pool.size();
    if (naa_rot[i] == 0) {
      cerr << "Fatal error: Position " << i + 1 << " has not available aa-rot"
           << endl;
      exit(1);
    }
    list[i] = new aa_rotamer[naa_rot[i]];
    total_rotamers += naa_rot[i];
    for (int ii = 0; ii < pool.size(); ii++) {
      list[i][ii].init(pool[ii]);
    }
    pool.clear();
  }
  return total_rotamers;
}

void FY_double2int(const double &phi, const double &psi, int &iF, int &iY) {
  if (phi == INF) {
    iF = FY_NBIN;
  } else {
    iF = static_cast<int>((phi / rpi - FY_MIN) / FY_BIN);
    // cout << iF << " " << phi/rpi-FY_MIN << " ";
    double dF = phi / rpi - iF * FY_BIN - FY_MIN;
    if (dF > FY_BIN / 2.0)
      iF++;
    if (iF >= FY_NBIN)
      iF -= FY_NBIN;
    // cout << iF << endl;
  }
  if (psi == INF) {
    iY = FY_NBIN;
  } else {
    iY = static_cast<int>((psi / rpi - FY_MIN) / FY_BIN);
    // cout << iY << " " << psi/rpi-FY_MIN << " ";
    double dY = psi / rpi - iY * FY_BIN - FY_MIN;
    if (dY > FY_BIN / 2.0)
      iY++;
    if (iY >= FY_NBIN)
      iY -= FY_NBIN;
    // cout << iY << endl;
  }
}

int gprotein::build_AA_Rotamer_List(BBDep_RotLib &rotlib, BBDep_AA &aalib,
                                    DesignTable &dt,
                                    BBDep_rotRecord **nativeRotamers,
                                    const double p_AA, const double p_ROT) {
  // use p=1% as the cutoff for amino acid
  double cutoff_E_AA = -log(p_AA);
  double cutoff_E_ROT = -log(p_ROT);
  int iF, iY;
  vector<aa_rotamer> pool;
  int total_rotamers = 0;
  for (int i = 0; i < length; i++) {
    FY_double2int(phi[i], psi[i], iF, iY);
    // if(phi[i]==INF) iF = FY_NBIN;
    // else iF=static_cast<int>((phi[i]/rpi-FY_MIN)/FY_BIN);
    // if(psi[i]==INF) iY = FY_NBIN;
    // else iY=static_cast<int>((psi[i]/rpi-FY_MIN)/FY_BIN);
    if ((iF == FY_NBIN) && (iY == FY_NBIN)) {
      cerr << "fatal error: phi/psi are not realistic at position" << i + 1
           << endl;
      exit(1);
    }

    TableEntry &tmp_entry = dt.getDT_Entry(i);
    int theAA;
    switch (dt.getDT_type(i)) {
    case ALLAA:
    case POLAR:
    case HYDPH:
    case AROMA:
    case PIKAA:
      for (int ia = 0; ia < tmp_entry.getNAvailableAA(); ia++) {
        theAA = tmp_entry.getAvailableAA(ia);
        if (aalib.getE_STAT(iF, iY, theAA) < cutoff_E_AA) {
          if ((theAA != GLY) && (theAA != ALA)) {
            int nRotamers = rotlib.get_resRotLib(iF, iY, theAA)->getNRotamers();
            for (int irot = 0; irot < nRotamers; irot++) {
              BBDep_rotRecord *theRotamer =
                  rotlib.get_resRotLib(iF, iY, theAA)->getRotamer(irot);
              if (theRotamer->get_ESTAT() < cutoff_E_ROT) {
                pool.push_back(aa_rotamer(theAA, theRotamer));
              }
            }
          } else {
            pool.push_back(aa_rotamer(theAA, NULL));
          }
        }
      }
      // native rot will be fored if no available aa,rot
      if (pool.size() <= 0) {
        cerr << "Position:" << i + 1 << " has no available rots" << endl;
        cerr << "use the native aa,rot" << endl;
        pool.push_back(aa_rotamer(getResidue(i)->getType(), nativeRotamers[i]));
      }
      break;
    case NATAA:
      theAA = getResidue(i)->getType();
      if (aalib.getE_STAT(iF, iY, theAA) < cutoff_E_AA) {
        if ((theAA != GLY) && (theAA != ALA)) {
          int nRotamers = rotlib.get_resRotLib(iF, iY, theAA)->getNRotamers();
          for (int irot = 0; irot < nRotamers; irot++) {
            BBDep_rotRecord *theRotamer =
                rotlib.get_resRotLib(iF, iY, theAA)->getRotamer(irot);
            if (theRotamer->get_ESTAT() < cutoff_E_ROT) {
              pool.push_back(aa_rotamer(theAA, theRotamer));
            }
          }
        } else {
          pool.push_back(aa_rotamer(theAA, NULL));
        }
      }
      // native rot will be fored is no available aa,rot
      if (pool.size() <= 0) {
        cerr << "Position:" << i + 1 << " has no available rots" << endl;
        cerr << "use the native rot" << endl;
        if (theAA != GLY && theAA != ALA) {
          if (nativeRotamers[i]->get_ESTAT() == INF)
            nativeRotamers[i]->setMinimumP();
          pool.push_back(
              aa_rotamer(getResidue(i)->getType(), nativeRotamers[i]));
        } else {
          pool.push_back(aa_rotamer(theAA, NULL));
        }
      } else { // Should include native rotamer into the list
        bool found = false;
        for (int ia = 0; ia < pool.size(); ia++) {
          if (pool[ia].getRotRecord() == nativeRotamers[i]) {
            found = true;
            break;
          }
        }
        if (!found) {
          cerr << "include native rotamer for position: " << i + 1 << endl;
          if (nativeRotamers[i]->get_ESTAT() == INF)
            nativeRotamers[i]->setMinimumP();
          pool.push_back(
              aa_rotamer(getResidue(i)->getType(), nativeRotamers[i]));
        }
      }

      break;
    case NAROT:
      pool.push_back(aa_rotamer(getResidue(i)->getType(), nativeRotamers[i]));
      break;
    case FIXNR:
      break;
    }

    naa_rot[i] = pool.size();
    if (naa_rot[i] > 0) {
      list[i] = new aa_rotamer[naa_rot[i]];
      total_rotamers += naa_rot[i];
      for (int ii = 0; ii < pool.size(); ii++) {
        list[i][ii].init(pool[ii]);
      }
    }
    pool.clear();
  }

  return total_rotamers;
}

int gprotein::build_AA_Rotamer_List(BBDep_RotLib &rotlib, BBDep_AA &aalib,
                                    DesignTable &dt, const double p_AA,
                                    const double p_ROT) {
  // use p=1% as the cutoff for amino acid
  double cutoff_E_AA = -log(p_AA);
  double cutoff_E_ROT = -log(p_ROT);
  int iF, iY;
  vector<aa_rotamer> pool;
  int total_rotamers = 0;
  for (int i = 0; i < length; i++) {
    FY_double2int(phi[i], psi[i], iF, iY);
    // if(phi[i]==INF) iF = FY_NBIN;
    // else iF=static_cast<int>((phi[i]/rpi-FY_MIN)/FY_BIN);
    // if(psi[i]==INF) iY = FY_NBIN;
    // else iY=static_cast<int>((psi[i]/rpi-FY_MIN)/FY_BIN);
    if ((iF == FY_NBIN) && (iY == FY_NBIN)) {
      cerr << "fatal error: phi/psi are not realistic at position" << i + 1
           << endl;
      exit(1);
    }

    TableEntry &tmp_entry = dt.getDT_Entry(i);
    int theAA;
    switch (dt.getDT_type(i)) {
    case ALLAA:
    case POLAR:
    case HYDPH:
    case AROMA:
    case PIKAA:
      for (int ia = 0; ia < tmp_entry.getNAvailableAA(); ia++) {
        theAA = tmp_entry.getAvailableAA(ia);
        if (aalib.getE_STAT(iF, iY, theAA) < cutoff_E_AA) {
          if ((theAA != GLY) && (theAA != ALA)) {
            int nRotamers = rotlib.get_resRotLib(iF, iY, theAA)->getNRotamers();
            for (int irot = 0; irot < nRotamers; irot++) {
              BBDep_rotRecord *theRotamer =
                  rotlib.get_resRotLib(iF, iY, theAA)->getRotamer(irot);
              if (theRotamer->get_ESTAT() < cutoff_E_ROT) {
                pool.push_back(aa_rotamer(theAA, theRotamer));
              }
            }
          } else {
            pool.push_back(aa_rotamer(theAA, NULL));
          }
        }
      }
      // native rot will be fored if no available aa,rot
      if (pool.size() <= 0) {
        cerr << "Position:" << i + 1 << " has no available rots" << endl;
        cerr << "Use the native aa" << endl;
        theAA = getResidue(i)->getType();
        if (aalib.getE_STAT(iF, iY, theAA) < cutoff_E_AA) {
          if ((theAA != GLY) && (theAA != ALA)) {
            int nRotamers = rotlib.get_resRotLib(iF, iY, theAA)->getNRotamers();
            for (int irot = 0; irot < nRotamers; irot++) {
              BBDep_rotRecord *theRotamer =
                  rotlib.get_resRotLib(iF, iY, theAA)->getRotamer(irot);
              if (theRotamer->get_ESTAT() < cutoff_E_ROT) {
                pool.push_back(aa_rotamer(theAA, theRotamer));
              }
            }
          } else {
            pool.push_back(aa_rotamer(theAA, NULL));
          }
        }
        // cerr << "use the native aa,rot" << endl;
        // pool.push_back(aa_rotamer(getResidue(i)->getType(),nativeRotamers[i]));
        if (pool.size() <= 0) {
          cerr << "fatal error" << endl;
          exit(1);
        }
      }
      break;
    case NATAA:
      theAA = getResidue(i)->getType();
      if (aalib.getE_STAT(iF, iY, theAA) < cutoff_E_AA) {
        if ((theAA != GLY) && (theAA != ALA)) {
          int nRotamers = rotlib.get_resRotLib(iF, iY, theAA)->getNRotamers();
          for (int irot = 0; irot < nRotamers; irot++) {
            BBDep_rotRecord *theRotamer =
                rotlib.get_resRotLib(iF, iY, theAA)->getRotamer(irot);
            if (theRotamer->get_ESTAT() < cutoff_E_ROT) {
              pool.push_back(aa_rotamer(theAA, theRotamer));
            }
          }
        } else {
          pool.push_back(aa_rotamer(theAA, NULL));
        }
      }
      // native rot will be fored is no available aa,rot
      if (pool.size() <= 0) {
        cerr << "Position:" << i + 1 << " has no available rots" << endl;
        // cerr << "use the native rot" << endl;
        // pool.push_back(aa_rotamer(getResidue(i)->getType(),nativeRotamers[i]));
        exit(1);
      }
      break;
    case NAROT:
      cerr << "Position:" << i + 1 << " has no NATIVE aa,rot information"
           << endl;
      exit(1);
      // pool.push_back(aa_rotamer(getResidue(i)->getType(),nativeRotamers[i]));
      break;
    case FIXNR:
      break;
    }

    naa_rot[i] = pool.size();
    if (naa_rot[i] > 0) {
      list[i] = new aa_rotamer[naa_rot[i]];
      total_rotamers += naa_rot[i];
      for (int ii = 0; ii < pool.size(); ii++) {
        list[i][ii].init(pool[ii]);
      }
    }
    pool.clear();
  }

  return total_rotamers;
}

void gprotein::updateList(BBDep_RotLib &rotlib, const int &iaa, const int &newF,
                          const int &newY) {
  for (int i = 0; i < naa_rot[iaa]; i++) {
    list[iaa][i].updateRotamer(rotlib, newF, newY);
  }
}

void gprotein::updateList(BBDep_RotLib &rotlib, const int &iaa, const int &newF,
                          const int &newY, BBDep_rotRecord *nativeRot) {
  double ERot_min = INF;
  int inative = -1;
  for (int i = 0; i < naa_rot[iaa]; i++) {
    if (nativeRot) {
      if ((nativeRot->getIChi()[0] ==
           list[iaa][i].getRotRecord()->getIChi()[0]) &&
          (nativeRot->getIChi()[1] ==
           list[iaa][i].getRotRecord()->getIChi()[1]) &&
          (nativeRot->getIChi()[2] ==
           list[iaa][i].getRotRecord()->getIChi()[2]) &&
          (nativeRot->getIChi()[3] ==
           list[iaa][i].getRotRecord()->getIChi()[3])) {
        inative = i;
      }
    }

    list[iaa][i].updateRotamer(rotlib, newF, newY);
    if (list[iaa][i].getRotRecord() &&
        list[iaa][i].getRotRecord()->get_ESTAT() < ERot_min) {
      ERot_min = list[iaa][i].getRotRecord()->get_ESTAT();
    }
  }
  if (nativeRot && naa_rot[iaa] > 0) {
    if (inative < 0) {
      cerr << "fatal error: native-like rotamer lost " << iaa << " "
           << nativeRot->getIChi()[0] << " " << nativeRot->getIChi()[1] << " "
           << nativeRot->getIChi()[2] << " " << nativeRot->getIChi()[3] << endl;
      exit(1);
    }
    if (ERot_min < list[iaa][inative].getRotRecord()->get_ESTAT()) {
      list[iaa][inative].getRotRecord()->setE_STAT(ERot_min);
    }
  }
}

int gprotein::rebuild_AA_Rotamer_List() {
  int ntotal = 0;
  for (int i = 0; i < length; i++) {
    vector<aa_rotamer> pool;
    for (int j = 0; j < naa_rot[i]; j++) {
      if (list[i][j].isAvailable()) {
        pool.push_back(list[i][j]);
      }
    }
    naa_rot[i] = pool.size();
    ntotal += naa_rot[i];
    delete[] list[i];
    list[i] = new aa_rotamer[naa_rot[i]];
    for (int j = 0; j < pool.size(); j++) {
      list[i][j].init(pool[j]);
    }
    pool.clear();
  }
  return ntotal;
}

int gprotein::rebuild_AA_Rotamer_List(const double &vdwr_cutoff) {
  int ntotal = 0;
  for (int i = 0; i < length; i++) {
    if (naa_rot[i] > 0) {
      vector<aa_rotamer> pool;
      for (int j = 0; j < naa_rot[i]; j++) {
        if (list[i][j].isAvailable()) {
          if (list[i][j].getMIN_BB_VDWR() < vdwr_cutoff ||
              list[i][j].getType() == PHE || list[i][j].getType() == TYR ||
              list[i][j].getType() == TRP) {
            pool.push_back(list[i][j]);
          }
        }
      }
      if (pool.size() > 0) {
        naa_rot[i] = pool.size();
        ntotal += naa_rot[i];
        delete[] list[i];
        list[i] = new aa_rotamer[naa_rot[i]];
        for (int j = 0; j < pool.size(); j++) {
          list[i][j].init(pool[j]);
        }
      } else {
        cerr << "error in rebuild AA:" << i + 1 << endl;
        exit(1);
      }
      pool.clear();
    }
  }
  return ntotal;
}

void gprotein::clear_Rotamer_List() {
  for (int i = 0; i < length; i++) {
    if (list[i])
      delete[] list[i];
  }
}

void gprotein::build_AA_Rotamer_List(BBDep_RotLib &rotlib,
                                     vector<int> &missing_idx) {
  int ia, iF, iY;
  vector<aa_rotamer> pool;
  for (int i = 0; i < length; i++) {
    FY_double2int(phi[i], psi[i], iF, iY);
    // if(phi[i]==INF) iF = FY_NBIN;
    // else iF=static_cast<int>((phi[i]/rpi-FY_MIN)/FY_BIN);
    // if(psi[i]==INF) iY = FY_NBIN;
    // else iY=static_cast<int>((psi[i]/rpi-FY_MIN)/FY_BIN);
    if ((iF == FY_NBIN) && (iY == FY_NBIN)) {
      cerr << "fatal error: phi/psi are not realistic at position" << i + 1
           << endl;
      exit(1);
    }

    int tobe_designed = 0;
    ia = getResidue(i)->getType();
    for (int j = 0; j < missing_idx.size(); j++) {
      if (i == missing_idx[j]) {
        tobe_designed = 1;
        // cout << i << " " << ia << endl;
        break;
      }
    }
    if (tobe_designed) {
      if (ia != GLY && ia != ALA) {
        int nRotamers = rotlib.get_resRotLib(iF, iY, ia)->getNRotamers();
        for (int irot = 0; irot < nRotamers; irot++) {
          BBDep_rotRecord *theRotamer =
              rotlib.get_resRotLib(iF, iY, ia)->getRotamer(irot);
          if (theRotamer->get_ESTAT() < INF) {
            pool.push_back(aa_rotamer(ia, theRotamer));
          }
        }
        if (pool.size() > 0) {
          // cout << pool.size() << endl;
          naa_rot[i] = pool.size();
          list[i] = new aa_rotamer[naa_rot[i]];
          for (int ii = 0; ii < pool.size(); ii++) {
            list[i][ii].init(pool[ii]);
          }
        } else {
          cout << "no available rotamer" << endl;
        }
      } else {
        naa_rot[i] = 1;
        list[i] = new aa_rotamer[naa_rot[i]];
        list[i][0].init(ia, NULL);
      }
    } else {
      naa_rot[i] = 0;
    }
    pool.clear();
  }
}

void gprotein::writePDB(ostream &out) {
  writePDB(out, 1, ' ');
  out << "TER\n";
}

void gprotein::writePDB(ostream &out, int startRes, char chainID) {
  char line[81];
  char chain = chainID;
  int iatom = 1;
  for (int i = 0; i < length; i++) {
    // print the heavy atoms first
    for (int j = 0; j < getResidue(i)->getNA(); j++) {
      sprintf(&line[0],
              "ATOM  %5d  %-4.4s%3.3s %1c%4d    %8.3lf%8.3lf%8.3lf  1.00  0.00 "
              "          %c  ",
              iatom++, getResidue(i)->getAtom(j)->getName(),
              getResidue(i)->getName(), chain, i + startRes,
              getResidue(i)->getAtom(j)->getR()->getX(),
              getResidue(i)->getAtom(j)->getR()->getY(),
              getResidue(i)->getAtom(j)->getR()->getZ(),
              getResidue(i)->getAtom(j)->getName()[0]);
      out << line << "\n";
    }
    // print the polar protons
    for (int j = 0; j < getResidue(i)->getNA(); j++) {
      for (int k = 0; k < getResidue(i)->getAtom(j)->n_attachedProtons; k++) {
        char atomname[10];
        int len;
        if (len =
                strlen(
                    getResidue(i)->getAtom(j)->attachedProtons[k].getName()) ==
                4) {
          atomname[0] =
              getResidue(i)->getAtom(j)->attachedProtons[k].getName()[3];
          for (int c = 0; c < 3; c++)
            atomname[c + 1] =
                getResidue(i)->getAtom(j)->attachedProtons[k].getName()[c];
          atomname[4] = '\0';
          sprintf(&line[0],
                  "ATOM  %5d %-4.4s %3.3s %1c%4d    %8.3lf%8.3lf%8.3lf  1.00  "
                  "0.00           %c  ",
                  iatom++, atomname, getResidue(i)->getName(), chain,
                  i + startRes,
                  getResidue(i)->getAtom(j)->attachedProtons[k].getR()->getX(),
                  getResidue(i)->getAtom(j)->attachedProtons[k].getR()->getY(),
                  getResidue(i)->getAtom(j)->attachedProtons[k].getR()->getZ(),
                  getResidue(i)->getAtom(j)->attachedProtons[k].getName()[0]);
        } else {
          sprintf(&line[0],
                  "ATOM  %5d  %-4.4s%3.3s %1c%4d    %8.3lf%8.3lf%8.3lf  1.00  "
                  "0.00           %c  ",
                  iatom++,
                  getResidue(i)->getAtom(j)->attachedProtons[k].getName(),
                  getResidue(i)->getName(), chain, i + startRes,
                  getResidue(i)->getAtom(j)->attachedProtons[k].getR()->getX(),
                  getResidue(i)->getAtom(j)->attachedProtons[k].getR()->getY(),
                  getResidue(i)->getAtom(j)->attachedProtons[k].getR()->getZ(),
                  getResidue(i)->getAtom(j)->attachedProtons[k].getName()[0]);
        }
        out << line << "\n";
      }
    }
  }
}

void gprotein::rotateTheta13(int ires, const double &dtheta) {
  if (ires > 0 && ires < (size() - 1)) {
    getDimondRotateMatrix(getResidue(ires - 1)->getCA()->r,
                          getResidue(ires + 1)->getCA()->r, dtheta, getM());
    // C,O of Residue ires-1
    atom *theAtom = getResidue(ires - 1)->getC();
    vec tmp = theAtom->r - getResidue(ires + 1)->getCA()->r;
    tmp.EularRotate(getM());
    theAtom->r = tmp + getResidue(ires + 1)->getCA()->r;

    theAtom = getResidue(ires - 1)->getO();
    tmp = theAtom->r - getResidue(ires + 1)->getCA()->r;
    tmp.EularRotate(getM());
    theAtom->r = tmp + getResidue(ires + 1)->getCA()->r;

    // N,H of Residue ires+1;
    theAtom = getResidue(ires + 1)->getN();
    tmp = theAtom->r - getResidue(ires + 1)->getCA()->r;
    tmp.EularRotate(getM());
    theAtom->r = tmp + getResidue(ires + 1)->getCA()->r;
    theAtom = theAtom->attachedProtons;
    if (getResidue(ires + 1)->getType() == PRO)
      theAtom = getResidue(ires + 1)->getAtom(6); // PROLINE, the CD Atom
    if (theAtom) {
      tmp = theAtom->r - getResidue(ires + 1)->getCA()->r;
      tmp.EularRotate(getM());
      theAtom->r = tmp + getResidue(ires + 1)->getCA()->r;
    }
    if (getResidue(ires + 1)->getType() == PRO)
      getResidue(ires + 1)->updateRotamers();

    // Residue ires
    for (int i = 0; i < getResidue(ires)->getNA(); i++) {
      theAtom = getResidue(ires)->getAtom(i);
      tmp = theAtom->r - getResidue(ires + 1)->getCA()->r;
      tmp.EularRotate(getM());
      theAtom->r = tmp + getResidue(ires + 1)->getCA()->r;
      for (int j = 0; j < getResidue(ires)->getAtom(i)->n_attachedProtons;
           j++) {
        theAtom = getResidue(ires)->getAtom(i)->attachedProtons + j;
        tmp = theAtom->r - getResidue(ires + 1)->getCA()->r;
        tmp.EularRotate(getM());
        theAtom->r = tmp + getResidue(ires + 1)->getCA()->r;
      }
    }
  }
}

void gprotein::rotateTheta12(int ires, const double &dtheta) {
  if (ires > 0) {
    getDimondRotateMatrix(getResidue(ires - 1)->getCA()->r,
                          getResidue(ires)->getCA()->r, dtheta, getM());
    // C,O of Residue ires-1
    atom *theAtom = getResidue(ires - 1)->getC();
    vec tmp = theAtom->r - getResidue(ires)->getCA()->r;
    tmp.EularRotate(getM());
    theAtom->r = tmp + getResidue(ires)->getCA()->r;

    theAtom = getResidue(ires - 1)->getO();
    tmp = theAtom->r - getResidue(ires)->getCA()->r;
    tmp.EularRotate(getM());
    theAtom->r = tmp + getResidue(ires)->getCA()->r;

    // N,H of Residue ires
    theAtom = getResidue(ires)->getN();
    tmp = theAtom->r - getResidue(ires)->getCA()->r;
    tmp.EularRotate(getM());
    theAtom->r = tmp + getResidue(ires)->getCA()->r;
    theAtom = theAtom->attachedProtons;
    if (getResidue(ires)->getType() == PRO)
      theAtom = getResidue(ires)->getAtom(6); // PROLINE, the CD Atom
    if (theAtom) {
      tmp = theAtom->r - getResidue(ires)->getCA()->r;
      tmp.EularRotate(getM());
      theAtom->r = tmp + getResidue(ires)->getCA()->r;
    }
    if (getResidue(ires)->getType() == PRO)
      getResidue(ires)->updateRotamers();
  }
}

/*Rotate around the axis of CAs of residues ires+1 and ires.
  dtheta is in RADIAN. */
void gprotein::rotateTheta32(int ires, const double &dtheta) {
  if (ires > 0) {
    getDimondRotateMatrix(getResidue(ires + 1)->getCA()->r,
                          getResidue(ires)->getCA()->r, dtheta, getM());
    // C,O of Residue ires
    atom *theAtom = getResidue(ires)->getC();
    vec tmp = theAtom->r - getResidue(ires)->getCA()->r;
    tmp.EularRotate(getM());
    theAtom->r = tmp + getResidue(ires)->getCA()->r;

    theAtom = getResidue(ires)->getO();
    tmp = theAtom->r - getResidue(ires)->getCA()->r;
    tmp.EularRotate(getM());
    theAtom->r = tmp + getResidue(ires)->getCA()->r;

    // N,H of Residue ires+1
    theAtom = getResidue(ires + 1)->getN();
    tmp = theAtom->r - getResidue(ires)->getCA()->r;
    tmp.EularRotate(getM());
    theAtom->r = tmp + getResidue(ires)->getCA()->r;
    theAtom = theAtom->attachedProtons;
    if (getResidue(ires + 1)->getType() == PRO)
      theAtom = getResidue(ires + 1)->getAtom(6); // PROLINE, the CD Atom
    if (theAtom) {
      tmp = theAtom->r - getResidue(ires)->getCA()->r;
      tmp.EularRotate(getM());
      theAtom->r = tmp + getResidue(ires)->getCA()->r;
    }
    if (getResidue(ires + 1)->getType() == PRO)
      getResidue(ires + 1)->updateRotamers();
  }
}

void gprotein::idealize() {
  gprotein *ps = this;
  // Update the rotamers first.
  for (int ires = 0; ires < ps->size(); ires++) {
    residue *theRes = ps->getResidue(ires);
    theRes->updateRotamers();
  }
  // Re-assign backbone coordinates from templates.
  residue *preRes = NULL;
  for (int ires = 0; ires < ps->size(); ires++) {
    residue *theRes = ps->getResidue(ires);
    vec xbar, ybar, zbar; // Vectors to align residue backbones.
    if (!preRes)          // For the first residue.
    {
      vec N = theRes->getN()->r;
      vec CA = theRes->getCA()->r;
      vec C = theRes->getC()->r;
      xbar = (CA - N).norm();
      zbar = xbar ^ (C - CA);
      zbar.Normalize();
      ybar = zbar ^ xbar;
      theRes->initTemplate();
      theRes->rotate(xbar, ybar, zbar);
      theRes->shift(CA);
      theRes->updateBBH();
    } else { // Not the first residue in a chain.
      vec nextNShift = preRes->nextN(xbar, ybar, zbar);
      theRes->initTemplate();
      theRes->rotate(xbar, ybar, zbar);
      vec shiftCA = nextNShift - theRes->getN()->r; // shift need for CA.
      theRes->shift(shiftCA);
      theRes->updateBBH(*preRes);
    }
    theRes->updatePolarSH();
    // I don't updateRotamers here so that the old chi angles are still
    // store in theRes->chis[] array and will be used to restore the sidechain.
    // theRes->updateRotamers();
    // Prepare for next residue.
    preRes = theRes;
  }
  ps->updateFY_UV();

  // Copy the dihedral angles.
  for (int ires = 0; ires < ps->size(); ires++) {
    double phi, psi, omega;
    ps->getFY(ires, phi, psi);
    omega = ps->getOmiga(ires);
    // Rotate phi angle.
    if (ires > 0) {
      ps->phi[ires] = PI;
      ps->rotatePhi(ires, phi);
    }
    // Rotate phs angle.
    {
      ps->psi[ires] = PI;
      ps->rotatePsi(ires, psi);
    }
    // Rotate omega angle.
    if (ires < ps->size() - 1) {
      ps->omiga[ires] = PI;
      ps->rotateOmiga(ires, omega);
    }
  }

  // Copy side chains dihedrals from old comformation
  // (chi angles are store in array chi[] in each residue) .
  for (int ires = 0; ires < ps->size(); ires++) {
    residue *theRes = ps->getResidue(ires);
    double oldChis[4];
    for (int iChi = 0; iChi < theRes->getNChi(); iChi++)
      oldChis[iChi] = theRes->getChi(iChi);
    theRes->updateRotamers();
    for (int iChi = 0; iChi < theRes->getNChi(); iChi++)
      theRes->rotateChi(iChi, oldChis[iChi], ps->getM());
  }
}

bool gprotein::updateFY(int i, int &newF, int &newY) {
  int iF, iY;
  FY_double2int(phi[i], psi[i], iF, iY);
  if (i > 0) {
    phi[i] = getDihedralAngle(
        gres[i - 1].getResidue()->getC()->r, gres[i].getResidue()->getN()->r,
        gres[i].getResidue()->getCA()->r, gres[i].getResidue()->getC()->r);
  } else {
    newF = iF;
  }
  if (i < size() - 1) {
    psi[i] = getDihedralAngle(
        gres[i].getResidue()->getN()->r, gres[i].getResidue()->getCA()->r,
        gres[i].getResidue()->getC()->r, gres[i + 1].getResidue()->getN()->r);
  } else {
    newY = iY;
  }
  FY_double2int(phi[i], psi[i], newF, newY);
  if ((newF != iF) || (newY != iY)) {
    return true;
  }
  return false;
}

void gen_grid::updateCellAtom(atom *theAtom) {
  int ix, iy, iz, cindex, pcindex, nh;
  ix = static_cast<int>((theAtom->r.x - min.x) / cell_box.x);
  iy = static_cast<int>((theAtom->r.y - min.y) / cell_box.y);
  iz = static_cast<int>((theAtom->r.z - min.z) / cell_box.z);
  cindex = ix + (iy + iz * ncell[1]) * ncell[0];
  if (theAtom->cell_index != cindex) {
    pcindex = theAtom->cell_index;
    // drop from the current cell
    if (cell_atoms[pcindex] == theAtom) { // top
      cell_atoms[pcindex] = theAtom->cell_next;
      if (theAtom->cell_next)
        theAtom->cell_next->cell_prev = NULL;
    } else {
      theAtom->cell_prev->cell_next = theAtom->cell_next;
      if (theAtom->cell_next) {
        theAtom->cell_next->cell_prev = theAtom->cell_prev;
      }
    }
    theAtom->cell_index = cindex;
    // add to new cell
    if (!cell_atoms[cindex]) { // FIRST ONE
      cell_atoms[cindex] = theAtom;
      cell_atoms[cindex]->cell_prev = NULL;
      cell_atoms[cindex]->cell_next = NULL;
    } else { // put on the top
      cell_atoms[cindex]->cell_prev = theAtom;
      theAtom->cell_next = cell_atoms[cindex];
      cell_atoms[cindex] = theAtom;
      cell_atoms[cindex]->cell_prev = NULL;
    }
  }

  nh = theAtom->n_attachedProtons;
  for (int k = 0; k < nh; k++) {
    atom *theH = theAtom->attachedProtons + k;

    ix = static_cast<int>((theH->r.x - min.x) / cell_box.x);
    iy = static_cast<int>((theH->r.y - min.y) / cell_box.y);
    iz = static_cast<int>((theH->r.z - min.z) / cell_box.z);
    cindex = ix + (iy + iz * ncell[1]) * ncell[0];
    if (theH->cell_index != cindex) {
      pcindex = theH->cell_index;
      // drop from the current cell
      if (cell_atoms[pcindex] == theH) { // top
        cell_atoms[pcindex] = theH->cell_next;
        if (theH->cell_next)
          theH->cell_next->cell_prev = NULL;
      } else {
        theH->cell_prev->cell_next = theH->cell_next;
        if (theH->cell_next) {
          theH->cell_next->cell_prev = theH->cell_prev;
        }
      }
      theH->cell_index = cindex;
      // add to new cell
      if (!cell_atoms[cindex]) { // FIST ONE
        cell_atoms[cindex] = theH;
        cell_atoms[cindex]->cell_prev = NULL;
        cell_atoms[cindex]->cell_next = NULL;
      } else { // put from the beginning
        cell_atoms[cindex]->cell_prev = theH;
        theH->cell_next = cell_atoms[cindex];
        cell_atoms[cindex] = theH;
        cell_atoms[cindex]->cell_prev = NULL;
      }
    }
  }
}

gen_grid::gen_grid(protein &theProtein, gprotein &theGProtein, VDW &theVDW,
                   double mir = 9.0)
    : grid(theProtein, theVDW, mir), gp(theGProtein) {
  max_ir = mir;
  max_ir2 = max_ir * max_ir;
  p.getBoundary(min, max);
  // ext range changed to 4 times maximum dimension of protein
  vec tmp = (max - min);
  double tmp2 = tmp.x;
  tmp2 = tmp2 > tmp.y ? tmp2 : tmp.y;
  tmp2 = tmp2 > tmp.z ? tmp2 : tmp.z;
  double ext = 4 * tmp2;
  // double ext = 200.0;
  min -= vec(ext, ext, ext);
  max += vec(ext, ext, ext);

  // assign cell_box, ncell;
  constructCells();
  ncells = ncell[0] * ncell[1] * ncell[2];
  // cout << ncells << endl;

  delete[] cell_atoms;
  cell_atoms = new atom *[ncells];
  for (int i = 0; i < ncells; i++)
    cell_atoms[i] = NULL;

  // initialize the cells_atoms
  for (int i = 0; i < gp.size(); i++) {
    residue *theRes = gp.getResidue(i);
    atom *theAtom = NULL;
    atom *theH = NULL;
    int ix, iy, iz, cindex, nh;
    for (int j = 0; j < theRes->getNA(); j++) {
      theAtom = theRes->getAtom(j);
      ix = static_cast<int>((theAtom->r.x - min.x) / cell_box.x);
      iy = static_cast<int>((theAtom->r.y - min.y) / cell_box.y);
      iz = static_cast<int>((theAtom->r.z - min.z) / cell_box.z);
      theAtom->cell_index = cindex = ix + (iy + iz * ncell[1]) * ncell[0];
      // cout << cindex << " " << theAtom->getName() << endl;
      if (!cell_atoms[cindex]) { // FIST ONE
        cell_atoms[cindex] = theAtom;
        cell_atoms[cindex]->cell_prev = NULL;
        cell_atoms[cindex]->cell_next = NULL;
      } else { // put from the beginning
        cell_atoms[cindex]->cell_prev = theAtom;
        theAtom->cell_next = cell_atoms[cindex];
        cell_atoms[cindex] = theAtom;
        cell_atoms[cindex]->cell_prev = NULL;
      }
      nh = theAtom->n_attachedProtons;
      // cout << "Protons" << nh << endl;
      for (int k = 0; k < nh; k++) {
        theH = theAtom->attachedProtons + k;
        ix = static_cast<int>((theH->r.x - min.x) / cell_box.x);
        iy = static_cast<int>((theH->r.y - min.y) / cell_box.y);
        iz = static_cast<int>((theH->r.z - min.z) / cell_box.z);
        theH->cell_index = cindex = ix + (iy + iz * ncell[1]) * ncell[0];
        if (!cell_atoms[cindex]) { // FIST ONE
          cell_atoms[cindex] = theH;
          cell_atoms[cindex]->cell_prev = NULL;
          cell_atoms[cindex]->cell_next = NULL;
        } else { // put from the beginning
          cell_atoms[cindex]->cell_prev = theH;
          theH->cell_next = cell_atoms[cindex];
          cell_atoms[cindex] = theH;
          cell_atoms[cindex]->cell_prev = NULL;
        }
      }
    }
  }
}

void gen_grid::dropCellResidue(int iRes) {
  int ix, iy, iz, cindex, pcindex, nh;
  residue *theRes = gp.getResidue(iRes);
  for (int i = 0; i < theRes->getNA(); i++) {
    atom *theAtom = theRes->getAtom(i);
    cindex = theAtom->cell_index;
    // drop from the current cell
    if (cell_atoms[cindex] == theAtom) { // top
      cell_atoms[cindex] = theAtom->cell_next;
      if (theAtom->cell_next)
        theAtom->cell_next->cell_prev = NULL;
    } else {
      theAtom->cell_prev->cell_next = theAtom->cell_next;
      if (theAtom->cell_next) {
        theAtom->cell_next->cell_prev = theAtom->cell_prev;
      }
    }
    nh = theAtom->n_attachedProtons;
    for (int k = 0; k < nh; k++) {
      atom *theH = theAtom->attachedProtons + k;
      cindex = theH->cell_index;
      if (cell_atoms[cindex] == theH) { // top
        cell_atoms[cindex] = theH->cell_next;
        if (theH->cell_next)
          theH->cell_next->cell_prev = NULL;
      } else {
        theH->cell_prev->cell_next = theH->cell_next;
        if (theH->cell_next) {
          theH->cell_next->cell_prev = theH->cell_prev;
        }
      }
    }
  }
}

void gen_grid::insertCellResidue(int iRes) {
  int ix, iy, iz, cindex, pcindex, nh;
  residue *theRes = gp.getResidue(iRes);
  for (int i = 0; i < theRes->getNA(); i++) {
    atom *theAtom = theRes->getAtom(i);

    ix = static_cast<int>((theAtom->r.x - min.x) / cell_box.x);
    iy = static_cast<int>((theAtom->r.y - min.y) / cell_box.y);
    iz = static_cast<int>((theAtom->r.z - min.z) / cell_box.z);
    cindex = ix + (iy + iz * ncell[1]) * ncell[0];
    theAtom->cell_index = cindex;
    // add to new cell
    if (!cell_atoms[cindex]) { // FIRST ONE
      cell_atoms[cindex] = theAtom;
      cell_atoms[cindex]->cell_prev = NULL;
      cell_atoms[cindex]->cell_next = NULL;
    } else { // put on the top
      cell_atoms[cindex]->cell_prev = theAtom;
      theAtom->cell_next = cell_atoms[cindex];
      cell_atoms[cindex] = theAtom;
      cell_atoms[cindex]->cell_prev = NULL;
    }
    nh = theAtom->n_attachedProtons;
    for (int k = 0; k < nh; k++) {
      atom *theH = theAtom->attachedProtons + k;

      ix = static_cast<int>((theH->r.x - min.x) / cell_box.x);
      iy = static_cast<int>((theH->r.y - min.y) / cell_box.y);
      iz = static_cast<int>((theH->r.z - min.z) / cell_box.z);
      cindex = ix + (iy + iz * ncell[1]) * ncell[0];
      theH->cell_index = cindex;
      // add to new cell
      if (!cell_atoms[cindex]) { // FIST ONE
        cell_atoms[cindex] = theH;
        cell_atoms[cindex]->cell_prev = NULL;
        cell_atoms[cindex]->cell_next = NULL;
      } else { // put from the beginning
        cell_atoms[cindex]->cell_prev = theH;
        theH->cell_next = cell_atoms[cindex];
        cell_atoms[cindex] = theH;
        cell_atoms[cindex]->cell_prev = NULL;
      }
    }
  }
}

double gen_grid::getEVS_Atom(atom *theAtom, double &E_VDW_attr,
                             double &E_VDW_rep, double &E_SOLV) {
  static double coef = 2.0 / (4.0 * PI * sqrt(PI));
  int ti, tj;
  int ix, iy, iz, cindex, kx, ky, kz;
  double d2, d6, B, A, e, d, sig, eps, cutoff, slope, intercept;
  double SOLVi = 0, xij, xji, solv_ij;
  atom *pt;

  ti = theAtom->FF_t;
  ix = static_cast<int>((theAtom->r.x - min.x) / cell_box.x);
  iy = static_cast<int>((theAtom->r.y - min.y) / cell_box.y);
  iz = static_cast<int>((theAtom->r.z - min.z) / cell_box.z);

  for (int cz = iz - step; cz <= iz + step; cz++) {
    kz = cz;
    if (cz < 0)
      kz += ncell[2];
    else if (cz >= ncell[2])
      kz -= ncell[2];

    for (int cy = iy - step; cy <= iy + step; cy++) {
      ky = cy;
      if (cy < 0)
        ky += ncell[1];
      else if (cy >= ncell[1])
        ky -= ncell[1];

      for (int cx = ix - step; cx <= ix + step; cx++) {
        kx = cx;
        if (cx < 0)
          kx += ncell[0];
        else if (cx >= ncell[0])
          kx -= ncell[0];
        cindex = ncell[0] * (ncell[1] * kz + ky) + kx;

        pt = cell_atoms[cindex];
        while (pt) {
          if (pt != theAtom &&
              !pt->isH()) { // not the same atom and not protons. Proton is not
                            // accounted for VDW and SOLV
            d2 = (pt->r - theAtom->r).mod2();
            d = sqrt(d2);
            tj = pt->FF_t;
            if (d2 < max_ir2) { // within the cutoff distance
              // vdw
              if (!isInternal(theAtom,
                              pt)) { // vdw does not include the internal pairs
                vdw.getVDWdata(ti, tj, B, A, eps, sig, cutoff, slope,
                               intercept);
                if (B && A) {
                  d6 = d2 * d2 * d2;
                  e = (B / d6 - A) / d6;
                  if (e > 0) {
                    if (d < cutoff) {
                      e = slope * d + intercept;
                    }
                    E_VDW_rep += e;
                  } else {
                    E_VDW_attr += e;
                  }
                }

                // solvation
                if (d < (EEF1_SOLV[ti][CHARMM_RADIUS] +
                         EEF1_SOLV[tj][CHARMM_RADIUS]) *
                            0.846) {
                  d = (EEF1_SOLV[ti][CHARMM_RADIUS] +
                       EEF1_SOLV[tj][CHARMM_RADIUS]) *
                      0.846;
                  d2 = d * d;
                }
                xij = (d - EEF1_SOLV[ti][CHARMM_RADIUS]) /
                      EEF1_SOLV[ti][EEF1_LAMBDA];
                xji = (d - EEF1_SOLV[tj][CHARMM_RADIUS]) /
                      EEF1_SOLV[tj][EEF1_LAMBDA];

                solv_ij =
                    (EEF1_SOLV[ti][EEF1_DG_FREE] * EEF1_SOLV[tj][EEF1_VOL] *
                         exp(-xij * xij) / EEF1_SOLV[ti][EEF1_LAMBDA] +
                     EEF1_SOLV[tj][EEF1_DG_FREE] * EEF1_SOLV[ti][EEF1_VOL] *
                         exp(-xji * xji) / EEF1_SOLV[tj][EEF1_LAMBDA]) *
                    coef / d2;
                SOLVi -= solv_ij;
              }
            }
          }
          pt = pt->cell_next;
        }
      }
    }
  }
  E_SOLV += SOLVi;
}

int gen_grid::updateCellResidue(int iRes) {
  residue *theRes = gp.getResidue(iRes);
  // cout << theRes->getNA() << endl;
  for (int i = 0; i < theRes->getNA(); i++) {
    atom *theAtom = theRes->getAtom(i);
    updateCellAtom(theAtom);
  }
}

double gen_grid::getEVS_RES_SC(int iRes, double &E_VDW_attr, double &E_VDW_rep,
                               double &E_SOLV, double *VDWA_SS, double *VDWR_SS,
                               double *SOLV_SS, double *VDWA_SB,
                               double *VDWR_SB, double *SOLV_SB) {
  E_VDW_attr = E_VDW_rep = E_SOLV = 0;
  double vdw_rep, vdw_attr, solv;
  residue *theRes = gp.getResidue(iRes);
  if (theRes->getType() == GLY)
    return 0; // GLY does not have sidechain(a H)
  atom *theAtom = NULL;
  static double coef = 2.0 / (4.0 * PI * sqrt(PI));
  int ti, tj;
  int ix, iy, iz, cindex, kx, ky, kz;
  int min_ix, min_iy, min_iz;
  double d2, d6, B, A, e, d, sig, eps, cutoff, slope, intercept;
  double SOLVi = 0, xij, xji, solv_ij;
  atom *pt;
  for (int i = 0; i < theRes->getNA(); i++) {
    theAtom = theRes->getAtom(i);
    if (theAtom->isSC()) { // sidechain atoms include CB
      ti = theAtom->FF_t;
      ix = static_cast<int>((theAtom->r.x - min.x) / cell_box.x);
      iy = static_cast<int>((theAtom->r.y - min.y) / cell_box.y);
      iz = static_cast<int>((theAtom->r.z - min.z) / cell_box.z);

      // SOLVi=EEF1_SOLV[ti][EEF1_DG_REF];
      SOLVi = 0;
      for (int cz = iz - step; cz <= iz + step; cz++) {
        kz = cz;
        if (cz < 0)
          kz += ncell[2];
        else if (cz >= ncell[2])
          kz -= ncell[2];

        for (int cy = iy - step; cy <= iy + step; cy++) {
          ky = cy;
          if (cy < 0)
            ky += ncell[1];
          else if (cy >= ncell[1])
            ky -= ncell[1];

          for (int cx = ix - step; cx <= ix + step; cx++) {
            kx = cx;
            if (cx < 0)
              kx += ncell[0];
            else if (cx >= ncell[0])
              kx -= ncell[0];
            cindex = ncell[0] * (ncell[1] * kz + ky) + kx;

            pt = cell_atoms[cindex];
            while (pt) {
              if (pt != theAtom &&
                  !pt->isH()) { // not the same atom and not protons. Proton is
                                // not accounted for VDW and SOLV
                d2 = (pt->r - theAtom->r).mod2();
                d = sqrt(d2);
                tj = pt->FF_t;
                if (d2 < max_ir2) { // within the cutoff distance
                  // vdw
                  if (isIncluded(
                          theAtom,
                          pt)) { // vdw does not include the internal pairs
                    vdw.getVDWdata(ti, tj, B, A, eps, sig, cutoff, slope,
                                   intercept);
                    if (B && A) {
                      d6 = d2 * d2 * d2;
                      e = (B / d6 - A) / d6;
                      if (e > 0) {
                        if (d < cutoff) {
                          e = slope * d + intercept;
                        }
                        // e = (sig-d)/sig*10.0;
                        E_VDW_rep += e;
                        // cout << theAtom->resID+1 << " " << pt->resID+1 << " "
                        // << d << " " << sig << " " << e <<  endl;
                        if (VDWR_SS && pt->isSC()) {
                          VDWR_SS[pt->resID] += e;
                        } else if (VDWR_SB && !pt->isSC()) {
                          VDWR_SB[pt->resID] += e;
                        }
                      } else {
                        E_VDW_attr += e;
                        if (VDWA_SS && pt->isSC()) {
                          VDWA_SS[pt->resID] += e;
                        } else if (VDWA_SB && !pt->isSC()) {
                          VDWA_SB[pt->resID] += e;
                        }
                      }
                    }

                    if (d < (EEF1_SOLV[ti][CHARMM_RADIUS] +
                             EEF1_SOLV[tj][CHARMM_RADIUS]) *
                                0.846) {
                      d = (EEF1_SOLV[ti][CHARMM_RADIUS] +
                           EEF1_SOLV[tj][CHARMM_RADIUS]) *
                          0.846;
                      d2 = d * d;
                    }
                    xij = (d - EEF1_SOLV[ti][CHARMM_RADIUS]) /
                          EEF1_SOLV[ti][EEF1_LAMBDA];
                    xji = (d - EEF1_SOLV[tj][CHARMM_RADIUS]) /
                          EEF1_SOLV[tj][EEF1_LAMBDA];

                    solv_ij =
                        (EEF1_SOLV[ti][EEF1_DG_FREE] * EEF1_SOLV[tj][EEF1_VOL] *
                             exp(-xij * xij) / EEF1_SOLV[ti][EEF1_LAMBDA] +
                         EEF1_SOLV[tj][EEF1_DG_FREE] * EEF1_SOLV[ti][EEF1_VOL] *
                             exp(-xji * xji) / EEF1_SOLV[tj][EEF1_LAMBDA]) *
                        coef / d2;
                    SOLVi -= solv_ij;
                    if (SOLV_SS && pt->isSC()) {
                      SOLV_SS[pt->resID] -= solv_ij;
                    } else if (SOLV_SB && !pt->isSC()) {
                      SOLV_SB[pt->resID] -= solv_ij;
                    }
                  }
                }
              }
              pt = pt->cell_next;
            }
          }
        }
      }
      E_SOLV += SOLVi;
    }
  }
}

/*include the SALT-BRIDGE interaction, ~QiQj/(Rij*Rij) */
double gen_grid::getEVS_RES_SC(int iRes, double &E_VDW_attr, double &E_VDW_rep,
                               double &E_SOLV, double &E_SB, double *VDWA_SS,
                               double *VDWR_SS, double *SOLV_SS,
                               double *VDWA_SB, double *VDWR_SB,
                               double *SOLV_SB, double *SB_SS) {
  E_VDW_attr = E_VDW_rep = E_SOLV = E_SB = 0;
  double vdw_rep, vdw_attr, solv;
  residue *theRes = gp.getResidue(iRes);
  if (theRes->getType() == GLY)
    return 0; // GLY does not have sidechain(a H)
  atom *theAtom = NULL;
  static double coef = 2.0 / (4.0 * PI * sqrt(PI));
  int ti, tj;
  int ix, iy, iz, cindex, kx, ky, kz;
  int min_ix, min_iy, min_iz;
  double d2, d6, B, A, e, d, sig, eps, cutoff, slope, intercept;
  double SOLVi = 0, xij, xji, solv_ij;
  atom *pt;
  for (int i = 0; i < theRes->getNA(); i++) {
    theAtom = theRes->getAtom(i);
    if (theAtom->isSC()) { // sidechain atoms include CB
      ti = theAtom->FF_t;
      ix = static_cast<int>((theAtom->r.x - min.x) / cell_box.x);
      iy = static_cast<int>((theAtom->r.y - min.y) / cell_box.y);
      iz = static_cast<int>((theAtom->r.z - min.z) / cell_box.z);

      SOLVi = 0;
      for (int cz = iz - step; cz <= iz + step; cz++) {
        kz = cz;
        if (cz < 0)
          kz += ncell[2];
        else if (cz >= ncell[2])
          kz -= ncell[2];

        for (int cy = iy - step; cy <= iy + step; cy++) {
          ky = cy;
          if (cy < 0)
            ky += ncell[1];
          else if (cy >= ncell[1])
            ky -= ncell[1];

          for (int cx = ix - step; cx <= ix + step; cx++) {
            kx = cx;
            if (cx < 0)
              kx += ncell[0];
            else if (cx >= ncell[0])
              kx -= ncell[0];
            cindex = ncell[0] * (ncell[1] * kz + ky) + kx;

            pt = cell_atoms[cindex];
            while (pt) {
              if (pt != theAtom &&
                  !pt->isH()) { // not the same atom and not protons. Proton is
                                // not accounted for VDW and SOLV
                d2 = (pt->r - theAtom->r).mod2();
                d = sqrt(d2);
                tj = pt->FF_t;
                if (d2 < max_ir2) { // within the cutoff distance
                  // vdw
                  if (isIncluded(
                          theAtom,
                          pt)) { // vdw does not include the internal pairs
                    vdw.getVDWdata(ti, tj, B, A, eps, sig, cutoff, slope,
                                   intercept);
                    if (B && A) {
                      d6 = d2 * d2 * d2;
                      e = (B / d6 - A) / d6;
                      if (e > 0) {
                        if (d < cutoff) {
                          e = slope * d + intercept;
                        }
                        E_VDW_rep += e;
                        if (VDWR_SS && pt->isSC()) {
                          VDWR_SS[pt->resID] += e;
                        } else if (VDWR_SB && !pt->isSC()) {
                          VDWR_SB[pt->resID] += e;
                        }
                      } else {
                        E_VDW_attr += e;
                        if (VDWA_SS && pt->isSC()) {
                          VDWA_SS[pt->resID] += e;
                        } else if (VDWA_SB && !pt->isSC()) {
                          VDWA_SB[pt->resID] += e;
                        }
                      }
                    }

                    if (theAtom->icharge && pt->icharge) { // i-charged:)
                      if (d < sig) {
                        e = (theAtom->icharge * pt->icharge) / (sig * sig);
                      } else {
                        e = (theAtom->icharge * pt->icharge) / d2;
                      }
                      E_SB += e;
                      if (SB_SS) {
                        SB_SS[pt->resID] += e;
                      }
                    }

                    // solvation
                    if (d < (EEF1_SOLV[ti][CHARMM_RADIUS] +
                             EEF1_SOLV[tj][CHARMM_RADIUS]) *
                                0.846) {
                      d = (EEF1_SOLV[ti][CHARMM_RADIUS] +
                           EEF1_SOLV[tj][CHARMM_RADIUS]) *
                          0.846;
                      d2 = d * d;
                    }
                    xij = (d - EEF1_SOLV[ti][CHARMM_RADIUS]) /
                          EEF1_SOLV[ti][EEF1_LAMBDA];
                    xji = (d - EEF1_SOLV[tj][CHARMM_RADIUS]) /
                          EEF1_SOLV[tj][EEF1_LAMBDA];

                    solv_ij =
                        (EEF1_SOLV[ti][EEF1_DG_FREE] * EEF1_SOLV[tj][EEF1_VOL] *
                             exp(-xij * xij) / EEF1_SOLV[ti][EEF1_LAMBDA] +
                         EEF1_SOLV[tj][EEF1_DG_FREE] * EEF1_SOLV[ti][EEF1_VOL] *
                             exp(-xji * xji) / EEF1_SOLV[tj][EEF1_LAMBDA]) *
                        coef / d2;
                    SOLVi -= solv_ij;
                    if (SOLV_SS && pt->isSC()) {
                      SOLV_SS[pt->resID] -= solv_ij;
                    } else if (SOLV_SB && !pt->isSC()) {
                      SOLV_SB[pt->resID] -= solv_ij;
                    }
                  }
                }
              }
              pt = pt->cell_next;
            }
          }
        }
      }
      E_SOLV += SOLVi;
    }
  }
}

/******************************************************************************************
  To obtain the derivative of VDW and SOLV.
  E = E(r)

  dE/dChi = dE/dr*dr/dChi.
  dE/dr can be calculated easily.

  dr/dChi = 1/r*{[Vx(R_m-R_o)]*(R_m-R_s)}
  Here, V is the unit vector of the rotation axis.
  R_o is the position of rotation end atom (CB of Chi1)
  R_m is the position of the rotating atom and R_s is the position of the static
 atom. The distance r is between static and rotating atoms. r=|R_m - R_s|

  REMBER!!! Before call this function, updateChiUV should be called at first

  !!!!!The derive is in the unit of PER RADIAN!!!!!!
 *******************************************************************************************/

double gen_grid::getE_DEV_RES_SC(int iRes, double &E_VDW_attr,
                                 double &E_VDW_rep, double &E_SOLV,
                                 double *DEV_VDW_attr, double *DEV_VDW_rep,
                                 double *DEV_SOLV) {
  E_VDW_attr = E_VDW_rep = E_SOLV = 0;
  double vdw_rep, vdw_attr, solv;
  residue *theRes = gp.getResidue(iRes);
  if (theRes->getType() == GLY)
    return 0; // GLY does not have sidechain(a H)
  atom *theAtom = NULL;
  static double coef = 2.0 / (4.0 * PI * sqrt(PI));
  static double coef2 = 1.0 / (PI * sqrt(PI));
  int ti, tj;
  int ix, iy, iz, cindex, kx, ky, kz;
  int min_ix, min_iy, min_iz;
  double d2, d6, B, A, e, d, sig, eps, cutoff, slope, intercept;
  double SOLVi = 0, xij, xji;
  atom *pt;

  int nCHI_affected;
  int nCHI = theRes->getNChi();
  double de;
  for (int ich = 0; ich < nCHI + 1; ich++) {
    DEV_VDW_attr[ich] = 0;
    DEV_VDW_rep[ich] = 0;
    DEV_SOLV[ich] = 0;
  }

  for (int i = 0; i < theRes->getNA(); i++) {
    theAtom = theRes->getAtom(i);
    nCHI_affected = theRes->getNChi_AFFECTED(i);
    if (theAtom->isSC()) { // sidechain atoms include CB
      ti = theAtom->FF_t;
      ix = static_cast<int>((theAtom->r.x - min.x) / cell_box.x);
      iy = static_cast<int>((theAtom->r.y - min.y) / cell_box.y);
      iz = static_cast<int>((theAtom->r.z - min.z) / cell_box.z);

      // SOLVi=EEF1_SOLV[ti][EEF1_DG_REF];
      SOLVi = 0;
      for (int cz = iz - step; cz <= iz + step; cz++) {
        // min_iz = abs(cz-iz);
        // if(min_iz>0) min_iz--;
        kz = cz;
        if (cz < 0)
          kz += ncell[2];
        else if (cz >= ncell[2])
          kz -= ncell[2];

        for (int cy = iy - step; cy <= iy + step; cy++) {
          // min_iy = abs(cy-iy);
          // if(min_iy>0) min_iy--;
          ky = cy;
          if (cy < 0)
            ky += ncell[1];
          else if (cy >= ncell[1])
            ky -= ncell[1];

          for (int cx = ix - step; cx <= ix + step; cx++) {
            // min_ix = abs(cx-ix);
            // if(min_ix>0) min_ix--;
            // if(min_ix*min_ix+min_iy*min_iy+min_iz*min_iz>=step*step) break;
            kx = cx;
            if (cx < 0)
              kx += ncell[0];
            else if (cx >= ncell[0])
              kx -= ncell[0];
            cindex = ncell[0] * (ncell[1] * kz + ky) + kx;

            pt = cell_atoms[cindex];
            while (pt) {
              // not the same atom and not protons. Proton is not accounted for
              // VDW and SOLV
              if (pt != theAtom && !theAtom->isH() && !pt->isH()) {
                d2 = (pt->r - theAtom->r).mod2();
                d = sqrt(d2);
                tj = pt->FF_t;
                if (d2 < max_ir2) { // within the cutoff distance
                  // vdw
                  vec Rsm = theAtom->r - pt->r;
                  if (isIncluded(
                          theAtom,
                          pt)) { // vdw does not include the internal pairs
                    vdw.getVDWdata(ti, tj, B, A, eps, sig, cutoff, slope,
                                   intercept);
                    if (B && A) {
                      d6 = d2 * d2 * d2;
                      e = (B / d6 - A) / d6;
                      if (e > 0) {
                        // e = (sig-d)/sig*10.0;
                        if (d < cutoff) {
                          e = slope * d + intercept;
                          de = slope;
                        } else {
                          de = 6.0 * (A - 2.0 * B / d6) / (d6 * d);
                        }
                        E_VDW_rep += e;
                        for (int ich = 0; ich < nCHI_affected; ich++) {
                          DEV_VDW_rep[ich] +=
                              de *
                              ((*theRes->getChiUV(ich) ^
                                (theAtom->r - *theRes->getChiP2(ich))) *
                               Rsm) /
                              d;
                        }
                      } else {
                        E_VDW_attr += e;
                        de = 6.0 * (A - 2.0 * B / d6) / (d6 * d);
                        for (int ich = 0; ich < nCHI_affected; ich++) {
                          DEV_VDW_attr[ich] +=
                              de *
                              ((*theRes->getChiUV(ich) ^
                                (theAtom->r - *theRes->getChiP2(ich))) *
                               Rsm) /
                              d;
                        }
                      }
                    }

                    // solvation
                    int solv_cut = 0;
                    if (d < (EEF1_SOLV[ti][CHARMM_RADIUS] +
                             EEF1_SOLV[tj][CHARMM_RADIUS]) *
                                0.846) {
                      d = (EEF1_SOLV[ti][CHARMM_RADIUS] +
                           EEF1_SOLV[tj][CHARMM_RADIUS]) *
                          0.846;
                      d2 = d * d;
                      solv_cut = 1;
                    }
                    xij = (d - EEF1_SOLV[ti][CHARMM_RADIUS]) /
                          EEF1_SOLV[ti][EEF1_LAMBDA];
                    xji = (d - EEF1_SOLV[tj][CHARMM_RADIUS]) /
                          EEF1_SOLV[tj][EEF1_LAMBDA];
                    SOLVi -=
                        (EEF1_SOLV[ti][EEF1_DG_FREE] * EEF1_SOLV[tj][EEF1_VOL] *
                             exp(-xij * xij) / EEF1_SOLV[ti][EEF1_LAMBDA] +
                         EEF1_SOLV[tj][EEF1_DG_FREE] * EEF1_SOLV[ti][EEF1_VOL] *
                             exp(-xji * xji) / EEF1_SOLV[tj][EEF1_LAMBDA]) *
                        coef / d2;
                    if (!solv_cut) {
                      de = (EEF1_SOLV[ti][EEF1_DG_FREE] *
                                EEF1_SOLV[tj][EEF1_VOL] * exp(-xij * xij) *
                                (xij / EEF1_SOLV[ti][EEF1_LAMBDA] + 1.0 / d) /
                                EEF1_SOLV[ti][EEF1_LAMBDA] +
                            EEF1_SOLV[tj][EEF1_DG_FREE] *
                                EEF1_SOLV[ti][EEF1_VOL] * exp(-xji * xji) *
                                (xji / EEF1_SOLV[tj][EEF1_LAMBDA] + 1.0 / d) /
                                EEF1_SOLV[tj][EEF1_LAMBDA]) *
                           coef2 / d2;
                      for (int ich = 0; ich < nCHI_affected; ich++) {
                        DEV_SOLV[ich] +=
                            de *
                            ((*theRes->getChiUV(ich) ^
                              (theAtom->r - *theRes->getChiP2(ich))) *
                             Rsm) /
                            d;
                      }
                    }
                  }
                }
              }
              pt = pt->cell_next;
            }
          }
        }
      }
      E_SOLV += SOLVi;
    }
  }
}

double gen_grid::getE_DEV_RES_SC(int iRes, double &E_VDW_attr,
                                 double &E_VDW_rep, double &E_SOLV,
                                 double &E_SB, double *DEV_VDW_attr,
                                 double *DEV_VDW_rep, double *DEV_SOLV,
                                 double *DEV_SB) {
  E_VDW_attr = E_VDW_rep = E_SOLV = 0;
  double vdw_rep, vdw_attr, solv;
  residue *theRes = gp.getResidue(iRes);
  if (theRes->getType() == GLY)
    return 0; // GLY does not have sidechain(a H)
  atom *theAtom = NULL;
  static double coef = 2.0 / (4.0 * PI * sqrt(PI));
  static double coef2 = 1.0 / (PI * sqrt(PI));
  int ti, tj;
  int ix, iy, iz, cindex, kx, ky, kz;
  int min_ix, min_iy, min_iz;
  double d2, d6, B, A, e, d, sig, eps, cutoff, slope, intercept;
  double SOLVi = 0, xij, xji;
  atom *pt;

  int nCHI_affected;
  int nCHI = theRes->getNChi();
  double de;
  for (int ich = 0; ich < nCHI + 1; ich++) {
    DEV_VDW_attr[ich] = 0;
    DEV_VDW_rep[ich] = 0;
    DEV_SOLV[ich] = 0;
    DEV_SB[ich] = 0;
  }

  for (int i = 0; i < theRes->getNA(); i++) {
    theAtom = theRes->getAtom(i);
    nCHI_affected = theRes->getNChi_AFFECTED(i);
    if (theAtom->isSC()) { // sidechain atoms include CB
      ti = theAtom->FF_t;
      ix = static_cast<int>((theAtom->r.x - min.x) / cell_box.x);
      iy = static_cast<int>((theAtom->r.y - min.y) / cell_box.y);
      iz = static_cast<int>((theAtom->r.z - min.z) / cell_box.z);

      // SOLVi=EEF1_SOLV[ti][EEF1_DG_REF];
      SOLVi = 0;
      for (int cz = iz - step; cz <= iz + step; cz++) {
        // min_iz = abs(cz-iz);
        // if(min_iz>0) min_iz--;
        kz = cz;
        if (cz < 0)
          kz += ncell[2];
        else if (cz >= ncell[2])
          kz -= ncell[2];

        for (int cy = iy - step; cy <= iy + step; cy++) {
          // min_iy = abs(cy-iy);
          // if(min_iy>0) min_iy--;
          ky = cy;
          if (cy < 0)
            ky += ncell[1];
          else if (cy >= ncell[1])
            ky -= ncell[1];

          for (int cx = ix - step; cx <= ix + step; cx++) {
            // min_ix = abs(cx-ix);
            // if(min_ix>0) min_ix--;
            // if(min_ix*min_ix+min_iy*min_iy+min_iz*min_iz>=step*step) break;
            kx = cx;
            if (cx < 0)
              kx += ncell[0];
            else if (cx >= ncell[0])
              kx -= ncell[0];
            cindex = ncell[0] * (ncell[1] * kz + ky) + kx;

            pt = cell_atoms[cindex];
            while (pt) {
              // not the same atom and not protons. Proton is not accounted for
              // VDW and SOLV
              if (pt != theAtom && !theAtom->isH() && !pt->isH()) {
                d2 = (pt->r - theAtom->r).mod2();
                d = sqrt(d2);
                tj = pt->FF_t;
                if (d2 < max_ir2) { // within the cutoff distance
                  // vdw
                  vec Rsm = theAtom->r - pt->r;
                  if (isIncluded(
                          theAtom,
                          pt)) { // vdw does not include the internal pairs
                    vdw.getVDWdata(ti, tj, B, A, eps, sig, cutoff, slope,
                                   intercept);
                    if (B && A) {
                      d6 = d2 * d2 * d2;
                      e = (B / d6 - A) / d6;
                      if (e > 0) {
                        if (d < cutoff) {
                          e = slope * d + intercept;
                          de = slope;
                        } else {
                          de = 6.0 * (A - 2.0 * B / d6) / (d6 * d);
                        }
                        // e = (sig-d)/sig*10.0;
                        E_VDW_rep += e;
                        // de = -10.0/sig;
                        for (int ich = 0; ich < nCHI_affected; ich++) {
                          DEV_VDW_rep[ich] +=
                              de *
                              ((*theRes->getChiUV(ich) ^
                                (theAtom->r - *theRes->getChiP2(ich))) *
                               Rsm) /
                              d;
                        }
                      } else {
                        E_VDW_attr += e;
                        de = 6.0 * (A - 2.0 * B / d6) / (d6 * d);
                        for (int ich = 0; ich < nCHI_affected; ich++) {
                          DEV_VDW_attr[ich] +=
                              de *
                              ((*theRes->getChiUV(ich) ^
                                (theAtom->r - *theRes->getChiP2(ich))) *
                               Rsm) /
                              d;
                        }
                      }
                    }

                    if (theAtom->icharge && pt->icharge) {
                      if (d < sig) {
                        e = (theAtom->icharge * pt->icharge) / (sig * sig);
                        E_SB += e;
                      } else {
                        e = (theAtom->icharge * pt->icharge) / d2;
                        E_SB += e;
                        de = -2.0 * e / d;
                        for (int ich = 0; ich < nCHI_affected; ich++) {
                          DEV_SB[ich] +=
                              de *
                              ((*theRes->getChiUV(ich) ^
                                (theAtom->r - *theRes->getChiP2(ich))) *
                               Rsm) /
                              d;
                        }
                      }
                    }

                    // solvation
                    int solv_cut = 0;
                    if (d < (EEF1_SOLV[ti][CHARMM_RADIUS] +
                             EEF1_SOLV[tj][CHARMM_RADIUS]) *
                                0.846) {
                      d = (EEF1_SOLV[ti][CHARMM_RADIUS] +
                           EEF1_SOLV[tj][CHARMM_RADIUS]) *
                          0.846;
                      d2 = d * d;
                      solv_cut = 1;
                    }
                    xij = (d - EEF1_SOLV[ti][CHARMM_RADIUS]) /
                          EEF1_SOLV[ti][EEF1_LAMBDA];
                    xji = (d - EEF1_SOLV[tj][CHARMM_RADIUS]) /
                          EEF1_SOLV[tj][EEF1_LAMBDA];
                    SOLVi -=
                        (EEF1_SOLV[ti][EEF1_DG_FREE] * EEF1_SOLV[tj][EEF1_VOL] *
                             exp(-xij * xij) / EEF1_SOLV[ti][EEF1_LAMBDA] +
                         EEF1_SOLV[tj][EEF1_DG_FREE] * EEF1_SOLV[ti][EEF1_VOL] *
                             exp(-xji * xji) / EEF1_SOLV[tj][EEF1_LAMBDA]) *
                        coef / d2;
                    if (!solv_cut) {
                      de = (EEF1_SOLV[ti][EEF1_DG_FREE] *
                                EEF1_SOLV[tj][EEF1_VOL] * exp(-xij * xij) *
                                (xij / EEF1_SOLV[ti][EEF1_LAMBDA] + 1.0 / d) /
                                EEF1_SOLV[ti][EEF1_LAMBDA] +
                            EEF1_SOLV[tj][EEF1_DG_FREE] *
                                EEF1_SOLV[ti][EEF1_VOL] * exp(-xji * xji) *
                                (xji / EEF1_SOLV[tj][EEF1_LAMBDA] + 1.0 / d) /
                                EEF1_SOLV[tj][EEF1_LAMBDA]) *
                           coef2 / d2;
                      for (int ich = 0; ich < nCHI_affected; ich++) {
                        DEV_SOLV[ich] +=
                            de *
                            ((*theRes->getChiUV(ich) ^
                              (theAtom->r - *theRes->getChiP2(ich))) *
                             Rsm) /
                            d;
                      }
                    }
                  }
                }
              }
              pt = pt->cell_next;
            }
          }
        }
      }
      E_SOLV += SOLVi;
    }
  }
}

/***************************************************************
  obtain the VDW and SOLV for GLY and PRO with respect to a
  reference aa,

  N(NZ with on proton) --- CA(CA1 with one proton) --- CO

  for GLY, CA1--> CA2;
  for PRO, NZ --> NZNP;
 ****************************************************************/
void gen_grid::getEVS_SPECIAL_RES(int iRes, double &E_VDW_attr,
                                  double &E_VDW_rep, double &E_SOLV,
                                  double *VDW_A, double *VDW_R, double *SOLV) {
  E_VDW_attr = E_VDW_rep = E_SOLV = 0;
  residue *theRes = gp.getResidue(iRes);
  int theType = theRes->getType();

  static double coef = 2.0 / (4.0 * PI * sqrt(PI));

  if (theType == GLY) { // compare CA
    atom *theAtom = theRes->getCA();
    int ti, t_ref, tj;
    ti = theAtom->FF_t;
    t_ref = gp.getGResidue(iRes)->get_residue(ALA)->getCA()->FF_t;

    if (ti == t_ref)
      return;

    int ix, iy, iz, cindex, kx, ky, kz;
    ;
    ix = static_cast<int>((theAtom->r.x - min.x) / cell_box.x);
    iy = static_cast<int>((theAtom->r.y - min.y) / cell_box.y);
    iz = static_cast<int>((theAtom->r.z - min.z) / cell_box.z);

    // double SOLVi = EEF1_SOLV[ti][EEF1_DG_REF] -
    // EEF1_SOLV[t_ref][EEF1_DG_REF];
    double SOLVi = 0;
    // if(SOLV){
    //   SOLV[theAtom->resID] += SOLVi;
    // }
    double d2, d6, B, A, e, d, sig, eps, cutoff, slope, intercept;
    double xij, xji;
    double solv_ij;
    atom *pt;

    for (int cz = iz - step; cz <= iz + step; cz++) {
      kz = cz;
      if (cz < 0)
        kz += ncell[2];
      else if (cz >= ncell[2])
        kz -= ncell[2];

      for (int cy = iy - step; cy <= iy + step; cy++) {
        ky = cy;
        if (cy < 0)
          ky += ncell[1];
        else if (cy >= ncell[1])
          ky -= ncell[1];

        for (int cx = ix - step; cx <= ix + step; cx++) {
          kx = cx;
          if (cx < 0)
            kx += ncell[0];
          else if (cx >= ncell[0])
            kx -= ncell[0];

          cindex = ncell[0] * (ncell[1] * kz + ky) + kx;

          pt = cell_atoms[cindex];
          while (pt) {
            if (pt != theAtom && !theAtom->isH() &&
                !pt->isH()) { // not the same atom and not protons. Proton is
                              // not accounted for VDW and SOLV
              d2 = (pt->r - theAtom->r).mod2();
              d = sqrt(d2);
              tj = pt->FF_t;
              if (d2 < max_ir2) { // within the cutoff distance
                // vdw
                if (isIncluded(theAtom,
                               pt)) { // vdw does not include the internal pairs
                  vdw.getVDWdata(ti, tj, B, A, eps, sig, cutoff, slope,
                                 intercept);
                  if (B && A) {
                    d6 = d2 * d2 * d2;
                    e = (B / d6 - A) / d6;
                    if (e > 0) {
                      // e = (sig-d)/sig*10.0;
                      if (d < cutoff) {
                        e = slope * d + intercept;
                      }
                      E_VDW_rep += e;
                      if (VDW_R && pt->isSC()) {
                        VDW_R[pt->resID] += e;
                      }
                    } else {
                      E_VDW_attr += e;
                      if (VDW_A && pt->isSC()) {
                        VDW_A[pt->resID] += e;
                      }
                    }
                  }
                  // ref
                  vdw.getVDWdata(t_ref, tj, B, A, eps, sig, cutoff, slope,
                                 intercept);
                  if (B && A) {
                    d6 = d2 * d2 * d2;
                    e = (B / d6 - A) / d6;
                    if (e > 0) {
                      // e = (sig-d)/sig*10.0;
                      if (d < cutoff) {
                        e = slope * d + intercept;
                      }
                      E_VDW_rep -= e;
                      if (VDW_R && pt->isSC()) {
                        VDW_R[pt->resID] -= e;
                      }
                    } else {
                      E_VDW_attr -= e;
                      if (VDW_A && pt->isSC()) {
                        VDW_A[pt->resID] -= e;
                      }
                    }
                  }

                  // solvation
                  // if tow atoms are too close, asuming a cutoff energy
                  if (d < (EEF1_SOLV[ti][CHARMM_RADIUS] +
                           EEF1_SOLV[tj][CHARMM_RADIUS]) *
                              0.846) {
                    d = (EEF1_SOLV[ti][CHARMM_RADIUS] +
                         EEF1_SOLV[tj][CHARMM_RADIUS]) *
                        0.846;
                    d2 = d * d;
                  }
                  xij = (d - EEF1_SOLV[ti][CHARMM_RADIUS]) /
                        EEF1_SOLV[ti][EEF1_LAMBDA];
                  xji = (d - EEF1_SOLV[tj][CHARMM_RADIUS]) /
                        EEF1_SOLV[tj][EEF1_LAMBDA];
                  solv_ij =
                      (EEF1_SOLV[ti][EEF1_DG_FREE] * EEF1_SOLV[tj][EEF1_VOL] *
                           exp(-xij * xij) / EEF1_SOLV[ti][EEF1_LAMBDA] +
                       EEF1_SOLV[tj][EEF1_DG_FREE] * EEF1_SOLV[ti][EEF1_VOL] *
                           exp(-xji * xji) / EEF1_SOLV[tj][EEF1_LAMBDA]) *
                      coef / d2;
                  SOLVi -= solv_ij;
                  if (SOLV && pt->isSC()) {
                    SOLV[pt->resID] -= solv_ij;
                  }
                  // ref
                  if (d < (EEF1_SOLV[t_ref][CHARMM_RADIUS] +
                           EEF1_SOLV[tj][CHARMM_RADIUS]) *
                              0.846) {
                    d = (EEF1_SOLV[t_ref][CHARMM_RADIUS] +
                         EEF1_SOLV[tj][CHARMM_RADIUS]) *
                        0.846;
                    d2 = d * d;
                  }
                  xij = (d - EEF1_SOLV[t_ref][CHARMM_RADIUS]) /
                        EEF1_SOLV[t_ref][EEF1_LAMBDA];
                  xji = (d - EEF1_SOLV[tj][CHARMM_RADIUS]) /
                        EEF1_SOLV[tj][EEF1_LAMBDA];
                  solv_ij = (EEF1_SOLV[t_ref][EEF1_DG_FREE] *
                                 EEF1_SOLV[tj][EEF1_VOL] * exp(-xij * xij) /
                                 EEF1_SOLV[t_ref][EEF1_LAMBDA] +
                             EEF1_SOLV[tj][EEF1_DG_FREE] *
                                 EEF1_SOLV[t_ref][EEF1_VOL] * exp(-xji * xji) /
                                 EEF1_SOLV[tj][EEF1_LAMBDA]) *
                            coef / d2;
                  SOLVi += solv_ij;
                  if (SOLV && pt->isSC()) {
                    SOLV[pt->resID] += solv_ij;
                  }
                }
              }
            }
            pt = pt->cell_next;
          }
        }
      }
    }
    E_SOLV = SOLVi;
  } else if (theType == PRO) { // compare N
    atom *theAtom = theRes->getN();
    int ti, t_ref, tj;
    ti = theAtom->FF_t;
    t_ref = gp.getGResidue(iRes)->get_residue(ALA)->getN()->FF_t;

    if (ti == t_ref)
      return;

    int ix, iy, iz, cindex, kx, ky, kz;
    ;
    ix = static_cast<int>((theAtom->r.x - min.x) / cell_box.x);
    iy = static_cast<int>((theAtom->r.y - min.y) / cell_box.y);
    iz = static_cast<int>((theAtom->r.z - min.z) / cell_box.z);

    double SOLVi = 0;
    // double SOLVi = EEF1_SOLV[ti][EEF1_DG_REF] -
    // EEF1_SOLV[t_ref][EEF1_DG_REF]; if(SOLV){
    //   SOLV[theAtom->resID] += SOLVi;
    // }
    double d2, d6, B, A, e, d, sig, eps, cutoff, slope, intercept;
    double xij, xji;
    double solv_ij;
    atom *pt;

    for (int cz = iz - step; cz <= iz + step; cz++) {
      kz = cz;
      if (cz < 0)
        kz += ncell[2];
      else if (cz >= ncell[2])
        kz -= ncell[2];

      for (int cy = iy - step; cy <= iy + step; cy++) {
        ky = cy;
        if (cy < 0)
          ky += ncell[1];
        else if (cy >= ncell[1])
          ky -= ncell[1];

        for (int cx = ix - step; cx <= ix + step; cx++) {
          kx = cx;
          if (cx < 0)
            kx += ncell[0];
          else if (cx >= ncell[0])
            kx -= ncell[0];

          cindex = ncell[0] * (ncell[1] * kz + ky) + kx;

          pt = cell_atoms[cindex];
          while (pt) {
            if (pt != theAtom && !theAtom->isH() &&
                !pt->isH()) { // not the same atom and not protons. Proton is
                              // not accounted for VDW and SOLV
              d2 = (pt->r - theAtom->r).mod2();
              d = sqrt(d2);
              tj = pt->FF_t;
              if (d2 < max_ir2) { // within the cutoff distance
                // vdw
                if (isIncluded(theAtom,
                               pt)) { // vdw does not include the internal pairs
                  vdw.getVDWdata(ti, tj, B, A, eps, sig, cutoff, slope,
                                 intercept);
                  if (B && A) {
                    d6 = d2 * d2 * d2;
                    e = (B / d6 - A) / d6;
                    if (e > 0) {
                      // e = (sig-d)/sig*10.0;
                      if (d < cutoff) {
                        e = slope * d + intercept;
                      }
                      E_VDW_rep += e;
                      if (VDW_R && pt->isSC()) {
                        VDW_R[pt->resID] += e;
                      }
                    } else {
                      E_VDW_attr += e;
                      if (VDW_A && pt->isSC()) {
                        VDW_A[pt->resID] += e;
                      }
                    }
                  }
                  // ref
                  vdw.getVDWdata(t_ref, tj, B, A, eps, sig, cutoff, slope,
                                 intercept);
                  if (B && A) {
                    d6 = d2 * d2 * d2;
                    e = (B / d6 - A) / d6;
                    if (e > 0) {
                      // e = (sig-d)/sig*10.0;
                      if (d < cutoff) {
                        e = slope * d + intercept;
                      }
                      E_VDW_rep -= e;
                      if (VDW_R && pt->isSC()) {
                        VDW_R[pt->resID] -= e;
                      }
                    } else {
                      E_VDW_attr -= e;
                      if (VDW_A && pt->isSC()) {
                        VDW_A[pt->resID] -= e;
                      }
                    }
                  }

                  // solvation
                  // if tow atoms are too close, asuming a cutoff energy
                  if (d < (EEF1_SOLV[ti][CHARMM_RADIUS] +
                           EEF1_SOLV[tj][CHARMM_RADIUS]) *
                              0.846) {
                    d = (EEF1_SOLV[ti][CHARMM_RADIUS] +
                         EEF1_SOLV[tj][CHARMM_RADIUS]) *
                        0.846;
                    d2 = d * d;
                  }
                  xij = (d - EEF1_SOLV[ti][CHARMM_RADIUS]) /
                        EEF1_SOLV[ti][EEF1_LAMBDA];
                  xji = (d - EEF1_SOLV[tj][CHARMM_RADIUS]) /
                        EEF1_SOLV[tj][EEF1_LAMBDA];
                  solv_ij =
                      (EEF1_SOLV[ti][EEF1_DG_FREE] * EEF1_SOLV[tj][EEF1_VOL] *
                           exp(-xij * xij) / EEF1_SOLV[ti][EEF1_LAMBDA] +
                       EEF1_SOLV[tj][EEF1_DG_FREE] * EEF1_SOLV[ti][EEF1_VOL] *
                           exp(-xji * xji) / EEF1_SOLV[tj][EEF1_LAMBDA]) *
                      coef / d2;
                  SOLVi -= solv_ij;
                  if (SOLV && pt->isSC()) {
                    SOLV[pt->resID] -= solv_ij;
                  }
                  // ref
                  if (d < (EEF1_SOLV[t_ref][CHARMM_RADIUS] +
                           EEF1_SOLV[tj][CHARMM_RADIUS]) *
                              0.846) {
                    d = (EEF1_SOLV[t_ref][CHARMM_RADIUS] +
                         EEF1_SOLV[tj][CHARMM_RADIUS]) *
                        0.846;
                    d2 = d * d;
                  }
                  xij = (d - EEF1_SOLV[t_ref][CHARMM_RADIUS]) /
                        EEF1_SOLV[t_ref][EEF1_LAMBDA];
                  xji = (d - EEF1_SOLV[tj][CHARMM_RADIUS]) /
                        EEF1_SOLV[tj][EEF1_LAMBDA];
                  solv_ij = (EEF1_SOLV[t_ref][EEF1_DG_FREE] *
                                 EEF1_SOLV[tj][EEF1_VOL] * exp(-xij * xij) /
                                 EEF1_SOLV[t_ref][EEF1_LAMBDA] +
                             EEF1_SOLV[tj][EEF1_DG_FREE] *
                                 EEF1_SOLV[t_ref][EEF1_VOL] * exp(-xji * xji) /
                                 EEF1_SOLV[tj][EEF1_LAMBDA]) *
                            coef / d2;
                  SOLVi += solv_ij;
                  if (SOLV && pt->isSC()) {
                    SOLV[pt->resID] += solv_ij;
                  }
                }
              }
            }
            pt = pt->cell_next;
          }
        }
      }
    }
    E_SOLV = SOLVi;
  } else
    return;
}

double gen_grid::getE_DEV(double &E_VDWA, double &E_VDWR, double &E_SOLV,
                          double *DVDWA, double *DVDWR, double *DSOLV) {
  E_VDWA = E_VDWR = E_SOLV = 0;
  residue *theRes = NULL;
  atom *theAtom = NULL;
  static double coef = 2.0 / (4.0 * PI * sqrt(PI));
  static double coef2 = 1.0 / (PI * sqrt(PI));
  int ti, tj;
  int ix, iy, iz, cindex, kx, ky, kz;
  int min_ix, min_iy, min_iz;
  double d2, d6, B, A, e, d, sig, eps, cutoff, slope, intercept;
  double SOLVi = 0, xij, xji, de;
  atom *pt;

  for (int i = 0; i < gp.size(); i++) {
    DVDWA[2 * i] = DVDWA[2 * i + 1] = 0;
    DVDWR[2 * i] = DVDWR[2 * i + 1] = 0;
    DSOLV[2 * i] = DSOLV[2 * i + 1] = 0;
  }

  gp.updateFY_UV();
  for (int ires = 0; ires < gp.size(); ires++) {
    theRes = gp.getResidue(ires);
    for (int iatom = 0; iatom < theRes->getNA(); iatom++) {
      theAtom = theRes->getAtom(iatom);
      ti = theAtom->FF_t;
      ix = static_cast<int>((theAtom->r.x - min.x) / cell_box.x);
      iy = static_cast<int>((theAtom->r.y - min.y) / cell_box.y);
      iz = static_cast<int>((theAtom->r.z - min.z) / cell_box.z);
      SOLVi = 0;
      for (int cz = iz - step; cz <= iz + step; cz++) {
        kz = cz;
        if (cz < 0)
          kz += ncell[2];
        else if (cz >= ncell[2])
          kz -= ncell[2];

        for (int cy = iy - step; cy <= iy + step; cy++) {
          ky = cy;
          if (cy < 0)
            ky += ncell[1];
          else if (cy >= ncell[1])
            ky -= ncell[1];

          for (int cx = ix - step; cx <= ix + step; cx++) {
            kx = cx;
            if (cx < 0)
              kx += ncell[0];
            else if (cx >= ncell[0])
              kx -= ncell[0];

            cindex = ncell[0] * (ncell[1] * kz + ky) + kx;
            pt = cell_atoms[cindex];
            while (pt) {
              int jres = pt->resID;
              vec Rsm = pt->r - theAtom->r;
              if (jres > ires + 1 && !pt->isH()) {
                d2 = (pt->r - theAtom->r).mod2();
                d = sqrt(d2);
                tj = pt->FF_t;
                if (d2 < max_ir2) { // within the cutoff distance
                  vdw.getVDWdata(ti, tj, B, A, eps, sig, cutoff, slope,
                                 intercept);
                  if (B && A) {
                    d6 = d2 * d2 * d2;
                    e = (B / d6 - A) / d6;
                    if (e > 0) {
                      // e = (sig-d)/sig*10.0;
                      if (d < cutoff) {
                        e = slope * d + intercept;
                        de = slope;
                      } else {
                        de = 6.0 * (A - 2.0 * B / d6) / (d6 * d);
                      }
                      E_VDWR += e;
                      // de = -10.0/sig;
                      // to calculate derative
                      // the psi of ires,
                      if (theAtom->PDBtype != _CA_ && theAtom->PDBtype != _C_ &&
                          theAtom->PDBtype != _O_) {
                        DVDWR[2 * ires + 1] += de *
                                               ((gp.getPsi_UV(ires) ^
                                                 (pt->r - theRes->getC()->r)) *
                                                Rsm) /
                                               d;
                        // cout <<
                        // de*((gp.getPsi_UV(ires)^(pt->r-theRes->getC()->r))*Rsm)/d
                        // << endl;
                      }
                      //(phi,psi) from ires+1 to jres-1;
                      for (int ii = ires + 1; ii <= jres - 1; ii++) {
                        DVDWR[2 * ii] +=
                            de *
                            ((gp.getPhi_UV(ii) ^
                              (pt->r - gp.getResidue(ii)->getCA()->r)) *
                             Rsm) /
                            d;
                        DVDWR[2 * ii + 1] +=
                            de *
                            ((gp.getPsi_UV(ii) ^
                              (pt->r - gp.getResidue(ii)->getC()->r)) *
                             Rsm) /
                            d;
                      }
                      // the phi of jres,
                      if (pt->PDBtype != _N_ && pt->PDBtype != _CA_) {
                        DVDWR[2 * jres] +=
                            de *
                            ((gp.getPhi_UV(jres) ^
                              (pt->r - gp.getResidue(jres)->getCA()->r)) *
                             Rsm) /
                            d;
                      }
                      if (pt->PDBtype == _O_) {
                        DVDWR[2 * jres + 1] +=
                            de *
                            ((gp.getPsi_UV(jres) ^
                              (pt->r - gp.getResidue(jres)->getC()->r)) *
                             Rsm) /
                            d;
                      }
                    } else {
                      E_VDWA += e;
                      de = 6.0 * (A - 2.0 * B / d6) / (d6 * d);
                      // to calculate derative
                      // the psi of ires,
                      if (theAtom->PDBtype != _CA_ && theAtom->PDBtype != _C_ &&
                          theAtom->PDBtype != _O_) {
                        DVDWA[2 * ires + 1] += de *
                                               ((gp.getPsi_UV(ires) ^
                                                 (pt->r - theRes->getC()->r)) *
                                                Rsm) /
                                               d;
                      }
                      //(phi,psi) from ires+1 to jres-1;
                      for (int ii = ires + 1; ii <= jres - 1; ii++) {
                        DVDWA[2 * ii] +=
                            de *
                            ((gp.getPhi_UV(ii) ^
                              (pt->r - gp.getResidue(ii)->getCA()->r)) *
                             Rsm) /
                            d;
                        DVDWA[2 * ii + 1] +=
                            de *
                            ((gp.getPsi_UV(ii) ^
                              (pt->r - gp.getResidue(ii)->getC()->r)) *
                             Rsm) /
                            d;
                      }
                      // the phi of jres,
                      if (pt->PDBtype != _N_ && pt->PDBtype != _CA_) {
                        DVDWA[2 * jres] +=
                            de *
                            ((gp.getPhi_UV(jres) ^
                              (pt->r - gp.getResidue(jres)->getCA()->r)) *
                             Rsm) /
                            d;
                      }
                      if (pt->PDBtype == _O_) {
                        DVDWA[2 * jres + 1] +=
                            de *
                            ((gp.getPsi_UV(jres) ^
                              (pt->r - gp.getResidue(jres)->getC()->r)) *
                             Rsm) /
                            d;
                      }
                    }
                  }
                  // solvation
                  int solv_cut = 0;
                  if (d < (EEF1_SOLV[ti][CHARMM_RADIUS] +
                           EEF1_SOLV[tj][CHARMM_RADIUS]) *
                              0.846) {
                    d = (EEF1_SOLV[ti][CHARMM_RADIUS] +
                         EEF1_SOLV[tj][CHARMM_RADIUS]) *
                        0.846;
                    d2 = d * d;
                    solv_cut = 1;
                  }
                  xij = (d - EEF1_SOLV[ti][CHARMM_RADIUS]) /
                        EEF1_SOLV[ti][EEF1_LAMBDA];
                  xji = (d - EEF1_SOLV[tj][CHARMM_RADIUS]) /
                        EEF1_SOLV[tj][EEF1_LAMBDA];
                  SOLVi -=
                      (EEF1_SOLV[ti][EEF1_DG_FREE] * EEF1_SOLV[tj][EEF1_VOL] *
                           exp(-xij * xij) / EEF1_SOLV[ti][EEF1_LAMBDA] +
                       EEF1_SOLV[tj][EEF1_DG_FREE] * EEF1_SOLV[ti][EEF1_VOL] *
                           exp(-xji * xji) / EEF1_SOLV[tj][EEF1_LAMBDA]) *
                      coef / d2;
                  if (!solv_cut) {
                    de = (EEF1_SOLV[ti][EEF1_DG_FREE] *
                              EEF1_SOLV[tj][EEF1_VOL] * exp(-xij * xij) *
                              (xij / EEF1_SOLV[ti][EEF1_LAMBDA] + 1.0 / d) /
                              EEF1_SOLV[ti][EEF1_LAMBDA] +
                          EEF1_SOLV[tj][EEF1_DG_FREE] *
                              EEF1_SOLV[ti][EEF1_VOL] * exp(-xji * xji) *
                              (xji / EEF1_SOLV[tj][EEF1_LAMBDA] + 1.0 / d) /
                              EEF1_SOLV[tj][EEF1_LAMBDA]) *
                         coef2 / d2;

                    // to calculate derative
                    // the psi of ires,
                    if (theAtom->PDBtype != _CA_ && theAtom->PDBtype != _C_ &&
                        theAtom->PDBtype != _O_) {
                      DSOLV[2 * ires + 1] +=
                          de *
                          ((gp.getPsi_UV(ires) ^ (pt->r - theRes->getC()->r)) *
                           Rsm) /
                          d;
                    }
                    //(phi,psi) from ires+1 to jres-1;
                    for (int ii = ires + 1; ii <= jres - 1; ii++) {
                      DSOLV[2 * ii] +=
                          de *
                          ((gp.getPhi_UV(ii) ^
                            (pt->r - gp.getResidue(ii)->getCA()->r)) *
                           Rsm) /
                          d;
                      DSOLV[2 * ii + 1] +=
                          de *
                          ((gp.getPsi_UV(ii) ^
                            (pt->r - gp.getResidue(ii)->getC()->r)) *
                           Rsm) /
                          d;
                    }
                    // the phi of jres,
                    if (pt->PDBtype != _N_ && pt->PDBtype != _CA_) {
                      DSOLV[2 * jres] +=
                          de *
                          ((gp.getPhi_UV(jres) ^
                            (pt->r - gp.getResidue(jres)->getCA()->r)) *
                           Rsm) /
                          d;
                    }
                    if (pt->PDBtype == _O_) {
                      DSOLV[2 * jres + 1] +=
                          de *
                          ((gp.getPsi_UV(jres) ^
                            (pt->r - gp.getResidue(jres)->getC()->r)) *
                           Rsm) /
                          d;
                    }
                  }
                }
              } else if (jres == ires + 1 && !pt->isH()) {
                // at first do not account for the short range bb interactions
                // since, the short range backbone interaction determines the
                // phi,psi, we assume that the phi,psi angles do not have large
                // deviations.
                if ((theAtom->isSC() && (pt->isSC() || pt->PDBtype == _C_ ||
                                         pt->PDBtype == _O_)) ||
                    (pt->isSC() &&
                     (theAtom->isSC() || theAtom->PDBtype == _N_))) {
                  d2 = (pt->r - theAtom->r).mod2();
                  d = sqrt(d2);
                  tj = pt->FF_t;
                  if (d2 < max_ir2) { // within the cutoff distance
                    vdw.getVDWdata(ti, tj, B, A, eps, sig, cutoff, slope,
                                   intercept);
                    if (B && A) {
                      d6 = d2 * d2 * d2;
                      e = (B / d6 - A) / d6;
                      if (e > 0) {
                        if (d < cutoff) {
                          e = slope * d + intercept;
                          de = slope;
                        } else {
                          de = 6.0 * (A - 2.0 * B / d6) / (d6 * d);
                        }
                        // e = (sig-d)/sig*10.0;
                        E_VDWR += e;
                        // de = -10.0/sig;
                        // to calculate derative
                        // the psi of ires,
                        if (theAtom->PDBtype != _CA_ &&
                            theAtom->PDBtype != _C_ &&
                            theAtom->PDBtype != _O_) {
                          DVDWR[2 * ires + 1] +=
                              de *
                              ((gp.getPsi_UV(ires) ^
                                (pt->r - theRes->getC()->r)) *
                               Rsm) /
                              d;
                        }
                        // the phi of jres,
                        if (pt->PDBtype != _N_ && pt->PDBtype != _CA_) {
                          DVDWR[2 * jres] +=
                              de *
                              ((gp.getPhi_UV(jres) ^
                                (pt->r - gp.getResidue(jres)->getCA()->r)) *
                               Rsm) /
                              d;
                        }
                        if (pt->PDBtype == _O_) {
                          DVDWR[2 * jres + 1] +=
                              de *
                              ((gp.getPsi_UV(jres) ^
                                (pt->r - gp.getResidue(jres)->getC()->r)) *
                               Rsm) /
                              d;
                        }
                      } else {
                        E_VDWA += e;
                        de = 6.0 * (A - 2.0 * B / d6) / (d6 * d);
                        // to calculate derative
                        // the psi of ires,
                        if (theAtom->PDBtype != _CA_ &&
                            theAtom->PDBtype != _C_ &&
                            theAtom->PDBtype != _O_) {
                          DVDWA[2 * ires + 1] +=
                              de *
                              ((gp.getPsi_UV(ires) ^
                                (pt->r - theRes->getC()->r)) *
                               Rsm) /
                              d;
                        }
                        // the phi of jres,
                        if (pt->PDBtype != _N_ && pt->PDBtype != _CA_) {
                          DVDWA[2 * jres] +=
                              de *
                              ((gp.getPhi_UV(jres) ^
                                (pt->r - gp.getResidue(jres)->getCA()->r)) *
                               Rsm) /
                              d;
                        }
                        if (pt->PDBtype == _O_) {
                          DVDWA[2 * jres + 1] +=
                              de *
                              ((gp.getPsi_UV(jres) ^
                                (pt->r - gp.getResidue(jres)->getC()->r)) *
                               Rsm) /
                              d;
                        }
                      }
                    }
                    // solvation
                    int solv_cut = 0;
                    if (d < (EEF1_SOLV[ti][CHARMM_RADIUS] +
                             EEF1_SOLV[tj][CHARMM_RADIUS]) *
                                0.846) {
                      d = (EEF1_SOLV[ti][CHARMM_RADIUS] +
                           EEF1_SOLV[tj][CHARMM_RADIUS]) *
                          0.846;
                      d2 = d * d;
                      solv_cut = 1;
                    }
                    xij = (d - EEF1_SOLV[ti][CHARMM_RADIUS]) /
                          EEF1_SOLV[ti][EEF1_LAMBDA];
                    xji = (d - EEF1_SOLV[tj][CHARMM_RADIUS]) /
                          EEF1_SOLV[tj][EEF1_LAMBDA];
                    SOLVi -=
                        (EEF1_SOLV[ti][EEF1_DG_FREE] * EEF1_SOLV[tj][EEF1_VOL] *
                             exp(-xij * xij) / EEF1_SOLV[ti][EEF1_LAMBDA] +
                         EEF1_SOLV[tj][EEF1_DG_FREE] * EEF1_SOLV[ti][EEF1_VOL] *
                             exp(-xji * xji) / EEF1_SOLV[tj][EEF1_LAMBDA]) *
                        coef / d2;
                    if (!solv_cut) {
                      de = (EEF1_SOLV[ti][EEF1_DG_FREE] *
                                EEF1_SOLV[tj][EEF1_VOL] * exp(-xij * xij) *
                                (xij / EEF1_SOLV[ti][EEF1_LAMBDA] + 1.0 / d) /
                                EEF1_SOLV[ti][EEF1_LAMBDA] +
                            EEF1_SOLV[tj][EEF1_DG_FREE] *
                                EEF1_SOLV[ti][EEF1_VOL] * exp(-xji * xji) *
                                (xji / EEF1_SOLV[tj][EEF1_LAMBDA] + 1.0 / d) /
                                EEF1_SOLV[tj][EEF1_LAMBDA]) *
                           coef2 / d2;

                      // to calculate derative
                      // the psi of ires,
                      if (theAtom->PDBtype != _CA_ && theAtom->PDBtype != _C_ &&
                          theAtom->PDBtype != _O_) {
                        DSOLV[2 * ires + 1] += de *
                                               ((gp.getPsi_UV(ires) ^
                                                 (pt->r - theRes->getC()->r)) *
                                                Rsm) /
                                               d;
                      }
                      // the phi of jres,
                      if (pt->PDBtype != _N_ && pt->PDBtype != _CA_) {
                        DSOLV[2 * jres] +=
                            de *
                            ((gp.getPhi_UV(jres) ^
                              (pt->r - gp.getResidue(jres)->getCA()->r)) *
                             Rsm) /
                            d;
                      }
                      if (pt->PDBtype == _O_) {
                        DSOLV[2 * jres + 1] +=
                            de *
                            ((gp.getPsi_UV(jres) ^
                              (pt->r - gp.getResidue(jres)->getC()->r)) *
                             Rsm) /
                            d;
                      }
                    }
                  }
                }
              }
              pt = pt->cell_next;
            }
          } // x
        }   // y
      }     // z
      E_SOLV += SOLVi;
    }
  }

  // cout << E_VDWA << " "
  //      << E_VDWR << " "
  //      << E_SOLV << endl;
  // for(int i=0; i<gp.size(); i++){
  //   cout << DVDWA[2*i]   << " " << DVDWR[2*i]   << " " << DSOLV[2*i] << endl;
  //   cout << DVDWA[2*i+1] << " " << DVDWR[2*i+1] << " " << DSOLV[2*i+1] <<
  //   endl;
  // }
  // cout << endl;
}

/*************************************************************
  using the trick of calculating the derivates by Go and Abe.
 **************************************************************/
double gen_grid::getE_DEV_GOTRICK(double &E_VDWA, double &E_VDWR,
                                  double &E_SOLV, double *DVDWA, double *DVDWR,
                                  double *DSOLV) {
  E_VDWA = E_VDWR = E_SOLV = 0;
  residue *theRes = NULL;
  atom *theAtom = NULL;
  static double coef = 2.0 / (4.0 * PI * sqrt(PI));
  static double coef2 = 1.0 / (PI * sqrt(PI));
  int ti, tj;
  int ix, iy, iz, cindex, kx, ky, kz;
  int min_ix, min_iy, min_iz;
  double d2, d6, B, A, e, d, sig, eps, cutoff, slope, intercept;
  double SOLVi = 0, xij, xji, de;
  atom *pt;

  for (int i = 0; i < gp.size(); i++) {
    DVDWA[2 * i] = DVDWA[2 * i + 1] = 0;
    DVDWR[2 * i] = DVDWR[2 * i + 1] = 0;
    DSOLV[2 * i] = DSOLV[2 * i + 1] = 0;
  }
  gp.updateFY_UV();

  vec *F_VDWA = new vec[2 * gp.size()];
  vec *G_VDWA = new vec[2 * gp.size()];
  vec *F_VDWR = new vec[2 * gp.size()];
  vec *G_VDWR = new vec[2 * gp.size()];
  vec *F_SOLV = new vec[2 * gp.size()];
  vec *G_SOLV = new vec[2 * gp.size()];
  vector<atom *> tmp_va;
  for (int ires = gp.size() - 1; ires >= 0; ires--) {
    // psi: Oi, Ni+1, Hi+1, Cai+1
    F_VDWA[2 * ires + 1] = vec(0, 0, 0);
    G_VDWA[2 * ires + 1] = vec(0, 0, 0);
    F_VDWR[2 * ires + 1] = vec(0, 0, 0);
    G_VDWR[2 * ires + 1] = vec(0, 0, 0);
    F_SOLV[2 * ires + 1] = vec(0, 0, 0);
    G_SOLV[2 * ires + 1] = vec(0, 0, 0);
    tmp_va.clear();
    theRes = gp.getResidue(ires);
    tmp_va.push_back(theRes->getO());
    if (ires < gp.size() - 1) {
      theRes = gp.getResidue(ires + 1);
      tmp_va.push_back(theRes->getN());
      tmp_va.push_back(theRes->getCA());
    }
    for (int iatom = 0; iatom < tmp_va.size(); iatom++) {
      theAtom = tmp_va[iatom];
      int id = theAtom->resID;
      ti = theAtom->FF_t;
      ix = static_cast<int>((theAtom->r.x - min.x) / cell_box.x);
      iy = static_cast<int>((theAtom->r.y - min.y) / cell_box.y);
      iz = static_cast<int>((theAtom->r.z - min.z) / cell_box.z);
      SOLVi = 0;
      for (int cz = iz - step; cz <= iz + step; cz++) {
        kz = cz;
        if (cz < 0)
          kz += ncell[2];
        else if (cz >= ncell[2])
          kz -= ncell[2];

        for (int cy = iy - step; cy <= iy + step; cy++) {
          ky = cy;
          if (cy < 0)
            ky += ncell[1];
          else if (cy >= ncell[1])
            ky -= ncell[1];

          for (int cx = ix - step; cx <= ix + step; cx++) {
            kx = cx;
            if (cx < 0)
              kx += ncell[0];
            else if (cx >= ncell[0])
              kx -= ncell[0];

            cindex = ncell[0] * (ncell[1] * kz + ky) + kx;
            pt = cell_atoms[cindex];
            while (pt) {
              int jres = pt->resID;
              if (!pt->isH()) {
                if ((abs(id - jres) > 1) ||
                    (jres == id + 1 &&
                     ((theAtom->isSC() && (pt->isSC() || pt->PDBtype == _C_ ||
                                           pt->PDBtype == _O_)) ||
                      (pt->isSC() && (theAtom->PDBtype == _N_)))) ||
                    (jres == id - 1 &&
                     ((theAtom->isSC() && (pt->isSC() || pt->PDBtype == _N_)) ||
                      (pt->isSC() && (theAtom->PDBtype == _C_ ||
                                      theAtom->PDBtype == _O_))))) {
                  d2 = (pt->r - theAtom->r).mod2();
                  d = sqrt(d2);
                  tj = pt->FF_t;
                  if (d2 < max_ir2) { // within the cutoff distance
                    // vdw
                    vdw.getVDWdata(ti, tj, B, A, eps, sig, cutoff, slope,
                                   intercept);
                    if (B && A) {
                      d6 = d2 * d2 * d2;
                      e = (B / d6 - A) / d6;
                      if (e > 0) {
                        // e = (sig-d)/sig*10.0;
                        if (d < cutoff) {
                          e = slope * d + intercept;
                          de = slope;
                        } else {
                          de = 6.0 * (A - 2.0 * B / d6) / (d6 * d);
                        }
                        E_VDWR += e;
                        // de = -10.0/sig;
                        F_VDWR[2 * ires + 1] += de / d * (theAtom->r ^ pt->r);
                        G_VDWR[2 * ires + 1] += de / d * (theAtom->r - pt->r);
                      } else {
                        E_VDWA += e;
                        de = 6.0 * (A - 2.0 * B / d6) / (d6 * d);
                        F_VDWA[2 * ires + 1] += de / d * (theAtom->r ^ pt->r);
                        G_VDWA[2 * ires + 1] += de / d * (theAtom->r - pt->r);
                      }
                    }
                    // solv
                    int solv_cut = 0;
                    if (d < (EEF1_SOLV[ti][CHARMM_RADIUS] +
                             EEF1_SOLV[tj][CHARMM_RADIUS]) *
                                0.846) {
                      d = (EEF1_SOLV[ti][CHARMM_RADIUS] +
                           EEF1_SOLV[tj][CHARMM_RADIUS]) *
                          0.846;
                      d2 = d * d;
                      solv_cut = 1;
                    }
                    xij = (d - EEF1_SOLV[ti][CHARMM_RADIUS]) /
                          EEF1_SOLV[ti][EEF1_LAMBDA];
                    xji = (d - EEF1_SOLV[tj][CHARMM_RADIUS]) /
                          EEF1_SOLV[tj][EEF1_LAMBDA];
                    SOLVi -=
                        (EEF1_SOLV[ti][EEF1_DG_FREE] * EEF1_SOLV[tj][EEF1_VOL] *
                             exp(-xij * xij) / EEF1_SOLV[ti][EEF1_LAMBDA] +
                         EEF1_SOLV[tj][EEF1_DG_FREE] * EEF1_SOLV[ti][EEF1_VOL] *
                             exp(-xji * xji) / EEF1_SOLV[tj][EEF1_LAMBDA]) *
                        coef / d2;
                    if (!solv_cut) {
                      de = (EEF1_SOLV[ti][EEF1_DG_FREE] *
                                EEF1_SOLV[tj][EEF1_VOL] * exp(-xij * xij) *
                                (xij / EEF1_SOLV[ti][EEF1_LAMBDA] + 1.0 / d) /
                                EEF1_SOLV[ti][EEF1_LAMBDA] +
                            EEF1_SOLV[tj][EEF1_DG_FREE] *
                                EEF1_SOLV[ti][EEF1_VOL] * exp(-xji * xji) *
                                (xji / EEF1_SOLV[tj][EEF1_LAMBDA] + 1.0 / d) /
                                EEF1_SOLV[tj][EEF1_LAMBDA]) *
                           coef2 / d2;
                      F_SOLV[2 * ires + 1] += de / d * (theAtom->r ^ pt->r);
                      G_SOLV[2 * ires + 1] += de / d * (theAtom->r - pt->r);
                    }
                  }
                }
              }
              pt = pt->cell_next;
            }
          } // end x
        }   // end y
      }     // end z
      E_SOLV += SOLVi;
    }
    if (ires < gp.size() - 1) {
      F_VDWR[2 * ires + 1] += F_VDWR[2 * ires + 2];
      G_VDWR[2 * ires + 1] += G_VDWR[2 * ires + 2];
      F_VDWA[2 * ires + 1] += F_VDWA[2 * ires + 2];
      G_VDWA[2 * ires + 1] += G_VDWA[2 * ires + 2];
      F_SOLV[2 * ires + 1] += F_SOLV[2 * ires + 2];
      G_SOLV[2 * ires + 1] += G_SOLV[2 * ires + 2];
    }

    vec tmp_g = gp.getPsi_UV(ires) ^ gp.getResidue(ires)->getCA()->r;
    DVDWR[2 * ires + 1] = -(gp.getPsi_UV(ires) * F_VDWR[2 * ires + 1]) -
                          tmp_g * G_VDWR[2 * ires + 1];
    DVDWA[2 * ires + 1] = -(gp.getPsi_UV(ires) * F_VDWA[2 * ires + 1]) -
                          tmp_g * G_VDWA[2 * ires + 1];
    DSOLV[2 * ires + 1] = -(gp.getPsi_UV(ires) * F_SOLV[2 * ires + 1]) -
                          tmp_g * G_SOLV[2 * ires + 1];

    // phi: Ci, Ri
    F_VDWR[2 * ires] = vec(0, 0, 0);
    G_VDWR[2 * ires] = vec(0, 0, 0);
    F_VDWA[2 * ires] = vec(0, 0, 0);
    G_VDWA[2 * ires] = vec(0, 0, 0);
    F_SOLV[2 * ires] = vec(0, 0, 0);
    G_SOLV[2 * ires] = vec(0, 0, 0);
    tmp_va.clear();
    theRes = gp.getResidue(ires);
    for (int iatom = 0; iatom < theRes->getNA(); iatom++) {
      theAtom = theRes->getAtom(iatom);
      if (theAtom->isSC()) {
        tmp_va.push_back(theAtom);
      } else if (theAtom->PDBtype == _C_) {
        tmp_va.push_back(theAtom);
      }
    }
    for (int iatom = 0; iatom < tmp_va.size(); iatom++) {
      theAtom = tmp_va[iatom];
      int id = theAtom->resID;
      ti = theAtom->FF_t;
      ix = static_cast<int>((theAtom->r.x - min.x) / cell_box.x);
      iy = static_cast<int>((theAtom->r.y - min.y) / cell_box.y);
      iz = static_cast<int>((theAtom->r.z - min.z) / cell_box.z);
      SOLVi = 0;
      for (int cz = iz - step; cz <= iz + step; cz++) {
        kz = cz;
        if (cz < 0)
          kz += ncell[2];
        else if (cz >= ncell[2])
          kz -= ncell[2];

        for (int cy = iy - step; cy <= iy + step; cy++) {
          ky = cy;
          if (cy < 0)
            ky += ncell[1];
          else if (cy >= ncell[1])
            ky -= ncell[1];

          for (int cx = ix - step; cx <= ix + step; cx++) {
            kx = cx;
            if (cx < 0)
              kx += ncell[0];
            else if (cx >= ncell[0])
              kx -= ncell[0];

            cindex = ncell[0] * (ncell[1] * kz + ky) + kx;
            pt = cell_atoms[cindex];
            while (pt) {
              int jres = pt->resID;
              if (!pt->isH()) {
                if ((abs(id - jres) > 1) ||
                    (jres == id + 1 &&
                     ((theAtom->isSC() && (pt->isSC() || pt->PDBtype == _C_ ||
                                           pt->PDBtype == _O_)) ||
                      (pt->isSC() && (theAtom->PDBtype == _N_)))) ||
                    (jres == id - 1 &&
                     ((theAtom->isSC() && (pt->isSC() || pt->PDBtype == _N_)) ||
                      (pt->isSC() && (theAtom->PDBtype == _C_ ||
                                      theAtom->PDBtype == _O_))))) {
                  d2 = (pt->r - theAtom->r).mod2();
                  d = sqrt(d2);
                  tj = pt->FF_t;
                  if (d2 < max_ir2) { // within the cutoff distance
                    // vdw
                    vdw.getVDWdata(ti, tj, B, A, eps, sig, cutoff, slope,
                                   intercept);
                    if (B && A) {
                      d6 = d2 * d2 * d2;
                      e = (B / d6 - A) / d6;
                      if (e > 0) {
                        // e = (sig-d)/sig*10.0;
                        if (d < cutoff) {
                          e = slope * d + intercept;
                          de = slope;
                        } else {
                          de = 6.0 * (A - 2.0 * B / d6) / (d6 * d);
                        }
                        E_VDWR += e;
                        // de = -10.0/sig;
                        F_VDWR[2 * ires] += de / d * (theAtom->r ^ pt->r);
                        G_VDWR[2 * ires] += de / d * (theAtom->r - pt->r);
                      } else {
                        E_VDWA += e;
                        de = 6.0 * (A - 2.0 * B / d6) / (d6 * d);
                        F_VDWA[2 * ires] += de / d * (theAtom->r ^ pt->r);
                        G_VDWA[2 * ires] += de / d * (theAtom->r - pt->r);
                      }
                    }
                    // solv
                    int solv_cut = 0;
                    if (d < (EEF1_SOLV[ti][CHARMM_RADIUS] +
                             EEF1_SOLV[tj][CHARMM_RADIUS]) *
                                0.846) {
                      d = (EEF1_SOLV[ti][CHARMM_RADIUS] +
                           EEF1_SOLV[tj][CHARMM_RADIUS]) *
                          0.846;
                      d2 = d * d;
                      solv_cut = 1;
                    }
                    xij = (d - EEF1_SOLV[ti][CHARMM_RADIUS]) /
                          EEF1_SOLV[ti][EEF1_LAMBDA];
                    xji = (d - EEF1_SOLV[tj][CHARMM_RADIUS]) /
                          EEF1_SOLV[tj][EEF1_LAMBDA];
                    SOLVi -=
                        (EEF1_SOLV[ti][EEF1_DG_FREE] * EEF1_SOLV[tj][EEF1_VOL] *
                             exp(-xij * xij) / EEF1_SOLV[ti][EEF1_LAMBDA] +
                         EEF1_SOLV[tj][EEF1_DG_FREE] * EEF1_SOLV[ti][EEF1_VOL] *
                             exp(-xji * xji) / EEF1_SOLV[tj][EEF1_LAMBDA]) *
                        coef / d2;
                    if (!solv_cut) {
                      de = (EEF1_SOLV[ti][EEF1_DG_FREE] *
                                EEF1_SOLV[tj][EEF1_VOL] * exp(-xij * xij) *
                                (xij / EEF1_SOLV[ti][EEF1_LAMBDA] + 1.0 / d) /
                                EEF1_SOLV[ti][EEF1_LAMBDA] +
                            EEF1_SOLV[tj][EEF1_DG_FREE] *
                                EEF1_SOLV[ti][EEF1_VOL] * exp(-xji * xji) *
                                (xji / EEF1_SOLV[tj][EEF1_LAMBDA] + 1.0 / d) /
                                EEF1_SOLV[tj][EEF1_LAMBDA]) *
                           coef2 / d2;
                      F_SOLV[2 * ires] += de / d * (theAtom->r ^ pt->r);
                      G_SOLV[2 * ires] += de / d * (theAtom->r - pt->r);
                    }
                  }
                }
              }
              pt = pt->cell_next;
            }
          } // end x
        }   // end y
      }     // end z
      E_SOLV += SOLVi;
    }

    F_VDWR[2 * ires] += F_VDWR[2 * ires + 1];
    G_VDWR[2 * ires] += G_VDWR[2 * ires + 1];
    F_VDWA[2 * ires] += F_VDWA[2 * ires + 1];
    G_VDWA[2 * ires] += G_VDWA[2 * ires + 1];
    F_SOLV[2 * ires] += F_SOLV[2 * ires + 1];
    G_SOLV[2 * ires] += G_SOLV[2 * ires + 1];

    tmp_g = gp.getPhi_UV(ires) ^ gp.getResidue(ires)->getCA()->r;
    DVDWR[2 * ires] =
        -(gp.getPhi_UV(ires) * F_VDWR[2 * ires]) - tmp_g * G_VDWR[2 * ires];
    DVDWA[2 * ires] =
        -(gp.getPhi_UV(ires) * F_VDWA[2 * ires]) - tmp_g * G_VDWA[2 * ires];
    DSOLV[2 * ires] =
        -(gp.getPhi_UV(ires) * F_SOLV[2 * ires]) - tmp_g * G_SOLV[2 * ires];
  }

  delete[] F_VDWR;
  delete[] G_VDWR;
  delete[] F_VDWA;
  delete[] G_VDWA;
  delete[] F_SOLV;
  delete[] G_SOLV;
}

/*************************************************************
  using the trick of calculating the derivates by Go and Abe.

  loop only once
 **************************************************************/
double gen_grid::getE_DEV_GOTRICK1(double &E_VDWA, double &E_VDWR,
                                   double &E_SOLV, double *DVDWA, double *DVDWR,
                                   double *DSOLV) {
  E_VDWA = E_VDWR = E_SOLV = 0;
  residue *theRes = NULL;
  atom *theAtom = NULL;
  static double coef = 2.0 / (4.0 * PI * sqrt(PI));
  static double coef2 = 1.0 / (PI * sqrt(PI));
  int ti, tj;
  int ix, iy, iz, cindex, kx, ky, kz;
  int min_ix, min_iy, min_iz;
  double d2, d6, B, A, e, d, sig, eps, cutoff, slope, intercept;
  double SOLVi = 0, xij, xji, de;
  atom *pt;

  vec *F_VDWA = new vec[2 * gp.size()];
  vec *G_VDWA = new vec[2 * gp.size()];
  vec *F_VDWR = new vec[2 * gp.size()];
  vec *G_VDWR = new vec[2 * gp.size()];
  vec *F_SOLV = new vec[2 * gp.size()];
  vec *G_SOLV = new vec[2 * gp.size()];
  for (int i = 0; i < gp.size(); i++) {
    DVDWA[2 * i] = DVDWA[2 * i + 1] = 0;
    DVDWR[2 * i] = DVDWR[2 * i + 1] = 0;
    DSOLV[2 * i] = DSOLV[2 * i + 1] = 0;
    F_VDWA[2 * i + 1] = vec(0, 0, 0);
    G_VDWA[2 * i + 1] = vec(0, 0, 0);
    F_VDWR[2 * i + 1] = vec(0, 0, 0);
    G_VDWR[2 * i + 1] = vec(0, 0, 0);
    F_SOLV[2 * i + 1] = vec(0, 0, 0);
    G_SOLV[2 * i + 1] = vec(0, 0, 0);
    F_VDWR[2 * i] = vec(0, 0, 0);
    G_VDWR[2 * i] = vec(0, 0, 0);
    F_VDWA[2 * i] = vec(0, 0, 0);
    G_VDWA[2 * i] = vec(0, 0, 0);
    F_SOLV[2 * i] = vec(0, 0, 0);
    G_SOLV[2 * i] = vec(0, 0, 0);
  }
  gp.updateFY_UV();

  vector<atom *> tmp_va;
  vec f_comp, g_comp;
  for (int ires = gp.size() - 1; ires >= 0; ires--) {
    // psi: Oi, Ni+1, Hi+1, Cai+1
    tmp_va.clear();
    theRes = gp.getResidue(ires);
    tmp_va.push_back(theRes->getO());
    if (ires < gp.size() - 1) {
      theRes = gp.getResidue(ires + 1);
      tmp_va.push_back(theRes->getN());
      tmp_va.push_back(theRes->getCA());
    }
    for (int iatom = 0; iatom < tmp_va.size(); iatom++) {
      theAtom = tmp_va[iatom];
      int id = theAtom->resID;
      ti = theAtom->FF_t;
      ix = static_cast<int>((theAtom->r.x - min.x) / cell_box.x);
      iy = static_cast<int>((theAtom->r.y - min.y) / cell_box.y);
      iz = static_cast<int>((theAtom->r.z - min.z) / cell_box.z);
      SOLVi = 0;
      for (int cz = iz - step; cz <= iz + step; cz++) {
        kz = cz;
        if (cz < 0)
          kz += ncell[2];
        else if (cz >= ncell[2])
          kz -= ncell[2];

        for (int cy = iy - step; cy <= iy + step; cy++) {
          ky = cy;
          if (cy < 0)
            ky += ncell[1];
          else if (cy >= ncell[1])
            ky -= ncell[1];

          for (int cx = ix - step; cx <= ix + step; cx++) {
            kx = cx;
            if (cx < 0)
              kx += ncell[0];
            else if (cx >= ncell[0])
              kx -= ncell[0];

            cindex = ncell[0] * (ncell[1] * kz + ky) + kx;
            pt = cell_atoms[cindex];
            while (pt) {
              int jres = pt->resID;
              if (!pt->isH()) {
                if ((id > jres + 1) ||
                    (jres == id - 1 &&
                     ((theAtom->isSC() && (pt->isSC() || pt->PDBtype == _N_)) ||
                      (pt->isSC() && (theAtom->PDBtype == _C_ ||
                                      theAtom->PDBtype == _O_))))) {
                  d2 = (pt->r - theAtom->r).mod2();
                  d = sqrt(d2);
                  tj = pt->FF_t;
                  if (d2 < max_ir2) { // within the cutoff distance
                    // vdw
                    vdw.getVDWdata(ti, tj, B, A, eps, sig, cutoff, slope,
                                   intercept);
                    if (B && A) {
                      d6 = d2 * d2 * d2;
                      e = (B / d6 - A) / d6;
                      if (e > 0) {
                        // e = (sig-d)/sig*10.0;
                        if (d < cutoff) {
                          e = slope * d + intercept;
                          de = slope;
                        } else {
                          de = 6.0 * (A - 2.0 * B / d6) / (d6 * d);
                        }
                        E_VDWR += e;
                        // de = -10.0/sig;
                        F_VDWR[2 * ires + 1] += f_comp =
                            de / d * (theAtom->r ^ pt->r);
                        G_VDWR[2 * ires + 1] += g_comp =
                            de / d * (theAtom->r - pt->r);

                        if (pt->isSC() || pt->PDBtype == _C_) {
                          F_VDWR[2 * jres] -= f_comp;
                          G_VDWR[2 * jres] -= g_comp;
                        } else if ((pt->PDBtype == _N_ ||
                                    pt->PDBtype == _CA_) &&
                                   jres > 0) {
                          F_VDWR[2 * (jres - 1) + 1] -= f_comp;
                          G_VDWR[2 * (jres - 1) + 1] -= g_comp;
                        } else if (pt->PDBtype == _O_) {
                          F_VDWR[2 * jres + 1] -= f_comp;
                          G_VDWR[2 * jres + 1] -= g_comp;
                        }

                      } else {
                        E_VDWA += e;
                        de = 6.0 * (A - 2.0 * B / d6) / (d6 * d);
                        F_VDWA[2 * ires + 1] += f_comp =
                            de / d * (theAtom->r ^ pt->r);
                        G_VDWA[2 * ires + 1] += g_comp =
                            de / d * (theAtom->r - pt->r);

                        if (pt->isSC() || pt->PDBtype == _C_) {
                          F_VDWA[2 * jres] -= f_comp;
                          G_VDWA[2 * jres] -= g_comp;
                        } else if ((pt->PDBtype == _N_ ||
                                    pt->PDBtype == _CA_) &&
                                   jres > 0) {
                          F_VDWA[2 * (jres - 1) + 1] -= f_comp;
                          G_VDWA[2 * (jres - 1) + 1] -= g_comp;
                        } else if (pt->PDBtype == _O_) {
                          F_VDWA[2 * jres + 1] -= f_comp;
                          G_VDWA[2 * jres + 1] -= g_comp;
                        }
                      }
                    }
                    // solv
                    int solv_cut = 0;
                    if (d < (EEF1_SOLV[ti][CHARMM_RADIUS] +
                             EEF1_SOLV[tj][CHARMM_RADIUS]) *
                                0.846) {
                      d = (EEF1_SOLV[ti][CHARMM_RADIUS] +
                           EEF1_SOLV[tj][CHARMM_RADIUS]) *
                          0.846;
                      d2 = d * d;
                      solv_cut = 1;
                    }
                    xij = (d - EEF1_SOLV[ti][CHARMM_RADIUS]) /
                          EEF1_SOLV[ti][EEF1_LAMBDA];
                    xji = (d - EEF1_SOLV[tj][CHARMM_RADIUS]) /
                          EEF1_SOLV[tj][EEF1_LAMBDA];
                    SOLVi -=
                        (EEF1_SOLV[ti][EEF1_DG_FREE] * EEF1_SOLV[tj][EEF1_VOL] *
                             exp(-xij * xij) / EEF1_SOLV[ti][EEF1_LAMBDA] +
                         EEF1_SOLV[tj][EEF1_DG_FREE] * EEF1_SOLV[ti][EEF1_VOL] *
                             exp(-xji * xji) / EEF1_SOLV[tj][EEF1_LAMBDA]) *
                        coef / d2;
                    if (!solv_cut) {
                      de = (EEF1_SOLV[ti][EEF1_DG_FREE] *
                                EEF1_SOLV[tj][EEF1_VOL] * exp(-xij * xij) *
                                (xij / EEF1_SOLV[ti][EEF1_LAMBDA] + 1.0 / d) /
                                EEF1_SOLV[ti][EEF1_LAMBDA] +
                            EEF1_SOLV[tj][EEF1_DG_FREE] *
                                EEF1_SOLV[ti][EEF1_VOL] * exp(-xji * xji) *
                                (xji / EEF1_SOLV[tj][EEF1_LAMBDA] + 1.0 / d) /
                                EEF1_SOLV[tj][EEF1_LAMBDA]) *
                           coef2 / d2;
                      F_SOLV[2 * ires + 1] += f_comp =
                          de / d * (theAtom->r ^ pt->r);
                      G_SOLV[2 * ires + 1] += g_comp =
                          de / d * (theAtom->r - pt->r);
                      if (pt->isSC() || pt->PDBtype == _C_) {
                        F_SOLV[2 * jres] -= f_comp;
                        G_SOLV[2 * jres] -= g_comp;
                      } else if ((pt->PDBtype == _N_ || pt->PDBtype == _CA_) &&
                                 jres > 0) {
                        F_SOLV[2 * (jres - 1) + 1] -= f_comp;
                        G_SOLV[2 * (jres - 1) + 1] -= g_comp;
                      } else if (pt->PDBtype == _O_) {
                        F_SOLV[2 * jres + 1] -= f_comp;
                        G_SOLV[2 * jres + 1] -= g_comp;
                      }
                    }
                  }
                }
              }
              pt = pt->cell_next;
            }
          } // end x
        }   // end y
      }     // end z
      E_SOLV += SOLVi;
    }
    if (ires < gp.size() - 1) {
      F_VDWR[2 * ires + 1] += F_VDWR[2 * ires + 2];
      G_VDWR[2 * ires + 1] += G_VDWR[2 * ires + 2];
      F_VDWA[2 * ires + 1] += F_VDWA[2 * ires + 2];
      G_VDWA[2 * ires + 1] += G_VDWA[2 * ires + 2];
      F_SOLV[2 * ires + 1] += F_SOLV[2 * ires + 2];
      G_SOLV[2 * ires + 1] += G_SOLV[2 * ires + 2];
    }

    vec tmp_g = gp.getPsi_UV(ires) ^ gp.getResidue(ires)->getCA()->r;
    DVDWR[2 * ires + 1] = -(gp.getPsi_UV(ires) * F_VDWR[2 * ires + 1]) -
                          tmp_g * G_VDWR[2 * ires + 1];
    DVDWA[2 * ires + 1] = -(gp.getPsi_UV(ires) * F_VDWA[2 * ires + 1]) -
                          tmp_g * G_VDWA[2 * ires + 1];
    DSOLV[2 * ires + 1] = -(gp.getPsi_UV(ires) * F_SOLV[2 * ires + 1]) -
                          tmp_g * G_SOLV[2 * ires + 1];

    // phi: Ci, Ri
    tmp_va.clear();
    theRes = gp.getResidue(ires);
    for (int iatom = 0; iatom < theRes->getNA(); iatom++) {
      theAtom = theRes->getAtom(iatom);
      if (theAtom->isSC()) {
        tmp_va.push_back(theAtom);
      } else if (theAtom->PDBtype == _C_) {
        tmp_va.push_back(theAtom);
      }
    }
    for (int iatom = 0; iatom < tmp_va.size(); iatom++) {
      theAtom = tmp_va[iatom];
      int id = theAtom->resID;
      ti = theAtom->FF_t;
      ix = static_cast<int>((theAtom->r.x - min.x) / cell_box.x);
      iy = static_cast<int>((theAtom->r.y - min.y) / cell_box.y);
      iz = static_cast<int>((theAtom->r.z - min.z) / cell_box.z);
      SOLVi = 0;
      for (int cz = iz - step; cz <= iz + step; cz++) {
        kz = cz;
        if (cz < 0)
          kz += ncell[2];
        else if (cz >= ncell[2])
          kz -= ncell[2];

        for (int cy = iy - step; cy <= iy + step; cy++) {
          ky = cy;
          if (cy < 0)
            ky += ncell[1];
          else if (cy >= ncell[1])
            ky -= ncell[1];

          for (int cx = ix - step; cx <= ix + step; cx++) {
            kx = cx;
            if (cx < 0)
              kx += ncell[0];
            else if (cx >= ncell[0])
              kx -= ncell[0];

            cindex = ncell[0] * (ncell[1] * kz + ky) + kx;
            pt = cell_atoms[cindex];
            while (pt) {
              int jres = pt->resID;
              if (!pt->isH()) {
                if ((id > jres + 1) ||
                    (jres == id - 1 &&
                     ((theAtom->isSC() && (pt->isSC() || pt->PDBtype == _N_)) ||
                      (pt->isSC() && (theAtom->PDBtype == _C_ ||
                                      theAtom->PDBtype == _O_))))) {
                  d2 = (pt->r - theAtom->r).mod2();
                  d = sqrt(d2);
                  tj = pt->FF_t;
                  if (d2 < max_ir2) { // within the cutoff distance
                    // vdw
                    vdw.getVDWdata(ti, tj, B, A, eps, sig, cutoff, slope,
                                   intercept);
                    if (B && A) {
                      d6 = d2 * d2 * d2;
                      e = (B / d6 - A) / d6;
                      if (e > 0) {
                        // e = (sig-d)/sig*10.0;
                        if (d < cutoff) {
                          e = slope * d + intercept;
                          de = slope;
                        } else {
                          de = 6.0 * (A - 2.0 * B / d6) / (d6 * d);
                        }
                        E_VDWR += e;
                        // de = -10.0/sig;
                        F_VDWR[2 * ires] += f_comp =
                            de / d * (theAtom->r ^ pt->r);
                        G_VDWR[2 * ires] += g_comp =
                            de / d * (theAtom->r - pt->r);
                        if (pt->isSC() || pt->PDBtype == _C_) {
                          F_VDWR[2 * jres] -= f_comp;
                          G_VDWR[2 * jres] -= g_comp;
                        } else if ((pt->PDBtype == _N_ ||
                                    pt->PDBtype == _CA_) &&
                                   jres > 0) {
                          F_VDWR[2 * (jres - 1) + 1] -= f_comp;
                          G_VDWR[2 * (jres - 1) + 1] -= g_comp;
                        } else if (pt->PDBtype == _O_) {
                          F_VDWR[2 * jres + 1] -= f_comp;
                          G_VDWR[2 * jres + 1] -= g_comp;
                        }
                      } else {
                        E_VDWA += e;
                        de = 6.0 * (A - 2.0 * B / d6) / (d6 * d);
                        F_VDWA[2 * ires] += f_comp =
                            de / d * (theAtom->r ^ pt->r);
                        G_VDWA[2 * ires] += g_comp =
                            de / d * (theAtom->r - pt->r);
                        if (pt->isSC() || pt->PDBtype == _C_) {
                          F_VDWA[2 * jres] -= f_comp;
                          G_VDWA[2 * jres] -= g_comp;
                        } else if ((pt->PDBtype == _N_ ||
                                    pt->PDBtype == _CA_) &&
                                   jres > 0) {
                          F_VDWA[2 * (jres - 1) + 1] -= f_comp;
                          G_VDWA[2 * (jres - 1) + 1] -= g_comp;
                        } else if (pt->PDBtype == _O_) {
                          F_VDWA[2 * jres + 1] -= f_comp;
                          G_VDWA[2 * jres + 1] -= g_comp;
                        }
                      }
                    }
                    // solv
                    int solv_cut = 0;
                    if (d < (EEF1_SOLV[ti][CHARMM_RADIUS] +
                             EEF1_SOLV[tj][CHARMM_RADIUS]) *
                                0.846) {
                      d = (EEF1_SOLV[ti][CHARMM_RADIUS] +
                           EEF1_SOLV[tj][CHARMM_RADIUS]) *
                          0.846;
                      d2 = d * d;
                      solv_cut = 1;
                    }
                    xij = (d - EEF1_SOLV[ti][CHARMM_RADIUS]) /
                          EEF1_SOLV[ti][EEF1_LAMBDA];
                    xji = (d - EEF1_SOLV[tj][CHARMM_RADIUS]) /
                          EEF1_SOLV[tj][EEF1_LAMBDA];
                    SOLVi -=
                        (EEF1_SOLV[ti][EEF1_DG_FREE] * EEF1_SOLV[tj][EEF1_VOL] *
                             exp(-xij * xij) / EEF1_SOLV[ti][EEF1_LAMBDA] +
                         EEF1_SOLV[tj][EEF1_DG_FREE] * EEF1_SOLV[ti][EEF1_VOL] *
                             exp(-xji * xji) / EEF1_SOLV[tj][EEF1_LAMBDA]) *
                        coef / d2;
                    if (!solv_cut) {
                      de = (EEF1_SOLV[ti][EEF1_DG_FREE] *
                                EEF1_SOLV[tj][EEF1_VOL] * exp(-xij * xij) *
                                (xij / EEF1_SOLV[ti][EEF1_LAMBDA] + 1.0 / d) /
                                EEF1_SOLV[ti][EEF1_LAMBDA] +
                            EEF1_SOLV[tj][EEF1_DG_FREE] *
                                EEF1_SOLV[ti][EEF1_VOL] * exp(-xji * xji) *
                                (xji / EEF1_SOLV[tj][EEF1_LAMBDA] + 1.0 / d) /
                                EEF1_SOLV[tj][EEF1_LAMBDA]) *
                           coef2 / d2;
                      F_SOLV[2 * ires] += f_comp =
                          de / d * (theAtom->r ^ pt->r);
                      G_SOLV[2 * ires] += g_comp =
                          de / d * (theAtom->r - pt->r);
                      if (pt->isSC() || pt->PDBtype == _C_) {
                        F_SOLV[2 * jres] -= f_comp;
                        G_SOLV[2 * jres] -= g_comp;
                      } else if ((pt->PDBtype == _N_ || pt->PDBtype == _CA_) &&
                                 jres > 0) {
                        F_SOLV[2 * (jres - 1) + 1] -= f_comp;
                        G_SOLV[2 * (jres - 1) + 1] -= g_comp;
                      } else if (pt->PDBtype == _O_) {
                        F_SOLV[2 * jres + 1] -= f_comp;
                        G_SOLV[2 * jres + 1] -= g_comp;
                      }
                    }
                  }
                }
              }
              pt = pt->cell_next;
            }
          } // end x
        }   // end y
      }     // end z
      E_SOLV += SOLVi;
    }

    F_VDWR[2 * ires] += F_VDWR[2 * ires + 1];
    G_VDWR[2 * ires] += G_VDWR[2 * ires + 1];
    F_VDWA[2 * ires] += F_VDWA[2 * ires + 1];
    G_VDWA[2 * ires] += G_VDWA[2 * ires + 1];
    F_SOLV[2 * ires] += F_SOLV[2 * ires + 1];
    G_SOLV[2 * ires] += G_SOLV[2 * ires + 1];

    tmp_g = gp.getPhi_UV(ires) ^ gp.getResidue(ires)->getCA()->r;
    DVDWR[2 * ires] =
        -(gp.getPhi_UV(ires) * F_VDWR[2 * ires]) - tmp_g * G_VDWR[2 * ires];
    DVDWA[2 * ires] =
        -(gp.getPhi_UV(ires) * F_VDWA[2 * ires]) - tmp_g * G_VDWA[2 * ires];
    DSOLV[2 * ires] =
        -(gp.getPhi_UV(ires) * F_SOLV[2 * ires]) - tmp_g * G_SOLV[2 * ires];
  }

  delete[] F_VDWR;
  delete[] G_VDWR;
  delete[] F_VDWA;
  delete[] G_VDWA;
  delete[] F_SOLV;
  delete[] G_SOLV;
}

/*************************************************************
  using the trick of calculating the derivates by Go and Abe.
  including phi, psi, and omiga
 **************************************************************/
double gen_grid::getE_DEV_FYO_GOTRICK(double &E_VDWA, double &E_VDWR,
                                      double &E_SOLV, double *DVDWA,
                                      double *DVDWR, double *DSOLV) {
  E_VDWA = E_VDWR = E_SOLV = 0;
  residue *theRes = NULL;
  atom *theAtom = NULL;
  static double coef = 2.0 / (4.0 * PI * sqrt(PI));
  static double coef2 = 1.0 / (PI * sqrt(PI));
  int ti, tj;
  int ix, iy, iz, cindex, kx, ky, kz;
  int min_ix, min_iy, min_iz;
  double d2, d6, B, A, e, d, sig, eps, cutoff, slope, intercept;
  double SOLVi = 0, xij, xji, de;
  atom *pt;

  for (int i = 0; i < gp.size(); i++) {
    DVDWA[3 * i] = DVDWA[3 * i + 1] = DVDWA[3 * i + 2] = 0;
    DVDWR[3 * i] = DVDWR[3 * i + 1] = DVDWR[3 * i + 2] = 0;
    DSOLV[3 * i] = DSOLV[3 * i + 1] = DSOLV[3 * i + 2] = 0;
  }
  gp.updateFY_UV();

  vec *F_VDWA = new vec[3 * gp.size()];
  vec *G_VDWA = new vec[3 * gp.size()];
  vec *F_VDWR = new vec[3 * gp.size()];
  vec *G_VDWR = new vec[3 * gp.size()];
  vec *F_SOLV = new vec[3 * gp.size()];
  vec *G_SOLV = new vec[3 * gp.size()];
  vector<atom *> tmp_va;
  for (int ires = gp.size() - 1; ires >= 0; ires--) {
    // omiga: Hi+1, CAi+1(2)
    F_VDWA[3 * ires + 2] = vec(0, 0, 0);
    G_VDWA[3 * ires + 2] = vec(0, 0, 0);
    F_VDWR[3 * ires + 2] = vec(0, 0, 0);
    G_VDWR[3 * ires + 2] = vec(0, 0, 0);
    F_SOLV[3 * ires + 2] = vec(0, 0, 0);
    G_SOLV[3 * ires + 2] = vec(0, 0, 0);
    if (ires < gp.size() - 1) {
      theAtom = gp.getResidue(ires + 1)->getCA();
      int id = ires + 1;
      ti = theAtom->FF_t;
      ix = static_cast<int>((theAtom->r.x - min.x) / cell_box.x);
      iy = static_cast<int>((theAtom->r.y - min.y) / cell_box.y);
      iz = static_cast<int>((theAtom->r.z - min.z) / cell_box.z);
      SOLVi = 0;
      for (int cz = iz - step; cz <= iz + step; cz++) {
        kz = cz;
        if (cz < 0)
          kz += ncell[2];
        else if (cz >= ncell[2])
          kz -= ncell[2];

        for (int cy = iy - step; cy <= iy + step; cy++) {
          ky = cy;
          if (cy < 0)
            ky += ncell[1];
          else if (cy >= ncell[1])
            ky -= ncell[1];

          for (int cx = ix - step; cx <= ix + step; cx++) {
            kx = cx;
            if (cx < 0)
              kx += ncell[0];
            else if (cx >= ncell[0])
              kx -= ncell[0];

            cindex = ncell[0] * (ncell[1] * kz + ky) + kx;
            pt = cell_atoms[cindex];
            while (pt) {
              int jres = pt->resID;
              if (!pt->isH()) {
                if (abs(id - jres) > 1) {
                  d2 = (pt->r - theAtom->r).mod2();
                  d = sqrt(d2);
                  tj = pt->FF_t;
                  if (d2 < max_ir2) { // within the cutoff distance
                    // vdw
                    vdw.getVDWdata(ti, tj, B, A, eps, sig, cutoff, slope,
                                   intercept);
                    if (B && A) {
                      d6 = d2 * d2 * d2;
                      e = (B / d6 - A) / d6;
                      if (e > 0) {
                        // e = (sig-d)/sig*10.0;
                        if (d < cutoff) {
                          e = slope * d + intercept;
                          de = slope;
                        } else {
                          de = 6.0 * (A - 2.0 * B / d6) / (d6 * d);
                        }
                        E_VDWR += e;
                        // de = -10.0/sig;
                        F_VDWR[3 * ires + 2] += de / d * (theAtom->r ^ pt->r);
                        G_VDWR[3 * ires + 2] += de / d * (theAtom->r - pt->r);
                      } else {
                        E_VDWA += e;
                        de = 6.0 * (A - 2.0 * B / d6) / (d6 * d);
                        F_VDWA[3 * ires + 2] += de / d * (theAtom->r ^ pt->r);
                        G_VDWA[3 * ires + 2] += de / d * (theAtom->r - pt->r);
                      }
                    }
                    // solv
                    int solv_cut = 0;
                    if (d < (EEF1_SOLV[ti][CHARMM_RADIUS] +
                             EEF1_SOLV[tj][CHARMM_RADIUS]) *
                                0.846) {
                      d = (EEF1_SOLV[ti][CHARMM_RADIUS] +
                           EEF1_SOLV[tj][CHARMM_RADIUS]) *
                          0.846;
                      d2 = d * d;
                      solv_cut = 1;
                    }
                    xij = (d - EEF1_SOLV[ti][CHARMM_RADIUS]) /
                          EEF1_SOLV[ti][EEF1_LAMBDA];
                    xji = (d - EEF1_SOLV[tj][CHARMM_RADIUS]) /
                          EEF1_SOLV[tj][EEF1_LAMBDA];
                    SOLVi -=
                        (EEF1_SOLV[ti][EEF1_DG_FREE] * EEF1_SOLV[tj][EEF1_VOL] *
                             exp(-xij * xij) / EEF1_SOLV[ti][EEF1_LAMBDA] +
                         EEF1_SOLV[tj][EEF1_DG_FREE] * EEF1_SOLV[ti][EEF1_VOL] *
                             exp(-xji * xji) / EEF1_SOLV[tj][EEF1_LAMBDA]) *
                        coef / d2;
                    if (!solv_cut) {
                      de = (EEF1_SOLV[ti][EEF1_DG_FREE] *
                                EEF1_SOLV[tj][EEF1_VOL] * exp(-xij * xij) *
                                (xij / EEF1_SOLV[ti][EEF1_LAMBDA] + 1.0 / d) /
                                EEF1_SOLV[ti][EEF1_LAMBDA] +
                            EEF1_SOLV[tj][EEF1_DG_FREE] *
                                EEF1_SOLV[ti][EEF1_VOL] * exp(-xji * xji) *
                                (xji / EEF1_SOLV[tj][EEF1_LAMBDA] + 1.0 / d) /
                                EEF1_SOLV[tj][EEF1_LAMBDA]) *
                           coef2 / d2;
                      F_SOLV[3 * ires + 2] += de / d * (theAtom->r ^ pt->r);
                      G_SOLV[3 * ires + 2] += de / d * (theAtom->r - pt->r);
                    }
                  }
                }
              }
              pt = pt->cell_next;
            }
          } // end x
        }   // end y
      }     // end z
      E_SOLV += SOLVi;
      F_VDWR[3 * ires + 2] += F_VDWR[3 * ires + 3];
      G_VDWR[3 * ires + 2] += G_VDWR[3 * ires + 3];
      F_VDWA[3 * ires + 2] += F_VDWA[3 * ires + 3];
      G_VDWA[3 * ires + 2] += G_VDWA[3 * ires + 3];
      F_SOLV[3 * ires + 2] += F_SOLV[3 * ires + 3];
      G_SOLV[3 * ires + 2] += G_SOLV[3 * ires + 3];

      vec tmp_g = gp.getOmiga_UV(ires) ^ gp.getResidue(ires + 1)->getN()->r;
      DVDWR[3 * ires + 2] = -(gp.getOmiga_UV(ires) * F_VDWR[3 * ires + 2]) -
                            tmp_g * G_VDWR[3 * ires + 2];
      DVDWA[3 * ires + 2] = -(gp.getOmiga_UV(ires) * F_VDWA[3 * ires + 2]) -
                            tmp_g * G_VDWA[3 * ires + 2];
      DSOLV[3 * ires + 2] = -(gp.getOmiga_UV(ires) * F_SOLV[3 * ires + 2]) -
                            tmp_g * G_SOLV[3 * ires + 2];
    }

    // psi: Oi, Ni+1(1)
    F_VDWA[3 * ires + 1] = vec(0, 0, 0);
    G_VDWA[3 * ires + 1] = vec(0, 0, 0);
    F_VDWR[3 * ires + 1] = vec(0, 0, 0);
    G_VDWR[3 * ires + 1] = vec(0, 0, 0);
    F_SOLV[3 * ires + 1] = vec(0, 0, 0);
    G_SOLV[3 * ires + 1] = vec(0, 0, 0);
    tmp_va.clear();
    theRes = gp.getResidue(ires);
    tmp_va.push_back(theRes->getO());
    if (ires < gp.size() - 1) {
      theRes = gp.getResidue(ires + 1);
      tmp_va.push_back(theRes->getN());
    }
    for (int iatom = 0; iatom < tmp_va.size(); iatom++) {
      theAtom = tmp_va[iatom];
      int id = theAtom->resID;
      ti = theAtom->FF_t;
      ix = static_cast<int>((theAtom->r.x - min.x) / cell_box.x);
      iy = static_cast<int>((theAtom->r.y - min.y) / cell_box.y);
      iz = static_cast<int>((theAtom->r.z - min.z) / cell_box.z);
      SOLVi = 0;
      for (int cz = iz - step; cz <= iz + step; cz++) {
        kz = cz;
        if (cz < 0)
          kz += ncell[2];
        else if (cz >= ncell[2])
          kz -= ncell[2];

        for (int cy = iy - step; cy <= iy + step; cy++) {
          ky = cy;
          if (cy < 0)
            ky += ncell[1];
          else if (cy >= ncell[1])
            ky -= ncell[1];

          for (int cx = ix - step; cx <= ix + step; cx++) {
            kx = cx;
            if (cx < 0)
              kx += ncell[0];
            else if (cx >= ncell[0])
              kx -= ncell[0];

            cindex = ncell[0] * (ncell[1] * kz + ky) + kx;
            pt = cell_atoms[cindex];
            while (pt) {
              int jres = pt->resID;
              if (!pt->isH()) {
                if ((abs(id - jres) > 1) ||
                    (jres == id + 1 &&
                     ((pt->isSC() && (theAtom->PDBtype == _N_)))) ||
                    (jres == id - 1 &&
                     ((pt->isSC() && (theAtom->PDBtype == _O_))))) {
                  d2 = (pt->r - theAtom->r).mod2();
                  d = sqrt(d2);
                  tj = pt->FF_t;
                  if (d2 < max_ir2) { // within the cutoff distance
                    // vdw
                    vdw.getVDWdata(ti, tj, B, A, eps, sig, cutoff, slope,
                                   intercept);
                    if (B && A) {
                      d6 = d2 * d2 * d2;
                      e = (B / d6 - A) / d6;
                      if (e > 0) {
                        // e = (sig-d)/sig*10.0;
                        if (d < cutoff) {
                          e = slope * d + intercept;
                          de = slope;
                        } else {
                          de = 6.0 * (A - 2.0 * B / d6) / (d6 * d);
                        }
                        E_VDWR += e;
                        // de = -10.0/sig;
                        F_VDWR[3 * ires + 1] += de / d * (theAtom->r ^ pt->r);
                        G_VDWR[3 * ires + 1] += de / d * (theAtom->r - pt->r);
                      } else {
                        E_VDWA += e;
                        de = 6.0 * (A - 2.0 * B / d6) / (d6 * d);
                        F_VDWA[3 * ires + 1] += de / d * (theAtom->r ^ pt->r);
                        G_VDWA[3 * ires + 1] += de / d * (theAtom->r - pt->r);
                      }
                    }
                    // solv
                    int solv_cut = 0;
                    if (d < (EEF1_SOLV[ti][CHARMM_RADIUS] +
                             EEF1_SOLV[tj][CHARMM_RADIUS]) *
                                0.846) {
                      d = (EEF1_SOLV[ti][CHARMM_RADIUS] +
                           EEF1_SOLV[tj][CHARMM_RADIUS]) *
                          0.846;
                      d2 = d * d;
                      solv_cut = 1;
                    }
                    xij = (d - EEF1_SOLV[ti][CHARMM_RADIUS]) /
                          EEF1_SOLV[ti][EEF1_LAMBDA];
                    xji = (d - EEF1_SOLV[tj][CHARMM_RADIUS]) /
                          EEF1_SOLV[tj][EEF1_LAMBDA];
                    SOLVi -=
                        (EEF1_SOLV[ti][EEF1_DG_FREE] * EEF1_SOLV[tj][EEF1_VOL] *
                             exp(-xij * xij) / EEF1_SOLV[ti][EEF1_LAMBDA] +
                         EEF1_SOLV[tj][EEF1_DG_FREE] * EEF1_SOLV[ti][EEF1_VOL] *
                             exp(-xji * xji) / EEF1_SOLV[tj][EEF1_LAMBDA]) *
                        coef / d2;
                    if (!solv_cut) {
                      de = (EEF1_SOLV[ti][EEF1_DG_FREE] *
                                EEF1_SOLV[tj][EEF1_VOL] * exp(-xij * xij) *
                                (xij / EEF1_SOLV[ti][EEF1_LAMBDA] + 1.0 / d) /
                                EEF1_SOLV[ti][EEF1_LAMBDA] +
                            EEF1_SOLV[tj][EEF1_DG_FREE] *
                                EEF1_SOLV[ti][EEF1_VOL] * exp(-xji * xji) *
                                (xji / EEF1_SOLV[tj][EEF1_LAMBDA] + 1.0 / d) /
                                EEF1_SOLV[tj][EEF1_LAMBDA]) *
                           coef2 / d2;
                      F_SOLV[3 * ires + 1] += de / d * (theAtom->r ^ pt->r);
                      G_SOLV[3 * ires + 1] += de / d * (theAtom->r - pt->r);
                    }
                  }
                }
              }
              pt = pt->cell_next;
            }
          } // end x
        }   // end y
      }     // end z
      E_SOLV += SOLVi;
    }
    F_VDWR[3 * ires + 1] += F_VDWR[3 * ires + 2];
    G_VDWR[3 * ires + 1] += G_VDWR[3 * ires + 2];
    F_VDWA[3 * ires + 1] += F_VDWA[3 * ires + 2];
    G_VDWA[3 * ires + 1] += G_VDWA[3 * ires + 2];
    F_SOLV[3 * ires + 1] += F_SOLV[3 * ires + 2];
    G_SOLV[3 * ires + 1] += G_SOLV[3 * ires + 2];

    vec tmp_g = gp.getPsi_UV(ires) ^ gp.getResidue(ires)->getC()->r;
    DVDWR[3 * ires + 1] = -(gp.getPsi_UV(ires) * F_VDWR[3 * ires + 1]) -
                          tmp_g * G_VDWR[3 * ires + 1];
    DVDWA[3 * ires + 1] = -(gp.getPsi_UV(ires) * F_VDWA[3 * ires + 1]) -
                          tmp_g * G_VDWA[3 * ires + 1];
    DSOLV[3 * ires + 1] = -(gp.getPsi_UV(ires) * F_SOLV[3 * ires + 1]) -
                          tmp_g * G_SOLV[3 * ires + 1];

    // phi: Ci, Ri
    F_VDWR[3 * ires] = vec(0, 0, 0);
    G_VDWR[3 * ires] = vec(0, 0, 0);
    F_VDWA[3 * ires] = vec(0, 0, 0);
    G_VDWA[3 * ires] = vec(0, 0, 0);
    F_SOLV[3 * ires] = vec(0, 0, 0);
    G_SOLV[3 * ires] = vec(0, 0, 0);
    tmp_va.clear();
    theRes = gp.getResidue(ires);
    for (int iatom = 0; iatom < theRes->getNA(); iatom++) {
      theAtom = theRes->getAtom(iatom);
      if (theAtom->isSC()) {
        tmp_va.push_back(theAtom);
      } else if (theAtom->PDBtype == _C_) {
        tmp_va.push_back(theAtom);
      }
    }
    for (int iatom = 0; iatom < tmp_va.size(); iatom++) {
      theAtom = tmp_va[iatom];
      int id = theAtom->resID;
      ti = theAtom->FF_t;
      ix = static_cast<int>((theAtom->r.x - min.x) / cell_box.x);
      iy = static_cast<int>((theAtom->r.y - min.y) / cell_box.y);
      iz = static_cast<int>((theAtom->r.z - min.z) / cell_box.z);
      SOLVi = 0;
      for (int cz = iz - step; cz <= iz + step; cz++) {
        kz = cz;
        if (cz < 0)
          kz += ncell[2];
        else if (cz >= ncell[2])
          kz -= ncell[2];

        for (int cy = iy - step; cy <= iy + step; cy++) {
          ky = cy;
          if (cy < 0)
            ky += ncell[1];
          else if (cy >= ncell[1])
            ky -= ncell[1];

          for (int cx = ix - step; cx <= ix + step; cx++) {
            kx = cx;
            if (cx < 0)
              kx += ncell[0];
            else if (cx >= ncell[0])
              kx -= ncell[0];

            cindex = ncell[0] * (ncell[1] * kz + ky) + kx;
            pt = cell_atoms[cindex];
            while (pt) {
              int jres = pt->resID;
              if (!pt->isH()) {
                if ((abs(id - jres) > 1) ||
                    (jres == id + 1 &&
                     ((theAtom->isSC() && (pt->isSC() || pt->PDBtype == _C_ ||
                                           pt->PDBtype == _O_)))) ||
                    (jres == id - 1 &&
                     ((theAtom->isSC() && (pt->isSC() || pt->PDBtype == _N_)) ||
                      (pt->isSC() && (theAtom->PDBtype == _C_))))) {
                  d2 = (pt->r - theAtom->r).mod2();
                  d = sqrt(d2);
                  tj = pt->FF_t;
                  if (d2 < max_ir2) { // within the cutoff distance
                    // vdw
                    vdw.getVDWdata(ti, tj, B, A, eps, sig, cutoff, slope,
                                   intercept);
                    if (B && A) {
                      d6 = d2 * d2 * d2;
                      e = (B / d6 - A) / d6;
                      if (e > 0) {
                        // e = (sig-d)/sig*10.0;
                        if (d < cutoff) {
                          e = slope * d + intercept;
                          de = slope;
                        } else {
                          de = 6.0 * (A - 2.0 * B / d6) / (d6 * d);
                        }
                        E_VDWR += e;
                        // de = -10.0/sig;
                        F_VDWR[3 * ires] += de / d * (theAtom->r ^ pt->r);
                        G_VDWR[3 * ires] += de / d * (theAtom->r - pt->r);
                      } else {
                        E_VDWA += e;
                        de = 6.0 * (A - 2.0 * B / d6) / (d6 * d);
                        F_VDWA[3 * ires] += de / d * (theAtom->r ^ pt->r);
                        G_VDWA[3 * ires] += de / d * (theAtom->r - pt->r);
                      }
                    }
                    // solv
                    int solv_cut = 0;
                    if (d < (EEF1_SOLV[ti][CHARMM_RADIUS] +
                             EEF1_SOLV[tj][CHARMM_RADIUS]) *
                                0.846) {
                      d = (EEF1_SOLV[ti][CHARMM_RADIUS] +
                           EEF1_SOLV[tj][CHARMM_RADIUS]) *
                          0.846;
                      d2 = d * d;
                      solv_cut = 1;
                    }
                    xij = (d - EEF1_SOLV[ti][CHARMM_RADIUS]) /
                          EEF1_SOLV[ti][EEF1_LAMBDA];
                    xji = (d - EEF1_SOLV[tj][CHARMM_RADIUS]) /
                          EEF1_SOLV[tj][EEF1_LAMBDA];
                    SOLVi -=
                        (EEF1_SOLV[ti][EEF1_DG_FREE] * EEF1_SOLV[tj][EEF1_VOL] *
                             exp(-xij * xij) / EEF1_SOLV[ti][EEF1_LAMBDA] +
                         EEF1_SOLV[tj][EEF1_DG_FREE] * EEF1_SOLV[ti][EEF1_VOL] *
                             exp(-xji * xji) / EEF1_SOLV[tj][EEF1_LAMBDA]) *
                        coef / d2;
                    if (!solv_cut) {
                      de = (EEF1_SOLV[ti][EEF1_DG_FREE] *
                                EEF1_SOLV[tj][EEF1_VOL] * exp(-xij * xij) *
                                (xij / EEF1_SOLV[ti][EEF1_LAMBDA] + 1.0 / d) /
                                EEF1_SOLV[ti][EEF1_LAMBDA] +
                            EEF1_SOLV[tj][EEF1_DG_FREE] *
                                EEF1_SOLV[ti][EEF1_VOL] * exp(-xji * xji) *
                                (xji / EEF1_SOLV[tj][EEF1_LAMBDA] + 1.0 / d) /
                                EEF1_SOLV[tj][EEF1_LAMBDA]) *
                           coef2 / d2;
                      F_SOLV[3 * ires] += de / d * (theAtom->r ^ pt->r);
                      G_SOLV[3 * ires] += de / d * (theAtom->r - pt->r);
                    }
                  }
                }
              }
              pt = pt->cell_next;
            }
          } // end x
        }   // end y
      }     // end z
      E_SOLV += SOLVi;
    }

    F_VDWR[3 * ires] += F_VDWR[3 * ires + 1];
    G_VDWR[3 * ires] += G_VDWR[3 * ires + 1];
    F_VDWA[3 * ires] += F_VDWA[3 * ires + 1];
    G_VDWA[3 * ires] += G_VDWA[3 * ires + 1];
    F_SOLV[3 * ires] += F_SOLV[3 * ires + 1];
    G_SOLV[3 * ires] += G_SOLV[3 * ires + 1];

    tmp_g = gp.getPhi_UV(ires) ^ gp.getResidue(ires)->getCA()->r;
    DVDWR[3 * ires] =
        -(gp.getPhi_UV(ires) * F_VDWR[3 * ires]) - tmp_g * G_VDWR[3 * ires];
    DVDWA[3 * ires] =
        -(gp.getPhi_UV(ires) * F_VDWA[3 * ires]) - tmp_g * G_VDWA[3 * ires];
    DSOLV[3 * ires] =
        -(gp.getPhi_UV(ires) * F_SOLV[3 * ires]) - tmp_g * G_SOLV[3 * ires];
  }

  // cout << "Debug information" << endl;
  //   for(int i=0; i<gp.size(); i++){
  //   printf("R %i \t D= %f for F %f %f %f\n",
  //   i,DSOLV[3*i],F_SOLV[3*i].x,F_SOLV[3*i].y,F_SOLV[3*i].z); printf(" \t D=
  //   %f for F %f %f %f\n",
  //   DSOLV[3*i+1],F_SOLV[3*i+1].x,F_SOLV[3*i+1].y,F_SOLV[3*i+1].z); printf("
  //   \t D= %f for F %f %f %f\n",
  //   DSOLV[3*i+2],F_SOLV[3*i+2].x,F_SOLV[3*i+2].y,F_SOLV[3*i+2].z);
  //   }

  delete[] F_VDWR;
  delete[] G_VDWR;
  delete[] F_VDWA;
  delete[] G_VDWA;
  delete[] F_SOLV;
  delete[] G_SOLV;
}

/*calculate the F and G value of one atom in the Go trick
 * to be used in the getE_DEV_FYOChi_GOTRICK function*/
void gen_grid::getGF_GOTRICK(atom *theAtom, vec &F_VDWR, vec &F_VDWA,
                             vec &F_SOLV, vec &G_VDWR, vec &G_VDWA, vec &G_SOLV,
                             double &SOLVi, double &E_VDWR, double &E_VDWA) {

  F_VDWR = F_VDWA = F_SOLV = vec(0, 0, 0);
  G_VDWR = G_VDWA = G_SOLV = vec(0, 0, 0);
  SOLVi = E_VDWA = E_VDWR = 0;

  static double coef = 2.0 / (4.0 * PI * sqrt(PI));
  static double coef2 = 1.0 / (PI * sqrt(PI));
  int ti, tj;
  int ix, iy, iz, cindex, kx, ky, kz;
  //	int min_ix, min_iy, min_iz;
  double d2, d6, B, A, e, d, sig, eps, cutoff, slope, intercept;
  double xij, xji, de;
  atom *pt;
  int id = theAtom->resID;
  ti = theAtom->FF_t;
  ix = static_cast<int>((theAtom->r.x - min.x) / cell_box.x);
  iy = static_cast<int>((theAtom->r.y - min.y) / cell_box.y);
  iz = static_cast<int>((theAtom->r.z - min.z) / cell_box.z);
  for (int cz = iz - step; cz <= iz + step; cz++) {
    kz = cz;
    if (cz < 0)
      kz += ncell[2];
    else if (cz >= ncell[2])
      kz -= ncell[2];
    for (int cy = iy - step; cy <= iy + step; cy++) {
      ky = cy;
      if (cy < 0)
        ky += ncell[1];
      else if (cy >= ncell[1])
        ky -= ncell[1];
      for (int cx = ix - step; cx <= ix + step; cx++) {
        kx = cx;
        if (cx < 0)
          kx += ncell[0];
        else if (cx >= ncell[0])
          kx -= ncell[0];
        cindex = ncell[0] * (ncell[1] * kz + ky) + kx;
        pt = cell_atoms[cindex];
        while (pt) {
          int jres = pt->resID;
          if (!pt->isH()) {
            if (abs(id - jres) > 1 ||
                abs(id - jres) == 1 && theAtom->isSC() && pt->isSC() ||
                theAtom->isSC() && (jres == id + 1 && pt->PDBtype == _C_) ||
                theAtom->isSC() && (jres == id + 1 && pt->PDBtype == _O_) ||
                theAtom->isSC() && (jres == id - 1 && pt->PDBtype == _N_) ||
                pt->isSC() && (id == jres + 1 && theAtom->PDBtype == _C_) ||
                pt->isSC() && (id == jres + 1 && theAtom->PDBtype == _O_) ||
                pt->isSC() && (id == jres - 1 && theAtom->PDBtype == _N_))
            // VDW+S interaction within same residue are not counted
            // for neighbor residue, only SC-SC and
            //	SC-{BB atom outside the omega dihedral plan} are considered
            //	ie, R_i<->R_i+/-1, R_i<->C_i+1, R_i<->O_i+1, R_i<->N_i-1
            // VDW+S interaction for non-neighbor residues are counted
            {
              d2 = (pt->r - theAtom->r).mod2();
              d = sqrt(d2);
              tj = pt->FF_t;
              if (d2 < max_ir2) { // within the cutoff distance
                /**** van de Waals and solvation energy calculation *****/
                // vdw
                vdw.getVDWdata(ti, tj, B, A, eps, sig, cutoff, slope,
                               intercept);
                if (B && A) {
                  d6 = d2 * d2 * d2;
                  e = (B / d6 - A) / d6;
                  if (e > 0) {
                    // e = (sig-d)/sig*10.0;
                    if (d < cutoff) {
                      e = slope * d + intercept;
                      de = slope;
                    } else {
                      de = 6.0 * (A - 2.0 * B / d6) / (d6 * d);
                    }
                    E_VDWR += e;
                    // de = -10.0/sig;
                    F_VDWR += de / d * (theAtom->r ^ pt->r);
                    G_VDWR += de / d * (theAtom->r - pt->r);
                  } else {
                    E_VDWA += e;
                    de = 6.0 * (A - 2.0 * B / d6) / (d6 * d);
                    F_VDWA += de / d * (theAtom->r ^ pt->r);
                    G_VDWA += de / d * (theAtom->r - pt->r);
                  }
                }
                // solv
                int solv_cut = 0;
                if (d < (EEF1_SOLV[ti][CHARMM_RADIUS] +
                         EEF1_SOLV[tj][CHARMM_RADIUS]) *
                            0.846) {
                  d = (EEF1_SOLV[ti][CHARMM_RADIUS] +
                       EEF1_SOLV[tj][CHARMM_RADIUS]) *
                      0.846;
                  d2 = d * d;
                  solv_cut = 1;
                }
                xij = (d - EEF1_SOLV[ti][CHARMM_RADIUS]) /
                      EEF1_SOLV[ti][EEF1_LAMBDA];
                xji = (d - EEF1_SOLV[tj][CHARMM_RADIUS]) /
                      EEF1_SOLV[tj][EEF1_LAMBDA];
                SOLVi -=
                    (EEF1_SOLV[ti][EEF1_DG_FREE] * EEF1_SOLV[tj][EEF1_VOL] *
                         exp(-xij * xij) / EEF1_SOLV[ti][EEF1_LAMBDA] +
                     EEF1_SOLV[tj][EEF1_DG_FREE] * EEF1_SOLV[ti][EEF1_VOL] *
                         exp(-xji * xji) / EEF1_SOLV[tj][EEF1_LAMBDA]) *
                    coef / d2;
                if (!solv_cut) {
                  de = (EEF1_SOLV[ti][EEF1_DG_FREE] * EEF1_SOLV[tj][EEF1_VOL] *
                            exp(-xij * xij) *
                            (xij / EEF1_SOLV[ti][EEF1_LAMBDA] + 1.0 / d) /
                            EEF1_SOLV[ti][EEF1_LAMBDA] +
                        EEF1_SOLV[tj][EEF1_DG_FREE] * EEF1_SOLV[ti][EEF1_VOL] *
                            exp(-xji * xji) *
                            (xji / EEF1_SOLV[tj][EEF1_LAMBDA] + 1.0 / d) /
                            EEF1_SOLV[tj][EEF1_LAMBDA]) *
                       coef2 / d2;
                  F_SOLV += de / d * (theAtom->r ^ pt->r);
                  G_SOLV += de / d * (theAtom->r - pt->r);
                }
                /****end van de Waals and solvation energy calculation *****/

              } // end if within cutoff distance
            }   // end if abs(id-ires)>1
          }     // end if !pt is H
          pt = pt->cell_next;
        } // end while pt
      }   // end x
    }     // end y
  }       // end z
}

double gen_grid::getE_DEV_FYOChi_GOTRICK(double &E_VDWA, double &E_VDWR,
                                         double &E_SOLV, double *DVDWA,
                                         double *DVDWR, double *DSOLV,
                                         int nDihedral, int *INDEX_MAPPING) {
  E_VDWA = E_VDWR = E_SOLV = 0;
  residue *theRes = NULL;
  atom *theAtom = NULL;
  static double coef = 2.0 / (4.0 * PI * sqrt(PI));
  static double coef2 = 1.0 / (PI * sqrt(PI));
  int ti, tj;
  int ix, iy, iz, cindex, kx, ky, kz;
  int min_ix, min_iy, min_iz;
  double d2, d6, B, A, e, d, sig, eps, cutoff, slope, intercept;
  double SOLVi = 0, xij, xji, de;
  atom *pt;
  for (int i = 0; i < nDihedral; i++) {
    DVDWA[i] = 0;
    DVDWR[i] = 0;
    DSOLV[i] = 0;
  }
  gp.updateFY_UV();
  for (int i = 0; i < gp.size(); i++) {
    gp.getResidue(i)->updateChiUV();
  }
  vec *F_VDWA = new vec[nDihedral];
  vec *G_VDWA = new vec[nDihedral];
  vec *F_VDWR = new vec[nDihedral];
  vec *G_VDWR = new vec[nDihedral];
  vec *F_SOLV = new vec[nDihedral];
  vec *G_SOLV = new vec[nDihedral];
  vec f_VDWR(0, 0, 0), f_VDWA(0, 0, 0), f_SOLV(0, 0, 0);
  vec g_VDWR(0, 0, 0), g_VDWA(0, 0, 0), g_SOLV(0, 0, 0);
  double e_SOLV = 0, e_VDWA = 0, e_VDWR = 0;
  int iIndex, iIndexNext;
  for (int i = 0; i < nDihedral; i++) {
    F_VDWA[i] = vec(0, 0, 0);
    G_VDWA[i] = vec(0, 0, 0);
    F_VDWR[i] = vec(0, 0, 0);
    G_VDWR[i] = vec(0, 0, 0);
    F_SOLV[i] = vec(0, 0, 0);
    G_SOLV[i] = vec(0, 0, 0);
  }
  vector<atom *> tmp_va;
  for (int ires = gp.size() - 1; ires >= 0; ires--) {
    theRes = gp.getResidue(ires);
    int NChi = theRes->getNChi();
    int type = theRes->getType();

    // omiga: Hi+1, CAi+1(2)
    iIndex =
        INDEX_MAPPING[3 * ires + 2]; // the index of omega bond in the array
    if (ires < gp.size() - 1) {
      theAtom = gp.getResidue(ires + 1)->getCA();
      getGF_GOTRICK(theAtom, f_VDWR, f_VDWA, f_SOLV, g_VDWR, g_VDWA, g_SOLV,
                    e_SOLV, e_VDWR, e_VDWA);
      E_SOLV += e_SOLV;
      E_VDWR += e_VDWR;
      E_VDWA += e_VDWA;
      iIndexNext =
          INDEX_MAPPING[3 *
                        (ires + 1)]; // the next phi bond index; add recurrently
      F_VDWR[iIndex] += F_VDWR[iIndexNext] + f_VDWR;
      G_VDWR[iIndex] += G_VDWR[iIndexNext] + g_VDWR;
      F_VDWA[iIndex] += F_VDWA[iIndexNext] + f_VDWA;
      G_VDWA[iIndex] += G_VDWA[iIndexNext] + g_VDWA;
      F_SOLV[iIndex] += F_SOLV[iIndexNext] + f_SOLV;
      G_SOLV[iIndex] += G_SOLV[iIndexNext] + g_SOLV;
      vec tmp_g =
          gp.getOmiga_UV(ires) ^
          gp.getResidue(ires + 1)->getN()->r; // omega is C_i->N_i+1 bond
      DVDWR[iIndex] =
          -(gp.getOmiga_UV(ires) * F_VDWR[iIndex]) - tmp_g * G_VDWR[iIndex];
      DVDWA[iIndex] =
          -(gp.getOmiga_UV(ires) * F_VDWA[iIndex]) - tmp_g * G_VDWA[iIndex];
      DSOLV[iIndex] =
          -(gp.getOmiga_UV(ires) * F_SOLV[iIndex]) - tmp_g * G_SOLV[iIndex];
    } // end if not last residue

    // psi: Oi, Ni+1(1)
    iIndex = INDEX_MAPPING[3 * ires + 1]; // the index of psi bond in the array
    iIndexNext = iIndex + 1; // the next bond, omega in the same residue
    tmp_va.clear();
    tmp_va.push_back(theRes->getO());
    if (ires < gp.size() - 1) {
      tmp_va.push_back(gp.getResidue(ires + 1)->getN());
    }
    for (int iatom = 0; iatom < tmp_va.size(); iatom++) {
      theAtom = tmp_va[iatom];
      getGF_GOTRICK(theAtom, f_VDWR, f_VDWA, f_SOLV, g_VDWR, g_VDWA, g_SOLV,
                    e_SOLV, e_VDWR, e_VDWA);
      E_SOLV += e_SOLV;
      E_VDWR += e_VDWR;
      E_VDWA += e_VDWA;
      F_VDWR[iIndex] += f_VDWR;
      G_VDWR[iIndex] += g_VDWR;
      F_VDWA[iIndex] += f_VDWA;
      G_VDWA[iIndex] += g_VDWA;
      F_SOLV[iIndex] += f_SOLV;
      G_SOLV[iIndex] += g_SOLV;
    }
    // add recurrently
    F_VDWR[iIndex] += F_VDWR[iIndexNext];
    G_VDWR[iIndex] += G_VDWR[iIndexNext];
    F_VDWA[iIndex] += F_VDWA[iIndexNext];
    G_VDWA[iIndex] += G_VDWA[iIndexNext];
    F_SOLV[iIndex] += F_SOLV[iIndexNext];
    G_SOLV[iIndex] += G_SOLV[iIndexNext];
    vec tmp_g = gp.getPsi_UV(ires) ^ gp.getResidue(ires)->getC()->r; // CA_i-C_i
    DVDWR[iIndex] =
        -(gp.getPsi_UV(ires) * F_VDWR[iIndex]) - tmp_g * G_VDWR[iIndex];
    DVDWA[iIndex] =
        -(gp.getPsi_UV(ires) * F_VDWA[iIndex]) - tmp_g * G_VDWA[iIndex];
    DSOLV[iIndex] =
        -(gp.getPsi_UV(ires) * F_SOLV[iIndex]) - tmp_g * G_SOLV[iIndex];

    // phi: Ci, CB_i
    iIndex = INDEX_MAPPING[3 * ires]; // the index of phi bond in the array
    // for CB, the sidechain bond should be calculated first
    /*** this part is for backbone chain, C_i atom is included ***/
    iIndexNext = iIndex + 1; // the next bond, psi in the same residue
    tmp_va.clear();
    tmp_va.push_back(theRes->getC()); // C_i atom
    if (theAtom = theRes->getCB()) {  // CB_i atom
      tmp_va.push_back(theAtom);
    }
    for (int iatom = 0; iatom < tmp_va.size(); iatom++) {
      theAtom = tmp_va[iatom];
      getGF_GOTRICK(theAtom, f_VDWR, f_VDWA, f_SOLV, g_VDWR, g_VDWA, g_SOLV,
                    e_SOLV, e_VDWR, e_VDWA);
      E_SOLV += e_SOLV;
      E_VDWR += e_VDWR;
      E_VDWA += e_VDWA;
      F_VDWR[iIndex] += f_VDWR;
      G_VDWR[iIndex] += g_VDWR;
      F_VDWA[iIndex] += f_VDWA;
      G_VDWA[iIndex] += g_VDWA;
      F_SOLV[iIndex] += f_SOLV;
      G_SOLV[iIndex] += g_SOLV;
    }
    // add up from the BB branch
    F_VDWR[iIndex] += F_VDWR[iIndexNext];
    G_VDWR[iIndex] += G_VDWR[iIndexNext];
    F_VDWA[iIndex] += F_VDWA[iIndexNext];
    G_VDWA[iIndex] += G_VDWA[iIndexNext];
    F_SOLV[iIndex] += F_SOLV[iIndexNext];
    G_SOLV[iIndex] += G_SOLV[iIndexNext];
    /** now go to the side chain **/
    int *scAtomDivN = theRes->getscAtomDivN();
    atom ***scAtomDiv = theRes->getscAtomDiv();
    /*** the following operation only applies for side chain rotamers ***/
    // for the end dihedral of the side chain
    int j = NChi - 1;
    if (j >= 0) {
      int jIndex = INDEX_MAPPING[3 * ires] + 3 + j;
      for (int jj = 0; jj < scAtomDivN[j];
           jj++) { //  for each atoms belong to the bond rotation set
        theAtom = scAtomDiv[j][jj];
        getGF_GOTRICK(theAtom, f_VDWR, f_VDWA, f_SOLV, g_VDWR, g_VDWA, g_SOLV,
                      e_SOLV, e_VDWR, e_VDWA);
        E_SOLV += e_SOLV;
        E_VDWR += e_VDWR;
        E_VDWA += e_VDWA;
        F_VDWR[jIndex] += f_VDWR;
        G_VDWR[jIndex] += g_VDWR;
        F_VDWA[jIndex] += f_VDWA;
        G_VDWA[jIndex] += g_VDWA;
        F_SOLV[jIndex] += f_SOLV;
        G_SOLV[jIndex] += g_SOLV;
      }
      vec temp_eb = *theRes->getChiUV(j);
      vec tmp_g = temp_eb ^ (*theRes->getChiP2(j)); //
      DVDWR[jIndex] = -(temp_eb * F_VDWR[jIndex]) - tmp_g * G_VDWR[jIndex];
      DVDWA[jIndex] = -(temp_eb * F_VDWA[jIndex]) - tmp_g * G_VDWA[jIndex];
      DSOLV[jIndex] = -(temp_eb * F_SOLV[jIndex]) - tmp_g * G_SOLV[jIndex];
    }
    // otherwise the F and G from the next bond in the sidechain has to be added
    for (int j = NChi - 2; j >= 0;
         j--) { // for each side chain dihedral, except the last
      int jIndex = INDEX_MAPPING[3 * ires] + 3 + j;
      int jIndexNext = jIndex + 1;
      for (int jj = 0; jj < scAtomDivN[j];
           jj++) { //  for each atoms belong to the bond rotation set
        theAtom = scAtomDiv[j][jj];
        getGF_GOTRICK(theAtom, f_VDWR, f_VDWA, f_SOLV, g_VDWR, g_VDWA, g_SOLV,
                      e_SOLV, e_VDWR, e_VDWA);
        E_SOLV += e_SOLV;
        E_VDWR += e_VDWR;
        E_VDWA += e_VDWA;
        F_VDWR[jIndex] += f_VDWR;
        G_VDWR[jIndex] += g_VDWR;
        F_VDWA[jIndex] += f_VDWA;
        G_VDWA[jIndex] += g_VDWA;
        F_SOLV[jIndex] += f_SOLV;
        G_SOLV[jIndex] += g_SOLV;
      }
      F_VDWR[jIndex] += F_VDWR[jIndexNext];
      G_VDWR[jIndex] += G_VDWR[jIndexNext];
      F_VDWA[jIndex] += F_VDWA[jIndexNext];
      G_VDWA[jIndex] += G_VDWA[jIndexNext];
      F_SOLV[jIndex] += F_SOLV[jIndexNext];
      G_SOLV[jIndex] += G_SOLV[jIndexNext];
      vec temp_eb = *theRes->getChiUV(j);
      vec tmp_g = temp_eb ^ (*theRes->getChiP2(j)); //
      DVDWR[jIndex] = -(temp_eb * F_VDWR[jIndex]) - tmp_g * G_VDWR[jIndex];
      DVDWA[jIndex] = -(temp_eb * F_VDWA[jIndex]) - tmp_g * G_VDWA[jIndex];
      DSOLV[jIndex] = -(temp_eb * F_SOLV[jIndex]) - tmp_g * G_SOLV[jIndex];
    }
    // add the sidechain contribution the the phi rotation
    if (NChi > 0) {
      iIndexNext = iIndex + 3; // the next sidechain bond
      F_VDWR[iIndex] += F_VDWR[iIndexNext];
      G_VDWR[iIndex] += G_VDWR[iIndexNext];
      F_VDWA[iIndex] += F_VDWA[iIndexNext];
      G_VDWA[iIndex] += G_VDWA[iIndexNext];
      F_SOLV[iIndex] += F_SOLV[iIndexNext];
      G_SOLV[iIndex] += G_SOLV[iIndexNext];
    }
    //  calculation Derivatives with F and G
    tmp_g = gp.getPhi_UV(ires) ^ gp.getResidue(ires)->getCA()->r; //
    DVDWR[iIndex] =
        -(gp.getPhi_UV(ires) * F_VDWR[iIndex]) - tmp_g * G_VDWR[iIndex];
    DVDWA[iIndex] =
        -(gp.getPhi_UV(ires) * F_VDWA[iIndex]) - tmp_g * G_VDWA[iIndex];
    DSOLV[iIndex] =
        -(gp.getPhi_UV(ires) * F_SOLV[iIndex]) - tmp_g * G_SOLV[iIndex];
  } // end for each residue

  // cout << "Debug information" << endl;
  //   for(int i=0; i<gp.size(); i++){
  //   printf("R %i \t D= %f for F %f %f %f\n",
  //   i,DSOLV[3*i],F_SOLV[3*i].x,F_SOLV[3*i].y,F_SOLV[3*i].z); printf(" \t D=
  //   %f for F %f %f %f\n",
  //   DSOLV[3*i+1],F_SOLV[3*i+1].x,F_SOLV[3*i+1].y,F_SOLV[3*i+1].z); printf("
  //   \t D= %f for F %f %f %f\n",
  //   DSOLV[3*i+2],F_SOLV[3*i+2].x,F_SOLV[3*i+2].y,F_SOLV[3*i+2].z);
  //   }

  delete[] F_VDWR;
  delete[] G_VDWR;
  delete[] F_VDWA;
  delete[] G_VDWA;
  delete[] F_SOLV;
  delete[] G_SOLV;
}

double gen_grid::getE(double &E_VDWA, double &E_VDWR, double &E_SOLV) {
  E_VDWA = E_VDWR = E_SOLV = 0;
  residue *theRes = NULL;
  atom *theAtom = NULL;
  static double coef = 2.0 / (4.0 * PI * sqrt(PI));
  static double coef2 = 1.0 / (PI * sqrt(PI));
  int ti, tj;
  int ix, iy, iz, cindex, kx, ky, kz;
  int min_ix, min_iy, min_iz;
  double d2, d6, B, A, e, d, sig, eps, cutoff, slope, intercept;
  double SOLVi = 0, xij, xji;
  atom *pt;

  gp.updateFY_UV();
  for (int ires = 0; ires < gp.size(); ires++) {
    theRes = gp.getResidue(ires);
    for (int iatom = 0; iatom < theRes->getNA(); iatom++) {
      theAtom = theRes->getAtom(iatom);
      ti = theAtom->FF_t;
      ix = static_cast<int>((theAtom->r.x - min.x) / cell_box.x);
      iy = static_cast<int>((theAtom->r.y - min.y) / cell_box.y);
      iz = static_cast<int>((theAtom->r.z - min.z) / cell_box.z);
      SOLVi = 0;
      for (int cz = iz - step; cz <= iz + step; cz++) {
        kz = cz;
        if (cz < 0)
          kz += ncell[2];
        else if (cz >= ncell[2])
          kz -= ncell[2];

        for (int cy = iy - step; cy <= iy + step; cy++) {
          ky = cy;
          if (cy < 0)
            ky += ncell[1];
          else if (cy >= ncell[1])
            ky -= ncell[1];

          for (int cx = ix - step; cx <= ix + step; cx++) {
            kx = cx;
            if (cx < 0)
              kx += ncell[0];
            else if (cx >= ncell[0])
              kx -= ncell[0];

            cindex = ncell[0] * (ncell[1] * kz + ky) + kx;
            pt = cell_atoms[cindex];
            while (pt) {
              int jres = pt->resID;
              vec Rsm = pt->r - theAtom->r;
              if (jres > ires + 1 && !pt->isH()) {
                d2 = (pt->r - theAtom->r).mod2();
                d = sqrt(d2);
                tj = pt->FF_t;
                if (d2 < max_ir2) { // within the cutoff distance
                  vdw.getVDWdata(ti, tj, B, A, eps, sig, cutoff, slope,
                                 intercept);
                  if (B && A) {
                    d6 = d2 * d2 * d2;
                    e = (B / d6 - A) / d6;
                    if (e > 0) {
                      // e = (sig-d)/sig*10.0;
                      if (d < cutoff) {
                        e = slope * d + intercept;
                      }
                      E_VDWR += e;
                    } else {
                      E_VDWA += e;
                    }
                  }
                  // solvation
                  int solv_cut = 0;
                  if (d < (EEF1_SOLV[ti][CHARMM_RADIUS] +
                           EEF1_SOLV[tj][CHARMM_RADIUS]) *
                              0.846) {
                    d = (EEF1_SOLV[ti][CHARMM_RADIUS] +
                         EEF1_SOLV[tj][CHARMM_RADIUS]) *
                        0.846;
                    d2 = d * d;
                    solv_cut = 1;
                  }
                  xij = (d - EEF1_SOLV[ti][CHARMM_RADIUS]) /
                        EEF1_SOLV[ti][EEF1_LAMBDA];
                  xji = (d - EEF1_SOLV[tj][CHARMM_RADIUS]) /
                        EEF1_SOLV[tj][EEF1_LAMBDA];
                  SOLVi -=
                      (EEF1_SOLV[ti][EEF1_DG_FREE] * EEF1_SOLV[tj][EEF1_VOL] *
                           exp(-xij * xij) / EEF1_SOLV[ti][EEF1_LAMBDA] +
                       EEF1_SOLV[tj][EEF1_DG_FREE] * EEF1_SOLV[ti][EEF1_VOL] *
                           exp(-xji * xji) / EEF1_SOLV[tj][EEF1_LAMBDA]) *
                      coef / d2;
                }
              } else if (jres == ires + 1 && !pt->isH()) {
                // at first do not account for the short range bb interactions
                // since, the short range backbone interaction determines the
                // phi,psi, we assume that the phi,psi angles do not have large
                // deviations.
                if ((theAtom->isSC() && (pt->isSC() || pt->PDBtype == _C_ ||
                                         pt->PDBtype == _O_)) ||
                    (pt->isSC() &&
                     (theAtom->isSC() || theAtom->PDBtype == _N_))) {
                  d2 = (pt->r - theAtom->r).mod2();
                  d = sqrt(d2);
                  tj = pt->FF_t;
                  if (d2 < max_ir2) { // within the cutoff distance
                    vdw.getVDWdata(ti, tj, B, A, eps, sig, cutoff, slope,
                                   intercept);
                    if (B && A) {
                      d6 = d2 * d2 * d2;
                      e = (B / d6 - A) / d6;
                      if (e > 0) {
                        // e = (sig-d)/sig*10.0;
                        if (d < cutoff) {
                          e = slope * d + intercept;
                        }
                        E_VDWR += e;
                      } else {
                        E_VDWA += e;
                      }
                    }
                    // solvation
                    int solv_cut = 0;
                    if (d < (EEF1_SOLV[ti][CHARMM_RADIUS] +
                             EEF1_SOLV[tj][CHARMM_RADIUS]) *
                                0.846) {
                      d = (EEF1_SOLV[ti][CHARMM_RADIUS] +
                           EEF1_SOLV[tj][CHARMM_RADIUS]) *
                          0.846;
                      d2 = d * d;
                      solv_cut = 1;
                    }
                    xij = (d - EEF1_SOLV[ti][CHARMM_RADIUS]) /
                          EEF1_SOLV[ti][EEF1_LAMBDA];
                    xji = (d - EEF1_SOLV[tj][CHARMM_RADIUS]) /
                          EEF1_SOLV[tj][EEF1_LAMBDA];
                    SOLVi -=
                        (EEF1_SOLV[ti][EEF1_DG_FREE] * EEF1_SOLV[tj][EEF1_VOL] *
                             exp(-xij * xij) / EEF1_SOLV[ti][EEF1_LAMBDA] +
                         EEF1_SOLV[tj][EEF1_DG_FREE] * EEF1_SOLV[ti][EEF1_VOL] *
                             exp(-xji * xji) / EEF1_SOLV[tj][EEF1_LAMBDA]) *
                        coef / d2;
                  }
                }
              }
              pt = pt->cell_next;
            }
          } // x
        }   // y
      }     // z
      E_SOLV += SOLVi;
    }
  }
}

void gen_grid::getDE_BACKRUB(const int ires, double &E_VDWA, double &E_VDWR,
                             double &E_SOLV, double *DVDWA, double *DVDWR,
                             double *DSOLV, vec *backrub_uv) {
  E_VDWA = E_VDWR = E_SOLV = 0;
  atom *theAtom = NULL;
  residue *theRes = gp.getResidue(ires);
  static double coef = 2.0 / (4.0 * PI * sqrt(PI));
  static double coef2 = 1.0 / (PI * sqrt(PI));
  int ti, tj;
  int ix, iy, iz, cindex, kx, ky, kz;
  double d2, d6, B, A, e, d, sig, eps, cutoff, slope, intercept;
  double SOLVi = 0, xij, xji, de;
  atom *pt;
  // clear the array
  for (int i = 0; i < 3; i++) {
    DVDWA[i] = DVDWR[i] = DSOLV[i] = 0;
  }
  // for residue ires
  for (int i = 0; i < theRes->getNA(); i++) {
    theAtom = theRes->getAtom(i);
    ti = theAtom->FF_t;
    ix = static_cast<int>((theAtom->r.x - min.x) / cell_box.x);
    iy = static_cast<int>((theAtom->r.y - min.y) / cell_box.y);
    iz = static_cast<int>((theAtom->r.z - min.z) / cell_box.z);

    SOLVi = 0;
    for (int cz = iz - step; cz <= iz + step; cz++) {
      kz = cz;
      if (cz < 0)
        kz += ncell[2];
      else if (cz >= ncell[2])
        kz -= ncell[2];

      for (int cy = iy - step; cy <= iy + step; cy++) {
        ky = cy;
        if (cy < 0)
          ky += ncell[1];
        else if (cy >= ncell[1])
          ky -= ncell[1];

        for (int cx = ix - step; cx <= ix + step; cx++) {
          kx = cx;
          if (cx < 0)
            kx += ncell[0];
          else if (cx >= ncell[0])
            kx -= ncell[0];
          cindex = ncell[0] * (ncell[1] * kz + ky) + kx;

          pt = cell_atoms[cindex];
          while (pt) {
            if (pt != theAtom &&
                !pt->isH()) { // not the same atom and not protons. Proton is
                              // not accounted for VDW and SOLV
              d2 = (pt->r - theAtom->r).mod2();
              d = sqrt(d2);
              tj = pt->FF_t;
              if (d2 < max_ir2) { // within the cutoff distance
                // vdw
                vec Rsm = theAtom->r - pt->r;
                if (!isInternal(
                        theAtom,
                        pt)) { // vdw does not include the internal pairs
                  vdw.getVDWdata(ti, tj, B, A, eps, sig, cutoff, slope,
                                 intercept);
                  if (B && A) {
                    d6 = d2 * d2 * d2;
                    e = (B / d6 - A) / d6;
                    if (e > 0) {
                      if (d < cutoff) {
                        e = slope * d + intercept;
                        de = slope;
                      } else {
                        de = 6.0 * (A - 2.0 * B / d6) / (d6 * d);
                      }
                      E_VDWR += e;
                      // derive derivative
                      // for theta13
                      DVDWR[0] += de *
                                  ((backrub_uv[0] ^
                                    (theAtom->r -
                                     gp.getResidue(ires + 1)->getCA()->r)) *
                                   Rsm) /
                                  d;
                      if (!theAtom->isSC()) {
                        if (theAtom->PDBtype == _N_) { // theta12
                          DVDWR[1] += de *
                                      ((backrub_uv[1] ^
                                        (theAtom->r -
                                         gp.getResidue(ires)->getCA()->r)) *
                                       Rsm) /
                                      d;
                        } else if (theAtom->PDBtype == _C_ ||
                                   theAtom->PDBtype == _O_) { // theta32
                          DVDWR[2] += de *
                                      ((backrub_uv[2] ^
                                        (theAtom->r -
                                         gp.getResidue(ires)->getCA()->r)) *
                                       Rsm) /
                                      d;
                        }
                      }
                    }

                    else {
                      E_VDWA += e;
                      de = 6.0 * (A - 2.0 * B / d6) / (d6 * d);
                      // derive derivative
                      DVDWA[0] += de *
                                  ((backrub_uv[0] ^
                                    (theAtom->r -
                                     gp.getResidue(ires + 1)->getCA()->r)) *
                                   Rsm) /
                                  d;
                      if (!theAtom->isSC()) {
                        if (theAtom->PDBtype == _N_) { // theta12
                          DVDWA[1] += de *
                                      ((backrub_uv[1] ^
                                        (theAtom->r -
                                         gp.getResidue(ires)->getCA()->r)) *
                                       Rsm) /
                                      d;
                        } else if (theAtom->PDBtype == _C_ ||
                                   theAtom->PDBtype == _O_) { // theta32
                          DVDWA[2] += de *
                                      ((backrub_uv[2] ^
                                        (theAtom->r -
                                         gp.getResidue(ires)->getCA()->r)) *
                                       Rsm) /
                                      d;
                        }
                      }
                    }
                  }

                  // solvation
                  int solv_cut = 0;
                  if (d < (EEF1_SOLV[ti][CHARMM_RADIUS] +
                           EEF1_SOLV[tj][CHARMM_RADIUS]) *
                              0.846) {
                    d = (EEF1_SOLV[ti][CHARMM_RADIUS] +
                         EEF1_SOLV[tj][CHARMM_RADIUS]) *
                        0.846;
                    d2 = d * d;
                    solv_cut = 1;
                  }
                  xij = (d - EEF1_SOLV[ti][CHARMM_RADIUS]) /
                        EEF1_SOLV[ti][EEF1_LAMBDA];
                  xji = (d - EEF1_SOLV[tj][CHARMM_RADIUS]) /
                        EEF1_SOLV[tj][EEF1_LAMBDA];

                  double solv_ij =
                      (EEF1_SOLV[ti][EEF1_DG_FREE] * EEF1_SOLV[tj][EEF1_VOL] *
                           exp(-xij * xij) / EEF1_SOLV[ti][EEF1_LAMBDA] +
                       EEF1_SOLV[tj][EEF1_DG_FREE] * EEF1_SOLV[ti][EEF1_VOL] *
                           exp(-xji * xji) / EEF1_SOLV[tj][EEF1_LAMBDA]) *
                      coef / d2;
                  if (!solv_cut) {
                    de = (EEF1_SOLV[ti][EEF1_DG_FREE] *
                              EEF1_SOLV[tj][EEF1_VOL] * exp(-xij * xij) *
                              (xij / EEF1_SOLV[ti][EEF1_LAMBDA] + 1.0 / d) /
                              EEF1_SOLV[ti][EEF1_LAMBDA] +
                          EEF1_SOLV[tj][EEF1_DG_FREE] *
                              EEF1_SOLV[ti][EEF1_VOL] * exp(-xji * xji) *
                              (xji / EEF1_SOLV[tj][EEF1_LAMBDA] + 1.0 / d) /
                              EEF1_SOLV[tj][EEF1_LAMBDA]) *
                         coef2 / d2;
                    DSOLV[0] +=
                        de *
                        ((backrub_uv[0] ^
                          (theAtom->r - gp.getResidue(ires + 1)->getCA()->r)) *
                         Rsm) /
                        d;
                    if (!theAtom->isSC()) {
                      if (theAtom->PDBtype == _N_) { // theta12
                        DSOLV[1] +=
                            de *
                            ((backrub_uv[1] ^
                              (theAtom->r - gp.getResidue(ires)->getCA()->r)) *
                             Rsm) /
                            d;
                      } else if (theAtom->PDBtype == _C_ ||
                                 theAtom->PDBtype == _O_) { // theta32
                        DSOLV[2] +=
                            de *
                            ((backrub_uv[2] ^
                              (theAtom->r - gp.getResidue(ires)->getCA()->r)) *
                             Rsm) /
                            d;
                      }
                    }
                  }
                  SOLVi -= solv_ij;
                }
              }
            }
            pt = pt->cell_next;
          }
        }
      }
    }
    E_SOLV += SOLVi;
  }
  // for atom C,O of ires-1
  vector<atom *> tmp_atoms;
  tmp_atoms.push_back(gp.getResidue(ires - 1)->getC());
  tmp_atoms.push_back(gp.getResidue(ires - 1)->getO());
  for (int i = 0; i < tmp_atoms.size(); i++) {
    theAtom = tmp_atoms[i];
    ti = theAtom->FF_t;
    ix = static_cast<int>((theAtom->r.x - min.x) / cell_box.x);
    iy = static_cast<int>((theAtom->r.y - min.y) / cell_box.y);
    iz = static_cast<int>((theAtom->r.z - min.z) / cell_box.z);

    SOLVi = 0;
    for (int cz = iz - step; cz <= iz + step; cz++) {
      kz = cz;
      if (cz < 0)
        kz += ncell[2];
      else if (cz >= ncell[2])
        kz -= ncell[2];

      for (int cy = iy - step; cy <= iy + step; cy++) {
        ky = cy;
        if (cy < 0)
          ky += ncell[1];
        else if (cy >= ncell[1])
          ky -= ncell[1];

        for (int cx = ix - step; cx <= ix + step; cx++) {
          kx = cx;
          if (cx < 0)
            kx += ncell[0];
          else if (cx >= ncell[0])
            kx -= ncell[0];
          cindex = ncell[0] * (ncell[1] * kz + ky) + kx;

          pt = cell_atoms[cindex];
          while (pt) {
            if (pt != theAtom &&
                !pt->isH()) { // not the same atom and not protons. Proton is
                              // not accounted for VDW and SOLV
              d2 = (pt->r - theAtom->r).mod2();
              d = sqrt(d2);
              tj = pt->FF_t;
              if (d2 < max_ir2) { // within the cutoff distance
                // vdw
                vec Rsm = theAtom->r - pt->r;
                if (!isInternal(
                        theAtom,
                        pt)) { // vdw does not include the internal pairs
                  vdw.getVDWdata(ti, tj, B, A, eps, sig, cutoff, slope,
                                 intercept);
                  if (B && A) {
                    d6 = d2 * d2 * d2;
                    e = (B / d6 - A) / d6;
                    if (e > 0) {
                      if (d < cutoff) {
                        e = slope * d + intercept;
                        de = slope;
                      } else {
                        de = 6.0 * (A - 2.0 * B / d6) / (d6 * d);
                      }
                      E_VDWR += e;
                      // derive derivative
                      // for theta13
                      DVDWR[0] += de *
                                  ((backrub_uv[0] ^
                                    (theAtom->r -
                                     gp.getResidue(ires + 1)->getCA()->r)) *
                                   Rsm) /
                                  d;
                      DVDWR[1] +=
                          de *
                          ((backrub_uv[1] ^
                            (theAtom->r - gp.getResidue(ires)->getCA()->r)) *
                           Rsm) /
                          d;
                    }

                    else {
                      E_VDWA += e;
                      de = 6.0 * (A - 2.0 * B / d6) / (d6 * d);
                      // derive derivative
                      DVDWA[0] += de *
                                  ((backrub_uv[0] ^
                                    (theAtom->r -
                                     gp.getResidue(ires + 1)->getCA()->r)) *
                                   Rsm) /
                                  d;
                      DVDWA[1] +=
                          de *
                          ((backrub_uv[1] ^
                            (theAtom->r - gp.getResidue(ires)->getCA()->r)) *
                           Rsm) /
                          d;
                    }
                  }

                  // solvation
                  int solv_cut = 0;
                  if (d < (EEF1_SOLV[ti][CHARMM_RADIUS] +
                           EEF1_SOLV[tj][CHARMM_RADIUS]) *
                              0.846) {
                    d = (EEF1_SOLV[ti][CHARMM_RADIUS] +
                         EEF1_SOLV[tj][CHARMM_RADIUS]) *
                        0.846;
                    d2 = d * d;
                    solv_cut = 1;
                  }
                  xij = (d - EEF1_SOLV[ti][CHARMM_RADIUS]) /
                        EEF1_SOLV[ti][EEF1_LAMBDA];
                  xji = (d - EEF1_SOLV[tj][CHARMM_RADIUS]) /
                        EEF1_SOLV[tj][EEF1_LAMBDA];

                  double solv_ij =
                      (EEF1_SOLV[ti][EEF1_DG_FREE] * EEF1_SOLV[tj][EEF1_VOL] *
                           exp(-xij * xij) / EEF1_SOLV[ti][EEF1_LAMBDA] +
                       EEF1_SOLV[tj][EEF1_DG_FREE] * EEF1_SOLV[ti][EEF1_VOL] *
                           exp(-xji * xji) / EEF1_SOLV[tj][EEF1_LAMBDA]) *
                      coef / d2;
                  if (!solv_cut) {
                    de = (EEF1_SOLV[ti][EEF1_DG_FREE] *
                              EEF1_SOLV[tj][EEF1_VOL] * exp(-xij * xij) *
                              (xij / EEF1_SOLV[ti][EEF1_LAMBDA] + 1.0 / d) /
                              EEF1_SOLV[ti][EEF1_LAMBDA] +
                          EEF1_SOLV[tj][EEF1_DG_FREE] *
                              EEF1_SOLV[ti][EEF1_VOL] * exp(-xji * xji) *
                              (xji / EEF1_SOLV[tj][EEF1_LAMBDA] + 1.0 / d) /
                              EEF1_SOLV[tj][EEF1_LAMBDA]) *
                         coef2 / d2;
                    DSOLV[0] +=
                        de *
                        ((backrub_uv[0] ^
                          (theAtom->r - gp.getResidue(ires + 1)->getCA()->r)) *
                         Rsm) /
                        d;
                    DSOLV[1] +=
                        de *
                        ((backrub_uv[1] ^
                          (theAtom->r - gp.getResidue(ires)->getCA()->r)) *
                         Rsm) /
                        d;
                  }
                  SOLVi -= solv_ij;
                }
              }
            }
            pt = pt->cell_next;
          }
        }
      }
    }
    E_SOLV += SOLVi;
  }

  // for atom N of ires+1
  {
    theAtom = gp.getResidue(ires + 1)->getN();
    ti = theAtom->FF_t;
    ix = static_cast<int>((theAtom->r.x - min.x) / cell_box.x);
    iy = static_cast<int>((theAtom->r.y - min.y) / cell_box.y);
    iz = static_cast<int>((theAtom->r.z - min.z) / cell_box.z);

    SOLVi = 0;
    for (int cz = iz - step; cz <= iz + step; cz++) {
      kz = cz;
      if (cz < 0)
        kz += ncell[2];
      else if (cz >= ncell[2])
        kz -= ncell[2];

      for (int cy = iy - step; cy <= iy + step; cy++) {
        ky = cy;
        if (cy < 0)
          ky += ncell[1];
        else if (cy >= ncell[1])
          ky -= ncell[1];

        for (int cx = ix - step; cx <= ix + step; cx++) {
          kx = cx;
          if (cx < 0)
            kx += ncell[0];
          else if (cx >= ncell[0])
            kx -= ncell[0];
          cindex = ncell[0] * (ncell[1] * kz + ky) + kx;

          pt = cell_atoms[cindex];
          while (pt) {
            if (pt != theAtom &&
                !pt->isH()) { // not the same atom and not protons. Proton is
                              // not accounted for VDW and SOLV
              d2 = (pt->r - theAtom->r).mod2();
              d = sqrt(d2);
              tj = pt->FF_t;
              if (d2 < max_ir2) { // within the cutoff distance
                // vdw
                vec Rsm = theAtom->r - pt->r;
                if (!isInternal(
                        theAtom,
                        pt)) { // vdw does not include the internal pairs
                  vdw.getVDWdata(ti, tj, B, A, eps, sig, cutoff, slope,
                                 intercept);
                  if (B && A) {
                    d6 = d2 * d2 * d2;
                    e = (B / d6 - A) / d6;
                    if (e > 0) {
                      if (d < cutoff) {
                        e = slope * d + intercept;
                        de = slope;
                      } else {
                        de = 6.0 * (A - 2.0 * B / d6) / (d6 * d);
                      }
                      E_VDWR += e;
                      // derive derivative
                      // for theta13
                      DVDWR[0] += de *
                                  ((backrub_uv[0] ^
                                    (theAtom->r -
                                     gp.getResidue(ires + 1)->getCA()->r)) *
                                   Rsm) /
                                  d;
                      DVDWR[2] +=
                          de *
                          ((backrub_uv[2] ^
                            (theAtom->r - gp.getResidue(ires)->getCA()->r)) *
                           Rsm) /
                          d;
                    }

                    else {
                      E_VDWA += e;
                      de = 6.0 * (A - 2.0 * B / d6) / (d6 * d);
                      // derive derivative
                      DVDWA[0] += de *
                                  ((backrub_uv[0] ^
                                    (theAtom->r -
                                     gp.getResidue(ires + 1)->getCA()->r)) *
                                   Rsm) /
                                  d;
                      DVDWA[2] +=
                          de *
                          ((backrub_uv[2] ^
                            (theAtom->r - gp.getResidue(ires)->getCA()->r)) *
                           Rsm) /
                          d;
                    }
                  }

                  // solvation
                  int solv_cut = 0;
                  if (d < (EEF1_SOLV[ti][CHARMM_RADIUS] +
                           EEF1_SOLV[tj][CHARMM_RADIUS]) *
                              0.846) {
                    d = (EEF1_SOLV[ti][CHARMM_RADIUS] +
                         EEF1_SOLV[tj][CHARMM_RADIUS]) *
                        0.846;
                    d2 = d * d;
                    solv_cut = 1;
                  }
                  xij = (d - EEF1_SOLV[ti][CHARMM_RADIUS]) /
                        EEF1_SOLV[ti][EEF1_LAMBDA];
                  xji = (d - EEF1_SOLV[tj][CHARMM_RADIUS]) /
                        EEF1_SOLV[tj][EEF1_LAMBDA];

                  double solv_ij =
                      (EEF1_SOLV[ti][EEF1_DG_FREE] * EEF1_SOLV[tj][EEF1_VOL] *
                           exp(-xij * xij) / EEF1_SOLV[ti][EEF1_LAMBDA] +
                       EEF1_SOLV[tj][EEF1_DG_FREE] * EEF1_SOLV[ti][EEF1_VOL] *
                           exp(-xji * xji) / EEF1_SOLV[tj][EEF1_LAMBDA]) *
                      coef / d2;
                  if (!solv_cut) {
                    de = (EEF1_SOLV[ti][EEF1_DG_FREE] *
                              EEF1_SOLV[tj][EEF1_VOL] * exp(-xij * xij) *
                              (xij / EEF1_SOLV[ti][EEF1_LAMBDA] + 1.0 / d) /
                              EEF1_SOLV[ti][EEF1_LAMBDA] +
                          EEF1_SOLV[tj][EEF1_DG_FREE] *
                              EEF1_SOLV[ti][EEF1_VOL] * exp(-xji * xji) *
                              (xji / EEF1_SOLV[tj][EEF1_LAMBDA] + 1.0 / d) /
                              EEF1_SOLV[tj][EEF1_LAMBDA]) *
                         coef2 / d2;
                    DSOLV[0] +=
                        de *
                        ((backrub_uv[0] ^
                          (theAtom->r - gp.getResidue(ires + 1)->getCA()->r)) *
                         Rsm) /
                        d;
                    DSOLV[2] +=
                        de *
                        ((backrub_uv[2] ^
                          (theAtom->r - gp.getResidue(ires)->getCA()->r)) *
                         Rsm) /
                        d;
                  }
                  SOLVi -= solv_ij;
                }
              }
            }
            pt = pt->cell_next;
          }
        }
      }
    }
    E_SOLV += SOLVi;
  }
}

gen_fine_grid::gen_fine_grid(protein &theProtein, gprotein &theGProtein)
    : fine_grid(theProtein, 1), gp(theGProtein) {
  max_r = MAX_ROTH_R;
  max_r2 = MAX_ROTH_R2;

  p.getBoundary(min, max);
  vec tmp = (max - min);
  double tmp2 = tmp.x;
  tmp2 = tmp2 > tmp.y ? tmp2 : tmp.y;
  tmp2 = tmp2 > tmp.z ? tmp2 : tmp.z;
  double ext = 4 * tmp2;
  // double ext = 200.0;
  min -= vec(ext, ext, ext);
  max += vec(ext, ext, ext);
  // assign cell_bex, ncell;
  constructCells();
  ncells = ncell[0] * ncell[1] * ncell[2];

  cell_protons = new atom *[ncells];
  cell_acceptors = new atom *[ncells];
  for (int i = 0; i < ncells; i++) {
    cell_protons[i] = NULL;
    cell_acceptors[i] = NULL;
  }

  // put the protons and acceptors into the file list
  for (int i = 0; i < gp.size(); i++) {
    residue *theRes = gp.getResidue(i);
    atom *theAcceptor = NULL;
    atom *theH = NULL;
    int ix, iy, iz, cindex;
    // first the protons
    atom **array = theRes->getPH_ARRAY();
    for (int j = 0; j < theRes->getNPH(); j++) {
      theH = array[j];
      ix = static_cast<int>((theH->r.x - min.x) / cell_box.x);
      iy = static_cast<int>((theH->r.y - min.y) / cell_box.y);
      iz = static_cast<int>((theH->r.z - min.z) / cell_box.z);
      theH->fine_cell_index = cindex = ix + (iy + iz * ncell[1]) * ncell[0];
      if (!cell_protons[cindex]) { // FIRST
        cell_protons[cindex] = theH;
        cell_protons[cindex]->fine_cell_prev = NULL;
        cell_protons[cindex]->fine_cell_next = NULL;
      } else { // insert from top
        cell_protons[cindex]->fine_cell_prev = theH;
        theH->fine_cell_next = cell_protons[cindex];
        cell_protons[cindex] = theH;
        cell_protons[cindex]->fine_cell_prev = NULL;
      }
    }
    // then the acceptors
    array = theRes->getHBAcceptor_ARRAY();
    for (int j = 0; j < theRes->getNHBAcceptor(); j++) {
      theAcceptor = array[j];
      ix = static_cast<int>((theAcceptor->r.x - min.x) / cell_box.x);
      iy = static_cast<int>((theAcceptor->r.y - min.y) / cell_box.y);
      iz = static_cast<int>((theAcceptor->r.z - min.z) / cell_box.z);
      theAcceptor->fine_cell_index = cindex =
          ix + (iy + iz * ncell[1]) * ncell[0];
      if (!cell_acceptors[cindex]) { // FIRST
        cell_acceptors[cindex] = theAcceptor;
        cell_acceptors[cindex]->fine_cell_prev = NULL;
        cell_acceptors[cindex]->fine_cell_next = NULL;
      } else { // insert from top
        cell_acceptors[cindex]->fine_cell_prev = theAcceptor;
        theAcceptor->fine_cell_next = cell_acceptors[cindex];
        cell_acceptors[cindex] = theAcceptor;
        cell_acceptors[cindex]->fine_cell_prev = NULL;
      }
    }
  }
}

void gen_fine_grid::updateFineCellProton(atom *theH) {
  int ix, iy, iz, cindex, pcindex;
  pcindex = theH->fine_cell_index;
  ix = static_cast<int>((theH->r.x - min.x) / cell_box.x);
  iy = static_cast<int>((theH->r.y - min.y) / cell_box.y);
  iz = static_cast<int>((theH->r.z - min.z) / cell_box.z);
  cindex = ix + (iy + iz * ncell[1]) * ncell[0];
  // if(theH->resID==17){
  //  cout << theH->r.x << " " << theH->r.y << " " << theH->r.z << endl;
  //  cout << theH->resID << " " << cindex << " " <<  pcindex << endl;
  //}
  if (cindex != pcindex) { // cell changed
    // drop from the current fine-cell
    if (cell_protons[pcindex] == theH) { // top
      cell_protons[pcindex] = theH->fine_cell_next;
      if (theH->fine_cell_next)
        theH->fine_cell_next->fine_cell_prev = NULL;
    } else {
      theH->fine_cell_prev->fine_cell_next = theH->fine_cell_next;
      if (theH->fine_cell_next)
        theH->fine_cell_next->fine_cell_prev = theH->fine_cell_prev;
    }
    theH->fine_cell_index = cindex;
    // add to the new cell
    if (!cell_protons[cindex]) {
      cell_protons[cindex] = theH;
      cell_protons[cindex]->fine_cell_prev = NULL;
      cell_protons[cindex]->fine_cell_next = NULL;
    } else { // insert from top
      cell_protons[cindex]->fine_cell_prev = theH;
      theH->fine_cell_next = cell_protons[cindex];
      cell_protons[cindex] = theH;
      cell_protons[cindex]->fine_cell_prev = NULL;
    }
  }
}

void gen_fine_grid::updateFineCellAcceptor(atom *theAcceptor) {
  int ix, iy, iz, cindex, pcindex;
  pcindex = theAcceptor->fine_cell_index;
  ix = static_cast<int>((theAcceptor->r.x - min.x) / cell_box.x);
  iy = static_cast<int>((theAcceptor->r.y - min.y) / cell_box.y);
  iz = static_cast<int>((theAcceptor->r.z - min.z) / cell_box.z);
  cindex = ix + (iy + iz * ncell[1]) * ncell[0];
  if (cindex != pcindex) { // cell changed
    // drop from the current fine-cell
    if (cell_acceptors[pcindex] == theAcceptor) { // top
      cell_acceptors[pcindex] = theAcceptor->fine_cell_next;
      if (theAcceptor->fine_cell_next)
        theAcceptor->fine_cell_next->fine_cell_prev = NULL;
    } else {
      theAcceptor->fine_cell_prev->fine_cell_next = theAcceptor->fine_cell_next;
      if (theAcceptor->fine_cell_next)
        theAcceptor->fine_cell_next->fine_cell_prev =
            theAcceptor->fine_cell_prev;
    }
    theAcceptor->fine_cell_index = cindex;

    // cout << iRes << " " << cindex << endl;
    // add to the new cell
    if (!cell_acceptors[cindex]) { // FIRST
      // cout << "first" << iRes << endl;
      cell_acceptors[cindex] = theAcceptor;
      cell_acceptors[cindex]->fine_cell_prev = NULL;
      cell_acceptors[cindex]->fine_cell_next = NULL;
    } else { // insert from top
      // cout << "top" << iRes << endl;
      cell_acceptors[cindex]->fine_cell_prev = theAcceptor;
      theAcceptor->fine_cell_next = cell_acceptors[cindex];
      cell_acceptors[cindex] = theAcceptor;
      cell_acceptors[cindex]->fine_cell_prev = NULL;
    }
  }
}

void gen_fine_grid::updateFineCellResidue(int iRes) {
  int ix, iy, iz, cindex, pcindex;
  residue *theRes = gp.getResidue(iRes);

  // protons
  atom *theH = NULL;
  atom **array = theRes->getPH_ARRAY();
  for (int i = 0; i < theRes->getNPH(); i++) {
    theH = array[i];
    updateFineCellProton(theH);
  }

  // acceptors
  atom *theAcceptor = NULL;
  array = theRes->getHBAcceptor_ARRAY();
  for (int i = 0; i < theRes->getNHBAcceptor(); i++) {
    theAcceptor = array[i];
    updateFineCellAcceptor(theAcceptor);
  }
}

void gen_fine_grid::dropFineCellResidue(int iRes) {
  int ix, iy, iz, cindex;
  residue *theRes = gp.getResidue(iRes);
  // protons
  atom *theH = NULL;
  atom **array = theRes->getPH_ARRAY();
  for (int i = 0; i < theRes->getNPH(); i++) {
    theH = array[i];
    cindex = theH->fine_cell_index;
    // drop from the current fine-cell
    if (cell_protons[cindex] == theH) { // top
      cell_protons[cindex] = theH->fine_cell_next;
      if (theH->fine_cell_next)
        theH->fine_cell_next->fine_cell_prev = NULL;
    } else {
      theH->fine_cell_prev->fine_cell_next = theH->fine_cell_next;
      if (theH->fine_cell_next)
        theH->fine_cell_next->fine_cell_prev = theH->fine_cell_prev;
    }
  }
  // acceptors
  atom *theAcceptor = NULL;
  array = theRes->getHBAcceptor_ARRAY();
  for (int i = 0; i < theRes->getNHBAcceptor(); i++) {
    theAcceptor = array[i];
    cindex = theAcceptor->fine_cell_index;
    // drop from the current fine-cell
    if (cell_acceptors[cindex] == theAcceptor) { // top
      cell_acceptors[cindex] = theAcceptor->fine_cell_next;
      if (theAcceptor->fine_cell_next)
        theAcceptor->fine_cell_next->fine_cell_prev = NULL;
    } else {
      theAcceptor->fine_cell_prev->fine_cell_next = theAcceptor->fine_cell_next;
      if (theAcceptor->fine_cell_next)
        theAcceptor->fine_cell_next->fine_cell_prev =
            theAcceptor->fine_cell_prev;
    }
  }
}

void gen_fine_grid::insertFineCellResidue(int iRes) {
  int ix, iy, iz, cindex;
  residue *theRes = gp.getResidue(iRes);
  // protons
  atom *theH = NULL;
  atom **array = theRes->getPH_ARRAY();
  for (int i = 0; i < theRes->getNPH(); i++) {
    theH = array[i];
    ix = static_cast<int>((theH->r.x - min.x) / cell_box.x);
    iy = static_cast<int>((theH->r.y - min.y) / cell_box.y);
    iz = static_cast<int>((theH->r.z - min.z) / cell_box.z);
    theH->fine_cell_index = cindex = ix + (iy + iz * ncell[1]) * ncell[0];
    if (!cell_protons[cindex]) { // FIRST
      cell_protons[cindex] = theH;
      cell_protons[cindex]->fine_cell_prev = NULL;
      cell_protons[cindex]->fine_cell_next = NULL;
    } else { // insert from top
      cell_protons[cindex]->fine_cell_prev = theH;
      theH->fine_cell_next = cell_protons[cindex];
      cell_protons[cindex] = theH;
      cell_protons[cindex]->fine_cell_prev = NULL;
    }
  }
  // then the acceptors
  atom *theAcceptor = NULL;
  array = theRes->getHBAcceptor_ARRAY();
  for (int j = 0; j < theRes->getNHBAcceptor(); j++) {
    theAcceptor = array[j];
    ix = static_cast<int>((theAcceptor->r.x - min.x) / cell_box.x);
    iy = static_cast<int>((theAcceptor->r.y - min.y) / cell_box.y);
    iz = static_cast<int>((theAcceptor->r.z - min.z) / cell_box.z);
    theAcceptor->fine_cell_index = cindex =
        ix + (iy + iz * ncell[1]) * ncell[0];
    if (!cell_acceptors[cindex]) { // FIRST
      cell_acceptors[cindex] = theAcceptor;
      cell_acceptors[cindex]->fine_cell_prev = NULL;
      cell_acceptors[cindex]->fine_cell_next = NULL;
    } else { // insert from top
      cell_acceptors[cindex]->fine_cell_prev = theAcceptor;
      theAcceptor->fine_cell_next = cell_acceptors[cindex];
      cell_acceptors[cindex] = theAcceptor;
      cell_acceptors[cindex]->fine_cell_prev = NULL;
    }
  }
}

double getHBE(const int &dSEQ, const int &isSC, const ACCEPTOR_HYBRID_t &hybrid,
              const double &rHA, const double &cosXAH, const double &cosAHD) {
  HB_STAT_t hb_stat;
  HB_TERM_t rHA_var, XAH_var, AHD_var;
  int SC_interp = 0;
  int rHA_interp = 0;
  int XAH_interp = 0;
  int AHD_interp = 0;
  rHA_var = allRHA;
  AHD_var = shortAHD;
  XAH_var = shortXAH;
  if (isSC) {
    if (hybrid == RING)
      hb_stat = SC_RING;
    else if (hybrid == SP2)
      hb_stat = SC_SP2;
    else if (hybrid == SP3)
      hb_stat = SC_SP3;
    else {
      cerr << "error, incorrect hybrid state of acceptor atom" << endl;
      exit(1);
    }
    rHA_var = allRHA;
    if (rHA > RHA_INTERP_MAX) {
      AHD_var = longAHD;
      XAH_var = longXAH;
    } else if (rHA > RHA_INTERP_MIN) {
      SC_interp = 1;
    }
  } else {
    if (dSEQ > 4)
      hb_stat = BB_SHEET;
    else if (dSEQ == 4)
      hb_stat = BB_HELIX;
    else if (dSEQ == 3)
      hb_stat = BB_HELIX;
    else
      return 0;
  }
  rHA_interp = (rHA > RHA_INTERP_EDGE);
  XAH_interp = (cosXAH < ANGLE_INTERP_EDGE);
  AHD_interp = (cosAHD < ANGLE_INTERP_EDGE);
  double E_rHA1 = 0, E_XAH1 = 0, E_AHD1 = 0;
  double E_rHA2 = 0, E_XAH2 = 0, E_AHD2 = 0;
  int table_index = -1;
  // rHA ---> E(rHA)
  // printf("before: %10.8lf %10.8lf %10.8lf\n", rHA, cosXAH, cosAHD);
  if (rHA > maxr[hb_stat]) {
    E_rHA1 = 0;
  } else {
    table_index = POLY_TABLE[hb_stat][rHA_var];
    E_rHA1 = POLYS[table_index].getValue(rHA);
  }
  // AHD ---> E(AHD)
  table_index = POLY_TABLE[hb_stat][AHD_var];
  E_AHD1 = POLYS[table_index].getValue(cosAHD);
  // XAH ---> E(XAH)
  table_index = POLY_TABLE[hb_stat][XAH_var];
  E_XAH1 = POLYS[table_index].getValue(cosXAH);

  // printf("%10.8lf %10.8lf %10.8lf\n", E_rHA1, E_XAH1, E_AHD1);

  // evaluate the interpolations for RHA
  double frac;
  if (SC_interp) {
    // cout << "SC_INTER: " << endl;
    frac = (RHA_INTERP_MAX - rHA) / (RHA_INTERP_MAX - RHA_INTERP_MIN);
    AHD_var = longAHD;
    XAH_var = longXAH;
    // AHD ---> E(AHD)
    table_index = POLY_TABLE[hb_stat][AHD_var];
    E_AHD2 = POLYS[table_index].getValue(cosAHD);
    // XAH ---> E(XAH)
    table_index = POLY_TABLE[hb_stat][XAH_var];
    E_XAH2 = POLYS[table_index].getValue(cosXAH);
    E_AHD1 = frac * E_AHD1 + (1.0 - frac) * E_AHD2;
    E_XAH1 = frac * E_XAH1 + (1.0 - frac) * E_XAH2;
  } else if (rHA_interp) {
    // cout << "RHA_INTER: " << endl;
    if (isSC) {
      // cout << isSC << endl;
      frac = (MAX_R - rHA) / (MAX_R - RHA_INTERP_MAX);
    } else {
      frac = (MAX_R - rHA) / (MAX_R - RHA_INTERP_EDGE);
    }
    E_AHD2 = 0;
    E_XAH2 = 0;
    E_AHD1 = frac * E_AHD1 + (1.0 - frac) * E_AHD2;
    E_XAH1 = frac * E_XAH1 + (1.0 - frac) * E_XAH2;
  }
  // evaluate the interpolations for RHA
  if (XAH_interp) {
    // cout << "XHA_INTER: " << endl;
    E_rHA2 = 0;
    E_AHD2 = 0;
    frac = (MIN_CXAH - cosXAH) / (MIN_CXAH - ANGLE_INTERP_EDGE);
    E_rHA1 = frac * E_rHA1 + (1.0 - frac) * E_rHA2;
    E_AHD1 = frac * E_AHD1 + (1.0 - frac) * E_AHD2;
  }
  // evaluate the interpolations for RHA
  if (AHD_interp) {
    // cout << "AHD_INTER: " << endl;
    E_rHA2 = 0;
    E_XAH2 = 0;
    frac = (MIN_CAHD - cosAHD) / (MIN_CAHD - ANGLE_INTERP_EDGE);
    E_rHA1 = frac * E_rHA1 + (1.0 - frac) * E_rHA2;
    E_XAH1 = frac * E_XAH1 + (1.0 - frac) * E_XAH2;
  }
  // cout << E_rHA1 << " " << E_XAH1 << " " << E_AHD1 << endl;
  // cout << E_rHA1+E_XAH1+E_AHD1 << endl;
  // cout << E_rHA1 << " " << E_XAH1 << " " << E_AHD1 << endl;
  // printf("%10.8lf %10.8lf %10.8lf\n", rHA, cosXAH, cosAHD);
  // printf("%10.8lf %10.8lf %10.8lf\n", E_rHA1, E_XAH1, E_AHD1);
  return E_rHA1 + E_XAH1 + E_AHD1;
}

double getHBE(const int &dSEQ, const int &isSC, const ACCEPTOR_HYBRID_t &hybrid,
              const double &rHA, const double &cosXAH, const double &cosAHD,
              double &EHA, double &EXAH, double &EAHD) {
  HB_STAT_t hb_stat;
  HB_TERM_t rHA_var, XAH_var, AHD_var;
  int SC_interp = 0;
  int rHA_interp = 0;
  int XAH_interp = 0;
  int AHD_interp = 0;
  rHA_var = allRHA;
  AHD_var = shortAHD;
  XAH_var = shortXAH;
  if (isSC) {
    if (hybrid == RING)
      hb_stat = SC_RING;
    else if (hybrid == SP2)
      hb_stat = SC_SP2;
    else if (hybrid == SP3)
      hb_stat = SC_SP3;
    else {
      cerr << "error, incorrect hybrid state of acceptor atom" << endl;
      exit(1);
    }
    rHA_var = allRHA;
    if (rHA > RHA_INTERP_MAX) {
      AHD_var = longAHD;
      XAH_var = longXAH;
    } else if (rHA > RHA_INTERP_MIN) {
      SC_interp = 1;
    }
  } else {
    if (dSEQ > 4)
      hb_stat = BB_SHEET;
    else if (dSEQ == 4)
      hb_stat = BB_HELIX;
    else if (dSEQ == 3)
      hb_stat = BB_HELIX;
    else
      return 0;
  }
  rHA_interp = (rHA > RHA_INTERP_EDGE);
  XAH_interp = (cosXAH < ANGLE_INTERP_EDGE);
  AHD_interp = (cosAHD < ANGLE_INTERP_EDGE);
  double E_rHA1 = 0, E_XAH1 = 0, E_AHD1 = 0;
  double E_rHA2 = 0, E_XAH2 = 0, E_AHD2 = 0;
  int table_index = -1;
  // rHA ---> E(rHA)
  // printf("before: %10.8lf %10.8lf %10.8lf\n", rHA, cosXAH, cosAHD);
  if (rHA > maxr[hb_stat]) {
    E_rHA1 = 0;
  } else {
    table_index = POLY_TABLE[hb_stat][rHA_var];
    E_rHA1 = POLYS[table_index].getValue(rHA);
  }
  // AHD ---> E(AHD)
  table_index = POLY_TABLE[hb_stat][AHD_var];
  E_AHD1 = POLYS[table_index].getValue(cosAHD);
  // XAH ---> E(XAH)
  table_index = POLY_TABLE[hb_stat][XAH_var];
  E_XAH1 = POLYS[table_index].getValue(cosXAH);

  // printf("%10.8lf %10.8lf %10.8lf\n", E_rHA1, E_XAH1, E_AHD1);

  // evaluate the interpolations for RHA
  double frac;
  if (SC_interp) {
    // cout << "SC_INTER: " << endl;
    frac = (RHA_INTERP_MAX - rHA) / (RHA_INTERP_MAX - RHA_INTERP_MIN);
    AHD_var = longAHD;
    XAH_var = longXAH;
    // AHD ---> E(AHD)
    table_index = POLY_TABLE[hb_stat][AHD_var];
    E_AHD2 = POLYS[table_index].getValue(cosAHD);
    // XAH ---> E(XAH)
    table_index = POLY_TABLE[hb_stat][XAH_var];
    E_XAH2 = POLYS[table_index].getValue(cosXAH);
    E_AHD1 = frac * E_AHD1 + (1.0 - frac) * E_AHD2;
    E_XAH1 = frac * E_XAH1 + (1.0 - frac) * E_XAH2;
  } else if (rHA_interp) {
    // cout << "RHA_INTER: " << endl;
    if (isSC) {
      // cout << isSC << endl;
      frac = (MAX_R - rHA) / (MAX_R - RHA_INTERP_MAX);
    } else {
      frac = (MAX_R - rHA) / (MAX_R - RHA_INTERP_EDGE);
    }
    E_AHD2 = 0;
    E_XAH2 = 0;
    E_AHD1 = frac * E_AHD1 + (1.0 - frac) * E_AHD2;
    E_XAH1 = frac * E_XAH1 + (1.0 - frac) * E_XAH2;
  }
  // evaluate the interpolations for RHA
  if (XAH_interp) {
    // cout << "XHA_INTER: " << endl;
    E_rHA2 = 0;
    E_AHD2 = 0;
    frac = (MIN_CXAH - cosXAH) / (MIN_CXAH - ANGLE_INTERP_EDGE);
    E_rHA1 = frac * E_rHA1 + (1.0 - frac) * E_rHA2;
    E_AHD1 = frac * E_AHD1 + (1.0 - frac) * E_AHD2;
  }
  // evaluate the interpolations for RHA
  if (AHD_interp) {
    // cout << "AHD_INTER: " << endl;
    E_rHA2 = 0;
    E_XAH2 = 0;
    frac = (MIN_CAHD - cosAHD) / (MIN_CAHD - ANGLE_INTERP_EDGE);
    E_rHA1 = frac * E_rHA1 + (1.0 - frac) * E_rHA2;
    E_XAH1 = frac * E_XAH1 + (1.0 - frac) * E_XAH2;
  }
  // cout << E_rHA1 << " " << E_XAH1 << " " << E_AHD1 << endl;
  // cout << E_rHA1+E_XAH1+E_AHD1 << endl;
  // cout << E_rHA1 << " " << E_XAH1 << " " << E_AHD1 << endl;
  // printf("%10.8lf %10.8lf %10.8lf\n", rHA, cosXAH, cosAHD);
  // printf("%10.8lf %10.8lf %10.8lf\n", E_rHA1, E_XAH1, E_AHD1);
  EHA = E_rHA1;
  EXAH = E_XAH1;
  EAHD = E_AHD1;
  return E_rHA1 + E_XAH1 + E_AHD1;
}

void gen_fine_grid::getHBE_Proton(atom *theH, double &EHB_bb_bb,
                                  double &EHB_bb_sc, double &EHB_sc_sc) {
  hbond_acceptor *theHBA = NULL;
  atom *theD = NULL;
  atom *theA = NULL;
  atom *theX = NULL;
  atom *pt;
  int cindex;
  int ix, iy, iz, cx, cy, cz, kx, ky, kz;
  double rHA = 0;
  double cXAH = -1;
  double cAHD = -1;
  int dSEQ = 0;
  int isSC = 0;
  double E = 0;
  int bonded = 0;

  // cout << theH->r.x << " " << theH->r.y << " " << theH->r.z << endl;
  // cout << cindex << " " << theH->fine_cell_index << endl;

  if (!theH->isHBonded()) { // the proton does not have a HBond partner;
    bonded = 0;
    theD = theH->attachedHA; // get the Donar heavy atom
    // get the cooridiate of the FineGrid(cell)
    ix = static_cast<int>((theH->r.x - min.x) / cell_box.x);
    iy = static_cast<int>((theH->r.y - min.y) / cell_box.y);
    iz = static_cast<int>((theH->r.z - min.z) / cell_box.z);
    cindex = ix + (iy + iz * ncell[1]) * ncell[0];
    if (theH->fine_cell_index != cindex) { // precaucious
      cout << "getHBE_Proton: fatal error in the cell index of protons" << endl;
      cout << theH->r.x << " " << theH->r.y << " " << theH->r.z << endl;
      cout << cindex << " " << theH->fine_cell_index << endl;
      cout << theH->resID << endl;
      exit(1);
    }
    // search the possible acceptors and calculate the energy
    for (cz = iz - 1; cz <= iz + 1; cz++) {
      kz = cz;
      if (cz < 0)
        kz += ncell[2];
      else if (cz >= ncell[2])
        kz -= ncell[2];

      for (cy = iy - 1; cy <= iy + 1; cy++) {
        ky = cy;
        if (cy < 0)
          ky += ncell[1];
        else if (cy >= ncell[1])
          ky -= ncell[1];

        for (cx = ix - 1; cx <= ix + 1; cx++) {
          kx = cx;
          if (cx < 0)
            kx += ncell[0];
          else if (cx >= ncell[0])
            kx -= ncell[0];

          cindex = ncell[0] * (ncell[1] * kz + ky) + kx;
          pt = cell_acceptors[cindex];

          while (pt) {
            // calculate the hbond energy
            theHBA = pt->acceptor;
            // precaucious of the accpetor
            if (!theHBA) {
              cout << "fatal error: Not a HBacceptor!" << endl;
              exit(1);
            }

            if (!theHBA->isSaturated()) { // not saturated

              theA = theHBA->getA();
              theX = theHBA->getX();
              rHA = theH->r.getDist(theA->r);
              vec XA = theA->r - theX->r;
              XA.Normalize();
              vec AH = theH->r - theA->r;
              AH.Normalize();
              vec HD = theD->r - theH->r;
              HD.Normalize();
              cXAH = XA * AH;
              cAHD = AH * HD;
              dSEQ = abs(theH->resID - theA->resID);
              if (rHA < MAX_R && rHA > MIN_R && cXAH < MAX_CXAH &&
                  cXAH > MIN_CXAH && cAHD < MAX_CAHD && cAHD > MIN_CAHD &&
                  dSEQ > 0) {
                // cout << rHA << " " << cAHD << " " << cXAH <<
                // endl;
                if (!theA->isSC() && !theH->isSC()) { // bb-bb
                  isSC = 0;
                  if ((E = getHBE(dSEQ, isSC, theHBA->getHybrid(), rHA, cXAH,
                                  cAHD)) < 0) {
                    theHBA->addH(theH, E);
                    theH->hbonded_acceptor = theHBA;
                    EHB_bb_bb += E;
                    bonded = 1;

                    // cout << "BB-BB:" << E << " " <<
                    // theH->resID << " " << theA->resID <<
                    // endl;
                  }
                } else if (theH->isSC() && theA->isSC()) { // sc-sc
                  isSC = 1;
                  if ((E = getHBE(dSEQ, isSC, theHBA->getHybrid(), rHA, cXAH,
                                  cAHD)) < 0) {
                    theH->hbonded_acceptor = theHBA;
                    theHBA->addH(theH, E);
                    EHB_sc_sc += E;
                    bonded = 1;

                    // cout << "SC-SC:" << theH->resID << "
                    // " << theA->resID << endl;
                  }
                } else { // bb-sc
                  isSC = 1;
                  if ((E = getHBE(dSEQ, isSC, theHBA->getHybrid(), rHA, cXAH,
                                  cAHD)) < 0) {
                    // cout << rHA << " " << cXAH << " " <<
                    // cAHD << " " << E << endl;
                    // printf("%10.8lf %10.8lf %10.8lf
                    // %10.8lf\n", rHA, cXAH, cAHD, E);
                    theHBA->addH(theH, E);
                    theH->hbonded_acceptor = theHBA;
                    EHB_bb_sc += E;
                    bonded = 1;

                    // cout << "H-A: BB-SC:" << theH->resID
                    // << " " << theA->resID << endl;
                  }
                }
              }
            }
            if (bonded)
              break;
            pt = pt->fine_cell_next;
          }
          if (bonded)
            break;
        }
        if (bonded)
          break;
      }
      if (bonded)
        break;
    }
  } else {
    theHBA = theH->hbonded_acceptor;
    theA = theHBA->getA();
    E = theHBA->getE(theH);
    if (E == INF) {
      cerr << "fatal error in HBONDing pairs" << endl;
      exit(1);
    }
    if (theH->isSC() && theA->isSC())
      EHB_sc_sc += E;
    else if (!theA->isSC() && !theH->isSC()) {
      // cout << theA->resID << " " << theH->resID << endl;
      EHB_bb_bb += E;
    } else
      EHB_bb_sc += E;
  }
}

void gen_fine_grid::getHBE_Acceptor(atom *theA, double &EHB_bb_bb,
                                    double &EHB_bb_sc, double &EHB_sc_sc) {
  hbond_acceptor *theHBA = NULL;
  atom *theH = NULL;
  atom *theD = NULL;
  atom *theX = NULL;
  atom *pt;
  int cindex;
  int ix, iy, iz, cx, cy, cz, kx, ky, kz;
  double rHA = 0;
  double cXAH = -1;
  double cAHD = -1;
  int dSEQ = 0;
  int isSC = 0;
  double E = 0;
  int bonded = 0;
  theHBA = theA->acceptor;
  if (!theHBA) {
    cout << "error the accreptor is not right" << endl;
    exit(1);
  }
  double hbe_bb_bb, hbe_bb_sc, hbe_sc_sc;
  theHBA->currentHBE(hbe_bb_bb, hbe_bb_sc, hbe_sc_sc);
  EHB_bb_bb += hbe_bb_bb;
  EHB_bb_sc += hbe_bb_sc;
  EHB_sc_sc += hbe_sc_sc;
  if (!theHBA->isSaturated()) { // not saturated
    // search the pissible un-hbonded protons
    theX = theHBA->getX();
    ix = static_cast<int>((theA->r.x - min.x) / cell_box.x);
    iy = static_cast<int>((theA->r.y - min.y) / cell_box.y);
    iz = static_cast<int>((theA->r.z - min.z) / cell_box.z);
    cindex = ix + (iy + iz * ncell[1]) * ncell[0];
    if (theA->fine_cell_index != cindex) { // precautions
      cout << "error in the accpetor cell index" << endl;
      exit(1);
    }
    // search teh neighboring cells to find the possible protons for Hbond
    for (cz = iz - 1; cz <= iz + 1; cz++) {
      kz = cz;
      if (cz < 0)
        kz += ncell[2];
      else if (cz >= ncell[2])
        kz -= ncell[2];

      for (cy = iy - 1; cy <= iy + 1; cy++) {
        ky = cy;
        if (cy < 0)
          ky += ncell[1];
        else if (cy >= ncell[1])
          ky -= ncell[1];

        for (cx = ix - 1; cx <= ix + 1; cx++) {
          kx = cx;
          if (cx < 0)
            kx += ncell[0];
          else if (cx >= ncell[0])
            kx -= ncell[0];

          cindex = ncell[0] * (ncell[1] * kz + ky) + kx;
          pt = cell_protons[cindex];

          while (pt) {
            theH = pt;
            theD = theH->attachedHA; // get the Donar
            // precaution
            if (!theD) {
              cout << "fatal error: Not a proton!" << endl;
              exit(1);
            }

            if (!theH->isHBonded()) { // not hbonded
              rHA = theH->r.getDist(theA->r);
              vec XA = theA->r - theX->r;
              XA.Normalize();
              vec AH = theH->r - theA->r;
              AH.Normalize();
              vec HD = theD->r - theH->r;
              HD.Normalize();
              cXAH = XA * AH;
              cAHD = AH * HD;
              dSEQ = abs(theH->resID - theA->resID);
              if (rHA < MAX_R && rHA > MIN_R && cXAH < MAX_CXAH &&
                  cXAH > MIN_CXAH && cAHD < MAX_CAHD && cAHD > MIN_CAHD &&
                  dSEQ > 0) {
                if (!theA->isSC() && !theH->isSC()) { // bb-bb
                  isSC = 0;
                  if ((E = getHBE(dSEQ, isSC, theHBA->getHybrid(), rHA, cXAH,
                                  cAHD)) < 0) {
                    theHBA->addH(theH, E);
                    theH->hbonded_acceptor = theHBA;
                    EHB_bb_bb += E;
                    bonded = 1;
                    // cout << "BB-BB:" << E << "  " <<
                    // theH->resID << " " << theA->resID <<
                    // endl;
                  }
                } else if (theH->isSC() && theA->isSC()) { // sc-sc
                  isSC = 1;
                  if ((E = getHBE(dSEQ, isSC, theHBA->getHybrid(), rHA, cXAH,
                                  cAHD)) < 0) {
                    theH->hbonded_acceptor = theHBA;
                    theHBA->addH(theH, E);
                    EHB_sc_sc += E;
                    bonded = 1;

                    // cout << "SC-SC:" << theH->resID << "
                    // " << theA->resID << endl;
                  }
                } else { // bb-sc
                  isSC = 1;
                  if ((E = getHBE(dSEQ, isSC, theHBA->getHybrid(), rHA, cXAH,
                                  cAHD)) < 0) {
                    theHBA->addH(theH, E);
                    theH->hbonded_acceptor = theHBA;
                    EHB_bb_sc += E;
                    bonded = 1;
                    // cout << "A-H: BB-SC:" << theA->resID
                    // << " " << theH->resID << endl;
                  }
                }
              }
            }
            if (theHBA->isSaturated())
              break;
            pt = pt->fine_cell_next;
          }
          if (theHBA->isSaturated())
            break;
        }
        if (theHBA->isSaturated())
          break;
      }
      if (theHBA->isSaturated())
        break;
    }
  }
}

void gen_fine_grid::getHBE_RES(int iRes, double &EHB_bb_bb, double &EHB_bb_sc,
                               double &EHB_sc_sc) {
  residue *theRes = gp.getResidue(iRes);
  atom *theH;
  atom *theA;
  atom **array = theRes->getPH_ARRAY();
  for (int iH = 0; iH < theRes->getNPH(); iH++) {
    theH = array[iH];
    getHBE_Proton(theH, EHB_bb_bb, EHB_bb_sc, EHB_sc_sc);
  }
  // acceptors
  array = theRes->getHBAcceptor_ARRAY();
  for (int iA = 0; iA < theRes->getNHBAcceptor(); iA++) {
    theA = array[iA];
    getHBE_Acceptor(theA, EHB_bb_bb, EHB_bb_sc, EHB_sc_sc);
  }
}

void gen_fine_grid::clearHBE_Proton(atom *theH, double &EHB_bb_bb,
                                    double &EHB_bb_sc, double &EHB_sc_sc) {
  atom *theA = NULL;
  hbond_acceptor *theHBA = NULL;
  double E;
  if (theH->isHBonded()) {
    theHBA = theH->hbonded_acceptor;
    theH->clearHBAcceptor();
    theA = theHBA->getA();
    theHBA->clearH(theH, E);

    if (theA->isSC() && theH->isSC()) { // sc_sc
      EHB_sc_sc += E;
    } else if (!theA->isSC() && !theH->isSC()) {
      EHB_bb_bb += E;
    } else {
      EHB_bb_sc += E;
    }
  }
}

void gen_fine_grid::clearHBE_Acceptor(atom *theA, double &EHB_bb_bb,
                                      double &EHB_bb_sc, double &EHB_sc_sc) {
  hbond_acceptor *theHBA = NULL;
  theHBA = theA->acceptor;
  theHBA->clearHs(EHB_bb_bb, EHB_bb_sc, EHB_sc_sc);
}

void gen_fine_grid::clearHBE_RES(int iRes, double &EHB_bb_bb, double &EHB_bb_sc,
                                 double &EHB_sc_sc) {
  residue *theRes = gp.getResidue(iRes);
  atom **array = theRes->getPH_ARRAY();
  atom *theH;
  atom *theA;
  hbond_acceptor *theHBA = NULL;
  double E;
  EHB_bb_bb = EHB_bb_sc = EHB_sc_sc = 0;
  for (int iH = 0; iH < theRes->getNPH(); iH++) {
    theH = array[iH];
    if (theH->isHBonded()) {
      theHBA = theH->hbonded_acceptor;
      theH->clearHBAcceptor();
      theA = theHBA->getA();
      theHBA->clearH(theH, E);

      if (theA->isSC() && theH->isSC()) { // sc_sc
        EHB_sc_sc += E;
      } else if (!theA->isSC() && !theH->isSC()) {
        EHB_bb_bb += E;
      } else {
        EHB_bb_sc += E;
      }
    }
  }

  array = theRes->getHBAcceptor_ARRAY();
  for (int iA = 0; iA < theRes->getNHBAcceptor(); iA++) {
    theA = array[iA];
    theHBA = theA->acceptor;
    theHBA->clearHs(EHB_bb_bb, EHB_bb_sc, EHB_sc_sc);
  }
}

void gen_fine_grid::clearSC_HBE_RES(int iRes, double &EHB_bb_sc,
                                    double &EHB_sc_sc) {
  residue *theRes = gp.getResidue(iRes);
  atom **array = theRes->getPH_ARRAY();
  atom *theH;
  atom *theA;
  hbond_acceptor *theHBA = NULL;
  double E;
  EHB_bb_sc = EHB_sc_sc = 0;
  for (int iH = 0; iH < theRes->getNPH(); iH++) {
    theH = array[iH];
    if (theH->isHBonded()) {
      theHBA = theH->hbonded_acceptor;
      theA = theHBA->getA();
      if (theA->isSC() || theH->isSC()) { // not a bb-bb
        theH->clearHBAcceptor();
        theHBA->clearH(theH, E);
        if (theA->isSC() && theH->isSC()) { // sc_sc
          EHB_sc_sc += E;
        } else { // bb_sc
          EHB_bb_sc += E;
        }
      }
    }
  }

  array = theRes->getHBAcceptor_ARRAY();
  for (int iA = 0; iA < theRes->getNHBAcceptor(); iA++) {
    theA = array[iA];
    theHBA = theA->acceptor;
    atom **Hs = theHBA->getHs();
    double *Es = theHBA->getEHBs();
    for (int ii = 0; ii < theHBA->getMAX_HBONDS(); ii++) {
      if (Hs[ii]) {
        if (Hs[ii]->isSC() || theA->isSC()) {
          if (Hs[ii]->isSC() && theA->isSC()) // sc_sc
            EHB_sc_sc += Es[ii];
          else // bb_sc
            EHB_bb_sc += Es[ii];
          Hs[ii]->clearHBAcceptor();
          Hs[ii] = NULL;
        }
      }
    }
  }
}

/*********************************************************************
  Calculate the derivative of hydogen bond energy with respect to
  to the chi angles!!!

 In addition, in order to calculate the DERIVATIVE, it is better
 to clear the HBONDS of sc before call this function!!!!

 For a sc-sc, sc-bb, IT IS alway TRUE that all the rotation angles chi
will affect the energy.

**********************************************************************/
void gen_fine_grid::getHBE_RES_DEV(int iRes, double &EHB_bb_bb,
                                   double &EHB_bb_sc, double &EHB_sc_sc,
                                   double *DEV_bb_sc, double *DEV_sc_sc) {
  residue *theRes = gp.getResidue(iRes);
  hbond_acceptor *theHBA = NULL;
  atom *theH = NULL;
  atom *theD = NULL;
  atom *theA = NULL;
  atom *theX = NULL;
  atom *pt;
  int cindex;
  int ix, iy, iz, cx, cy, cz, kx, ky, kz;
  double rHA = 0;
  double cXAH = -1;
  double cAHD = -1;
  int dSEQ = 0;
  int isSC = 0;
  double E = 0;
  int bonded = 0;
  double dE_dr;
  double dE_dxH;
  double dE_dxD;

  // protons
  int type = theRes->getType();
  int nCHI = theRes->getNChi();
  // cout << nCHI << endl;

  int hasRotH = 0;
  if (type == SER || type == THR || type == TYR || type == LYS) {
    hasRotH = 1;
  }
  for (int ich = 0; ich < nCHI; ich++) {
    DEV_bb_sc[ich] = 0;
    DEV_sc_sc[ich] = 0;
  }
  if (hasRotH) {
    DEV_bb_sc[nCHI] = 0;
    DEV_sc_sc[nCHI] = 0;
  }

  atom **array = theRes->getPH_ARRAY();
  for (int iH = 0; iH < theRes->getNPH(); iH++) {
    theH = array[iH];
    if (!theH->isHBonded()) { // the proton does not have a HBond partner;
      bonded = 0;
      theD = theH->attachedHA; // get the Donar heavy atom

      // get the cooridiate of the FineGrid(cell)
      ix = static_cast<int>((theH->r.x - min.x) / cell_box.x);
      iy = static_cast<int>((theH->r.y - min.y) / cell_box.y);
      iz = static_cast<int>((theH->r.z - min.z) / cell_box.z);
      cindex = ix + (iy + iz * ncell[1]) * ncell[0];
      if (theH->fine_cell_index != cindex) { // precaucious
        cout << "HBE_RES_DEV: fatal error in the cell index of protons" << endl;
        cout << theH->resID << endl;
        exit(1);
      }
      // search the possible acceptors and calculate the energy
      for (cz = iz - 1; cz <= iz + 1; cz++) {
        kz = cz;
        if (cz < 0)
          kz += ncell[2];
        else if (cz >= ncell[2])
          kz -= ncell[2];

        for (cy = iy - 1; cy <= iy + 1; cy++) {
          ky = cy;
          if (cy < 0)
            ky += ncell[1];
          else if (cy >= ncell[1])
            ky -= ncell[1];

          for (cx = ix - 1; cx <= ix + 1; cx++) {
            kx = cx;
            if (cx < 0)
              kx += ncell[0];
            else if (cx >= ncell[0])
              kx -= ncell[0];

            cindex = ncell[0] * (ncell[1] * kz + ky) + kx;
            pt = cell_acceptors[cindex];

            while (pt) {
              // calculate the hbond energy
              theHBA = pt->acceptor;
              // precaucious of the accpetor
              if (!theHBA) {
                cout << "fatal error: Not a HBacceptor!" << endl;
                exit(1);
              }

              if (!theHBA->isSaturated()) { // not saturated

                theA = theHBA->getA();
                theX = theHBA->getX();
                rHA = theH->r.getDist(theA->r);
                vec XA = theA->r - theX->r;
                XA.Normalize();
                vec AH = theH->r - theA->r;
                AH.Normalize();
                vec HD = theD->r - theH->r;
                HD.Normalize();
                cXAH = XA * AH;
                cAHD = AH * HD;
                dSEQ = abs(theH->resID - theA->resID);
                if (rHA < MAX_R && rHA > MIN_R && cXAH < MAX_CXAH &&
                    cXAH > MIN_CXAH && cAHD < MAX_CAHD && cAHD > MIN_CAHD &&
                    dSEQ > 0) {
                  // cout << nCHI << "hi" << endl;
                  // cout << rHA << endl;
                  // cout << rHA << " " << cAHD << " " << cXAH
                  // << endl;
                  if (!theA->isSC() && !theH->isSC()) { // bb-bb
                    isSC = 0;
                    if ((E = getHBE(dSEQ, isSC, theHBA->getHybrid(), rHA, cXAH,
                                    cAHD)) < 0) {
                      theHBA->addH(theH, E);
                      theH->hbonded_acceptor = theHBA;
                      EHB_bb_bb += E;
                      bonded = 1;

                      // cout << "BB-BB:" << theH->resID
                      // << " " << theA->resID << endl;
                    }
                  } else if (theH->isSC() && theA->isSC()) { // sc-sc
                    isSC = 1;
                    if ((E = getHBE_DEV(dSEQ, isSC, theHBA->getHybrid(), rHA,
                                        cXAH, cAHD, dE_dr, dE_dxH, dE_dxD)) <
                        0) {
                      theH->hbonded_acceptor = theHBA;
                      theHBA->addH(theH, E);
                      EHB_sc_sc += E;
                      bonded = 1;

                      // cout << "H: sc-sc" << rHA << " "
                      // << nCHI << " " << dE_dr << " " <<
                      // dE_dxD << " " << dE_dxH << endl;
                      /*get the derivate! for all the
                        protons, all the chi angles
                        contribute*/
                      for (int ich = 0; ich < nCHI; ich++) {
                        /*rHA*/
                        vec R2 = theH->r - *theRes->getChiP2(ich);
                        vec VxR2_norm = *theRes->getChiUV(ich) ^ R2 / rHA;

                        // cout << rHA << endl;

                        DEV_sc_sc[ich] += dE_dr * (VxR2_norm * AH) * rHA;

                        DEV_sc_sc[ich] +=
                            dE_dxD * ((*theRes->getChiUV(ich) ^ HD) * AH +
                                      VxR2_norm * (HD - cAHD * AH));

                        // DEV_sc_sc[ich] += dE_dxD*(
                        // (*theRes->getChiUV[ich]^HD_fix)*AH_fix
                        // + 		   VxR2*HD_fix -
                        //		   (HD_fix*AH_fix)*(AH_fix*VxR2)/(rHA*rHA)
                        //)/(rHA*HD_fix.getMag());

                        DEV_sc_sc[ich] +=
                            dE_dxH * (VxR2_norm * (XA - cXAH * AH));

                        // DEV_sc_sc[ich] += dE_dxH*(
                        // VxR2*XA -
                        // cXAH*(VxR2*AH_fix)/rHA )/rHA;
                      }
                      if (hasRotH) { // ser,thr,tyr,LYS;
                                     // donar do not move
                        vec R2 = theH->r - theD->r;
                        vec VxR2_norm = *theRes->getChiUV(nCHI) ^ R2 / rHA;

                        DEV_sc_sc[nCHI] += dE_dr * (VxR2_norm * AH) * rHA;

                        DEV_sc_sc[nCHI] +=
                            dE_dxD * ((*theRes->getChiUV(nCHI) ^ HD) * AH +
                                      -(VxR2_norm * AH) * cAHD);

                        DEV_sc_sc[nCHI] +=
                            dE_dxH * (VxR2_norm * (XA - cXAH * AH));
                      }

                      // cout << "SC-SC:" << theH->resID
                      // << " " << theA->resID << endl;
                    }
                  } else { // bb-sc
                    isSC = 1;
                    if ((E = getHBE_DEV(dSEQ, isSC, theHBA->getHybrid(), rHA,
                                        cXAH, cAHD, dE_dr, dE_dxH, dE_dxD)) <
                        0) {
                      // cout << rHA << " " << cXAH << " "
                      // << cAHD << endl;
                      theHBA->addH(theH, E);
                      theH->hbonded_acceptor = theHBA;
                      EHB_bb_sc += E;
                      bonded = 1;

                      // cout << "H: bb_sc " << rHA << " "
                      // << cXAH<< " " << cAHD << " " <<
                      // dE_dr << " " << dE_dxD << " " <<
                      // dE_dxH << endl; only the Proton is
                      // in the sidechain. We do the
                      // calculation.
                      if (theH->isSC()) {
                        // cout << "H: bb_sc " << rHA <<
                        // " " << cXAH<< " " << cAHD << "
                        // " <<E <<" "<< dE_dr << " " <<
                        // dE_dxH << " " << dE_dxD <<
                        // endl; printf("%10.8lf %10.8lf
                        // %10.8lf %10.8lf %10.8lf
                        // %10.8lf %10.8lf\n", rHA, cXAH,
                        // cAHD, E, dE_dr, dE_dxH,
                        // dE_dxD);
                        /*get the derivate! for all the
                          protons, all the chi angles
                          contribute*/
                        for (int ich = 0; ich < nCHI; ich++) {
                          /*rHA*/
                          vec R2 = theH->r - *theRes->getChiP2(ich);
                          vec VxR2_norm = *theRes->getChiUV(ich) ^ R2 / rHA;

                          DEV_bb_sc[ich] += dE_dr * (VxR2_norm * AH) * rHA;

                          DEV_bb_sc[ich] +=
                              dE_dxD * ((*theRes->getChiUV(ich) ^ HD) * AH +
                                        VxR2_norm * (HD - cAHD * AH));

                          // DEV_sc_sc[ich] +=
                          // dE_dxD*(
                          // (*theRes->getChiUV[ich]^HD_fix)*AH_fix
                          // + 		   VxR2*HD_fix -
                          //		   (HD_fix*AH_fix)*(AH_fix*VxR2)/(rHA*rHA)
                          //)/(rHA*HD_fix.getMag());

                          DEV_bb_sc[ich] +=
                              dE_dxH * (VxR2_norm * (XA - cXAH * AH));

                          // DEV_sc_sc[ich] +=
                          // dE_dxH*( VxR2*XA -
                          // cXAH*(VxR2*AH_fix)/rHA
                          // )/rHA;
                        }
                        if (hasRotH) { // ser,thr,tyr,LYS;
                                       // donar do not
                                       // move
                          // vec R2 = theH->r -
                          // *theRes->getChiP2(nCHI);
                          vec R2 = theH->r - theD->r;
                          vec VxR2_norm = *theRes->getChiUV(nCHI) ^ R2 / rHA;

                          DEV_bb_sc[nCHI] += dE_dr * (VxR2_norm * AH) * rHA;

                          DEV_bb_sc[nCHI] +=
                              dE_dxD * ((*theRes->getChiUV(nCHI) ^ HD) * AH +
                                        -(VxR2_norm * AH) * cAHD);

                          DEV_bb_sc[nCHI] +=
                              dE_dxH * (VxR2_norm * (XA - cXAH * AH));
                        }
                      }
                      // cout << "BB-SC:" << theH->resID
                      // << " " << theA->resID << endl;
                    }
                  }
                }
              }
              if (bonded)
                break;
              pt = pt->fine_cell_next;
            }
            if (bonded)
              break;
          }
          if (bonded)
            break;
        }
        if (bonded)
          break;
      }
    } else {
      theHBA = theH->hbonded_acceptor;
      theA = theHBA->getA();
      E = theHBA->getE(theH);
      if (E == INF) {
        cerr << "fatal error in HBONDing pairs" << endl;
        exit(1);
      }
      if (theH->isSC() && theA->isSC())
        EHB_sc_sc += E;
      else if (!theA->isSC() && !theH->isSC())
        EHB_bb_bb += E;
      else
        EHB_bb_sc += E;
    }
  }
  // acceptors
  array = theRes->getHBAcceptor_ARRAY();
  double hbe_bb_bb, hbe_bb_sc, hbe_sc_sc;
  for (int iA = 0; iA < theRes->getNHBAcceptor(); iA++) {
    theA = array[iA];
    theHBA = theA->acceptor;
    if (!theHBA) {
      cout << "error the accreptor is not right" << endl;
      exit(1);
    }
    theHBA->currentHBE(hbe_bb_bb, hbe_bb_sc, hbe_sc_sc);
    EHB_bb_bb += hbe_bb_bb;
    EHB_bb_sc += hbe_bb_sc;
    EHB_sc_sc += hbe_sc_sc;
    if (!theHBA->isSaturated()) { // not saturated
      // search the pissible un-hbonded protons
      theX = theHBA->getX();
      ix = static_cast<int>((theA->r.x - min.x) / cell_box.x);
      iy = static_cast<int>((theA->r.y - min.y) / cell_box.y);
      iz = static_cast<int>((theA->r.z - min.z) / cell_box.z);
      cindex = ix + (iy + iz * ncell[1]) * ncell[0];
      if (theA->fine_cell_index != cindex) { // precautions
        cout << "error in the accpetor cell index" << endl;
        exit(1);
      }
      // search teh neighboring cells to find the possible protons for
      // Hbond
      for (cz = iz - 1; cz <= iz + 1; cz++) {
        kz = cz;
        if (cz < 0)
          kz += ncell[2];
        else if (cz >= ncell[2])
          kz -= ncell[2];

        for (cy = iy - 1; cy <= iy + 1; cy++) {
          ky = cy;
          if (cy < 0)
            ky += ncell[1];
          else if (cy >= ncell[1])
            ky -= ncell[1];

          for (cx = ix - 1; cx <= ix + 1; cx++) {
            kx = cx;
            if (cx < 0)
              kx += ncell[0];
            else if (cx >= ncell[0])
              kx -= ncell[0];

            cindex = ncell[0] * (ncell[1] * kz + ky) + kx;
            pt = cell_protons[cindex];

            while (pt) {
              theH = pt;
              theD = theH->attachedHA; // get the Donar
              // precaution
              if (!theD) {
                cout << "fatal error: Not a proton!" << endl;
                exit(1);
              }

              if (!theH->isHBonded()) { // not hbonded
                rHA = theH->r.getDist(theA->r);
                vec XA = theA->r - theX->r;
                XA.Normalize();
                vec AH = theH->r - theA->r;
                AH.Normalize();
                vec HD = theD->r - theH->r;
                HD.Normalize();
                cXAH = XA * AH;
                cAHD = AH * HD;
                dSEQ = abs(theH->resID - theA->resID);
                if (rHA < MAX_R && rHA > MIN_R && cXAH < MAX_CXAH &&
                    cXAH > MIN_CXAH && cAHD < MAX_CAHD && cAHD > MIN_CAHD &&
                    dSEQ > 0) {
                  // cout << nCHI << "hi" << endl;

                  if (!theA->isSC() && !theH->isSC()) { // bb-bb
                    isSC = 0;
                    if ((E = getHBE(dSEQ, isSC, theHBA->getHybrid(), rHA, cXAH,
                                    cAHD)) < 0) {
                      theHBA->addH(theH, E);
                      theH->hbonded_acceptor = theHBA;
                      EHB_bb_bb += E;
                      bonded = 1;
                      // cout << "BB-BB:" << theH->resID
                      // << " " << theA->resID << endl;
                    }
                  } else if (theH->isSC() && theA->isSC()) { // sc-sc
                    isSC = 1;
                    if ((E = getHBE_DEV(dSEQ, isSC, theHBA->getHybrid(), rHA,
                                        cXAH, cAHD, dE_dr, dE_dxH, dE_dxD)) <
                        0) {
                      theH->hbonded_acceptor = theHBA;
                      theHBA->addH(theH, E);
                      EHB_sc_sc += E;
                      bonded = 1;

                      // cout << "A: SC-SC:" <<
                      // theH->resID << " " << theA->resID
                      // << endl;

                      /*get the derivaste for all the
                       * acceptors*/
                      for (int ich = 0; ich < nCHI; ich++) {
                        vec R2 = theA->r - *theRes->getChiP2(ich);
                        vec VxR2_norm = *theRes->getChiUV(ich) ^ R2 / rHA;

                        DEV_sc_sc[ich] -= dE_dr * (VxR2_norm * AH) * rHA;

                        DEV_sc_sc[ich] -=
                            dE_dxD * (VxR2_norm * (HD - cAHD * AH));

                        DEV_sc_sc[ich] +=
                            dE_dxH * ((*theRes->getChiUV(ich) ^ XA) * AH -
                                      VxR2_norm * (XA - cXAH * AH));
                      }
                    }
                  } else { // bb-sc
                    isSC = 1;
                    if ((E = getHBE_DEV(dSEQ, isSC, theHBA->getHybrid(), rHA,
                                        cXAH, cAHD, dE_dr, dE_dxH, dE_dxD)) <
                        0) {
                      theHBA->addH(theH, E);
                      theH->hbonded_acceptor = theHBA;
                      EHB_bb_sc += E;
                      bonded = 1;

                      // cout << "A: BB-SC:" <<
                      // theH->resID << " " << theA->resID
                      // << endl;
                      if (theA->isSC()) {
                        for (int ich = 0; ich < nCHI; ich++) {
                          vec R2 = theA->r - *theRes->getChiP2(ich);
                          vec VxR2_norm = *theRes->getChiUV(ich) ^ R2 / rHA;

                          DEV_bb_sc[ich] -= dE_dr * (VxR2_norm * AH) * rHA;

                          DEV_bb_sc[ich] -=
                              dE_dxD * (VxR2_norm * (HD - cAHD * AH));

                          DEV_bb_sc[ich] +=
                              dE_dxH * ((*theRes->getChiUV(ich) ^ XA) * AH -
                                        VxR2_norm * (XA - cXAH * AH));
                        }
                      }
                    }
                  }
                }
              }
              if (theHBA->isSaturated())
                break;
              pt = pt->fine_cell_next;
            }
            if (theHBA->isSaturated())
              break;
          }
          if (theHBA->isSaturated())
            break;
        }
        if (theHBA->isSaturated())
          break;
      }
    }
  }
}

/*********************************************************************
To calculate the hydrogen bonds for the whole protein. To avoid over counting,
we only search the protons.
*********************************************************************/
void gen_fine_grid::getEHB_DEV(double &EHB_bb_bb, double &EHB_bb_sc,
                               double &EHB_sc_sc, double *DEV_bb_bb,
                               double *DEV_bb_sc, double *DEV_sc_sc) {
  residue *theRes = NULL;
  hbond_acceptor *theHBA = NULL;
  atom *theH = NULL;
  atom *theD = NULL;
  atom *theA = NULL;
  atom *theX = NULL;
  atom *pt;
  int cindex;
  int ix, iy, iz, cx, cy, cz, kx, ky, kz;
  double rHA = 0;
  double cXAH = -1;
  double cAHD = -1;
  int dSEQ = 0;
  int isSC = 0;
  double E = 0;
  int bonded = 0;
  double dE_dr;
  double dE_dxH;
  double dE_dxD;

  // protons
  // int type=theRes->getType();
  // int nCHI = theRes->getNChi();
  // cout << nCHI << endl;

  EHB_bb_bb = EHB_bb_sc = EHB_sc_sc = 0;
  int jres;
  for (int ires = 0; ires < gp.size(); ires++) {
    theRes = gp.getResidue(ires);
    atom **array = theRes->getPH_ARRAY();
    for (int iH = 0; iH < theRes->getNPH(); iH++) {
      // cout << ires << " " << theRes->getNPH() << endl;
      theH = array[iH];
      theD = theH->attachedHA;
      bonded = 0;
      if (!theH->isHBonded()) {
        // get the cooridiate of the FineGrid(cell)
        ix = static_cast<int>((theH->r.x - min.x) / cell_box.x);
        iy = static_cast<int>((theH->r.y - min.y) / cell_box.y);
        iz = static_cast<int>((theH->r.z - min.z) / cell_box.z);
        cindex = ix + (iy + iz * ncell[1]) * ncell[0];
        if (theH->fine_cell_index != cindex) { // precaucious
          cout << "getEHB_DEV: fatal error in the cell index of "
                  "protons"
               << endl;
          cout << theH->resID << endl;
          exit(1);
        }
        // search the possible acceptors and calculate the energy
        for (cz = iz - 1; cz <= iz + 1; cz++) {
          kz = cz;
          if (cz < 0)
            kz += ncell[2];
          else if (cz >= ncell[2])
            kz -= ncell[2];

          for (cy = iy - 1; cy <= iy + 1; cy++) {
            ky = cy;
            if (cy < 0)
              ky += ncell[1];
            else if (cy >= ncell[1])
              ky -= ncell[1];

            for (cx = ix - 1; cx <= ix + 1; cx++) {
              kx = cx;
              if (cx < 0)
                kx += ncell[0];
              else if (cx >= ncell[0])
                kx -= ncell[0];

              cindex = ncell[0] * (ncell[1] * kz + ky) + kx;
              pt = cell_acceptors[cindex];

              while (pt) {
                // calculate the hbond energy
                theHBA = pt->acceptor;
                jres = pt->resID;
                // precaucious of the accpetor
                if (!theHBA) {
                  cout << "fatal error: Not a HBacceptor!" << endl;
                  exit(1);
                }

                if (!theHBA->isSaturated()) { // not saturated
                  theA = theHBA->getA();
                  theX = theHBA->getX();
                  rHA = theH->r.getDist(theA->r);
                  vec XA = theA->r - theX->r;
                  XA.Normalize();
                  vec AH = theH->r - theA->r;
                  AH.Normalize();
                  vec HD = theD->r - theH->r;
                  HD.Normalize();
                  cXAH = XA * AH;
                  cAHD = AH * HD;
                  dSEQ = abs(theH->resID - theA->resID);
                  if (rHA < MAX_R && rHA > MIN_R && cXAH < MAX_CXAH &&
                      cXAH > MIN_CXAH && cAHD < MAX_CAHD && cAHD > MIN_CAHD &&
                      dSEQ > 0) {
                    if (!theA->isSC() && !theH->isSC()) { // bb-bb
                      isSC = 0;
                      if ((E = getHBE_DEV(dSEQ, isSC, theHBA->getHybrid(), rHA,
                                          cXAH, cAHD, dE_dr, dE_dxH, dE_dxD)) <
                          0) {
                        theHBA->addH(theH, E);
                        theH->hbonded_acceptor = theHBA;
                        EHB_bb_bb += E;
                        bonded = 1;
                        // cout << E << " " << ires << "
                        // " << jres << endl;
                        if (ires < jres) {
                          // psi of ires
                          vec R2;
                          vec VxR2_norm;
                          for (int ii = ires; ii <= jres; ii++) {
                            // phi of ii
                            R2 = theA->r - gp.getResidue(ii)->getCA()->r;
                            VxR2_norm = gp.getPhi_UV(ii) ^ R2 / rHA;

                            DEV_bb_bb[2 * ii] -= dE_dr * (VxR2_norm * AH) * rHA;
                            DEV_bb_bb[2 * ii] -=
                                dE_dxD * (VxR2_norm * (HD - cAHD * AH));
                            DEV_bb_bb[2 * ii] +=
                                dE_dxH * ((gp.getPhi_UV(ii) ^ XA) * AH -
                                          VxR2_norm * (XA - cXAH * AH));
                            // psi of ii
                            R2 = theA->r - gp.getResidue(ii)->getC()->r;
                            VxR2_norm = gp.getPsi_UV(ii) ^ R2 / rHA;

                            DEV_bb_bb[2 * ii + 1] -=
                                dE_dr * (VxR2_norm * AH) * rHA;
                            DEV_bb_bb[2 * ii + 1] -=
                                dE_dxD * (VxR2_norm * (HD - cAHD * AH));
                            DEV_bb_bb[2 * ii + 1] +=
                                dE_dxH * ((gp.getPsi_UV(ii) ^ XA) * AH -
                                          VxR2_norm * (XA - cXAH * AH));

                            // cout << dE_dr << " "
                            // << dE_dxD << " " <<
                            // dE_dxH << endl;
                          }
                        } else if (ires > jres) {
                          vec R2;
                          vec VxR2_norm;
                          for (int ii = jres + 1; ii <= ires - 1; ii++) {
                            // phi of ii
                            R2 = theH->r - gp.getResidue(ii)->getCA()->r;
                            VxR2_norm = gp.getPhi_UV(ii) ^ R2 / rHA;
                            DEV_bb_bb[2 * ii] += dE_dr * (VxR2_norm * AH) * rHA;
                            DEV_bb_bb[2 * ii] +=
                                dE_dxD * ((gp.getPhi_UV(ii) ^ HD) * AH +
                                          VxR2_norm * (HD - cAHD * AH));
                            DEV_bb_bb[2 * ii] +=
                                dE_dxH * (VxR2_norm * (XA - cXAH * AH));

                            // psi of ii
                            R2 = theH->r - gp.getResidue(ii)->getC()->r;
                            VxR2_norm = gp.getPsi_UV(ii) ^ R2 / rHA;
                            DEV_bb_bb[2 * ii + 1] +=
                                dE_dr * (VxR2_norm * AH) * rHA;
                            DEV_bb_bb[2 * ii + 1] +=
                                dE_dxD * ((gp.getPsi_UV(ii) ^ HD) * AH +
                                          VxR2_norm * (HD - cAHD * AH));
                            DEV_bb_bb[2 * ii + 1] +=
                                dE_dxH * (VxR2_norm * (XA - cXAH * AH));
                          }
                        }
                      }
                    } else if (theH->isSC() && theA->isSC()) { // sc-sc
                      isSC = 1;
                      if ((E = getHBE_DEV(dSEQ, isSC, theHBA->getHybrid(), rHA,
                                          cXAH, cAHD, dE_dr, dE_dxH, dE_dxD)) <
                          0) {
                        theH->hbonded_acceptor = theHBA;
                        theHBA->addH(theH, E);
                        EHB_sc_sc += E;
                        bonded = 1;

                        // cout << E << " " << ires << "
                        // " << jres << endl;
                        if (ires < jres) {
                          // psi of ires
                          vec R2 = theA->r - gp.getResidue(ires)->getC()->r;
                          vec VxR2_norm = gp.getPsi_UV(ires) ^ R2 / rHA;

                          DEV_sc_sc[2 * ires + 1] -=
                              dE_dr * (VxR2_norm * AH) * rHA;
                          DEV_sc_sc[2 * ires + 1] -=
                              dE_dxD * (VxR2_norm * (HD - cAHD * AH));
                          DEV_sc_sc[2 * ires + 1] +=
                              dE_dxH * ((gp.getPsi_UV(ires) ^ XA) * AH -
                                        VxR2_norm * (XA - cXAH * AH));
                          for (int ii = ires + 1; ii <= jres - 1; ii++) {
                            // phi of ii
                            R2 = theA->r - gp.getResidue(ii)->getCA()->r;
                            VxR2_norm = gp.getPhi_UV(ii) ^ R2 / rHA;

                            DEV_sc_sc[2 * ii] -= dE_dr * (VxR2_norm * AH) * rHA;
                            DEV_sc_sc[2 * ii] -=
                                dE_dxD * (VxR2_norm * (HD - cAHD * AH));
                            DEV_sc_sc[2 * ii] +=
                                dE_dxH * ((gp.getPhi_UV(ii) ^ XA) * AH -
                                          VxR2_norm * (XA - cXAH * AH));
                            // psi of ii
                            R2 = theA->r - gp.getResidue(ii)->getC()->r;
                            VxR2_norm = gp.getPsi_UV(ii) ^ R2 / rHA;

                            DEV_sc_sc[2 * ii + 1] -=
                                dE_dr * (VxR2_norm * AH) * rHA;
                            DEV_sc_sc[2 * ii + 1] -=
                                dE_dxD * (VxR2_norm * (HD - cAHD * AH));
                            DEV_sc_sc[2 * ii + 1] +=
                                dE_dxH * ((gp.getPsi_UV(ii) ^ XA) * AH -
                                          VxR2_norm * (XA - cXAH * AH));
                          }
                          // phi of jres
                          R2 = theA->r - gp.getResidue(jres)->getCA()->r;
                          VxR2_norm = gp.getPhi_UV(jres) ^ R2 / rHA;

                          DEV_sc_sc[2 * jres] -= dE_dr * (VxR2_norm * AH) * rHA;
                          DEV_sc_sc[2 * jres] -=
                              dE_dxD * (VxR2_norm * (HD - cAHD * AH));
                          DEV_sc_sc[2 * jres] +=
                              dE_dxH * ((gp.getPhi_UV(jres) ^ XA) * AH -
                                        VxR2_norm * (XA - cXAH * AH));

                        } else if (ires > jres) {
                          // psi of jres
                          vec R2 = theH->r - gp.getResidue(jres)->getC()->r;
                          vec VxR2_norm = gp.getPsi_UV(jres) ^ R2 / rHA;
                          DEV_sc_sc[2 * jres + 1] +=
                              dE_dr * (VxR2_norm * AH) * rHA;
                          DEV_sc_sc[2 * jres + 1] +=
                              dE_dxD * ((gp.getPsi_UV(jres) ^ HD) * AH +
                                        VxR2_norm * (HD - cAHD * AH));
                          DEV_sc_sc[2 * jres + 1] +=
                              dE_dxH * (VxR2_norm * (XA - cXAH * AH));
                          for (int ii = jres + 1; ii <= ires - 1; ii++) {
                            // phi of ii
                            R2 = theH->r - gp.getResidue(ii)->getCA()->r;
                            VxR2_norm = gp.getPhi_UV(ii) ^ R2 / rHA;
                            DEV_sc_sc[2 * ii] += dE_dr * (VxR2_norm * AH) * rHA;
                            DEV_sc_sc[2 * ii] +=
                                dE_dxD * ((gp.getPhi_UV(ii) ^ HD) * AH +
                                          VxR2_norm * (HD - cAHD * AH));
                            DEV_sc_sc[2 * ii] +=
                                dE_dxH * (VxR2_norm * (XA - cXAH * AH));

                            // psi of ii
                            R2 = theH->r - gp.getResidue(ii)->getC()->r;
                            VxR2_norm = gp.getPsi_UV(ii) ^ R2 / rHA;
                            DEV_sc_sc[2 * ii + 1] +=
                                dE_dr * (VxR2_norm * AH) * rHA;
                            DEV_sc_sc[2 * ii + 1] +=
                                dE_dxD * ((gp.getPsi_UV(ii) ^ HD) * AH +
                                          VxR2_norm * (HD - cAHD * AH));
                            DEV_sc_sc[2 * ii + 1] +=
                                dE_dxH * (VxR2_norm * (XA - cXAH * AH));
                          }
                          // phi of ires
                          R2 = theH->r - gp.getResidue(ires)->getCA()->r;
                          VxR2_norm = gp.getPhi_UV(ires) ^ R2 / rHA;
                          DEV_sc_sc[2 * ires] += dE_dr * (VxR2_norm * AH) * rHA;
                          DEV_sc_sc[2 * ires] +=
                              dE_dxD * ((gp.getPhi_UV(ires) ^ HD) * AH +
                                        VxR2_norm * (HD - cAHD * AH));
                          DEV_sc_sc[2 * ires] +=
                              dE_dxH * (VxR2_norm * (XA - cXAH * AH));
                        }
                      }
                    } else { // bb-sc
                      isSC = 1;
                      if ((E = getHBE_DEV(dSEQ, isSC, theHBA->getHybrid(), rHA,
                                          cXAH, cAHD, dE_dr, dE_dxH, dE_dxD)) <
                          0) {
                        // cout << rHA << " " << cXAH <<
                        // " " << cAHD << endl;
                        theHBA->addH(theH, E);
                        theH->hbonded_acceptor = theHBA;
                        EHB_bb_sc += E;
                        bonded = 1;

                        // cout << E << " " << ires << "
                        // " << jres << endl;
                        if (ires < jres) {
                          vec R2;
                          vec VxR2_norm;
                          // phi of ires
                          if (!theH->isSC()) {
                            R2 = theA->r - gp.getResidue(ires)->getCA()->r;
                            VxR2_norm = gp.getPhi_UV(ires) ^ R2 / rHA;

                            DEV_bb_sc[2 * ires] -=
                                dE_dr * (VxR2_norm * AH) * rHA;
                            DEV_bb_sc[2 * ires] -=
                                dE_dxD * (VxR2_norm * (HD - cAHD * AH));
                            DEV_bb_sc[2 * ires] +=
                                dE_dxH * ((gp.getPhi_UV(ires) ^ XA) * AH -
                                          VxR2_norm * (XA - cXAH * AH));
                          }
                          // psi of ires
                          R2 = theA->r - gp.getResidue(ires)->getC()->r;
                          VxR2_norm = gp.getPsi_UV(ires) ^ R2 / rHA;

                          DEV_bb_sc[2 * ires + 1] -=
                              dE_dr * (VxR2_norm * AH) * rHA;
                          DEV_bb_sc[2 * ires + 1] -=
                              dE_dxD * (VxR2_norm * (HD - cAHD * AH));
                          DEV_bb_sc[2 * ires + 1] +=
                              dE_dxH * ((gp.getPsi_UV(ires) ^ XA) * AH -
                                        VxR2_norm * (XA - cXAH * AH));

                          for (int ii = ires + 1; ii <= jres - 1; ii++) {
                            // phi of ii
                            R2 = theA->r - gp.getResidue(ii)->getCA()->r;
                            VxR2_norm = gp.getPhi_UV(ii) ^ R2 / rHA;

                            DEV_bb_sc[2 * ii] -= dE_dr * (VxR2_norm * AH) * rHA;
                            DEV_bb_sc[2 * ii] -=
                                dE_dxD * (VxR2_norm * (HD - cAHD * AH));
                            DEV_bb_sc[2 * ii] +=
                                dE_dxH * ((gp.getPhi_UV(ii) ^ XA) * AH -
                                          VxR2_norm * (XA - cXAH * AH));
                            // psi of ii
                            R2 = theA->r - gp.getResidue(ii)->getC()->r;
                            VxR2_norm = gp.getPsi_UV(ii) ^ R2 / rHA;

                            DEV_bb_sc[2 * ii + 1] -=
                                dE_dr * (VxR2_norm * AH) * rHA;
                            DEV_bb_sc[2 * ii + 1] -=
                                dE_dxD * (VxR2_norm * (HD - cAHD * AH));
                            DEV_bb_sc[2 * ii + 1] +=
                                dE_dxH * ((gp.getPsi_UV(ii) ^ XA) * AH -
                                          VxR2_norm * (XA - cXAH * AH));
                          }
                          // phi of jres
                          R2 = theA->r - gp.getResidue(jres)->getCA()->r;
                          VxR2_norm = gp.getPhi_UV(jres) ^ R2 / rHA;

                          DEV_bb_sc[2 * jres] -= dE_dr * (VxR2_norm * AH) * rHA;
                          DEV_bb_sc[2 * jres] -=
                              dE_dxD * (VxR2_norm * (HD - cAHD * AH));
                          DEV_bb_sc[2 * jres] +=
                              dE_dxH * ((gp.getPhi_UV(jres) ^ XA) * AH -
                                        VxR2_norm * (XA - cXAH * AH));
                          // psi of jres
                          if (theH->isSC()) {
                            R2 = theA->r - gp.getResidue(jres)->getC()->r;
                            VxR2_norm = gp.getPsi_UV(jres) ^ R2 / rHA;

                            DEV_bb_sc[2 * jres + 1] -=
                                dE_dr * (VxR2_norm * AH) * rHA;
                            DEV_bb_sc[2 * jres + 1] -=
                                dE_dxD * (VxR2_norm * (HD - cAHD * AH));
                            DEV_bb_sc[2 * jres + 1] +=
                                dE_dxH * ((gp.getPsi_UV(jres) ^ XA) * AH -
                                          VxR2_norm * (XA - cXAH * AH));
                          }
                        } else if (ires > jres) {
                          // psi of jres
                          vec R2;
                          vec VxR2_norm;

                          if (!theH->isSC()) {
                            R2 = theH->r - gp.getResidue(jres)->getC()->r;
                            VxR2_norm = gp.getPsi_UV(jres) ^ R2 / rHA;

                            DEV_bb_sc[2 * jres + 1] +=
                                dE_dr * (VxR2_norm * AH) * rHA;
                            DEV_bb_sc[2 * jres + 1] +=
                                dE_dxD * ((gp.getPsi_UV(jres) ^ HD) * AH +
                                          VxR2_norm * (HD - cAHD * AH));
                            DEV_bb_sc[2 * jres + 1] +=
                                dE_dxH * (VxR2_norm * (XA - cXAH * AH));
                          }
                          for (int ii = jres + 1; ii <= ires - 1; ii++) {
                            // phi of ii
                            R2 = theH->r - gp.getResidue(ii)->getCA()->r;
                            VxR2_norm = gp.getPhi_UV(ii) ^ R2 / rHA;
                            DEV_bb_sc[2 * ii] += dE_dr * (VxR2_norm * AH) * rHA;
                            DEV_bb_sc[2 * ii] +=
                                dE_dxD * ((gp.getPhi_UV(ii) ^ HD) * AH +
                                          VxR2_norm * (HD - cAHD * AH));
                            DEV_bb_sc[2 * ii] +=
                                dE_dxH * (VxR2_norm * (XA - cXAH * AH));

                            // psi of ii
                            R2 = theH->r - gp.getResidue(ii)->getC()->r;
                            VxR2_norm = gp.getPsi_UV(ii) ^ R2 / rHA;
                            DEV_bb_sc[2 * ii + 1] +=
                                dE_dr * (VxR2_norm * AH) * rHA;
                            DEV_bb_sc[2 * ii + 1] +=
                                dE_dxD * ((gp.getPsi_UV(ii) ^ HD) * AH +
                                          VxR2_norm * (HD - cAHD * AH));
                            DEV_bb_sc[2 * ii + 1] +=
                                dE_dxH * (VxR2_norm * (XA - cXAH * AH));
                          }
                          // phi of ires
                          if (theH->isSC()) {
                            R2 = theH->r - gp.getResidue(ires)->getCA()->r;
                            VxR2_norm = gp.getPhi_UV(ires) ^ R2 / rHA;
                            DEV_bb_sc[2 * ires] +=
                                dE_dr * (VxR2_norm * AH) * rHA;
                            DEV_bb_sc[2 * ires] +=
                                dE_dxD * ((gp.getPhi_UV(ires) ^ HD) * AH +
                                          VxR2_norm * (HD - cAHD * AH));
                            DEV_bb_sc[2 * ires] +=
                                dE_dxH * (VxR2_norm * (XA - cXAH * AH));
                          }
                        }
                      }
                    }
                  }
                } // end of saturation
                else {
                  // cout << "saturated" << endl;
                }

                if (bonded)
                  break;
                pt = pt->fine_cell_next;
              }
              if (bonded)
                break;
            }
            if (bonded)
              break;
          }
          if (bonded)
            break;
        }
      } else {
        cout << "it should be removed" << endl;
        exit(1);
      }
    }
  }

  // for(int i=0; i<gp.size(); i++){
  //  cout << i << " " <<  DEV_bb_bb[2*i] << " " << DEV_bb_sc[2*i] << " " <<
  //  DEV_sc_sc[2*i] << endl; cout << i << " " << DEV_bb_bb[2*i+1] << " " <<
  //  DEV_bb_sc[2*i+1] << " " << DEV_sc_sc[2*i+1] << endl;
  //}
}

/*********************************************************************
To calculate the hydrogen bonds for the whole protein. To avoid over counting,
we only search the protons.
*********************************************************************/
void gen_fine_grid::getEHB_FYO_DEV(double &EHB_bb_bb, double &EHB_bb_sc,
                                   double &EHB_sc_sc, double *DEV_bb_bb,
                                   double *DEV_bb_sc, double *DEV_sc_sc) {
  residue *theRes = NULL;
  hbond_acceptor *theHBA = NULL;
  atom *theH = NULL;
  atom *theD = NULL;
  atom *theA = NULL;
  atom *theX = NULL;
  atom *pt;
  int cindex;
  int ix, iy, iz, cx, cy, cz, kx, ky, kz;
  double rHA = 0;
  double cXAH = -1;
  double cAHD = -1;
  int dSEQ = 0;
  int isSC = 0;
  double E = 0;
  int bonded = 0;
  double dE_dr;
  double dE_dxH;
  double dE_dxD;

  // protons
  // int type=theRes->getType();
  // int nCHI = theRes->getNChi();
  // cout << nCHI << endl;

  EHB_bb_bb = EHB_bb_sc = EHB_sc_sc = 0;
  for (int i = 0; i < gp.size(); i++) {
    DEV_bb_bb[3 * i] = DEV_bb_bb[3 * i + 1] = DEV_bb_bb[3 * i + 2] = 0;
    DEV_bb_sc[3 * i] = DEV_bb_sc[3 * i + 1] = DEV_bb_sc[3 * i + 2] = 0;
    DEV_sc_sc[3 * i] = DEV_sc_sc[3 * i + 1] = DEV_sc_sc[3 * i + 2] = 0;
  }
  int jres;
  for (int ires = 0; ires < gp.size(); ires++) {
    theRes = gp.getResidue(ires);
    atom **array = theRes->getPH_ARRAY();
    for (int iH = 0; iH < theRes->getNPH(); iH++) {
      // cout << ires << " " << theRes->getNPH() << endl;
      theH = array[iH];
      theD = theH->attachedHA;
      bonded = 0;
      if (!theH->isHBonded()) {
        // get the cooridiate of the FineGrid(cell)
        ix = static_cast<int>((theH->r.x - min.x) / cell_box.x);
        iy = static_cast<int>((theH->r.y - min.y) / cell_box.y);
        iz = static_cast<int>((theH->r.z - min.z) / cell_box.z);
        cindex = ix + (iy + iz * ncell[1]) * ncell[0];
        if (theH->fine_cell_index != cindex) { // precaucious
          cout << "fatal error in the cell index of protons" << endl;
          cout << theH->resID << endl;
          exit(1);
        }
        // search the possible acceptors and calculate the energy
        for (cz = iz - 1; cz <= iz + 1; cz++) {
          kz = cz;
          if (cz < 0)
            kz += ncell[2];
          else if (cz >= ncell[2])
            kz -= ncell[2];

          for (cy = iy - 1; cy <= iy + 1; cy++) {
            ky = cy;
            if (cy < 0)
              ky += ncell[1];
            else if (cy >= ncell[1])
              ky -= ncell[1];

            for (cx = ix - 1; cx <= ix + 1; cx++) {
              kx = cx;
              if (cx < 0)
                kx += ncell[0];
              else if (cx >= ncell[0])
                kx -= ncell[0];

              cindex = ncell[0] * (ncell[1] * kz + ky) + kx;
              pt = cell_acceptors[cindex];

              while (pt) {
                // calculate the hbond energy
                theHBA = pt->acceptor;
                jres = pt->resID;
                // precaucious of the accpetor
                if (!theHBA) {
                  cout << "fatal error: Not a HBacceptor!" << endl;
                  exit(1);
                }

                if (!theHBA->isSaturated()) { // not saturated
                  theA = theHBA->getA();
                  theX = theHBA->getX();
                  rHA = theH->r.getDist(theA->r);
                  vec XA = theA->r - theX->r;
                  XA.Normalize();
                  vec AH = theH->r - theA->r;
                  AH.Normalize();
                  vec HD = theD->r - theH->r;
                  HD.Normalize();
                  cXAH = XA * AH;
                  cAHD = AH * HD;
                  dSEQ = abs(theH->resID - theA->resID);
                  if (rHA < MAX_R && rHA > MIN_R && cXAH < MAX_CXAH &&
                      cXAH > MIN_CXAH && cAHD < MAX_CAHD && cAHD > MIN_CAHD &&
                      dSEQ > 0) {
                    if (!theA->isSC() && !theH->isSC()) { // bb-bb
                      isSC = 0;
                      if ((E = getHBE_DEV(dSEQ, isSC, theHBA->getHybrid(), rHA,
                                          cXAH, cAHD, dE_dr, dE_dxH, dE_dxD)) <
                          0) {
                        theHBA->addH(theH, E);
                        theH->hbonded_acceptor = theHBA;
                        EHB_bb_bb += E;
                        bonded = 1;
                        // cout <<"FYO: " <<  E << " "
                        // << ires << " " << jres <<
                        // endl;
                        if (ires < jres) {
                          // psi of ires
                          vec R2;
                          vec VxR2_norm;
                          for (int ii = ires; ii <= jres; ii++) {
                            // phi of ii
                            R2 = theA->r - gp.getResidue(ii)->getCA()->r;
                            VxR2_norm = gp.getPhi_UV(ii) ^ R2 / rHA;
                            DEV_bb_bb[3 * ii] -= dE_dr * (VxR2_norm * AH) * rHA;
                            DEV_bb_bb[3 * ii] -=
                                dE_dxD * (VxR2_norm * (HD - cAHD * AH));
                            DEV_bb_bb[3 * ii] +=
                                dE_dxH * ((gp.getPhi_UV(ii) ^ XA) * AH -
                                          VxR2_norm * (XA - cXAH * AH));
                            // psi of ii
                            R2 = theA->r - gp.getResidue(ii)->getC()->r;
                            VxR2_norm = gp.getPsi_UV(ii) ^ R2 / rHA;

                            DEV_bb_bb[3 * ii + 1] -=
                                dE_dr * (VxR2_norm * AH) * rHA;
                            DEV_bb_bb[3 * ii + 1] -=
                                dE_dxD * (VxR2_norm * (HD - cAHD * AH));
                            DEV_bb_bb[3 * ii + 1] +=
                                dE_dxH * ((gp.getPsi_UV(ii) ^ XA) * AH -
                                          VxR2_norm * (XA - cXAH * AH));

                            // omiga of ii (<jres)
                            if (ii < jres) {
                              R2 = theA->r - gp.getResidue(ii + 1)->getN()->r;
                              VxR2_norm = gp.getOmiga_UV(ii) ^ R2 / rHA;

                              DEV_bb_bb[3 * ii + 2] -=
                                  dE_dr * (VxR2_norm * AH) * rHA;
                              DEV_bb_bb[3 * ii + 2] -=
                                  dE_dxD * (VxR2_norm * (HD - cAHD * AH));
                              DEV_bb_bb[3 * ii + 2] +=
                                  dE_dxH * ((gp.getOmiga_UV(ii) ^ XA) * AH -
                                            VxR2_norm * (XA - cXAH * AH));
                            }

                            // cout << dE_dr << " "
                            // << dE_dxD << " " <<
                            // dE_dxH << endl;
                          }
                        } else if (ires > jres) {
                          vec R2;
                          vec VxR2_norm;
                          // omiga of jres
                          R2 = theH->r - gp.getResidue(jres + 1)->getN()->r;
                          VxR2_norm = gp.getOmiga_UV(jres) ^ R2 / rHA;
                          DEV_bb_bb[3 * jres + 2] +=
                              dE_dr * (VxR2_norm * AH) * rHA;
                          DEV_bb_bb[3 * jres + 2] +=
                              dE_dxD * ((gp.getOmiga_UV(jres) ^ HD) * AH +
                                        VxR2_norm * (HD - cAHD * AH));
                          DEV_bb_bb[3 * jres + 2] +=
                              dE_dxH * (VxR2_norm * (XA - cXAH * AH));

                          for (int ii = jres + 1; ii <= ires - 1; ii++) {
                            // phi of ii
                            R2 = theH->r - gp.getResidue(ii)->getCA()->r;
                            VxR2_norm = gp.getPhi_UV(ii) ^ R2 / rHA;
                            DEV_bb_bb[3 * ii] += dE_dr * (VxR2_norm * AH) * rHA;
                            DEV_bb_bb[3 * ii] +=
                                dE_dxD * ((gp.getPhi_UV(ii) ^ HD) * AH +
                                          VxR2_norm * (HD - cAHD * AH));
                            DEV_bb_bb[3 * ii] +=
                                dE_dxH * (VxR2_norm * (XA - cXAH * AH));

                            // psi of ii
                            R2 = theH->r - gp.getResidue(ii)->getC()->r;
                            VxR2_norm = gp.getPsi_UV(ii) ^ R2 / rHA;
                            DEV_bb_bb[3 * ii + 1] +=
                                dE_dr * (VxR2_norm * AH) * rHA;
                            DEV_bb_bb[3 * ii + 1] +=
                                dE_dxD * ((gp.getPsi_UV(ii) ^ HD) * AH +
                                          VxR2_norm * (HD - cAHD * AH));
                            DEV_bb_bb[3 * ii + 1] +=
                                dE_dxH * (VxR2_norm * (XA - cXAH * AH));

                            // omiga of ii
                            R2 = theH->r - gp.getResidue(ii + 1)->getN()->r;
                            VxR2_norm = gp.getOmiga_UV(ii) ^ R2 / rHA;
                            DEV_bb_bb[3 * ii + 2] +=
                                dE_dr * (VxR2_norm * AH) * rHA;
                            DEV_bb_bb[3 * ii + 2] +=
                                dE_dxD * ((gp.getOmiga_UV(ii) ^ HD) * AH +
                                          VxR2_norm * (HD - cAHD * AH));
                            DEV_bb_bb[3 * ii + 2] +=
                                dE_dxH * (VxR2_norm * (XA - cXAH * AH));
                          }
                        }
                      }
                    } else if (theH->isSC() && theA->isSC()) { // sc-sc
                      isSC = 1;
                      if ((E = getHBE_DEV(dSEQ, isSC, theHBA->getHybrid(), rHA,
                                          cXAH, cAHD, dE_dr, dE_dxH, dE_dxD)) <
                          0) {
                        theH->hbonded_acceptor = theHBA;
                        theHBA->addH(theH, E);
                        EHB_sc_sc += E;
                        bonded = 1;

                        // cout << E << " " << ires << "
                        // " << jres << endl;
                        if (ires < jres) {
                          // psi of ires
                          vec R2 = theA->r - gp.getResidue(ires)->getC()->r;
                          vec VxR2_norm = gp.getPsi_UV(ires) ^ R2 / rHA;

                          DEV_sc_sc[3 * ires + 1] -=
                              dE_dr * (VxR2_norm * AH) * rHA;
                          DEV_sc_sc[3 * ires + 1] -=
                              dE_dxD * (VxR2_norm * (HD - cAHD * AH));
                          DEV_sc_sc[3 * ires + 1] +=
                              dE_dxH * ((gp.getPsi_UV(ires) ^ XA) * AH -
                                        VxR2_norm * (XA - cXAH * AH));
                          // omiga of ires
                          R2 = theA->r - gp.getResidue(ires + 1)->getN()->r;
                          VxR2_norm = gp.getOmiga_UV(ires) ^ R2 / rHA;

                          DEV_sc_sc[3 * ires + 2] -=
                              dE_dr * (VxR2_norm * AH) * rHA;
                          DEV_sc_sc[3 * ires + 2] -=
                              dE_dxD * (VxR2_norm * (HD - cAHD * AH));
                          DEV_sc_sc[3 * ires + 2] +=
                              dE_dxH * ((gp.getOmiga_UV(ires) ^ XA) * AH -
                                        VxR2_norm * (XA - cXAH * AH));

                          for (int ii = ires + 1; ii <= jres - 1; ii++) {
                            // phi of ii
                            R2 = theA->r - gp.getResidue(ii)->getCA()->r;
                            VxR2_norm = gp.getPhi_UV(ii) ^ R2 / rHA;

                            DEV_sc_sc[3 * ii] -= dE_dr * (VxR2_norm * AH) * rHA;
                            DEV_sc_sc[3 * ii] -=
                                dE_dxD * (VxR2_norm * (HD - cAHD * AH));
                            DEV_sc_sc[3 * ii] +=
                                dE_dxH * ((gp.getPhi_UV(ii) ^ XA) * AH -
                                          VxR2_norm * (XA - cXAH * AH));
                            // psi of ii
                            R2 = theA->r - gp.getResidue(ii)->getC()->r;
                            VxR2_norm = gp.getPsi_UV(ii) ^ R2 / rHA;

                            DEV_sc_sc[3 * ii + 1] -=
                                dE_dr * (VxR2_norm * AH) * rHA;
                            DEV_sc_sc[3 * ii + 1] -=
                                dE_dxD * (VxR2_norm * (HD - cAHD * AH));
                            DEV_sc_sc[3 * ii + 1] +=
                                dE_dxH * ((gp.getPsi_UV(ii) ^ XA) * AH -
                                          VxR2_norm * (XA - cXAH * AH));
                            // omiga of ii
                            R2 = theA->r - gp.getResidue(ii + 1)->getN()->r;
                            VxR2_norm = gp.getOmiga_UV(ii) ^ R2 / rHA;

                            DEV_sc_sc[3 * ii + 2] -=
                                dE_dr * (VxR2_norm * AH) * rHA;
                            DEV_sc_sc[3 * ii + 2] -=
                                dE_dxD * (VxR2_norm * (HD - cAHD * AH));
                            DEV_sc_sc[3 * ii + 2] +=
                                dE_dxH * ((gp.getOmiga_UV(ii) ^ XA) * AH -
                                          VxR2_norm * (XA - cXAH * AH));
                          }
                          // phi of jres
                          R2 = theA->r - gp.getResidue(jres)->getCA()->r;
                          VxR2_norm = gp.getPhi_UV(jres) ^ R2 / rHA;

                          DEV_sc_sc[3 * jres] -= dE_dr * (VxR2_norm * AH) * rHA;
                          DEV_sc_sc[3 * jres] -=
                              dE_dxD * (VxR2_norm * (HD - cAHD * AH));
                          DEV_sc_sc[3 * jres] +=
                              dE_dxH * ((gp.getPhi_UV(jres) ^ XA) * AH -
                                        VxR2_norm * (XA - cXAH * AH));

                        } else if (ires > jres) {
                          // psi of jres
                          vec R2 = theH->r - gp.getResidue(jres)->getC()->r;
                          vec VxR2_norm = gp.getPsi_UV(jres) ^ R2 / rHA;
                          DEV_sc_sc[3 * jres + 1] +=
                              dE_dr * (VxR2_norm * AH) * rHA;
                          DEV_sc_sc[3 * jres + 1] +=
                              dE_dxD * ((gp.getPsi_UV(jres) ^ HD) * AH +
                                        VxR2_norm * (HD - cAHD * AH));
                          DEV_sc_sc[3 * jres + 1] +=
                              dE_dxH * (VxR2_norm * (XA - cXAH * AH));
                          // omiga of jres
                          R2 = theH->r - gp.getResidue(jres + 1)->getN()->r;
                          VxR2_norm = gp.getOmiga_UV(jres) ^ R2 / rHA;
                          DEV_sc_sc[3 * jres + 2] +=
                              dE_dr * (VxR2_norm * AH) * rHA;
                          DEV_sc_sc[3 * jres + 2] +=
                              dE_dxD * ((gp.getOmiga_UV(jres) ^ HD) * AH +
                                        VxR2_norm * (HD - cAHD * AH));
                          DEV_sc_sc[3 * jres + 2] +=
                              dE_dxH * (VxR2_norm * (XA - cXAH * AH));

                          for (int ii = jres + 1; ii <= ires - 1; ii++) {
                            // phi of ii
                            R2 = theH->r - gp.getResidue(ii)->getCA()->r;
                            VxR2_norm = gp.getPhi_UV(ii) ^ R2 / rHA;
                            DEV_sc_sc[3 * ii] += dE_dr * (VxR2_norm * AH) * rHA;
                            DEV_sc_sc[3 * ii] +=
                                dE_dxD * ((gp.getPhi_UV(ii) ^ HD) * AH +
                                          VxR2_norm * (HD - cAHD * AH));
                            DEV_sc_sc[3 * ii] +=
                                dE_dxH * (VxR2_norm * (XA - cXAH * AH));

                            // psi of ii
                            R2 = theH->r - gp.getResidue(ii)->getC()->r;
                            VxR2_norm = gp.getPsi_UV(ii) ^ R2 / rHA;
                            DEV_sc_sc[3 * ii + 1] +=
                                dE_dr * (VxR2_norm * AH) * rHA;
                            DEV_sc_sc[3 * ii + 1] +=
                                dE_dxD * ((gp.getPsi_UV(ii) ^ HD) * AH +
                                          VxR2_norm * (HD - cAHD * AH));
                            DEV_sc_sc[3 * ii + 1] +=
                                dE_dxH * (VxR2_norm * (XA - cXAH * AH));

                            // omiga of ii
                            R2 = theH->r - gp.getResidue(ii + 1)->getN()->r;
                            VxR2_norm = gp.getOmiga_UV(ii) ^ R2 / rHA;
                            DEV_sc_sc[3 * ii + 2] +=
                                dE_dr * (VxR2_norm * AH) * rHA;
                            DEV_sc_sc[3 * ii + 2] +=
                                dE_dxD * ((gp.getOmiga_UV(ii) ^ HD) * AH +
                                          VxR2_norm * (HD - cAHD * AH));
                            DEV_sc_sc[3 * ii + 2] +=
                                dE_dxH * (VxR2_norm * (XA - cXAH * AH));
                          }
                          // phi of ires
                          R2 = theH->r - gp.getResidue(ires)->getCA()->r;
                          VxR2_norm = gp.getPhi_UV(ires) ^ R2 / rHA;
                          DEV_sc_sc[3 * ires] += dE_dr * (VxR2_norm * AH) * rHA;
                          DEV_sc_sc[3 * ires] +=
                              dE_dxD * ((gp.getPhi_UV(ires) ^ HD) * AH +
                                        VxR2_norm * (HD - cAHD * AH));
                          DEV_sc_sc[3 * ires] +=
                              dE_dxH * (VxR2_norm * (XA - cXAH * AH));
                        }
                      }
                    } else { // bb-sc
                      isSC = 1;
                      if ((E = getHBE_DEV(dSEQ, isSC, theHBA->getHybrid(), rHA,
                                          cXAH, cAHD, dE_dr, dE_dxH, dE_dxD)) <
                          0) {
                        // cout << rHA << " " << cXAH <<
                        // " " << cAHD << endl;
                        theHBA->addH(theH, E);
                        theH->hbonded_acceptor = theHBA;
                        EHB_bb_sc += E;
                        bonded = 1;

                        // cout << E << " " << ires << "
                        // " << jres << endl;
                        if (ires < jres) {
                          vec R2;
                          vec VxR2_norm;
                          // phi of ires
                          if (!theH->isSC()) {
                            R2 = theA->r - gp.getResidue(ires)->getCA()->r;
                            VxR2_norm = gp.getPhi_UV(ires) ^ R2 / rHA;

                            DEV_bb_sc[3 * ires] -=
                                dE_dr * (VxR2_norm * AH) * rHA;
                            DEV_bb_sc[3 * ires] -=
                                dE_dxD * (VxR2_norm * (HD - cAHD * AH));
                            DEV_bb_sc[3 * ires] +=
                                dE_dxH * ((gp.getPhi_UV(ires) ^ XA) * AH -
                                          VxR2_norm * (XA - cXAH * AH));
                          }
                          // psi of ires
                          R2 = theA->r - gp.getResidue(ires)->getC()->r;
                          VxR2_norm = gp.getPsi_UV(ires) ^ R2 / rHA;

                          DEV_bb_sc[3 * ires + 1] -=
                              dE_dr * (VxR2_norm * AH) * rHA;
                          DEV_bb_sc[3 * ires + 1] -=
                              dE_dxD * (VxR2_norm * (HD - cAHD * AH));
                          DEV_bb_sc[3 * ires + 1] +=
                              dE_dxH * ((gp.getPsi_UV(ires) ^ XA) * AH -
                                        VxR2_norm * (XA - cXAH * AH));
                          // omiga of ires
                          R2 = theA->r - gp.getResidue(ires + 1)->getN()->r;
                          VxR2_norm = gp.getOmiga_UV(ires) ^ R2 / rHA;

                          DEV_bb_sc[3 * ires + 2] -=
                              dE_dr * (VxR2_norm * AH) * rHA;
                          DEV_bb_sc[3 * ires + 2] -=
                              dE_dxD * (VxR2_norm * (HD - cAHD * AH));
                          DEV_bb_sc[3 * ires + 2] +=
                              dE_dxH * ((gp.getOmiga_UV(ires) ^ XA) * AH -
                                        VxR2_norm * (XA - cXAH * AH));

                          for (int ii = ires + 1; ii <= jres - 1; ii++) {
                            // phi of ii
                            R2 = theA->r - gp.getResidue(ii)->getCA()->r;
                            VxR2_norm = gp.getPhi_UV(ii) ^ R2 / rHA;

                            DEV_bb_sc[3 * ii] -= dE_dr * (VxR2_norm * AH) * rHA;
                            DEV_bb_sc[3 * ii] -=
                                dE_dxD * (VxR2_norm * (HD - cAHD * AH));
                            DEV_bb_sc[3 * ii] +=
                                dE_dxH * ((gp.getPhi_UV(ii) ^ XA) * AH -
                                          VxR2_norm * (XA - cXAH * AH));
                            // psi of ii
                            R2 = theA->r - gp.getResidue(ii)->getC()->r;
                            VxR2_norm = gp.getPsi_UV(ii) ^ R2 / rHA;

                            DEV_bb_sc[3 * ii + 1] -=
                                dE_dr * (VxR2_norm * AH) * rHA;
                            DEV_bb_sc[3 * ii + 1] -=
                                dE_dxD * (VxR2_norm * (HD - cAHD * AH));
                            DEV_bb_sc[3 * ii + 1] +=
                                dE_dxH * ((gp.getPsi_UV(ii) ^ XA) * AH -
                                          VxR2_norm * (XA - cXAH * AH));
                            // omiga of ii
                            R2 = theA->r - gp.getResidue(ii + 1)->getN()->r;
                            VxR2_norm = gp.getOmiga_UV(ii) ^ R2 / rHA;

                            DEV_bb_sc[3 * ii + 2] -=
                                dE_dr * (VxR2_norm * AH) * rHA;
                            DEV_bb_sc[3 * ii + 2] -=
                                dE_dxD * (VxR2_norm * (HD - cAHD * AH));
                            DEV_bb_sc[3 * ii + 2] +=
                                dE_dxH * ((gp.getOmiga_UV(ii) ^ XA) * AH -
                                          VxR2_norm * (XA - cXAH * AH));
                          }
                          // phi of jres
                          R2 = theA->r - gp.getResidue(jres)->getCA()->r;
                          VxR2_norm = gp.getPhi_UV(jres) ^ R2 / rHA;

                          DEV_bb_sc[3 * jres] -= dE_dr * (VxR2_norm * AH) * rHA;
                          DEV_bb_sc[3 * jres] -=
                              dE_dxD * (VxR2_norm * (HD - cAHD * AH));
                          DEV_bb_sc[3 * jres] +=
                              dE_dxH * ((gp.getPhi_UV(jres) ^ XA) * AH -
                                        VxR2_norm * (XA - cXAH * AH));
                          // psi of jres
                          if (theH->isSC()) {
                            R2 = theA->r - gp.getResidue(jres)->getC()->r;
                            VxR2_norm = gp.getPsi_UV(jres) ^ R2 / rHA;

                            DEV_bb_sc[3 * jres + 1] -=
                                dE_dr * (VxR2_norm * AH) * rHA;
                            DEV_bb_sc[3 * jres + 1] -=
                                dE_dxD * (VxR2_norm * (HD - cAHD * AH));
                            DEV_bb_sc[3 * jres + 1] +=
                                dE_dxH * ((gp.getPsi_UV(jres) ^ XA) * AH -
                                          VxR2_norm * (XA - cXAH * AH));
                          }
                        } else if (ires > jres) {
                          // psi of jres
                          vec R2;
                          vec VxR2_norm;

                          if (!theH->isSC()) {
                            R2 = theH->r - gp.getResidue(jres)->getC()->r;
                            VxR2_norm = gp.getPsi_UV(jres) ^ R2 / rHA;

                            DEV_bb_sc[3 * jres + 1] +=
                                dE_dr * (VxR2_norm * AH) * rHA;
                            DEV_bb_sc[3 * jres + 1] +=
                                dE_dxD * ((gp.getPsi_UV(jres) ^ HD) * AH +
                                          VxR2_norm * (HD - cAHD * AH));
                            // DEV_sc_sc[3*jres+1]
                            // += dE_dxH
                            // *(VxR2_norm*(XA-cXAH*AH));
                            DEV_bb_sc[3 * jres + 1] +=
                                dE_dxH * (VxR2_norm * (XA - cXAH * AH));
                          }
                          // omiga of jres
                          R2 = theH->r - gp.getResidue(jres + 1)->getN()->r;
                          VxR2_norm = gp.getOmiga_UV(jres) ^ R2 / rHA;

                          DEV_bb_sc[3 * jres + 2] +=
                              dE_dr * (VxR2_norm * AH) * rHA;
                          DEV_bb_sc[3 * jres + 2] +=
                              dE_dxD * ((gp.getOmiga_UV(jres) ^ HD) * AH +
                                        VxR2_norm * (HD - cAHD * AH));
                          // DEV_sc_sc[3*jres+2] +=
                          // dE_dxH
                          // *(VxR2_norm*(XA-cXAH*AH));
                          DEV_bb_sc[3 * jres + 2] +=
                              dE_dxH * (VxR2_norm * (XA - cXAH * AH));

                          for (int ii = jres + 1; ii <= ires - 1; ii++) {
                            // phi of ii
                            R2 = theH->r - gp.getResidue(ii)->getCA()->r;
                            VxR2_norm = gp.getPhi_UV(ii) ^ R2 / rHA;
                            DEV_bb_sc[3 * ii] += dE_dr * (VxR2_norm * AH) * rHA;
                            DEV_bb_sc[3 * ii] +=
                                dE_dxD * ((gp.getPhi_UV(ii) ^ HD) * AH +
                                          VxR2_norm * (HD - cAHD * AH));
                            DEV_bb_sc[3 * ii] +=
                                dE_dxH * (VxR2_norm * (XA - cXAH * AH));

                            // psi of ii
                            R2 = theH->r - gp.getResidue(ii)->getC()->r;
                            VxR2_norm = gp.getPsi_UV(ii) ^ R2 / rHA;
                            DEV_bb_sc[3 * ii + 1] +=
                                dE_dr * (VxR2_norm * AH) * rHA;
                            DEV_bb_sc[3 * ii + 1] +=
                                dE_dxD * ((gp.getPsi_UV(ii) ^ HD) * AH +
                                          VxR2_norm * (HD - cAHD * AH));
                            DEV_bb_sc[3 * ii + 1] +=
                                dE_dxH * (VxR2_norm * (XA - cXAH * AH));

                            // omiga of ii
                            R2 = theH->r - gp.getResidue(ii + 1)->getN()->r;
                            VxR2_norm = gp.getOmiga_UV(ii) ^ R2 / rHA;
                            DEV_bb_sc[3 * ii + 2] +=
                                dE_dr * (VxR2_norm * AH) * rHA;
                            DEV_bb_sc[3 * ii + 2] +=
                                dE_dxD * ((gp.getOmiga_UV(ii) ^ HD) * AH +
                                          VxR2_norm * (HD - cAHD * AH));
                            DEV_bb_sc[3 * ii + 2] +=
                                dE_dxH * (VxR2_norm * (XA - cXAH * AH));
                          }
                          // phi of ires
                          if (theH->isSC()) {
                            R2 = theH->r - gp.getResidue(ires)->getCA()->r;
                            VxR2_norm = gp.getPhi_UV(ires) ^ R2 / rHA;
                            DEV_bb_sc[3 * ires] +=
                                dE_dr * (VxR2_norm * AH) * rHA;
                            DEV_bb_sc[3 * ires] +=
                                dE_dxD * ((gp.getPhi_UV(ires) ^ HD) * AH +
                                          VxR2_norm * (HD - cAHD * AH));
                            DEV_bb_sc[3 * ires] +=
                                dE_dxH * (VxR2_norm * (XA - cXAH * AH));
                          }
                        }
                      }
                    }
                  }
                } // end of saturation
                else {
                  // cout << "saturated" << endl;
                }

                if (bonded)
                  break;
                pt = pt->fine_cell_next;
              }
              if (bonded)
                break;
            }
            if (bonded)
              break;
          }
          if (bonded)
            break;
        }
      } else {
        cout << "it should be removed" << endl;
        exit(1);
      }
    }
  }

  // for(int i=0; i<gp.size(); i++){
  //  cout << i << " " <<  DEV_bb_bb[2*i] << " " << DEV_bb_sc[2*i] << " " <<
  //  DEV_sc_sc[2*i] << endl; cout << i << " " << DEV_bb_bb[2*i+1] << " " <<
  //  DEV_bb_sc[2*i+1] << " " << DEV_sc_sc[2*i+1] << endl;
  //}
}

void gen_fine_grid::clearHBE() {
  double ehb_bb = 0, ehb_bs = 0, ehb_ss = 0;
  //  for(int i=0;i<ncell[0];i++){
  //  for(int j=0;j<ncell[1];j++){
  //  for(int k=0;k<ncell[2];k++){
  //  	int cindex= i+(j+k*ncell[1])*ncell[0];
  //  	if (atom* pt=cell_protons[cindex]){
  //  		cout << "cindex = "<< cindex << ":";
  //  		do{
  //  			cout << pt->resID << " "<< pt->r.x <<" "<<
  //  pt->r.y << " "<< pt->r.z; 			pt = pt->fine_cell_next;
  //  } while(pt); cout << endl;
  //  	}
  //  }
  //  }
  //  } //debug
  for (int i = 0; i < gp.size(); i++)
    clearHBE_RES(i, ehb_bb, ehb_bs, ehb_ss);
}

void gen_fine_grid::getEHB(double &EHB_bb_bb, double &EHB_bb_sc,
                           double &EHB_sc_sc) {
  residue *theRes = NULL;
  hbond_acceptor *theHBA = NULL;
  atom *theH = NULL;
  atom *theD = NULL;
  atom *theA = NULL;
  atom *theX = NULL;
  atom *pt;
  int cindex;
  int ix, iy, iz, cx, cy, cz, kx, ky, kz;
  double rHA = 0;
  double cXAH = -1;
  double cAHD = -1;
  int dSEQ = 0;
  int isSC = 0;
  double E = 0;
  int bonded = 0;
  double dE_dr;
  double dE_dxH;
  double dE_dxD;

  // protons
  // int type=theRes->getType();
  // int nCHI = theRes->getNChi();
  // cout << nCHI << endl;

  EHB_bb_bb = EHB_bb_sc = EHB_sc_sc = 0;
  int jres;
  for (int ires = 0; ires < gp.size(); ires++) {
    theRes = gp.getResidue(ires);
    atom **array = theRes->getPH_ARRAY();
    for (int iH = 0; iH < theRes->getNPH(); iH++) {
      // cout << ires << " " << theRes->getNPH() << endl;
      theH = array[iH];
      theD = theH->attachedHA;
      bonded = 0;
      if (!theH->isHBonded()) {
        // get the cooridiate of the FineGrid(cell)
        ix = static_cast<int>((theH->r.x - min.x) / cell_box.x);
        iy = static_cast<int>((theH->r.y - min.y) / cell_box.y);
        iz = static_cast<int>((theH->r.z - min.z) / cell_box.z);
        cindex = ix + (iy + iz * ncell[1]) * ncell[0];
        if (theH->fine_cell_index != cindex) { // precaucious
          cout << "getEHB: fatal error in the cell index of protons" << endl;
          exit(1);
        }
        // search the possible acceptors and calculate the energy
        for (cz = iz - 1; cz <= iz + 1; cz++) {
          kz = cz;
          if (cz < 0)
            kz += ncell[2];
          else if (cz >= ncell[2])
            kz -= ncell[2];

          for (cy = iy - 1; cy <= iy + 1; cy++) {
            ky = cy;
            if (cy < 0)
              ky += ncell[1];
            else if (cy >= ncell[1])
              ky -= ncell[1];

            for (cx = ix - 1; cx <= ix + 1; cx++) {
              kx = cx;
              if (cx < 0)
                kx += ncell[0];
              else if (cx >= ncell[0])
                kx -= ncell[0];

              cindex = ncell[0] * (ncell[1] * kz + ky) + kx;
              pt = cell_acceptors[cindex];

              while (pt) {
                // calculate the hbond energy
                theHBA = pt->acceptor;
                jres = pt->resID;
                // precaucious of the accpetor
                if (!theHBA) {
                  cout << "fatal error: Not a HBacceptor!" << endl;
                  exit(1);
                }

                if (!theHBA->isSaturated()) { // not saturated
                  theA = theHBA->getA();
                  theX = theHBA->getX();
                  rHA = theH->r.getDist(theA->r);
                  vec XA = theA->r - theX->r;
                  XA.Normalize();
                  vec AH = theH->r - theA->r;
                  AH.Normalize();
                  vec HD = theD->r - theH->r;
                  HD.Normalize();
                  cXAH = XA * AH;
                  cAHD = AH * HD;
                  dSEQ = abs(theH->resID - theA->resID);
                  if (rHA < MAX_R && rHA > MIN_R && cXAH < MAX_CXAH &&
                      cXAH > MIN_CXAH && cAHD < MAX_CAHD && cAHD > MIN_CAHD &&
                      dSEQ > 0) {
                    if (!theA->isSC() && !theH->isSC()) { // bb-bb
                      isSC = 0;
                      if ((E = getHBE(dSEQ, isSC, theHBA->getHybrid(), rHA,
                                      cXAH, cAHD)) < 0) {
                        theHBA->addH(theH, E);
                        theH->hbonded_acceptor = theHBA;
                        EHB_bb_bb += E;
                        bonded = 1;
                      }
                    } else if (theH->isSC() && theA->isSC()) { // sc-sc
                      isSC = 1;
                      if ((E = getHBE_DEV(dSEQ, isSC, theHBA->getHybrid(), rHA,
                                          cXAH, cAHD, dE_dr, dE_dxH, dE_dxD)) <
                          0) {
                        theH->hbonded_acceptor = theHBA;
                        theHBA->addH(theH, E);
                        EHB_sc_sc += E;
                        bonded = 1;
                      }
                    } else { // bb-sc
                      isSC = 1;
                      if ((E = getHBE_DEV(dSEQ, isSC, theHBA->getHybrid(), rHA,
                                          cXAH, cAHD, dE_dr, dE_dxH, dE_dxD)) <
                          0) {
                        // cout << rHA << " " << cXAH <<
                        // " " << cAHD << endl;
                        theHBA->addH(theH, E);
                        theH->hbonded_acceptor = theHBA;
                        EHB_bb_sc += E;
                        bonded = 1;
                      }
                    }
                  }
                } // end of saturation
                else {
                  // cout << "saturated" << endl;
                }

                if (bonded)
                  break;
                pt = pt->fine_cell_next;
              }
              if (bonded)
                break;
            }
            if (bonded)
              break;
          }
          if (bonded)
            break;
        }
      } else {
        cout << "it should be removed" << endl;
        exit(1);
      }
    }
  }

  // for(int i=0; i<gp.size(); i++){
  //  cout << i << " " <<  DEV_bb_bb[2*i] << " " << DEV_bb_sc[2*i] << " " <<
  //  DEV_sc_sc[2*i] << endl; cout << i << " " << DEV_bb_bb[2*i+1] << " " <<
  //  DEV_bb_sc[2*i+1] << " " << DEV_sc_sc[2*i+1] << endl;
  //}
}

void gen_fine_grid::getHH_EVR_RES_SC(int iRes, double &VDW_REP,
                                     double *VDWR_SS_1D, double *VDWR_SB_1D) {
  residue *theRes = gp.getResidue(iRes);
  atom *theH = NULL;
  atom *nxtH = NULL;
  int cindex;
  int ix, iy, iz, cx, cy, cz, kx, ky, kz;
  double d2, d, e;

  atom **array = theRes->getPH_ARRAY();
  for (int iH = 0; iH < theRes->getNPH(); iH++) {
    theH = array[iH];
    if (theH->isSC()) { // only consider the sidechain proton
      ix = static_cast<int>((theH->r.x - min.x) / cell_box.x);
      iy = static_cast<int>((theH->r.y - min.y) / cell_box.y);
      iz = static_cast<int>((theH->r.z - min.z) / cell_box.z);
      cindex = ix + (iy + iz * ncell[1]) * ncell[0];
      if (theH->fine_cell_index != cindex) { // precaucious
        cout << "fatal error in the cell index of protons" << endl;
        exit(1);
      }
      // search the possible acceptors and calculate the energy
      for (cz = iz - 1; cz <= iz + 1; cz++) {
        kz = cz;
        if (cz < 0)
          kz += ncell[2];
        else if (cz >= ncell[2])
          kz -= ncell[2];

        for (cy = iy - 1; cy <= iy + 1; cy++) {
          ky = cy;
          if (cy < 0)
            ky += ncell[1];
          else if (cy >= ncell[1])
            ky -= ncell[1];

          for (cx = ix - 1; cx <= ix + 1; cx++) {
            kx = cx;
            if (cx < 0)
              kx += ncell[0];
            else if (cx >= ncell[0])
              kx -= ncell[0];

            cindex = ncell[0] * (ncell[1] * kz + ky) + kx;
            nxtH = cell_protons[cindex];
            while (nxtH) {
              if (abs(theH->resID - nxtH->resID) >= 2) {
                d2 = (theH->r - nxtH->r).mod2();
                d = sqrt(d2);
                if (d2 < 4.0) { // within repulsion range
                  e = 10.0 - 5.0 * d;
                  VDW_REP += e;
                  if (VDWR_SS_1D && nxtH->isSC()) {
                    VDWR_SS_1D[nxtH->resID] += e;
                  } else if (VDWR_SB_1D && !nxtH->isSC()) {
                    VDWR_SB_1D[nxtH->resID] += e;
                  }
                }
              }
              nxtH = nxtH->fine_cell_next;
            }
          } // end x
        }   // end y
      }     // end z
    }       // end SC-proton IF
  }         // end of proton-array loop
}

void gen_fine_grid::getHH_EVR_DEV_SC(int iRes, double &VDWR, double *DEV_VDWR) {
  residue *theRes = gp.getResidue(iRes);
  atom *theH = NULL;
  atom *nxtH = NULL;
  int cindex;
  int ix, iy, iz, cx, cy, cz, kx, ky, kz;
  double d2, d, e, de;

  // protons
  int type = theRes->getType();
  int nCHI = theRes->getNChi();

  int hasRotH = 0;
  if (type == SER || type == THR || type == TYR || type == LYS) {
    hasRotH = 1;
  }

  atom **array = theRes->getPH_ARRAY();
  for (int iH = 0; iH < theRes->getNPH(); iH++) {
    theH = array[iH];
    if (theH->isSC()) { // only consider the sidechain proton
      ix = static_cast<int>((theH->r.x - min.x) / cell_box.x);
      iy = static_cast<int>((theH->r.y - min.y) / cell_box.y);
      iz = static_cast<int>((theH->r.z - min.z) / cell_box.z);
      cindex = ix + (iy + iz * ncell[1]) * ncell[0];
      if (theH->fine_cell_index != cindex) { // precaucious
        cout << "fatal error in the cell index of protons" << endl;
        exit(1);
      }
      // search the possible acceptors and calculate the energy
      for (cz = iz - 1; cz <= iz + 1; cz++) {
        kz = cz;
        if (cz < 0)
          kz += ncell[2];
        else if (cz >= ncell[2])
          kz -= ncell[2];

        for (cy = iy - 1; cy <= iy + 1; cy++) {
          ky = cy;
          if (cy < 0)
            ky += ncell[1];
          else if (cy >= ncell[1])
            ky -= ncell[1];

          for (cx = ix - 1; cx <= ix + 1; cx++) {
            kx = cx;
            if (cx < 0)
              kx += ncell[0];
            else if (cx >= ncell[0])
              kx -= ncell[0];

            cindex = ncell[0] * (ncell[1] * kz + ky) + kx;
            nxtH = cell_protons[cindex];
            while (nxtH) {
              if (abs(theH->resID - nxtH->resID) >= 2) {
                d2 = (theH->r - nxtH->r).mod2();
                d = sqrt(d2);
                if (d2 < 4.0) { // within repulsion range
                  e = 10.0 - 5.0 * d;
                  VDWR += e;
                  de = -5.0;
                  /*for all side-chain polar protons,
                    the affected rotamers are all-included*/
                  vec Rsm = theH->r - nxtH->r;
                  for (int ich = 0; ich < nCHI; ich++) {
                    DEV_VDWR[ich] += de *
                                     ((*theRes->getChiUV(ich) ^
                                       (theH->r - *theRes->getChiP2(ich))) *
                                      Rsm) /
                                     d;
                    ;
                  }
                  /*rotatable Protons are a little bit
                   * different*/
                  if (hasRotH) {
                    DEV_VDWR[nCHI] += de *
                                      ((*theRes->getChiUV(nCHI) ^
                                        (theH->r - *theRes->getChiP2(nCHI))) *
                                       Rsm) /
                                      d;
                    ;
                  }
                }
              }
              nxtH = nxtH->fine_cell_next;
            }
          } // end x
        }   // end y
      }     // end z
    }       // end SC-proton IF
  }         // end of proton-array loop
}

double getHBE_DEV(const int &dSEQ, const int &isSC,
                  const ACCEPTOR_HYBRID_t &hybrid, const double &rHA,
                  const double &cosXAH, const double &cosAHD, double &dE_dr,
                  double &dE_dxH, double &dE_dxD) {
  HB_STAT_t hb_stat;
  HB_TERM_t rHA_var, XAH_var, AHD_var;
  int SC_interp = 0;
  int rHA_interp = 0;
  int XAH_interp = 0;
  int AHD_interp = 0;
  rHA_var = allRHA;
  AHD_var = shortAHD;
  XAH_var = shortXAH;

  if (isSC) {
    if (hybrid == RING)
      hb_stat = SC_RING;
    else if (hybrid == SP2)
      hb_stat = SC_SP2;
    else if (hybrid == SP3)
      hb_stat = SC_SP3;
    else {
      cerr << "error, incorrect hybrid state of acceptor atom" << endl;
      exit(1);
    }
    rHA_var = allRHA;
    if (rHA > RHA_INTERP_MAX) {
      AHD_var = longAHD;
      XAH_var = longXAH;
    } else if (rHA > RHA_INTERP_MIN) {
      SC_interp = 1;
    }
  } else {
    if (dSEQ > 4)
      hb_stat = BB_SHEET;
    else if (dSEQ == 4)
      hb_stat = BB_HELIX;
    else if (dSEQ == 3)
      hb_stat = BB_HELIX;
    else
      return 0;
  }

  rHA_interp = (rHA > RHA_INTERP_EDGE);
  XAH_interp = (cosXAH < ANGLE_INTERP_EDGE);
  AHD_interp = (cosAHD < ANGLE_INTERP_EDGE);
  double E_rHA1 = 0, E_XAH1 = 0, E_AHD1 = 0;
  double E_rHA2 = 0, E_XAH2 = 0, E_AHD2 = 0;
  int table_index = -1;
  double dEr_dr1, dExD_dxD1, dExH_dxH1;
  double dExD_dxD2, dExH_dxH2;
  double dEr_dxH, dExD_dxH;
  double dEr_dxD, dExH_dxD;
  double dExH_dr, dExD_dr;
  if (rHA > maxr[hb_stat]) {
    E_rHA1 = 0;
    dEr_dr1 = 0;
  } else {
    table_index = POLY_TABLE[hb_stat][rHA_var];
    E_rHA1 = POLYS[table_index].getValue(rHA);
    dEr_dr1 = POLYS[table_index].getDerivative(rHA);
  }

  // AHD ---> E(AHD)
  table_index = POLY_TABLE[hb_stat][AHD_var];
  E_AHD1 = POLYS[table_index].getValue(cosAHD);
  dExD_dxD1 = POLYS[table_index].getDerivative(cosAHD);

  // XAH ---> E(XAH)
  table_index = POLY_TABLE[hb_stat][XAH_var];
  E_XAH1 = POLYS[table_index].getValue(cosXAH);
  dExH_dxH1 = POLYS[table_index].getDerivative(cosXAH);

  // evaluate the interpolations for RHA
  dExH_dr = 0;
  dExD_dr = 0;
  double frac;
  if (SC_interp) {
    // cout << "SC_INTER" << endl;
    frac = (RHA_INTERP_MAX - rHA) / (RHA_INTERP_MAX - RHA_INTERP_MIN);
    AHD_var = longAHD;
    XAH_var = longXAH;
    // AHD ---> E(AHD)
    table_index = POLY_TABLE[hb_stat][AHD_var];
    E_AHD2 = POLYS[table_index].getValue(cosAHD);
    dExD_dxD2 = POLYS[table_index].getDerivative(cosAHD);
    // XAH ---> E(XAH)
    table_index = POLY_TABLE[hb_stat][XAH_var];
    E_XAH2 = POLYS[table_index].getValue(cosXAH);
    dExH_dxH2 = POLYS[table_index].getDerivative(cosXAH);
    dExD_dr = (E_AHD2 - E_AHD1) / (RHA_INTERP_MAX - RHA_INTERP_MIN);
    dExH_dr = (E_XAH2 - E_XAH1) / (RHA_INTERP_MAX - RHA_INTERP_MIN);
    E_AHD1 = frac * E_AHD1 + (1.0 - frac) * E_AHD2;
    E_XAH1 = frac * E_XAH1 + (1.0 - frac) * E_XAH2;
    dExD_dxD1 = frac * dExD_dxD1 + (1.0 - frac) * dExD_dxD2;
    dExH_dxH1 = frac * dExH_dxH1 + (1.0 - frac) * dExH_dxH2;
  } else if (rHA_interp) {
    // cout << "RHA_INTER" << endl;
    double denominator;
    if (isSC) {
      // cout << isSC << endl;
      frac = (MAX_R - rHA) / (MAX_R - RHA_INTERP_MAX);
      denominator = 1.0 / (MAX_R - RHA_INTERP_MAX);
    } else {
      frac = (MAX_R - rHA) / (MAX_R - RHA_INTERP_EDGE);
      denominator = 1.0 / (MAX_R - RHA_INTERP_EDGE);
    }
    E_AHD2 = 0;
    E_XAH2 = 0;
    dExD_dxD2 = 0;
    dExH_dxH2 = 0;
    dExD_dr = (E_AHD2 - E_AHD1) * denominator;
    dExH_dr = (E_XAH2 - E_XAH1) * denominator;
    // cout << MAX_R - rHA << endl;
    E_AHD1 = frac * E_AHD1 + (1.0 - frac) * E_AHD2;
    E_XAH1 = frac * E_XAH1 + (1.0 - frac) * E_XAH2;

    dExD_dxD1 = frac * dExD_dxD1 + (1.0 - frac) * dExD_dxD2;
    dExH_dxH1 = frac * dExH_dxH1 + (1.0 - frac) * dExH_dxH2;
  }

  // evaluate the interpolations for RHA
  dEr_dxD = 0;
  dExH_dxD = 0;
  if (AHD_interp) {
    // cout << "AHD_INTER" << endl;
    E_rHA2 = 0;
    E_XAH2 = 0;
    double dEr_dr2 = 0;
    double dExH_dxH2 = 0;
    double dExH_dr2 = 0;
    frac = (MIN_CAHD - cosAHD) / (MIN_CAHD - ANGLE_INTERP_EDGE);
    dEr_dxD = (E_rHA2 - E_rHA1) / (MIN_CAHD - ANGLE_INTERP_EDGE);
    dExH_dxD = (E_XAH2 - E_XAH1) / (MIN_CAHD - ANGLE_INTERP_EDGE);

    E_rHA1 = frac * E_rHA1 + (1.0 - frac) * E_rHA2;
    E_XAH1 = frac * E_XAH1 + (1.0 - frac) * E_XAH2;

    dExH_dr = frac * dExH_dr + (1.0 - frac) * dExH_dr2;

    dEr_dr1 = frac * dEr_dr1 + (1.0 - frac) * dEr_dr2;
    dExH_dxH1 = frac * dExH_dxH1 + (1.0 - frac) * dExH_dxH2;
  }
  // evaluate the interpolations for RHA
  dEr_dxH = 0;
  dExD_dxH = 0;
  if (XAH_interp) {
    E_rHA2 = 0;
    E_AHD2 = 0;
    double dEr_dr2 = 0;
    double dEr_dxD2 = 0;
    double dExD_dxD2 = 0;
    double dExD_dr2 = 0;
    frac = (MIN_CXAH - cosXAH) / (MIN_CXAH - ANGLE_INTERP_EDGE);
    dEr_dxH = (E_rHA2 - E_rHA1) / (MIN_CXAH - ANGLE_INTERP_EDGE);
    dExD_dxH = (E_AHD2 - E_AHD1) / (MIN_CXAH - ANGLE_INTERP_EDGE);
    E_rHA1 = frac * E_rHA1 + (1.0 - frac) * E_rHA2;
    E_AHD1 = frac * E_AHD1 + (1.0 - frac) * E_AHD2;

    dEr_dr1 = frac * dEr_dr1 + (1.0 - frac) * dEr_dr2;
    dEr_dxD = frac * dEr_dxD + (1.0 - frac) * dEr_dxD2;

    dExD_dxD1 = frac * dExD_dxD1 + (1.0 - frac) * dExD_dxD2;
    dExD_dr = frac * dExD_dr + (1.0 - frac) * dExD_dr2;
  }

  dE_dr = dEr_dr1 + dExH_dr + dExD_dr;
  dE_dxD = dExD_dxD1 + dEr_dxD + dExH_dxD;
  dE_dxH = dExH_dxH1 + dEr_dxH + dExD_dxH;
  return E_rHA1 + E_XAH1 + E_AHD1;
}

void gen_fine_grid::getEHB_DEV_BACKRUB(const int ires, double &EHB_bb_bb,
                                       double &EHB_bb_sc, double &EHB_sc_sc,
                                       double *DEV_bb_bb, double *DEV_bb_sc,
                                       double *DEV_sc_sc, vec *backrub_uv) {
  residue *theRes = gp.getResidue(ires);
  hbond_acceptor *theHBA = NULL;
  atom *theH = NULL;
  atom *theD = NULL;
  atom *theA = NULL;
  atom *theX = NULL;
  atom *pt;
  int cindex;
  int ix, iy, iz, cx, cy, cz, kx, ky, kz;
  double rHA = 0;
  double cXAH = -1;
  double cAHD = -1;
  int dSEQ = 0;
  int isSC = 0;
  double E = 0;
  int bonded = 0;
  double dE_dr;
  double dE_dxH;
  double dE_dxD;
  for (int i = 0; i < 3; i++) {
    DEV_bb_bb[i] = 0;
    DEV_bb_sc[i] = 0;
    DEV_sc_sc[i] = 0;
  }
  EHB_bb_bb = EHB_bb_sc = EHB_sc_sc = 0;

  double hbe_bb_bb, hbe_bb_sc, hbe_sc_sc;
  { // sidechain of ires
    atom **array = theRes->getPH_ARRAY();
    for (int iH = 0; iH < theRes->getNPH(); iH++) {
      theH = array[iH];
      if (theH->isSC()) {
        if (!theH->isHBonded()) { // the proton does not have a HBond
                                  // partner;
          bonded = 0;
          theD = theH->attachedHA; // get the Donar heavy atom

          // get the cooridiate of the FineGrid(cell)
          ix = static_cast<int>((theH->r.x - min.x) / cell_box.x);
          iy = static_cast<int>((theH->r.y - min.y) / cell_box.y);
          iz = static_cast<int>((theH->r.z - min.z) / cell_box.z);
          cindex = ix + (iy + iz * ncell[1]) * ncell[0];
          if (theH->fine_cell_index != cindex) { // precaucious
            cout << "fatal error in the cell index of protons" << endl;
            exit(1);
          }
          // search the possible acceptors and calculate the energy
          for (cz = iz - 1; cz <= iz + 1; cz++) {
            kz = cz;
            if (cz < 0)
              kz += ncell[2];
            else if (cz >= ncell[2])
              kz -= ncell[2];

            for (cy = iy - 1; cy <= iy + 1; cy++) {
              ky = cy;
              if (cy < 0)
                ky += ncell[1];
              else if (cy >= ncell[1])
                ky -= ncell[1];

              for (cx = ix - 1; cx <= ix + 1; cx++) {
                kx = cx;
                if (cx < 0)
                  kx += ncell[0];
                else if (cx >= ncell[0])
                  kx -= ncell[0];

                cindex = ncell[0] * (ncell[1] * kz + ky) + kx;
                pt = cell_acceptors[cindex];

                while (pt) {
                  // calculate the hbond energy
                  theHBA = pt->acceptor;
                  // precaucious of the accpetor
                  if (!theHBA) {
                    cout << "fatal error: Not a HBacceptor!" << endl;
                    exit(1);
                  }

                  if (!theHBA->isSaturated()) { // not saturated

                    theA = theHBA->getA();
                    theX = theHBA->getX();
                    rHA = theH->r.getDist(theA->r);
                    vec XA = theA->r - theX->r;
                    XA.Normalize();
                    vec AH = theH->r - theA->r;
                    AH.Normalize();
                    vec HD = theD->r - theH->r;
                    HD.Normalize();
                    cXAH = XA * AH;
                    cAHD = AH * HD;
                    dSEQ = abs(theH->resID - theA->resID);
                    if (rHA < MAX_R && rHA > MIN_R && cXAH < MAX_CXAH &&
                        cXAH > MIN_CXAH && cAHD < MAX_CAHD && cAHD > MIN_CAHD &&
                        dSEQ > 0) {
                      if (theA->isSC()) { // sc-sc
                        isSC = 1;
                        if ((E = getHBE_DEV(dSEQ, isSC, theHBA->getHybrid(),
                                            rHA, cXAH, cAHD, dE_dr, dE_dxH,
                                            dE_dxD)) < 0) {
                          theH->hbonded_acceptor = theHBA;
                          theHBA->addH(theH, E);
                          EHB_sc_sc += E;
                          bonded = 1;

                          { // for theta13
                            vec R2 =
                                theH->r - gp.getResidue(ires + 1)->getCA()->r;
                            vec VxR2_norm = backrub_uv[0] ^ R2 / rHA;
                            DEV_sc_sc[0] += dE_dr * (VxR2_norm * AH) * rHA;
                            DEV_sc_sc[0] +=
                                dE_dxD * ((backrub_uv[0] ^ HD) * AH +
                                          VxR2_norm * (HD - cAHD * AH));
                            DEV_sc_sc[0] +=
                                dE_dxH * (VxR2_norm * (XA - cXAH * AH));
                          }
                        }
                      } else { // bb-sc
                        isSC = 1;
                        if ((E = getHBE_DEV(dSEQ, isSC, theHBA->getHybrid(),
                                            rHA, cXAH, cAHD, dE_dr, dE_dxH,
                                            dE_dxD)) < 0) {
                          theHBA->addH(theH, E);
                          theH->hbonded_acceptor = theHBA;
                          EHB_bb_sc += E;
                          bonded = 1;

                          { // for theta13
                            vec R2 =
                                theH->r - gp.getResidue(ires + 1)->getCA()->r;
                            vec VxR2_norm = backrub_uv[0] ^ R2 / rHA;
                            DEV_bb_sc[0] += dE_dr * (VxR2_norm * AH) * rHA;
                            DEV_bb_sc[0] +=
                                dE_dxD * ((backrub_uv[0] ^ HD) * AH +
                                          VxR2_norm * (HD - cAHD * AH));
                            DEV_bb_sc[0] +=
                                dE_dxH * (VxR2_norm * (XA - cXAH * AH));
                          }
                        }
                      }
                    }
                  }
                  if (bonded)
                    break;
                  pt = pt->fine_cell_next;
                }
                if (bonded)
                  break;
              }
              if (bonded)
                break;
            }
            if (bonded)
              break;
          }
        } else {
          theHBA = theH->hbonded_acceptor;
          theA = theHBA->getA();
          E = theHBA->getE(theH);
          if (E == INF) {
            cerr << "fatal error in HBONDing pairs" << endl;
            exit(1);
          }
          if (theH->isSC() && theA->isSC())
            EHB_sc_sc += E;
          else if (!theA->isSC() && !theH->isSC())
            EHB_bb_bb += E;
          else
            EHB_bb_sc += E;
        }
      }
    }
    // acceptors
    array = theRes->getHBAcceptor_ARRAY();
    for (int iA = 0; iA < theRes->getNHBAcceptor(); iA++) {
      theA = array[iA];
      if (theA->isSC()) {
        theHBA = theA->acceptor;
        if (!theHBA) {
          cout << "error the accreptor is not right" << endl;
          exit(1);
        }
        theHBA->currentHBE(hbe_bb_bb, hbe_bb_sc, hbe_sc_sc);
        EHB_bb_bb += hbe_bb_bb;
        EHB_bb_sc += hbe_bb_sc;
        EHB_sc_sc += hbe_sc_sc;
        if (!theHBA->isSaturated()) { // not saturated
          // search the pissible un-hbonded protons
          theX = theHBA->getX();
          ix = static_cast<int>((theA->r.x - min.x) / cell_box.x);
          iy = static_cast<int>((theA->r.y - min.y) / cell_box.y);
          iz = static_cast<int>((theA->r.z - min.z) / cell_box.z);
          cindex = ix + (iy + iz * ncell[1]) * ncell[0];
          if (theA->fine_cell_index != cindex) { // precautions
            cout << "error in the accpetor cell index" << endl;
            exit(1);
          }
          // search teh neighboring cells to find the possible protons
          // for Hbond
          for (cz = iz - 1; cz <= iz + 1; cz++) {
            kz = cz;
            if (cz < 0)
              kz += ncell[2];
            else if (cz >= ncell[2])
              kz -= ncell[2];

            for (cy = iy - 1; cy <= iy + 1; cy++) {
              ky = cy;
              if (cy < 0)
                ky += ncell[1];
              else if (cy >= ncell[1])
                ky -= ncell[1];

              for (cx = ix - 1; cx <= ix + 1; cx++) {
                kx = cx;
                if (cx < 0)
                  kx += ncell[0];
                else if (cx >= ncell[0])
                  kx -= ncell[0];

                cindex = ncell[0] * (ncell[1] * kz + ky) + kx;
                pt = cell_protons[cindex];

                while (pt) {
                  theH = pt;
                  theD = theH->attachedHA; // get the Donar
                  // precaution
                  if (!theD) {
                    cout << "fatal error: Not a proton!" << endl;
                    exit(1);
                  }

                  if (!theH->isHBonded()) { // not hbonded
                    rHA = theH->r.getDist(theA->r);
                    vec XA = theA->r - theX->r;
                    XA.Normalize();
                    vec AH = theH->r - theA->r;
                    AH.Normalize();
                    vec HD = theD->r - theH->r;
                    HD.Normalize();
                    cXAH = XA * AH;
                    cAHD = AH * HD;
                    dSEQ = abs(theH->resID - theA->resID);
                    if (rHA < MAX_R && rHA > MIN_R && cXAH < MAX_CXAH &&
                        cXAH > MIN_CXAH && cAHD < MAX_CAHD && cAHD > MIN_CAHD &&
                        dSEQ > 0) {
                      if (theH->isSC()) { // sc-sc
                        isSC = 1;
                        if ((E = getHBE_DEV(dSEQ, isSC, theHBA->getHybrid(),
                                            rHA, cXAH, cAHD, dE_dr, dE_dxH,
                                            dE_dxD)) < 0) {
                          theH->hbonded_acceptor = theHBA;
                          theHBA->addH(theH, E);
                          EHB_sc_sc += E;
                          bonded = 1;

                          /*get the derivaste for all
                           * the acceptors*/
                          { // theta13
                            vec R2 =
                                theA->r - gp.getResidue(ires + 1)->getCA()->r;
                            vec VxR2_norm = backrub_uv[0] ^ R2 / rHA;

                            DEV_sc_sc[0] -= dE_dr * (VxR2_norm * AH) * rHA;
                            DEV_sc_sc[0] -=
                                dE_dxD * (VxR2_norm * (HD - cAHD * AH));
                            DEV_sc_sc[0] +=
                                dE_dxH * ((backrub_uv[0] ^ XA) * AH -
                                          VxR2_norm * (XA - cXAH * AH));
                          }
                        }
                      } else { // bb-sc
                        isSC = 1;
                        if ((E = getHBE_DEV(dSEQ, isSC, theHBA->getHybrid(),
                                            rHA, cXAH, cAHD, dE_dr, dE_dxH,
                                            dE_dxD)) < 0) {
                          theHBA->addH(theH, E);
                          theH->hbonded_acceptor = theHBA;
                          EHB_bb_sc += E;
                          bonded = 1;

                          { // theta13
                            vec R2 =
                                theA->r - gp.getResidue(ires + 1)->getCA()->r;
                            vec VxR2_norm = backrub_uv[0] ^ R2 / rHA;

                            DEV_bb_sc[0] -= dE_dr * (VxR2_norm * AH) * rHA;
                            DEV_bb_sc[0] -=
                                dE_dxD * (VxR2_norm * (HD - cAHD * AH));
                            DEV_bb_sc[0] +=
                                dE_dxH * ((backrub_uv[0] ^ XA) * AH -
                                          VxR2_norm * (XA - cXAH * AH));
                          }
                        }
                      }
                    }
                  }
                  if (theHBA->isSaturated())
                    break;
                  pt = pt->fine_cell_next;
                }
                if (theHBA->isSaturated())
                  break;
              }
              if (theHBA->isSaturated())
                break;
            }
            if (theHBA->isSaturated())
              break;
          }
        }
      }
    }
  } // end of SideChain of iRES

  { // O of ires-1; H of ires
    if (theH = theRes->getHN()) {
      if (!theH->isHBonded()) { // the proton does not have a HBond
                                // partner;
        bonded = 0;
        theD = theH->attachedHA; // get the Donar heavy atom

        // get the cooridiate of the FineGrid(cell)
        ix = static_cast<int>((theH->r.x - min.x) / cell_box.x);
        iy = static_cast<int>((theH->r.y - min.y) / cell_box.y);
        iz = static_cast<int>((theH->r.z - min.z) / cell_box.z);
        cindex = ix + (iy + iz * ncell[1]) * ncell[0];
        if (theH->fine_cell_index != cindex) { // precaucious
          cout << "fatal error in the cell index of protons" << endl;
          exit(1);
        }
        // search the possible acceptors and calculate the energy
        for (cz = iz - 1; cz <= iz + 1; cz++) {
          kz = cz;
          if (cz < 0)
            kz += ncell[2];
          else if (cz >= ncell[2])
            kz -= ncell[2];

          for (cy = iy - 1; cy <= iy + 1; cy++) {
            ky = cy;
            if (cy < 0)
              ky += ncell[1];
            else if (cy >= ncell[1])
              ky -= ncell[1];

            for (cx = ix - 1; cx <= ix + 1; cx++) {
              kx = cx;
              if (cx < 0)
                kx += ncell[0];
              else if (cx >= ncell[0])
                kx -= ncell[0];

              cindex = ncell[0] * (ncell[1] * kz + ky) + kx;
              pt = cell_acceptors[cindex];

              while (pt) {
                // calculate the hbond energy
                theHBA = pt->acceptor;
                // precaucious of the accpetor
                if (!theHBA) {
                  cout << "fatal error: Not a HBacceptor!" << endl;
                  exit(1);
                }

                if (!theHBA->isSaturated()) { // not saturated

                  theA = theHBA->getA();
                  theX = theHBA->getX();
                  rHA = theH->r.getDist(theA->r);
                  vec XA = theA->r - theX->r;
                  XA.Normalize();
                  vec AH = theH->r - theA->r;
                  AH.Normalize();
                  vec HD = theD->r - theH->r;
                  HD.Normalize();
                  cXAH = XA * AH;
                  cAHD = AH * HD;
                  dSEQ = abs(theH->resID - theA->resID);
                  if (rHA < MAX_R && rHA > MIN_R && cXAH < MAX_CXAH &&
                      cXAH > MIN_CXAH && cAHD < MAX_CAHD && cAHD > MIN_CAHD &&
                      dSEQ > 0) {
                    if (theA->isSC()) { // bb-sc
                      isSC = 1;
                      if ((E = getHBE_DEV(dSEQ, isSC, theHBA->getHybrid(), rHA,
                                          cXAH, cAHD, dE_dr, dE_dxH, dE_dxD)) <
                          0) {
                        theH->hbonded_acceptor = theHBA;
                        theHBA->addH(theH, E);
                        EHB_bb_sc += E;
                        bonded = 1;

                        {
                          // for theta13
                          vec R2 =
                              theH->r - gp.getResidue(ires + 1)->getCA()->r;
                          vec VxR2_norm = backrub_uv[0] ^ R2 / rHA;
                          DEV_bb_sc[0] += dE_dr * (VxR2_norm * AH) * rHA;
                          DEV_bb_sc[0] +=
                              dE_dxD * ((backrub_uv[0] ^ HD) * AH +
                                        VxR2_norm * (HD - cAHD * AH));
                          DEV_bb_sc[0] +=
                              dE_dxH * (VxR2_norm * (XA - cXAH * AH));
                          // for theta12
                          R2 = theH->r - gp.getResidue(ires)->getCA()->r;
                          VxR2_norm = backrub_uv[1] ^ R2 / rHA;
                          DEV_bb_sc[1] += dE_dr * (VxR2_norm * AH) * rHA;
                          DEV_bb_sc[1] +=
                              dE_dxD * ((backrub_uv[1] ^ HD) * AH +
                                        VxR2_norm * (HD - cAHD * AH));
                          DEV_bb_sc[1] +=
                              dE_dxH * (VxR2_norm * (XA - cXAH * AH));
                        }
                      }
                    } else { // bb-bb
                      isSC = 0;
                      if ((E = getHBE_DEV(dSEQ, isSC, theHBA->getHybrid(), rHA,
                                          cXAH, cAHD, dE_dr, dE_dxH, dE_dxD)) <
                          0) {
                        theHBA->addH(theH, E);
                        theH->hbonded_acceptor = theHBA;
                        EHB_bb_bb += E;
                        bonded = 1;

                        {
                          // for theta13
                          vec R2 =
                              theH->r - gp.getResidue(ires + 1)->getCA()->r;
                          vec VxR2_norm = backrub_uv[0] ^ R2 / rHA;
                          DEV_bb_bb[0] += dE_dr * (VxR2_norm * AH) * rHA;
                          DEV_bb_bb[0] +=
                              dE_dxD * ((backrub_uv[0] ^ HD) * AH +
                                        VxR2_norm * (HD - cAHD * AH));
                          DEV_bb_bb[0] +=
                              dE_dxH * (VxR2_norm * (XA - cXAH * AH));
                          // for theta12
                          R2 = theH->r - gp.getResidue(ires)->getCA()->r;
                          VxR2_norm = backrub_uv[1] ^ R2 / rHA;
                          DEV_bb_bb[1] += dE_dr * (VxR2_norm * AH) * rHA;
                          DEV_bb_bb[1] +=
                              dE_dxD * ((backrub_uv[1] ^ HD) * AH +
                                        VxR2_norm * (HD - cAHD * AH));
                          DEV_bb_bb[1] +=
                              dE_dxH * (VxR2_norm * (XA - cXAH * AH));
                        }
                      }
                    }
                  }
                }
                if (bonded)
                  break;
                pt = pt->fine_cell_next;
              }
              if (bonded)
                break;
            }
            if (bonded)
              break;
          }
          if (bonded)
            break;
        }
      } else {
        theHBA = theH->hbonded_acceptor;
        theA = theHBA->getA();
        E = theHBA->getE(theH);
        if (E == INF) {
          cerr << "fatal error in HBONDing pairs" << endl;
          exit(1);
        }
        if (theH->isSC() && theA->isSC())
          EHB_sc_sc += E;
        else if (!theA->isSC() && !theH->isSC())
          EHB_bb_bb += E;
        else
          EHB_bb_sc += E;
      }

    } // theH of ires
    { // Atom O of ires-1
      theA = gp.getResidue(ires - 1)->getO();
      theHBA = theA->acceptor;
      if (!theHBA) {
        cout << "error the accreptor is not right" << endl;
        exit(1);
      }
      theHBA->currentHBE(hbe_bb_bb, hbe_bb_sc, hbe_sc_sc);
      EHB_bb_bb += hbe_bb_bb;
      EHB_bb_sc += hbe_bb_sc;
      EHB_sc_sc += hbe_sc_sc;
      if (!theHBA->isSaturated()) { // not saturated
        // search the pissible un-hbonded protons
        theX = theHBA->getX();
        ix = static_cast<int>((theA->r.x - min.x) / cell_box.x);
        iy = static_cast<int>((theA->r.y - min.y) / cell_box.y);
        iz = static_cast<int>((theA->r.z - min.z) / cell_box.z);
        cindex = ix + (iy + iz * ncell[1]) * ncell[0];
        if (theA->fine_cell_index != cindex) { // precautions
          cout << "error in the accpetor cell index" << endl;
          exit(1);
        }
        // search teh neighboring cells to find the possible protons for
        // Hbond
        for (cz = iz - 1; cz <= iz + 1; cz++) {
          kz = cz;
          if (cz < 0)
            kz += ncell[2];
          else if (cz >= ncell[2])
            kz -= ncell[2];

          for (cy = iy - 1; cy <= iy + 1; cy++) {
            ky = cy;
            if (cy < 0)
              ky += ncell[1];
            else if (cy >= ncell[1])
              ky -= ncell[1];

            for (cx = ix - 1; cx <= ix + 1; cx++) {
              kx = cx;
              if (cx < 0)
                kx += ncell[0];
              else if (cx >= ncell[0])
                kx -= ncell[0];

              cindex = ncell[0] * (ncell[1] * kz + ky) + kx;
              pt = cell_protons[cindex];

              while (pt) {
                theH = pt;
                theD = theH->attachedHA; // get the Donar
                // precaution
                if (!theD) {
                  cout << "fatal error: Not a proton!" << endl;
                  exit(1);
                }

                if (!theH->isHBonded()) { // not hbonded
                  rHA = theH->r.getDist(theA->r);
                  vec XA = theA->r - theX->r;
                  XA.Normalize();
                  vec AH = theH->r - theA->r;
                  AH.Normalize();
                  vec HD = theD->r - theH->r;
                  HD.Normalize();
                  cXAH = XA * AH;
                  cAHD = AH * HD;
                  dSEQ = abs(theH->resID - theA->resID);
                  if (rHA < MAX_R && rHA > MIN_R && cXAH < MAX_CXAH &&
                      cXAH > MIN_CXAH && cAHD < MAX_CAHD && cAHD > MIN_CAHD &&
                      dSEQ > 0) {
                    if (theH->isSC()) { // bb-sc
                      isSC = 1;
                      if ((E = getHBE_DEV(dSEQ, isSC, theHBA->getHybrid(), rHA,
                                          cXAH, cAHD, dE_dr, dE_dxH, dE_dxD)) <
                          0) {
                        theH->hbonded_acceptor = theHBA;
                        theHBA->addH(theH, E);
                        EHB_bb_sc += E;
                        bonded = 1;

                        /*get the derivaste for all the
                         * acceptors*/
                        {
                          // theta13
                          vec R2 =
                              theA->r - gp.getResidue(ires + 1)->getCA()->r;
                          vec VxR2_norm = backrub_uv[0] ^ R2 / rHA;
                          DEV_bb_sc[0] -= dE_dr * (VxR2_norm * AH) * rHA;
                          DEV_bb_sc[0] -=
                              dE_dxD * (VxR2_norm * (HD - cAHD * AH));
                          DEV_bb_sc[0] +=
                              dE_dxH * ((backrub_uv[0] ^ XA) * AH -
                                        VxR2_norm * (XA - cXAH * AH));
                          // for theta 12
                          R2 = theA->r - gp.getResidue(ires)->getCA()->r;
                          VxR2_norm = backrub_uv[1] ^ R2 / rHA;
                          DEV_bb_sc[1] -= dE_dr * (VxR2_norm * AH) * rHA;
                          DEV_bb_sc[1] -=
                              dE_dxD * (VxR2_norm * (HD - cAHD * AH));
                          DEV_bb_sc[1] +=
                              dE_dxH * ((backrub_uv[1] ^ XA) * AH -
                                        VxR2_norm * (XA - cXAH * AH));
                        }
                      }
                    } else { // bb-bb
                      isSC = 0;
                      if ((E = getHBE_DEV(dSEQ, isSC, theHBA->getHybrid(), rHA,
                                          cXAH, cAHD, dE_dr, dE_dxH, dE_dxD)) <
                          0) {
                        theHBA->addH(theH, E);
                        theH->hbonded_acceptor = theHBA;
                        EHB_bb_bb += E;
                        bonded = 1;

                        { // theta13
                          vec R2 =
                              theA->r - gp.getResidue(ires + 1)->getCA()->r;
                          vec VxR2_norm = backrub_uv[0] ^ R2 / rHA;
                          DEV_bb_bb[0] -= dE_dr * (VxR2_norm * AH) * rHA;
                          DEV_bb_bb[0] -=
                              dE_dxD * (VxR2_norm * (HD - cAHD * AH));
                          DEV_bb_bb[0] +=
                              dE_dxH * ((backrub_uv[0] ^ XA) * AH -
                                        VxR2_norm * (XA - cXAH * AH));
                          // theta12
                          R2 = theA->r - gp.getResidue(ires)->getCA()->r;
                          VxR2_norm = backrub_uv[1] ^ R2 / rHA;
                          DEV_bb_bb[1] -= dE_dr * (VxR2_norm * AH) * rHA;
                          DEV_bb_bb[1] -=
                              dE_dxD * (VxR2_norm * (HD - cAHD * AH));
                          DEV_bb_bb[1] +=
                              dE_dxH * ((backrub_uv[1] ^ XA) * AH -
                                        VxR2_norm * (XA - cXAH * AH));
                        }
                      }
                    }
                  }
                }
                if (theHBA->isSaturated())
                  break;
                pt = pt->fine_cell_next;
              }
              if (theHBA->isSaturated())
                break;
            }
            if (theHBA->isSaturated())
              break;
          }
          if (theHBA->isSaturated())
            break;
        }
      }
    } // theO of ires-1
  }

  { // O of ires; H of ires+1
    if (theH = gp.getResidue(ires + 1)->getHN()) {
      if (!theH->isHBonded()) { // the proton does not have a HBond
                                // partner;
        bonded = 0;
        theD = theH->attachedHA; // get the Donar heavy atom

        // get the cooridiate of the FineGrid(cell)
        ix = static_cast<int>((theH->r.x - min.x) / cell_box.x);
        iy = static_cast<int>((theH->r.y - min.y) / cell_box.y);
        iz = static_cast<int>((theH->r.z - min.z) / cell_box.z);
        cindex = ix + (iy + iz * ncell[1]) * ncell[0];
        if (theH->fine_cell_index != cindex) { // precaucious
          cout << "fatal error in the cell index of protons" << endl;
          exit(1);
        }
        // search the possible acceptors and calculate the energy
        for (cz = iz - 1; cz <= iz + 1; cz++) {
          kz = cz;
          if (cz < 0)
            kz += ncell[2];
          else if (cz >= ncell[2])
            kz -= ncell[2];

          for (cy = iy - 1; cy <= iy + 1; cy++) {
            ky = cy;
            if (cy < 0)
              ky += ncell[1];
            else if (cy >= ncell[1])
              ky -= ncell[1];

            for (cx = ix - 1; cx <= ix + 1; cx++) {
              kx = cx;
              if (cx < 0)
                kx += ncell[0];
              else if (cx >= ncell[0])
                kx -= ncell[0];

              cindex = ncell[0] * (ncell[1] * kz + ky) + kx;
              pt = cell_acceptors[cindex];

              while (pt) {
                // calculate the hbond energy
                theHBA = pt->acceptor;
                // precaucious of the accpetor
                if (!theHBA) {
                  cout << "fatal error: Not a HBacceptor!" << endl;
                  exit(1);
                }

                if (!theHBA->isSaturated()) { // not saturated

                  theA = theHBA->getA();
                  theX = theHBA->getX();
                  rHA = theH->r.getDist(theA->r);
                  vec XA = theA->r - theX->r;
                  XA.Normalize();
                  vec AH = theH->r - theA->r;
                  AH.Normalize();
                  vec HD = theD->r - theH->r;
                  HD.Normalize();
                  cXAH = XA * AH;
                  cAHD = AH * HD;
                  dSEQ = abs(theH->resID - theA->resID);
                  if (rHA < MAX_R && rHA > MIN_R && cXAH < MAX_CXAH &&
                      cXAH > MIN_CXAH && cAHD < MAX_CAHD && cAHD > MIN_CAHD &&
                      dSEQ > 0) {
                    if (theA->isSC()) { // bb-sc
                      isSC = 1;
                      if ((E = getHBE_DEV(dSEQ, isSC, theHBA->getHybrid(), rHA,
                                          cXAH, cAHD, dE_dr, dE_dxH, dE_dxD)) <
                          0) {
                        theH->hbonded_acceptor = theHBA;
                        theHBA->addH(theH, E);
                        EHB_bb_sc += E;
                        bonded = 1;

                        {
                          // for theta13
                          vec R2 =
                              theH->r - gp.getResidue(ires + 1)->getCA()->r;
                          vec VxR2_norm = backrub_uv[0] ^ R2 / rHA;
                          DEV_bb_sc[0] += dE_dr * (VxR2_norm * AH) * rHA;
                          DEV_bb_sc[0] +=
                              dE_dxD * ((backrub_uv[0] ^ HD) * AH +
                                        VxR2_norm * (HD - cAHD * AH));
                          DEV_bb_sc[0] +=
                              dE_dxH * (VxR2_norm * (XA - cXAH * AH));
                          // for theta12
                          R2 = theH->r - gp.getResidue(ires)->getCA()->r;
                          VxR2_norm = backrub_uv[2] ^ R2 / rHA;
                          DEV_bb_sc[2] += dE_dr * (VxR2_norm * AH) * rHA;
                          DEV_bb_sc[2] +=
                              dE_dxD * ((backrub_uv[2] ^ HD) * AH +
                                        VxR2_norm * (HD - cAHD * AH));
                          DEV_bb_sc[2] +=
                              dE_dxH * (VxR2_norm * (XA - cXAH * AH));
                        }
                      }
                    } else { // bb-bb
                      isSC = 0;
                      if ((E = getHBE_DEV(dSEQ, isSC, theHBA->getHybrid(), rHA,
                                          cXAH, cAHD, dE_dr, dE_dxH, dE_dxD)) <
                          0) {
                        theHBA->addH(theH, E);
                        theH->hbonded_acceptor = theHBA;
                        EHB_bb_bb += E;
                        bonded = 1;

                        {
                          // for theta13
                          vec R2 =
                              theH->r - gp.getResidue(ires + 1)->getCA()->r;
                          vec VxR2_norm = backrub_uv[0] ^ R2 / rHA;
                          DEV_bb_bb[0] += dE_dr * (VxR2_norm * AH) * rHA;
                          DEV_bb_bb[0] +=
                              dE_dxD * ((backrub_uv[0] ^ HD) * AH +
                                        VxR2_norm * (HD - cAHD * AH));
                          DEV_bb_bb[0] +=
                              dE_dxH * (VxR2_norm * (XA - cXAH * AH));
                          // for theta12
                          R2 = theH->r - gp.getResidue(ires)->getCA()->r;
                          VxR2_norm = backrub_uv[2] ^ R2 / rHA;
                          DEV_bb_bb[2] += dE_dr * (VxR2_norm * AH) * rHA;
                          DEV_bb_bb[2] +=
                              dE_dxD * ((backrub_uv[2] ^ HD) * AH +
                                        VxR2_norm * (HD - cAHD * AH));
                          DEV_bb_bb[2] +=
                              dE_dxH * (VxR2_norm * (XA - cXAH * AH));
                        }
                      }
                    }
                  }
                }
                if (bonded)
                  break;
                pt = pt->fine_cell_next;
              }
              if (bonded)
                break;
            }
            if (bonded)
              break;
          }
          if (bonded)
            break;
        }
      } else {
        theHBA = theH->hbonded_acceptor;
        theA = theHBA->getA();
        E = theHBA->getE(theH);
        if (E == INF) {
          cerr << "fatal error in HBONDing pairs" << endl;
          exit(1);
        }
        if (theH->isSC() && theA->isSC())
          EHB_sc_sc += E;
        else if (!theA->isSC() && !theH->isSC())
          EHB_bb_bb += E;
        else
          EHB_bb_sc += E;
      }

    } // theH of ires
    {
      theA = gp.getResidue(ires)->getO();
      theHBA = theA->acceptor;
      if (!theHBA) {
        cout << "error the accreptor is not right" << endl;
        exit(1);
      }
      theHBA->currentHBE(hbe_bb_bb, hbe_bb_sc, hbe_sc_sc);
      EHB_bb_bb += hbe_bb_bb;
      EHB_bb_sc += hbe_bb_sc;
      EHB_sc_sc += hbe_sc_sc;
      if (!theHBA->isSaturated()) { // not saturated
        // search the pissible un-hbonded protons
        theX = theHBA->getX();
        ix = static_cast<int>((theA->r.x - min.x) / cell_box.x);
        iy = static_cast<int>((theA->r.y - min.y) / cell_box.y);
        iz = static_cast<int>((theA->r.z - min.z) / cell_box.z);
        cindex = ix + (iy + iz * ncell[1]) * ncell[0];
        if (theA->fine_cell_index != cindex) { // precautions
          cout << "error in the accpetor cell index" << endl;
          exit(1);
        }
        // search teh neighboring cells to find the possible protons for
        // Hbond
        for (cz = iz - 1; cz <= iz + 1; cz++) {
          kz = cz;
          if (cz < 0)
            kz += ncell[2];
          else if (cz >= ncell[2])
            kz -= ncell[2];

          for (cy = iy - 1; cy <= iy + 1; cy++) {
            ky = cy;
            if (cy < 0)
              ky += ncell[1];
            else if (cy >= ncell[1])
              ky -= ncell[1];

            for (cx = ix - 1; cx <= ix + 1; cx++) {
              kx = cx;
              if (cx < 0)
                kx += ncell[0];
              else if (cx >= ncell[0])
                kx -= ncell[0];

              cindex = ncell[0] * (ncell[1] * kz + ky) + kx;
              pt = cell_protons[cindex];

              while (pt) {
                theH = pt;
                theD = theH->attachedHA; // get the Donar
                // precaution
                if (!theD) {
                  cout << "fatal error: Not a proton!" << endl;
                  exit(1);
                }

                if (!theH->isHBonded()) { // not hbonded
                  rHA = theH->r.getDist(theA->r);
                  vec XA = theA->r - theX->r;
                  XA.Normalize();
                  vec AH = theH->r - theA->r;
                  AH.Normalize();
                  vec HD = theD->r - theH->r;
                  HD.Normalize();
                  cXAH = XA * AH;
                  cAHD = AH * HD;
                  dSEQ = abs(theH->resID - theA->resID);
                  if (rHA < MAX_R && rHA > MIN_R && cXAH < MAX_CXAH &&
                      cXAH > MIN_CXAH && cAHD < MAX_CAHD && cAHD > MIN_CAHD &&
                      dSEQ > 0) {
                    if (theH->isSC()) { // bb-sc
                      isSC = 1;
                      if ((E = getHBE_DEV(dSEQ, isSC, theHBA->getHybrid(), rHA,
                                          cXAH, cAHD, dE_dr, dE_dxH, dE_dxD)) <
                          0) {
                        theH->hbonded_acceptor = theHBA;
                        theHBA->addH(theH, E);
                        EHB_bb_sc += E;
                        bonded = 1;

                        /*get the derivaste for all the
                         * acceptors*/
                        {
                          // theta13
                          vec R2 =
                              theA->r - gp.getResidue(ires + 1)->getCA()->r;
                          vec VxR2_norm = backrub_uv[0] ^ R2 / rHA;
                          DEV_bb_sc[0] -= dE_dr * (VxR2_norm * AH) * rHA;
                          DEV_bb_sc[0] -=
                              dE_dxD * (VxR2_norm * (HD - cAHD * AH));
                          DEV_bb_sc[0] +=
                              dE_dxH * ((backrub_uv[0] ^ XA) * AH -
                                        VxR2_norm * (XA - cXAH * AH));
                          // for theta 12
                          R2 = theA->r - gp.getResidue(ires)->getCA()->r;
                          VxR2_norm = backrub_uv[2] ^ R2 / rHA;
                          DEV_bb_sc[2] -= dE_dr * (VxR2_norm * AH) * rHA;
                          DEV_bb_sc[2] -=
                              dE_dxD * (VxR2_norm * (HD - cAHD * AH));
                          DEV_bb_sc[2] +=
                              dE_dxH * ((backrub_uv[2] ^ XA) * AH -
                                        VxR2_norm * (XA - cXAH * AH));
                        }
                      }
                    } else { // bb-bb
                      isSC = 0;
                      if ((E = getHBE_DEV(dSEQ, isSC, theHBA->getHybrid(), rHA,
                                          cXAH, cAHD, dE_dr, dE_dxH, dE_dxD)) <
                          0) {
                        theHBA->addH(theH, E);
                        theH->hbonded_acceptor = theHBA;
                        EHB_bb_bb += E;
                        bonded = 1;

                        { // theta13
                          vec R2 =
                              theA->r - gp.getResidue(ires + 1)->getCA()->r;
                          vec VxR2_norm = backrub_uv[0] ^ R2 / rHA;
                          DEV_bb_bb[0] -= dE_dr * (VxR2_norm * AH) * rHA;
                          DEV_bb_bb[0] -=
                              dE_dxD * (VxR2_norm * (HD - cAHD * AH));
                          DEV_bb_bb[0] +=
                              dE_dxH * ((backrub_uv[0] ^ XA) * AH -
                                        VxR2_norm * (XA - cXAH * AH));
                          // theta12
                          R2 = theA->r - gp.getResidue(ires)->getCA()->r;
                          VxR2_norm = backrub_uv[2] ^ R2 / rHA;
                          DEV_bb_bb[2] -= dE_dr * (VxR2_norm * AH) * rHA;
                          DEV_bb_bb[2] -=
                              dE_dxD * (VxR2_norm * (HD - cAHD * AH));
                          DEV_bb_bb[2] +=
                              dE_dxH * ((backrub_uv[2] ^ XA) * AH -
                                        VxR2_norm * (XA - cXAH * AH));
                        }
                      }
                    }
                  }
                }
                if (theHBA->isSaturated())
                  break;
                pt = pt->fine_cell_next;
              }
              if (theHBA->isSaturated())
                break;
            }
            if (theHBA->isSaturated())
              break;
          }
          if (theHBA->isSaturated())
            break;
        }
      }
    } // theO of ires-1
  }
}

int topology::fftid(const string &fftype) const {
  for (int i = 0; i < ntype; i++) {
    if (mass_table[i].FFname == fftype) {
      return mass_table[i].FFTindex;
    }
  }
  return -1;
}

void topology::addTopology(const string &name) {
  strcpy(file, name.c_str());
  ifstream in(file, ios::in);
  if (!in) {
    cerr << "Failed to open topology file " << file << endl;
    exit(1);
  }
  string line;
  char buf[1000];
  char ch_tmp[1000];
  double f_dummy;
  int i_dummy;
  vector<string> pool;
  vector<string> bond_pool;
  vector<string> dihedral_pool;
  vector<string> improper_pool;
  // read untill the MASS records
  while (getline(in, line)) {
    if (line.length() > 4) {
      if (line.substr(0, 4) == "MASS") {
        pool.push_back(line);
        break;
      }
    }
  }
  // read the mass records
  do {
    if (in.eof() || !getline(in, line)) {
      break;
    }
    if (line.substr(0, 4) == "MASS")
      pool.push_back(line);
  } while (line.substr(0, 3) != "END");
  // loop all the new mass_types, and only add those not in the current table
  for (int i = 0; i < pool.size(); i++) {
    i_dummy = pool[i].length() - 4;
    pool[i].substr(4, i_dummy).copy(buf, i_dummy);
    buf[i_dummy] = '\0';
    sscanf(buf, "%s %lf", ch_tmp, &f_dummy);
    bool exist = false;
    for (int j = 0; j < mass_table.size(); j++) {
      if (mass_table[j].FFname == ch_tmp) {
        exist = true;
        if (mass_table[j].atm_mass != f_dummy) {
          cerr << "pre-existing type:" << ch_tmp << " has different mass"
               << endl;
          exit(1);
        }
        break;
      }
    }
    if (!exist) {
      mass_table.push_back(TopAtomType());
      mass_table.back().init(ch_tmp, f_dummy);
      mass_table.back().FFTindex = int(mass_table.size()) - 1;
    }
  }
  ntype = mass_table.size();
  // rewind
  // cout << line << " " << in.eof() << endl;
  if (in.eof()) {
    in.close();
    in.open(file, ios::in);
  }

  int start = 0;
  int pt = 0;
  char resiname[100];
  char pdbtype[100];
  char fftype[100];
  char atom_i[100], atom_j[100], atom_k[100], atom_l[100];
  double q, x, y, z;
  int fft;
  int is_rna;
  // construct the residue topology
  while (getline(in, line)) {
    line.copy(buf, line.length());
    buf[line.length()] = '\0';
    if (!start) { // start a new records?
      sscanf(buf, "%s %s", ch_tmp, resiname);
      bool is_rna = jnc::string_starts_with(ch_tmp, "RNA");
      if (jnc::string_istarts_with(ch_tmp, "RESI") || is_rna) {
        start = 1;
        pool.clear();
        bond_pool.clear();
        dihedral_pool.clear();
        improper_pool.clear();
        // cout << resiname << endl;
      }
    } else {
      sscanf(buf, "%s", ch_tmp);
      if (!strncasecmp(ch_tmp, "END", 3)) {
        // search whether it was there already
        bool exist = false;
        for (int i = 0; i < resi_table.size(); i++) {
          if (jnc::string_iequals(resi_table[i]->PDBname, resiname)) {
            // if (!strcasecmp(resi_table[i]->getName(), resiname)) {
            cerr << "Residue " << resiname << " is defined previously" << endl;
            exist = true;
            break;
          }
        }
        if (!exist) {
          start = 0;
          // start to construct the read residue infor
          TopResidue *res = new TopResidue(resiname, pool.size());
          // cout << "RESIDUE NAME:"<<resiname << endl;
          for (int i = 0; i < pool.size(); i++) {
            int nchar = pool[i].length();
            pool[i].copy(buf, nchar);
            buf[nchar] = '\0';
            for (int j = 0; j < nchar; j++) { // substitute '=' as ' '
              if (buf[j] == '=')
                buf[j] = ' ';
            }
            if (is_rna) {
              sscanf(buf, "%s %s %s %s %s %lf %s (%lf %lf %lf)", ch_tmp,
                     pdbtype, ch_tmp, fftype, ch_tmp, &q, ch_tmp, &x, &y, &z);
            } else {
              sscanf(buf, "%s %s %s %s %s %lf", ch_tmp, pdbtype, ch_tmp, fftype,
                     ch_tmp, &q);
            }
            // initialize the residue topology
            if ((fft = fftid(fftype)) != -1) {
              res->getAtom(i)->init(pdbtype, fftid(fftype), q);
              if (is_rna) {
                res->getAtom(i)->ic = {x, y, z};
                // res->getAtom(i)->setIC(x, y, z);
                // cout << resiname << " " << pdbtype << " " << fftid(fftype) <<
                // " " << x << " " << y << " " << z << endl;
              }
              // cout << pdbtype << " " << fftype << " " << q << endl;
            } else {
              cerr << "topology file:" << file << " has a mistake" << endl;
              cerr << fftype << endl;
              exit(2);
            }
          }
          // construct the bond_list
          res->nbond = bond_pool.size();
          res->bond = new int *[res->nbond];
          res->bond[0] = new int[2 * res->nbond];
          for (int i = 1; i < res->nbond; i++)
            res->bond[i] = res->bond[i - 1] + 2;
          for (int i = 0; i < bond_pool.size(); i++) {
            int nchar = bond_pool[i].length();
            bond_pool[i].copy(buf, nchar);
            buf[nchar] = '\0';
            sscanf(buf, "%s %s %s", ch_tmp, atom_i, atom_j);
            res->bond[i][0] = res->getIndex(atom_i);
            res->bond[i][1] = res->getIndex(atom_j);
            // cout<< res->bond[i][0] << " " << res->bond[i][1] << endl;
          }
          // construct the dihedral_list
          res->ndihedral = dihedral_pool.size();
          if (res->ndihedral != 0) {
            res->dihedral = new int *[res->ndihedral];
            res->dihedral[0] = new int[4 * res->ndihedral];
            for (int i = 1; i < res->ndihedral; i++)
              res->dihedral[i] = res->dihedral[i - 1] + 4;
            for (int i = 0; i < dihedral_pool.size(); i++) {
              int nchar = dihedral_pool[i].length();
              dihedral_pool[i].copy(buf, nchar);
              buf[nchar] = '\0';
              sscanf(buf, "%s %s %s %s %s", ch_tmp, atom_i, atom_j, atom_k,
                     atom_l);
              res->dihedral[i][0] = res->getIndex(atom_i);
              res->dihedral[i][1] = res->getIndex(atom_j);
              res->dihedral[i][2] = res->getIndex(atom_k);
              res->dihedral[i][3] = res->getIndex(atom_l);
              // cout << res->dihedral[i][0] << " "
              //     << res->dihedral[i][1] << " "
              //     << res->dihedral[i][2] << " "
              //     << res->dihedral[i][3] << endl;
            }
          }
          // construct the impromer dihedral_list
          res->nimproper = improper_pool.size();
          if (res->nimproper != 0) {
            res->improper = new int *[res->nimproper];
            res->improper[0] = new int[4 * res->nimproper];
            for (int i = 1; i < res->nimproper; i++)
              res->improper[i] = res->improper[i - 1] + 4;
            for (int i = 0; i < improper_pool.size(); i++) {
              int nchar = improper_pool[i].length();
              improper_pool[i].copy(buf, nchar);
              buf[nchar] = '\0';
              sscanf(buf, "%s %s %s %s %s", ch_tmp, atom_i, atom_j, atom_k,
                     atom_l);
              res->improper[i][0] = res->getIndex(atom_i);
              res->improper[i][1] = res->getIndex(atom_j);
              res->improper[i][2] = res->getIndex(atom_k);
              res->improper[i][3] = res->getIndex(atom_l);
              // cout << res->improper[i][0] << " "
              //     << res->improper[i][1] << " "
              //     << res->improper[i][2] << " "
              //     << res->improper[i][3] << endl;
            }
          }
          resi_table.push_back(res);
        }
      } else if (!strncasecmp(ch_tmp, "ATOM", 4)) { // atom records
        pool.push_back(line);
      } else if (!strncasecmp(ch_tmp, "BOND", 4)) { // bond infor
        // to upper
        line = jnc::string_upper_c(line);
        int start = line.find("BOND");
        int next;
        if (start != -1) {
          while ((next = line.substr(start + 4).find("BOND")) != -1) {
            bond_pool.push_back(line.substr(start, next + 4));
            // cout << line.substr(start, next+4) << endl;
            start += next + 4;
          }
          bond_pool.push_back(line.substr(start));
          // cout << line.substr(start) << endl;
        } else {
          cout << "fatal error: BOND" << endl;
          exit(1);
        }
      } else if (!strncasecmp(ch_tmp, "DIHE", 4)) { // dihedral
        line = jnc::string_upper_c(line);
        int start = line.find("DIHE");
        int next;
        if (start != -1) {
          while ((next = line.substr(start + 4).find("DIHE")) != -1) {
            dihedral_pool.push_back(line.substr(start, next + 4));
            // cout << line.substr(start, next+4) << endl;
            start += next + 4;
          }
          dihedral_pool.push_back(line.substr(start));
          // cout << line.substr(start) << endl;
        } else {
          cout << "fatal error: Dihedral" << endl;
          exit(1);
        }
      } else if (!strncasecmp(ch_tmp, "IMPR", 4)) { // improper dihedral
        line = jnc::string_upper_c(line);
        int start = line.find("IMPR");
        int next;
        if (start != -1) {
          while ((next = line.substr(start + 4).find("IMPR")) != -1) {
            improper_pool.push_back(line.substr(start, next + 4));
            // cout << line.substr(start, next+4) << endl;
            start += next + 4;
          }
          improper_pool.push_back(line.substr(start));
          // cout << line.substr(start) << endl;
        } else {
          cout << "fatal error: improper" << endl;
          exit(1);
        }
      }
    }
    // clear the tmperary char array
    line[0] = '\0';
    ch_tmp[0] = '\0';
  }
  in.close();
}

int topology::getAtomInfor(const string &res_name, const string &atm_name,
                           int &type, double &mass, double &charge) const {
  int res_index = -1;
  for (int i = 0; i < resi_table.size(); i++) {
    if (resi_table[i]->PDBname == res_name) {
      // if (!strcmp(resi_table[i]->PDBname, res_name)) {  // correct resname
      res_index = i;
      break;
    }
  }
  if (res_index == -1) {
    cerr << "No match of residue name" << res_name << " : " << atm_name << endl;
    return -1;
  }
  for (int j = 0; j < resi_table[res_index]->size(); j++) {
    if (resi_table[res_index]->getAtom(j)->PDBname ==
        atm_name) { // corr atm name
      type = resi_table[res_index]->getAtom(j)->FFtype;
      mass = mass_table[type].atm_mass;
      charge = resi_table[res_index]->getAtom(j)->charge;
      return 1;
    }
  }
  string aname;
  aname = atm_name;
  int len = (int)aname.size();
  if (aname[len - 1] >= '0' &&
      aname[len - 1] <= '9') { // delete the last number digit
    aname.pop_back();
    // aname[len - 1] = '\0';
    int nmatch = 0;
    int atom_index = -1;
    for (int j = 0; j < resi_table[res_index]->size(); j++) {
      //            cerr << resi_table[res_index]->getAtom(j)->getName() <<
      //            ' ' << aname << endl;
      if (resi_table[res_index]->getAtom(j)->PDBname ==
          aname) { // corr atm name
        nmatch++;
        atom_index = j;
      }
    }
    if (nmatch == 0) {
      cerr << "warning: No match of atom '" << atm_name << "' in '" << res_name
           << "'" << endl;
      return 0;
    } else if (nmatch > 1) {
      cerr << "Multiple match" << endl;
      return 0;
    } else {
      type = resi_table[res_index]->getAtom(atom_index)->FFtype;
      mass = mass_table[type].atm_mass;
      charge = resi_table[res_index]->getAtom(atom_index)->charge;
      return 1;
    }
  } else
    return 0;
}

BBDep_AA::BBDep_AA() {
  fy_aa = new double **[FY_NBIN + 1];
  fy_aa[0] = new double *[(FY_NBIN + 1) * (FY_NBIN + 1)];
  for (int i = 1; i < FY_NBIN + 1; i++)
    fy_aa[i] = fy_aa[i - 1] + FY_NBIN + 1;
  double *pt = new double[(FY_NBIN + 1) * (FY_NBIN + 1) * 20];
  for (int i = 0; i < (FY_NBIN + 1) * (FY_NBIN + 1) * 20; i++)
    pt[i] = 0;
  for (int i = 0; i < FY_NBIN + 1; i++)
    for (int j = 0; j < FY_NBIN + 1; j++)
      fy_aa[i][j] = pt + ((FY_NBIN + 1) * i + j) * 20;

  e_stat = new double **[FY_NBIN + 1];
  e_stat[0] = new double *[(FY_NBIN + 1) * (FY_NBIN + 1)];
  for (int i = 1; i < FY_NBIN + 1; i++)
    e_stat[i] = e_stat[i - 1] + FY_NBIN + 1;
  double *f_pt = new double[(FY_NBIN + 1) * (FY_NBIN + 1) * 20];
  for (int i = 0; i < (FY_NBIN + 1) * (FY_NBIN + 1) * 20; i++)
    f_pt[i] = INF;
  for (int i = 0; i < FY_NBIN + 1; i++)
    for (int j = 0; j < FY_NBIN + 1; j++)
      e_stat[i][j] = f_pt + ((FY_NBIN + 1) * i + j) * 20;
}

BBDep_AA::~BBDep_AA() {
  delete[] fy_aa[0][0];
  delete[] fy_aa[0];
  delete[] fy_aa;

  delete[] e_stat[0][0];
  delete[] e_stat[0];
  delete[] e_stat;
}

/*****************************************************
  To take into consider the N- or C-terminal residues,
  we use the summed phi or psi distribution
  [FY_NBIN][ipsi][ia]---N-terminal
  [iphi][FY_NBIN][ia]---C-terminal
******************************************************/
void BBDep_AA::updateCAP() {
  /*N-terminal*/
  for (int iY = 0; iY < FY_NBIN; iY++) {
    for (int ia = 0; ia < 20; ia++) {
      double n_psi_ia = 0;
      for (int iF = 0; iF < FY_NBIN; iF++)
        n_psi_ia += fy_aa[iF][iY][ia];
      fy_aa[FY_NBIN][iY][ia] = n_psi_ia;
    }
  }
  /*C-terminal*/
  for (int iF = 0; iF < FY_NBIN; iF++) {
    for (int ia = 0; ia < 20; ia++) {
      double n_phi_ia = 0;
      for (int iY = 0; iY < FY_NBIN; iY++)
        n_phi_ia += fy_aa[iF][iY][ia];
      fy_aa[iF][FY_NBIN][ia] = n_phi_ia;
    }
  }
}

void BBDep_AA::readFile(const std::string &file) {
  ifstream in(file.c_str(), ios::in);
  char buf[100];
  char name[10];
  double nrecords;
  int id;
  double phi, psi;
  int iF, iY;
  while (in.getline(buf, 100)) {
    sscanf(buf, "%lf %lf %s %lf", &phi, &psi, name, &nrecords);
    id = residueid(name);
    iF = static_cast<int>((phi - FY_MIN) / FY_BIN);
    iY = static_cast<int>((psi - FY_MIN) / FY_BIN);
    // cout << phi << " " << psi  << " " << name << " " << nrecords << endl;
    fy_aa[iF][iY][id] = nrecords;
  }
  updateCAP();
  in.close();
}

void BBDep_AA::updateStatisticalEnergy() {
  for (int iF = 0; iF < FY_NBIN; iF++) {
    for (int iY = 0; iY < FY_NBIN; iY++) {
      double n_total = 0;
      for (int ia = 0; ia < 20; ia++) {
        n_total += fy_aa[iF][iY][ia];
      }
      if (n_total > 0) {
        for (int ia = 0; ia < 20; ia++) {
          double p = static_cast<double>(fy_aa[iF][iY][ia]) /
                     static_cast<double>(n_total);
          if (p > 0) {
            e_stat[iF][iY][ia] = -log(p);
            // cout << iF << " " << iY << " " << ia << " " << fy_aa[iF][iY][ia]
            // << " " << e_stat[iF][iY][ia] << endl;
          }
        }
      }
    }
  }

  // N-terminal residues
  for (int iY = 0; iY < FY_NBIN; iY++) {
    double n_psi = 0;
    for (int ia = 0; ia < 20; ia++)
      n_psi += fy_aa[FY_NBIN][iY][ia];
    if (n_psi > 0)
      for (int ia = 0; ia < 20; ia++) {
        double p = static_cast<double>(fy_aa[FY_NBIN][iY][ia]) / n_psi;
        if (p > 0)
          e_stat[FY_NBIN][iY][ia] = -log(p);
      }
  }
  // C-terminal residues
  for (int iF = 0; iF < FY_NBIN; iF++) {
    double n_phi = 0;
    for (int ia = 0; ia < 20; ia++)
      n_phi += fy_aa[iF][FY_NBIN][ia];
    if (n_phi > 0)
      for (int ia = 0; ia < 20; ia++) {
        double p = static_cast<double>(fy_aa[iF][FY_NBIN][ia]) / n_phi;
        if (p > 0)
          e_stat[iF][FY_NBIN][ia] = -log(p);
      }
  }
}

BBDep_RotLib::BBDep_RotLib() {
  records = new BBDep_resRotLib **[FY_NBIN + 1];
  records[0] = new BBDep_resRotLib *[(FY_NBIN + 1) * (FY_NBIN + 1)];
  for (int i = 1; i < FY_NBIN + 1; i++)
    records[i] = records[i - 1] + FY_NBIN + 1;
  for (int i = 0; i < (FY_NBIN + 1) * (FY_NBIN + 1); i++) {
    records[0][i] = new BBDep_resRotLib[20];
  }
}

BBDep_RotLib::~BBDep_RotLib() {
  for (int i = 0; i < (FY_NBIN + 1) * (FY_NBIN + 1); i++)
    delete[] records[0][i];
  delete[] records[0];
  delete[] records;
}

void BBDep_RotLib::readFile(const std::string &rotFile) {
  ifstream in(rotFile.c_str(), ios::in);

  vector<string> rot_pool;
  string line;
  char buf[1000];

  char name[10];
  int id = -1;

  int phi, psi;
  int prev_phi = -180, prev_psi = -180;
  while (getline(in, line)) {
    line.copy(buf, line.length());
    buf[line.length()] = '\0';
    sscanf(buf, "%s %d %d", name, &phi, &psi);
    if (phi == prev_phi && psi == prev_psi) {
      id = residueid(name);
      rot_pool.push_back(line);
    } else {
      set_record(prev_phi, prev_psi, id, rot_pool);
      prev_phi = phi;
      prev_psi = psi;
      rot_pool.clear();
      id = residueid(name);
      rot_pool.push_back(line);
    }
  }
  set_record(phi, psi, id, rot_pool);
  updateCAP();
  in.close();

  /*check the ordering*/
  for (int iF = 0; iF < FY_NBIN; iF++) {
    for (int iY = 0; iY < FY_NBIN; iY++) {
      for (int iaa = 0; iaa < 20; iaa++) {
        if (iaa != ALA && iaa != GLY)
          records[iF][iY][iaa].ordering();
      }
    }
  }
  for (int iFY = 0; iFY < FY_NBIN; iFY++) {
    for (int iaa = 0; iaa < 20; iaa++) {
      if (iaa != ALA && iaa != GLY) {
        records[FY_NBIN][iFY][iaa].ordering(); // N-terminal
        records[iFY][FY_NBIN][iaa].ordering(); // C-terminal
      }
    }
  }
}

void BBDep_RotLib::set_record(int phi, int psi, int id, vector<string> &pool) {
  int iphi = static_cast<int>((phi - FY_MIN) / FY_BIN);
  int ipsi = static_cast<int>((psi - FY_MIN) / FY_BIN);
  iphi %= FY_NBIN;
  ipsi %= FY_NBIN;

  char buf[1000];
  char name[10];
  int f, y;
  int nrecords;
  int prev_nrecords;
  int ichi1, ichi2, ichi3, ichi4;
  double p;
  double achi1, achi2, achi3, achi4;
  double schi1, schi2, schi3, schi4;

  BBDep_rotRecord *rotlib = records[iphi][ipsi][id].getRotamers();
  int nrot = records[iphi][ipsi][id].getNRotamers();
  int nrec = records[iphi][ipsi][id].getNRecords();

  // number of records
  pool[0].copy(buf, pool[0].length());
  buf[pool[0].length()] = '\0';
  sscanf(buf, "%s %d %d %d", name, &f, &y, &nrecords);
  prev_nrecords = nrecords;
  if (rotlib) { // already assigned
    if (nrec != nrecords) {
      cout << "fatal error!" << nrec << " " << nrecords << endl;
      exit(2);
    } else {
      return;
    }
  }

  records[iphi][ipsi][id].init(id, pool.size());
  rotlib = records[iphi][ipsi][id].getRotamers();
  records[iphi][ipsi][id].setNRecords(nrecords);
  for (int i = 0; i < pool.size(); i++) {
    pool[i].copy(buf, pool[i].length());
    buf[pool[i].length()] = '\0';
    sscanf(buf, "%s %d %d %d %d %d %d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf",
           name, &f, &y, &nrecords, &ichi1, &ichi2, &ichi3, &ichi4, &p, &achi1,
           &achi2, &achi3, &achi4, &schi1, &schi2, &schi3, &schi4);
    rotlib[i].setIChi(ichi1, ichi2, ichi3, ichi4);
    rotlib[i].setChi(achi1, schi1, achi2, schi2, achi3, schi3, achi4, schi4);
    rotlib[i].setP(p);
    if (nrecords != prev_nrecords) {
      cout << "fatal error in the bbdep rotamer library" << endl;
      exit(1);
    }
  }
}

/*in principle the term p(a|phi,psi) can be determined from a larger PDB set
  Right now, we only use these outlined from the Dunbrac library*/
void BBDep_RotLib::updateStatisticalEnergy() {
  // int N_aa[20];
  // for(int i=0; i<20; i++){
  //   N_aa[i]=0;
  //   for(int iF=0; iF<FY_NBIN; iF++){
  //     for(int iY=0; iY<FY_NBIN; iY++){
  //	N_aa[i]+=records[iF][iY][i].getNRecords();
  //     }
  //   }
  //   cout << residuename(i) << " " << N_aa[i] << endl;
  // }

  for (int iF = 0; iF < FY_NBIN; iF++) {
    for (int iY = 0; iY < FY_NBIN; iY++) {
      for (int iaa = 0; iaa < 20; iaa++) {
        // double p = static_cast<double>(records[iF][iY][iaa].getNRecords());
        // p /= static_cast<double>(N_aa[iaa]);
        records[iF][iY][iaa].setE_STAT();
      }
    }
  }
  // N- and C- terminals
  for (int iFY = 0; iFY < FY_NBIN; iFY++) {
    for (int iaa = 0; iaa < 20; iaa++) {
      if (iaa != ALA && iaa != GLY) {
        records[FY_NBIN][iFY][iaa].setE_STAT(); // N-terminal
        records[iFY][FY_NBIN][iaa].setE_STAT(); // C-terminal
      }
    }
  }
}

/*In the case of capping, the averaged distribution is assumed
  For the N-terminal residues (x20), phi is not defined but psi is defined.
  for each psi-bins, we average all the phis
  The same operation for C-terminal capping, where is the average is over psi.*/
void BBDep_RotLib::updateCAP() {
  /*N-teminal*/
  for (int ia = 0; ia < 20; ia++) {
    if (ia != GLY && ia != ALA) {
      int nRot = records[0][0][ia].getNRotamers();
      for (int iY = 0; iY < FY_NBIN; iY++) {
        records[FY_NBIN][iY][ia].init(ia, nRot);
        for (int irot = 0; irot < nRot; irot++) {
          double N_psi_chi = 0;
          int N_psi = 0;
          BBDep_rotRecord *rpt = records[FY_NBIN][iY][ia].getRotamer(irot);
          double chi[4];
          double prev_chi[4];
          double sig[4];
          for (int k = 0; k < 4; k++) {
            prev_chi[k] = -3600;
            chi[k] = sig[k] = 0;
          }
          const double *aveChi = NULL;
          const double *sigChi = NULL;
          const unsigned char *iChi = NULL;
          double p = 0;
          for (int iF = 0; iF < FY_NBIN; iF++) { // average over phi
            aveChi = records[iF][iY][ia].getRotamer(irot)->getAveChi();
            sigChi = records[iF][iY][ia].getRotamer(irot)->getSigChi();
            iChi = records[iF][iY][ia].getRotamer(irot)->getIChi();
            p = records[iF][iY][ia].getRotamer(irot)->getP();
            int n_phi_psi = records[iF][iY][ia].getNRecords();
            N_psi += n_phi_psi;
            N_psi_chi += n_phi_psi * p;
            for (int k = 0; k < 4; k++) {
              double angle = aveChi[k];
              if (prev_chi[k] < -360)
                prev_chi[k] = aveChi[k]; // initialization
              if (aveChi[k] - prev_chi[k] > 180)
                angle -= 360;
              if (aveChi[k] - prev_chi[k] < -180)
                angle += 360;

              chi[k] += n_phi_psi * p * angle;
              sig[k] += n_phi_psi * p * sigChi[k];

              if (N_psi > 0)
                prev_chi[k] = chi[k] / N_psi;
              else
                prev_chi[k] = chi[k];
            }
          }

          double p_psi_chi = p;
          if (N_psi > 0 && N_psi_chi > 0) {
            p_psi_chi = N_psi_chi / N_psi;
            for (int k = 0; k < 4; k++) {
              chi[k] /= N_psi_chi;
              sig[k] /= N_psi_chi;

              if (chi[k] > 180)
                chi[k] -= 360;
              if (chi[k] < -180)
                chi[k] += 360;
            }
          } else {
            for (int k = 0; k < 4; k++) {
              chi[k] = aveChi[k];
              sig[k] = sigChi[k];
            }
          }
          rpt->setP(p_psi_chi);
          rpt->setIChi(iChi[0], iChi[1], iChi[2], iChi[3]);
          rpt->setChi(chi[0], sig[0], chi[1], sig[1], chi[2], sig[2], chi[3],
                      sig[3]);
          // cout << static_cast<int>(iChi[0]) << " " << chi[0] << " " << sig[0]
          // << " "
          //      << static_cast<int>(iChi[1]) << " " << chi[1] << " " << sig[1]
          //      << " "
          //      << static_cast<int>(iChi[2]) << " " << chi[2] << " " << sig[2]
          //      << " "
          //      << static_cast<int>(iChi[3]) << " " << chi[3] << " " << sig[3]
          //      << endl;
        }
        // double p=0;
        // for(int irot=0; irot<nRot; irot++){
        //   p+= records[FY_NBIN][iY][ia].getRotamer(irot)->getP();
        // }
        // cout << p << endl;
      }
    }
  }

  /*C-teminal*/
  for (int ia = 0; ia < 20; ia++) {
    if (ia != GLY && ia != ALA) {
      int nRot = records[0][0][ia].getNRotamers();
      for (int iF = 0; iF < FY_NBIN; iF++) {
        records[iF][FY_NBIN][ia].init(ia, nRot);
        for (int irot = 0; irot < nRot; irot++) {
          double N_phi_chi = 0;
          int N_phi = 0;
          BBDep_rotRecord *rpt = records[iF][FY_NBIN][ia].getRotamer(irot);
          double chi[4];
          double prev_chi[4];
          double sig[4];
          for (int k = 0; k < 4; k++) {
            prev_chi[k] = -3600;
            chi[k] = sig[k] = 0;
          }
          const double *aveChi = NULL;
          const double *sigChi = NULL;
          const unsigned char *iChi = NULL;
          double p = 0;
          for (int iY = 0; iY < FY_NBIN; iY++) { // average over psi
            aveChi = records[iF][iY][ia].getRotamer(irot)->getAveChi();
            sigChi = records[iF][iY][ia].getRotamer(irot)->getSigChi();
            iChi = records[iF][iY][ia].getRotamer(irot)->getIChi();
            p = records[iF][iY][ia].getRotamer(irot)->getP();
            int n_phi_psi = records[iF][iY][ia].getNRecords();
            N_phi += n_phi_psi;
            N_phi_chi += n_phi_psi * p;
            for (int k = 0; k < 4; k++) {
              double angle = aveChi[k];
              if (prev_chi[k] < -360)
                prev_chi[k] = aveChi[k]; // initialization
              if (aveChi[k] - prev_chi[k] > 180)
                angle -= 360;
              if (aveChi[k] - prev_chi[k] < -180)
                angle += 360;

              chi[k] += n_phi_psi * p * angle;
              sig[k] += n_phi_psi * p * sigChi[k];

              if (N_phi > 0)
                prev_chi[k] = chi[k] / N_phi;
              else
                prev_chi[k] = chi[k];
            }
          }

          double p_phi_chi = p;
          if (N_phi > 0 && N_phi_chi > 0) {
            p_phi_chi = N_phi_chi / N_phi;
            for (int k = 0; k < 4; k++) {
              chi[k] /= N_phi_chi;
              sig[k] /= N_phi_chi;

              if (chi[k] > 180)
                chi[k] -= 360;
              if (chi[k] < -180)
                chi[k] += 360;
            }
          } else {
            for (int k = 0; k < 4; k++) {
              chi[k] = aveChi[k];
              sig[k] = sigChi[k];
            }
          }
          rpt->setP(p_phi_chi);
          rpt->setIChi(iChi[0], iChi[1], iChi[2], iChi[3]);
          rpt->setChi(chi[0], sig[0], chi[1], sig[1], chi[2], sig[2], chi[3],
                      sig[3]);
          // cout << static_cast<int>(iChi[0]) << " " << chi[0] << " " << sig[0]
          // << " "
          //      << static_cast<int>(iChi[1]) << " " << chi[1] << " " << sig[1]
          //      << " "
          //      << static_cast<int>(iChi[2]) << " " << chi[2] << " " << sig[2]
          //      << " "
          //      << static_cast<int>(iChi[3]) << " " << chi[3] << " " << sig[3]
          //      << endl;
        }
        // double p=0;
        // for(int irot=0; irot<nRot; irot++){
        //   p+= records[FY_NBIN][iY][ia].getRotamer(irot)->getP();
        // }
        // cout << p << endl;
      }
    }
  }
}

BBDep_rotRecord *BBDep_RotLib::newRotamer(BBDep_rotRecord *oldRotamer,
                                          const int &type, const int &newF,
                                          const int &newY) {
  BBDep_resRotLib *resLib = get_resRotLib(newF, newY, type);
  return resLib->getRotamer(oldRotamer->getIChi());
}

/**
 * Allocate memory for 2D Array.
 */
template <class T> T **initArray2D(int size) {
  T **tmp = new T *[size];
  tmp[0] = new T[size * size];
  for (int i = 1; i < size; i++) {
    tmp[i] = tmp[i - 1] + size;
  }
  return tmp;
}

void medusa_score() {
  char *paramDir;
  char *ipdb;
  int longRange, shortRange, midRange, reducedVDW;

  // opt_getE(opt.argc, opt.argv, paramDir, ipdb, longRange, shortRange,
  // midRange, reducedVDW);

  double mir = 9.0;
  // if (longRange) {
  //   mir = 9.0;
  // } else if (shortRange) {
  //   mir = 5.5;
  // } else if (midRange) {
  //   mir = 6.5;
  // }

  topology top(jnc::string_format("%s/cedutop.pro", paramDir));

  /*backbone dependent rotamer library*/
  BBDep_RotLib bbdep_rotlib;
  bbdep_rotlib.readFile(jnc::string_format("%s/bbdep02.May.lib", paramDir));
  bbdep_rotlib.updateStatisticalEnergy();

  /*backbone dependent amino acid distribution*/
  BBDep_AA bbdep_aa;
  bbdep_aa.readFile(jnc::string_format("%s/PDB30.PHI-PSI-AA-DIST", paramDir));
  bbdep_aa.updateStatisticalEnergy();

  /*nonbonded-vdw*/
  VDW vdw;
  vdw.init(jnc::string_format("%s/medupar.nonb", paramDir), top);
  //  vdw.ced2charmm();//change to CHARMM19 force field, to match EEF1
  if (reducedVDW) {
    vdw.shrink(0.95); // use a smaller VDW distance
  }
  vdw.HBond_RFix(top); // to avoid the disfavoration of HBonding(HBA-D)

  /*weight-coef of energy terms*/
  WEIGHT_COEF weight;
  if (!weight.readFile(jnc::string_format("%s/full-vdw.coef", paramDir))) {
    cerr << "can not open the weight file";
    exit(1);
  }

  /*initialize the protein*/
  protein p("TEST", ipdb, _PDB_);
  p.updateTopology(top);
  p.setResidueIndex();

  /*generalized protein*/
  gprotein gp(p);
  gp.updateTopology(top);
  //  gp.writePDB(cout);
  /*create a rotamer array*/
  BBDep_rotRecord **rotamers = new BBDep_rotRecord *[gp.size()];

  for (int i = 0; i < gp.size(); i++)
    rotamers[i] = NULL;
  /*initially, reshuffle the aa-rot sequece space*/
  double phi, psi;
  int iF, iY;
  /*update the ICHS for each residue*/
  for (int i = 0; i https://www.machinelearningplus.com/python/parallel-processing-python/iF = FY_NBIN;
    // else iF = static_cast<int>((phi/rpi-FY_MIN)/FY_BIN);
    // if(psi==INF) iY = FY_NBIN;
    // else iY = static_cast<int>((psi/rpi-FY_MIN)/FY_BIN);

    gp.getGResidue(i)->updateICHIs_ALL(bbdep_rotlib, phi, psi);

    int the_type = gp.getResidue(i)->getType();
    BBDep_resRotLib *res_rotlib = bbdep_rotlib.get_resRotLib(iF, iY, the_type);

    BBDep_rotRecord *temp = res_rotlib->getRotamer(
        gp.getResidue(i)->getNChi(), gp.getResidue(i)->getIChis());

    rotamers[i] = temp;
  }

  /*load the protein into the grids, for the calculation of VDW,SOLV,HBONDS*/
  // double mir = 9.0;
  // double mir = 5.5;
  gen_grid grid(p, gp, vdw, mir);
  gen_fine_grid fgrid(p, gp);

  double VDW_ATTR, VDW_REP, SOLV;
  double EHB_bb_bb, EHB_bb_sc, EHB_sc_sc;
  double E_F_Y_CHI;
  double E_F_Y_AA;

  // at first construct all the hydrogen bonds
  EHB_bb_bb = EHB_bb_sc = EHB_sc_sc = 0;
  for (int i = 0; i < gp.size(); i++) {
    fgrid.getHBE_RES(i, EHB_bb_bb, EHB_bb_sc, EHB_sc_sc);
  }
  // cout << EHB_bb_bb <<" " << EHB_bb_sc << " " << EHB_sc_sc << endl;

  // also calculate VDW/SOLV/ESB and fill the array. The energies of
  // VDW/SOLV/ESB is additive. this is used to save the calculation time
  double **VDWA_SS_2D = initArray2D<double>(gp.size());
  double **VDWR_SS_2D = initArray2D<double>(gp.size());
  double **SOLV_SS_2D = initArray2D<double>(gp.size());
  double **ESB_SS_2D = initArray2D<double>(gp.size());
  double **VDWA_SB_2D = initArray2D<double>(gp.size());
  double **VDWR_SB_2D = initArray2D<double>(gp.size());
  double **SOLV_SB_2D =https://www.machinelearningplus.com/python/parallel-processing-python/= 0;
    }
  }

  double *VDWA_SS_1D = new double[gp.size()];
  double *VDWR_SS_1D = new double[gp.size()];
  double *SOLV_SS_1D = new double[gp.size()];
  double *ESB_SS_1D = new double[gp.size()];
  double *VDWA_SB_1D = new double[gp.size()];
  double *VDWR_SB_1D = new double[gp.size()];
  double *SOLV_SB_1D = new double[gp.size()];
  for (int i = 0; i < gp.size(); i++) {
    VDWA_SS_1D[i] = 0;
    VDWR_SS_1D[i] = 0;
    SOLV_SS_1D[i] = 0;
    ESB_SS_1D[i] = 0;
    VDWA_SB_1D[i] = 0;
    VDWR_SB_1D[i] = 0;
    SOLV_SB_1D[i] = 0;
  }

  /*calculate the initial energy*/
  for (int i = 0; i < gp.size(); i++) {
    double vdw_a, vdw_r, solv, sb;
    int t = gp.getResidue(i)->getType();
    for (int k = 0; k < gp.size(); k++) {
      VDWA_SS_1D[k] = VDWR_SS_1D[k] = SOLV_SS_1D[k] = 0;
      VDWA_SB_1D[k] = VDWR_SB_1D[k] = SOLV_SB_1D[k] = 0;
    }
    /*
       grid.getEVS_RES_SC(i, vdw_a, vdw_r, solv, sb,
       VDWA_SS_1D, VDWR_SS_1D, SOLV_SS_1D,
       VDWA_SB_1D, VDWR_SB_1D, SOLV_SB_1D, ESB_SS_1D);
       */
    grid.getEVS_RES_SC(i, vdw_a, vdw_r, solv, VDWA_SS_1D, VDWR_SS_1D,
                       SOLV_SS_1D, VDWA_SB_1D, VDWR_SB_1D, SOLV_SB_1D);
    fgrid.getHH_EVR_RES_SC(i, vdw_r, VDWR_SS_1D, VDWR_SB_1D);

    for (int j = i + 1; j < gp.size(); j++) {
      VDWA_SS_2D[i][j] = VDWA_SS_2D[j][i] += VDWA_SS_1D[j];
      VDWR_SS_2D[i][j] = VDWR_SS_2D[j][i] += VDWR_SS_1D[j];
      SOLV_SS_2D[i][j] = SOLV_SS_2D[j][i] += SOLV_SS_1D[j];
    }

    for (int j = 0; j < gphttps://www.machinelearningplus.com/python/parallel-processing-python/_SS_1D[i];
  }

  // END OF THE CALCULATION OF INITIAL ENERGY
  double evdw_a = 0, evdw_r = 0, esolv = 0, esb = 0;
  double efy_aa = 0;
  double efy_chi = 0;
  int the_aa_type;
  double e_aa_ref = 0;
  for (int iaa = 0; iaa < gp.size(); iaa++) {
    esolv += SOLV_SS_2D[iaa][iaa];
    for (int jaa = iaa + 1; jaa < gp.size(); jaa++) {
      evdw_a += VDWA_SS_2D[iaa][jaa];
      evdw_r += VDWR_SS_2D[iaa][jaa];
      esolv += SOLV_SS_2D[iaa][jaa];
      // esb    += ESB_SS_2D[iaa][jaa];
    }
    for (int jaa = 0; jaa < gp.size(); jaa++) {
      evdw_a += VDWA_SB_2D[iaa][jaa];
      evdw_r += VDWR_SB_2D[iaa][jaa];
      esolv += SOLV_SB_2D[iaa][jaa];
      cout << iaa + 1 << " " << jaa + 1 << " " << VDWR_SB_2D[iaa][jaa] << endl;
    }
    gp.getFY(iaa, phi, psi);
    FY_double2int(phi, psi, iF, iY);
    the_aa_type = gp.getResidue(iaa)->getType();

    if (the_aa_type != GLY && the_aa_type != ALA)
      efy_chi += rotamers[iaa]->get_ESTAT();
    efy_aa += bbdep_aa.getE_STAT(iF, iY, the_aa_type);
    e_aa_ref += weight.getWeight()[CYS_W + the_aa_type];
  }
  double ehb_bb = 0, ehb_sb = 0, ehb_ss = 0;
  for (int iaa = 0; iaa < gp.size(); iaa++) {
    fgrid.getHBE_RES(iaa, ehb_bb, ehb_sb, ehb_ss);
  }
  ehb_bb /= 2.0;
  ehb_sb /= 2.0;
  ehb_ss /= 2.0;

  double the_E = evdw_a + evdw_r * weight.getWeight()[VDW_REP_W] +
                 esolv * weight.getWeight()[SOLV_EEF1_W]
                 //+ esb*weight.getWeight()[SB_W]
                 + ehb_bb * weight.getWeight()[HB_BB_BB_W] +
                 ehb_sb * weight.getWeight()[HB_BB_SC_W] +
                 ehb_ss * weight.getWeight()[HB_SC_SC_W] +
                 efy_aa * weight.getWeight()[FY_AA_W] +
                 efy_chi * weight.getWeight()[FY_AA_CHI_W] + e_aa_ref;

  cout << the_E << " " << evdw_a << " " << evdw_r << " " << esolv << " "
       << ehb_bb << " " << ehb_sb << " " << ehb_ss << " " << efy_aa << " "
       << efy_chi << " " << e_aa_ref << endl;
}


} // namespace medusa
