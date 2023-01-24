#include "ms.h"
#include "jnc.h"

#include <fstream>
#include <iostream>
#include <regex>
#include <sstream>

namespace medusa {

void VDWpair_init(double iB = 0, double iA = 0, double iB1_4 = 0, double iA1_4 = 0) {
  B = iB;
  A = iA;
  B1_4 = iB1_4;
  A1_4 = iA1_4;

  if (A && B) {
    epsilon = (A * A) / (4.0 * B);
    sigma = pow(B / A, 1.0 / 6.0);
    rmin = sigma * pow(2.0, 1.0 / 6.0);

    d_cutoff = sigma * cutoff_ratio;
    slope = -24 * epsilon * (2 * pow(cutoff_ratio, -13.0) - pow(cutoff_ratio, -7.0)) / sigma;
    y_intercept = 4.0 * epsilon * (pow(cutoff_ratio, -12.0) - pow(cutoff_ratio, -6.0)) - slope * d_cutoff;
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
    slope1_4 = -24 * epsilon1_4 * (2 * pow(cutoff_ratio, -13.0) - pow(cutoff_ratio, -7.0)) / sigma1_4;
    y_intercept1_4 = 4.0 * epsilon1_4 * (pow(cutoff_ratio, -12.0) - pow(cutoff_ratio, -6.0)) - slope1_4 * d_cutoff1_4;
  } else {
    sigma1_4 = INF;
    epsilon1_4 = 0;
    rmin1_4 = INF;

    d_cutoff1_4 = INF;
    slope1_4 = INF;
    y_intercept1_4 = INF;
  }
}

void vdw_init(VDW &vdw, const std::string &file, const Topology &top) {
  using namespace std;

  nType = top.getNType();
  const vector<TopAtomType> &atm = top.getMassTable();
  // add initialization of the solvation parameters
  // //initialized to unphysical value
  EEF1_SOLV.assign(nType, vector<double>());
  for (int i = 0; i < nType; i++)
    EEF1_SOLV[i].assign(N_EEF1_t, -INF);

  // end initialization of solvation parameters
  pairs.assign(nType, vector<VDWpair>());
  for (int i = 0; i < nType; i++) {
    pairs[i].assign(nType, VDWpair());
  }

  ifstream in(file, ios::in);
  char buf[1000];

  // construct a temperary array which stores
  // status_value[0/1], eps, sigma, eps(1_4), sigma(1_4);
  // the first array stores whether the current type has
  // been assgined the eps and sigma values
  E_S.assign(nType, vector<double>());
  for (int i = 0; i < nType; i++) {
    E_S[i].assign(5, 0);
  }

  double eps, sigma, eps1_4, sigma1_4;
  double B, A, B1_4, A1_4;
  char ch_dummy[10];
  char atm_i[10];
  char atm_j[10];
  int i_dummy;
  double f_dummy;
  string line;
  while (getline(in, line)) {
    line.copy(buf, line.length());
    buf[line.length()] = '\0';
    // get rid of the comments
    for (int i = 0; i < line.length(); i++) {
      if (buf[i] == '!') {
        buf[i] = '\0';
        break;
      }
    }
    sscanf(buf, "%s", ch_dummy);
    if (!strncasecmp(ch_dummy, "NONB", 4)) { // nonbond
      i_dummy = strlen(buf);
      // take away the '(' and ')'
      for (int i = 0; i < i_dummy; i++) {
        if (buf[i] == '(' || buf[i] == ')')
          buf[i] = ' ';
      }
      // read the records
      sscanf(buf, "%s %s %s %lf %lf %lf %lf", ch_dummy, ch_dummy, atm_i, &eps, &sigma, &eps1_4, &sigma1_4);
      i_dummy = strlen(atm_i);
      int wildcard = 0;
      if (atm_i[i_dummy - 1] == '*') { // a wild card comparison
        i_dummy -= 1;
        wildcard = 1;
      }
      // search the mass type table and found the type
      bool matchflag;
      for (int i = 0; i < nType; i++) {
        if (wildcard)
          matchflag = !strncasecmp(atm[i].getName(), atm_i, i_dummy);
        else
          matchflag = !strcasecmp(atm[i].getName(), atm_i);
        if (matchflag) { // a match!
          if (E_S[i][0] == 1) {
            cerr << "previously assigned EPS, SIGMA values" << atm_i << endl;
          }
          E_S[i][0] = 1;
          E_S[i][1] = eps;
          E_S[i][2] = sigma;
          E_S[i][3] = eps1_4;
          E_S[i][4] = sigma1_4;
          // construct the table
          for (int j = 0; j < nType; j++) {
            if (E_S[j][0]) { // an already EPS_SIGMA assinged type
              // assign the VDW interaction values for I and J types
              eps = sqrt(E_S[i][1] * E_S[j][1]);
              sigma = (E_S[i][2] + E_S[j][2]) / 2.0;
              eps1_4 = sqrt(E_S[i][3] * E_S[j][3]);
              sigma1_4 = (E_S[i][4] + E_S[j][4]) / 2.0;
              f_dummy = pow(sigma, 6);
              B = 4.0 * f_dummy * eps;
              A = B * f_dummy;
              pairs[i][j].setBA(A, B);
              f_dummy = pow(sigma1_4, 6);
              B = 4.0 * f_dummy * eps1_4;
              A = B * f_dummy;
              pairs[i][j].setBA1_4(A, B);
              if (i != j)
                pairs[j][i] = pairs[i][j];
            }
          }
          if (!wildcard)
            break;
        }
      }
    } else if (!strncasecmp(ch_dummy, "SOLV", 4)) { // solvation
      double vol, dg_ref, dg_free, lambda, charmm_r;
      i_dummy = strlen(buf);
      // take away the '(' and ')'
      for (int i = 0; i < i_dummy; i++) {
        if (buf[i] == '(' || buf[i] == ')')
          buf[i] = ' ';
      }
      // read the records
      sscanf(buf, "%s %s %s %lf %lf %lf %lf %lf", ch_dummy, ch_dummy, atm_i, &vol, &dg_ref, &dg_free, &lambda,
             &charmm_r);
      i_dummy = strlen(atm_i);
      int wildcard = 0;
      if (atm_i[i_dummy - 1] == '*') { // a wild card comparison
        i_dummy -= 1;
        wildcard = 1;
      }
      // search the mass type table and found the type
      bool matchflag;
      for (int i = 0; i < nType; i++) {
        if (wildcard)
          matchflag = !strncasecmp(atm[i].getName(), atm_i, i_dummy);
        else
          matchflag = !strcasecmp(atm[i].getName(), atm_i);
        if (matchflag) { // a match!
          EEF1_SOLV[i][EEF1_VOL] = vol;
          EEF1_SOLV[i][EEF1_DG_REF] = dg_ref;
          EEF1_SOLV[i][EEF1_DG_FREE] = dg_free;
          EEF1_SOLV[i][EEF1_LAMBDA] = lambda;
          EEF1_SOLV[i][CHARMM_RADIUS] = charmm_r;
          // construct the table
          if (!wildcard)
            break;
        }
      }
    } else if (!strncasecmp(ch_dummy, "NBFIX",
                            5)) { // the specific VDW interaction
      i_dummy = strlen(buf);
      // take away the '(' and ')'
      for (int i = 0; i < i_dummy; i++) {
        if (buf[i] == '(' || buf[i] == ')')
          buf[i] = ' ';
      }
      sscanf(buf, "%s %s %s %s %s %lf %lf %lf %lf", ch_dummy, ch_dummy, atm_i, ch_dummy, atm_j, &B, &A, &B1_4, &A1_4);
      int wildcard_i = 0;
      int wildcard_j = 0;
      int natm_i = strlen(atm_i);
      // is atom i has a wild card
      if (atm_i[natm_i - 1] == '*') {
        wildcard_i = 1;
        natm_i -= 1;
      }
      int natm_j = strlen(atm_j);
      // is atom j has a wildcard
      if (atm_j[natm_j - 1] == '*') {
        wildcard_j = 1;
        natm_j -= 1;
      }
      for (int i = 0; i < nType; i++) {
        if (!strncasecmp(atm[i].getName(), atm_i, natm_i)) { // match of atom i
          // cout << "atom i " << atm[i].getName() << " " << atm_i << " " <<
          // endl;
          for (int j = 0; j < nType; j++) {
            if (!strncasecmp(atm[j].getName(), atm_j,
                             natm_j)) { // match of atom j
              // cout << "atom j " << atm[j].getName() << " " << atm_j << " " <<
              // endl;
              pairs[i][j].setBA(B, A);
              pairs[i][j].setBA1_4(B1_4, A1_4);
              pairs[j][i] = pairs[i][j];
              if (!wildcard_j)
                break;
            }
          }
          if (!wildcard_i)
            break;
        }
      }
    }
    ch_dummy[0] = '\0';
    buf[0] = '\0';
  }
  in.close();
}

void vdw_fix_hbond(Topology &top) {
  // Fix the RMin between the HBond Donar and the HBond Acceptor
  // to favor the hydrogen bonds
  // Donar list:
  // NCK, NCR, NZ, NZNQ, OZH, OW
  // Acceptor list:
  // NR, OC, OZ, OW, OZH
  using namespace std;
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
      cout << "unkown donar type " << Donar[i] << ",VDW radius adjustment fails" << endl;
      exit(1);
    }
  }
  vector<int> acceptors;
  for (int i = 0; i < NAcceptor; i++) {
    fft = top.fftid(Acceptor[i]);
    if (fft >= 0)
      acceptors.push_back(fft);
    else {
      cout << "unkown acceptor type " << Acceptor[i] << ", VDW radius adjustment fails" << endl;
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

void initialize_energies(Env &e) {
  for (int i = 0; i < gp.size(); i++) {
    double vdw_a, vdw_r, solv, sb;
    int t = gp.getResidue(i)->getType();
    for (int k = 0; k < gp.size(); k++) {
      e.vdwa_ss_1d[k] = e.vdwr_ss_1d[k] = e.solv_ss_1d[k] = 0;
      e.vdwa_sb_1d[k] = e.vdwr_sb_1d[k] = e.solv_sb_1d[k] = 0;
    }
    grid.getEVS_RES_SC(i, vdw_a, vdw_r, solv, e.vdwa_ss_1d, e.vdwr_ss_1d, e.solv_ss_1d, e.vdwa_sb_1d, e.vdwr_sb_1d,
                       e.solv_sb_1d);
    fgrid.getHH_EVR_RES_SC(i, vdw_r, e.vdwr_ss_1d, e.vdwr_sb_1d);

    for (int j = i + 1; j < gp.size(); j++) {
      e.vdwa_ss_2d[i][j] = e.vdwa_ss_2d[j][i] += e.vdwa_ss_1d[j];
      e.vdwr_ss_2d[i][j] = e.vdwr_ss_2d[j][i] += e.vdwr_ss_1d[j];
      e.solv_ss_2d[i][j] = e.solv_ss_2d[j][i] += e.solv_ss_1d[j];
    }

    for (int j = 0; j < gp.size(); j++) {
      e.vdwa_sb_2d[i][j] += e.vdwa_sb_1d[j];
      e.vdwr_sb_2d[i][j] += e.vdwr_sb_1d[j];
      e.solv_sb_2d[i][j] += e.solv_sb_1d[j];
    }
    e.solv_ss_2d[i][i] = e.solv_ss_1d[i];
  }
}

void calculate_energies(Env &e) {
  using namespace std;
  double evdw_a = 0, evdw_r = 0, esolv = 0, esb = 0;
  double efy_aa = 0;
  double efy_chi = 0;
  int the_aa_type;
  double e_aa_ref = 0;
  for (int iaa = 0; iaa < gp.size(); iaa++) {
    esolv += e.solv_ss_2d[iaa][iaa];
    for (int jaa = iaa + 1; jaa < gp.size(); jaa++) {
      evdw_a += e.vdwa_ss_2d[iaa][jaa];
      evdw_r += e.vdwr_ss_2d[iaa][jaa];
      esolv += e.solv_ss_2d[iaa][jaa];
    }
    for (int jaa = 0; jaa < gp.size(); jaa++) {
      evdw_a += e.vdwa_sb_2d[iaa][jaa];
      evdw_r += e.vdwr_sb_2d[iaa][jaa];
      esolv += e.solv_sb_2d[iaa][jaa];
      cout << iaa + 1 << " " << jaa + 1 << " " << e.vdwr_sb_2d[iaa][jaa] << endl;
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
}

void build_grid(Env &e, Conf &conf) {
  gen_grid grid(p, gp, vdw, mir);
  gen_fine_grid fgrid(p, gp);
}

void construct_hydrogen_bonds(Env &e) {
  EHB_bb_bb = EHB_bb_sc = EHB_sc_sc = 0;
  for (int i = 0; i < gp.size(); i++) {
    fgrid.getHBE_RES(i, EHB_bb_bb, EHB_bb_sc, EHB_sc_sc);
  }
}

void update_ichs(Env &e) {
  /*initially, reshuffle the aa-rot sequece space*/
  double phi, psi;
  int iF, iY;
  /*update the ICHS for each residue*/
  for (int i = 0; i < gp.size(); i++) {
    gp.getFY(i, phi, psi);
    FY_double2int(phi, psi, iF, iY);

    gp.getGResidue(i)->updateICHIs_ALL(bbdep_rotlib, phi, psi);

    int the_type = gp.getResidue(i)->getType();
    BBDep_resRotLib *res_rotlib = bbdep_rotlib.get_resRotLib(iF, iY, the_type);

    BBDep_rotRecord *temp = res_rotlib->getRotamer(gp.getResidue(i)->getNChi(), gp.getResidue(i)->getIChis());

    rotamers[i] = temp;
  }
}

void create_rotamers(Env &e) {
  BBDep_rotRecord **rotamers = new BBDep_rotRecord *[gp.size()];
  for (int i = 0; i < gp.size(); i++) {
    rotamers[i] = NULL;
  }
}

void read_protein(Env &e) {
  /*initialize the protein*/
  protein p("TEST", ipdb, _PDB_);
  p.updateTopology(top);
  p.setResidueIndex();

  /*generalized protein*/
  gprotein gp(p);
  gp.updateTopology(top);
  //  gp.writePDB(cout);
}

void BBDep_RotLib_init(BBDep_RotLib &lib, const std::string &rotFile) {
  using namespace std;
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

void BBDep_RotLib_update_statistical_energy(BBDep_RotLib &lib) {
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

void bbdep_aa_init(BBDep_AA &aa, const std::string &file) {
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

void bbdep_aa_update_statistical_energy(BBDep_AA &aa) {
  for (int iF = 0; iF < FY_NBIN; iF++) {
    for (int iY = 0; iY < FY_NBIN; iY++) {
      double n_total = 0;
      for (int ia = 0; ia < 20; ia++) {
        n_total += fy_aa[iF][iY][ia];
      }
      if (n_total > 0) {
        for (int ia = 0; ia < 20; ia++) {
          double p = static_cast<double>(fy_aa[iF][iY][ia]) / static_cast<double>(n_total);
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

std::map<int, double> foo() {
  Env env;
  Config conf;

  std::string paramDir = "";

  std::map<int, double> ms; // Medusa Score
  std::vector<int> energy_categories{MS_VDW_ATTR, MS_VDW_REP,  MS_SOLV_EEF1, MS_HB_BB_BB,
                                     MS_HB_BB_SC, MS_HB_SC_SC, MS_FY_AA,     MS_FY_AA_CHI};

  Toplogy(conf.top, jnc::string_format("%s/cedutop.pro", paramDir));

  BBDep_RotLib_init(conf.bbdep_rotlib, jnc::string_format("%s/bbdep02.May.lib", paramDir));
  BBDep_RotLib_update_statistical_energy(conf.bbdep_rotlib);

  BBDep_AA bbdep_aa;
  bbdep_aa_init(bbdep_aa, jnc::string_format("%s/PDB30.PHI-PSI-AA-DIST", paramDir));
  bbdep_aa_update_statistical_energy(bbdep_aa);

  vdw_init(conf.vdw, jnc::string_format("%s/medupar.nonb", paramDir), conf.top);
  vdw_fix_hbond(conf.top);

  read_protein(env);

  create_rotamers(env);

  update_ichs(env);

  build_grid(env, conf);

  construct_hydrogen_bonds(env);

  initialize_energies(env);

  calculate_energies(env);

  for (auto ec : energy_categories) {
    ms[ec] *= conf.weights[ec];
  }
  return ms;
}

} // namespace medusa