#pragma once

#include "jnc.h"

#include <array>
#include <map>
#include <string>
#include <vector>

namespace medusa {

using vec = std::array<double, 3>;

const static double INF = 1.0e+38;
const static double PI = 4.0 * atan(1.0);

enum {
  MS_VDW_ATTR,
  MS_VDW_REP,
  MS_SOLV_EEF1,
  MS_SB,
  MS_HB_BB_BB,
  MS_HB_BB_SC,
  MS_HB_SC_SC,
  MS_FY_AA,
  MS_FY_AA_CHI,
  MS_CYS,
  MS_MET,
  MS_PHE,
  MS_ILE,
  MS_LEU,
  MS_VAL,
  MS_TRP,
  MS_TYR,
  MS_ALA,
  MS_GLY,
  MS_THR,
  MS_SER,
  MS_GLN,
  MS_ASN,
  MS_GLU,
  MS_ASP,
  MS_HIS,
  MS_ARG,
  MS_LYS,
  MS_PRO,
  MSATER_HB,
  MS_CONS
};

using Array1d = std::vector<double>;
using Array2d = std::vector<std::vector<double>>;

struct VDWpair {
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
};

struct VDW {
  int nType = 0;
  std::vector<std::vector<VDWpair>> pairs;
  std::vector<std::vector<double>> EEF1_SOLV;
  std::vector<std::vector<double>> E_S;
};

struct BBDep_RotLib {
  //
};

struct BBDep_AA {
  //
};

struct Config {
  std::map<int, double> weights{
      {MS_VDW_ATTR, 1.0000}, {MS_VDW_REP, 0.8748},  {MS_SOLV_EEF1, 0.7302},
      {MS_SB, 0.0000},       {MS_HB_BB_BB, 2.0000}, {MS_HB_BB_SC, 1.6325},
      {MS_HB_SC_SC, 1.8357}, {MS_FY_AA, 0.6306},    {MS_FY_AA_CHI, 0.6157},
      {MS_CYS, 10.0000},     {MS_MET, -0.4924},     {MS_PHE, 2.7415},
      {MS_ILE, 0.5196},      {MS_LEU, 0.5669},      {MS_VAL, -0.0652},
      {MS_TRP, 4.3255},      {MS_TYR, 2.4628},      {MS_ALA, -1.0662},
      {MS_GLY, -1.4106},     {MS_THR, -0.8646},     {MS_SER, -1.1980},
      {MS_GLN, -0.2150},     {MS_ASN, -1.1531},     {MS_GLU, -2.1617},
      {MS_ASP, -2.4339},     {MS_HIS, 0.8818},      {MS_ARG, -0.0425},
      {MS_LYS, -1.3347},     {MS_PRO, 0.9400},      {MS_WATER_HB, 1.0},
      {MS_CONS, 1.0}};

  double mir = 9.0;
  VDW vdw;
  Topology top;
  BBDep_RotLib bbdep_rotlib;
  BBDep_AA bbdep_aa;
};

struct Env {
  Array1d vdwa_ss_1d;
  Array1d vdwa_sb_1d;
  Array1d vdwr_ss_1d;
  Array1d vdwr_sb_1d;
  Array1d solv_ss_1d;
  Array1d solv_sb_1d;
  Array2d vdwa_ss_2d;
  Array2d vdwa_sb_2d;
  Array2d vdwr_ss_2d;
  Array2d vdwr_sb_2d;
  Array2d solv_ss_2d;
  Array2d solv_sb_2d;
};

} // namespace medusa